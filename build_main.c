/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* build_main.c
 * Build an hmm from a starting alignment. 
 * 
 * Sun Dec 12 09:36:41 1993
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "squid.h" 
#include "states.h"
#include "externs.h" 
#include "version.h"

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "dhp:r:vwA:B:C:D:E:MP:RW:"

static char usage[]  = "\
Usage: hmmb [-options] <hmmfile output> <alignment file>\n\
   where options are:\n\
   -d        : maximum discrimination (shepherd) model\n\
   -h        : print brief help on usage\n\
   -p <pfile>: use prior from <pfile>\n\
   -r <rfile>: read random model from <rfile>\n\
   -v        : weighted maximum likelihood (Voronoi rule)\n\
   -w        : weighted maximum likelihood (Sonnhammer rule)\n\
\n\
   Experimental options:\n\
   -A <prior>   : set architecture prior to <prior> (0 to 1.0)\n\
   -B <pamwt>   : set the weight on a -P PAM prior (default 20.0)\n\
   -C <cthresh> : set -M, -d convergence threshold (0.0001)\n\
   -D <damp>    : set damp factor (0..1) for shepherd rules\n\
   -E <epsilon> : set learning rate for maximin rule (0.1 or less)\n\
   -M           : maximin (gradient) shepherd model\n\
   -P <pam>     : use PAM-based prior, using PAM matrix in <pam>\n\
   -R           : ragged alignment: ignore exterior gaps\n\
   -W <wtfile>  : re-save alignment (including final weights from\n\
                    -w or -d options) in <wtfile>\n";

static char banner[] = "hmmb: hidden Markov model construction from alignment";

int
main(int argc, char **argv)
{ 
  struct hmm_struc *hmm;        /* calculated hidden Markov model   */ 
  struct shmm_s    *shmm;       /* HMM in search form               */
  struct prior_s   *prior;	/* prior HMM model configuration    */
  char **aseqs;                 /* aligned sequences                */
  int    nseqs;			/* number of aseqs                  */ 
  AINFO  ainfo;			/* additional alignment information */
  char  *seqfile;               /* sequence file                    */
  int    format;		/* format of sequence file          */
  char  *hmmfile; 		/* OUTPUT: learned hmm              */ 
  FILE  *hmmfp;			/* OUTPUT: pointer to open hmmfile  */
  int    idx;			/* counter for sequences            */
  struct trace_s **tr;          /* fake tracebacks of model/aseq alignments */
  double convergence_thresh;	/* convergence threshold for shepherd rules */
  double converge_criterion;
  double weighted_average;
  float  score;
  float  tot_score;
  int    iteration;
  int    worstidx;
  float  worstscore;		/* worst current score */
  float  bestscore;		/* best current score  */
  float  sqsum;			/* for calculating std. dev. of scores */
  float  stddev;
  int    k, x;
  float *wt;			/* weights array */
  float  randomseq[MAXABET];	/* random background model */
  char  *randomfile;		

  float  mpri;			/* architecture prior, prob. of a match node */
  int    save_binary;
  int    ragged_alignment;
  char  *priorfile;
  char  *pamfile;
  float  pamwt;
  char  *wtfile;
  FILE  *wfp;
  double damp_factor;		/* damping factor for shepherd rules */
  double epsilon;               /* learning rate for maximin */
  double maximin_window;        /* 0-1, fraction of score window to update with in maximin */
  enum strategy_e { NORMAL, MINUSQ, MAXIMIN, ERIK_WEIGHT, VORONOI } model_style;
  
  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  save_binary        = TRUE;
  model_style        = NORMAL;
  damp_factor        = 0.99;
  convergence_thresh = 0.001;
  wtfile             = NULL;
  epsilon            = 0.1;
  maximin_window     = 0.1;
  priorfile          = NULL;
  pamfile            = NULL;
  ragged_alignment   = FALSE;
  mpri               = 0.85;
  pamwt              = 20.;
  randomfile         = NULL;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'd': model_style        = MINUSQ;       break;
    case 'v': model_style        = VORONOI;      break;
    case 'w': model_style        = ERIK_WEIGHT;  break;
    case 'M': model_style        = MAXIMIN;      break;

    case 'p': priorfile          = optarg;       break;
    case 'r': randomfile         = optarg;       break;

    case 'A': mpri               = atof(optarg); break;
    case 'B': pamwt              = atof(optarg); break; 
    case 'C': convergence_thresh = atof(optarg); break;
    case 'D': damp_factor        = atof(optarg); break;
    case 'E': epsilon            = atof(optarg); break;
    case 'P': pamfile            = optarg;       break;
    case 'R': ragged_alignment   = TRUE;         break;
    case 'W': wtfile             = optarg;       break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", 
	     banner, RELEASE, RELEASEDATE, usage); 
      exit(0);

    default:
      Die(usage);
    }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  
  hmmfile = argv[argc-2];
  seqfile = argv[argc-1]; 

  sre_srandom(time(0));

  /*********************************************** 
   * Get sequence data
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Alignment file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

				/* read the training seqs from file */
  if (! ReadAlignment(seqfile, format, &aseqs, &nseqs, &ainfo))
    Die("Failed to read aligned sequence file %s", seqfile);
  
  for (idx = 0; idx < nseqs; idx++)
    s2upper(aseqs[idx]);

  if (! DetermineAlphabet(aseqs, nseqs))
    Die("Failed to determine alphabet used in seqs in %s", seqfile);

  if (randomfile == NULL) DefaultRandomModel(randomseq);
  else                    ReadRandomModel(randomfile, randomseq);

  if (priorfile == NULL)  DefaultSimplePrior(&prior);
  else                    ReadPrior(priorfile, &prior);
    
  /* Ad hoc PAM-based prior, extension on a simple Dirichlet prior.
   */
  if (pamfile != NULL)
    {
      FILE *pamfp;
      int **pam;
      float scale;
      
      if (Alphabet_type != kAmino)
	Die("-P option for PAM matrix use is only valid for protein sequences");
      if (prior->strategy != PRI_SIMPLE)
	Die("You may use a custom or a PAM matrix-based prior, but not both!");
      if ((pamfp = fopen(pamfile, "r")) == NULL &&
	  (pamfp = EnvFileOpen(pamfile, "BLASTMAT")) == NULL)
	Die("Failed to open PAM scoring matrix file %s", pamfile);
      if (! ParsePAMFile(pamfp, &pam, &scale))
	Die("Failed to parse PAM scoring matrix file %s", pamfile);
      ToPAMPrior(pam, scale, pamwt, prior);

      fclose(pamfp);
      Free2DArray(pam,27);
    }

  if (ragged_alignment)
    BlockRaggedEdgedAlignment(aseqs, nseqs, ainfo.alen);


  /*********************************************** 
   * Create the model 
   ***********************************************/

  printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("Training alignment:                %s\n", seqfile);
  printf("Number of sequences:               %d\n", nseqs);
  printf("Model output to:                   %s\n", hmmfile);
  printf("Model construction strategy:       ");
  switch (model_style) {
  case NORMAL:        puts("Max likelihood");                 break;
  case ERIK_WEIGHT:   puts("Weighted ML (Sonnhammer rule)");  break;
  case VORONOI:       puts("Weighted ML (Voronoi rule)");     break;
  case MINUSQ:        puts("Max discrimination");             break;
  case MAXIMIN:       puts("Max discrimination (maximin)");   break;
  default: Die("No such strategy");
  }
  if (priorfile != NULL)
  printf("Customized prior read from:        %s\n", priorfile);
  printf("Prior strategy:                    ");
  switch(prior->strategy) {
  case PRI_SIMPLE: puts("simple Dirichlet");                break;
  case PRI_PAM:    
    puts("ad hoc PAM matrix based mixture"); 
    printf("PAM matrix used for prior:         %s\n",   pamfile); 
    printf("Weight on PAM prior:               %.1f\n", pamwt);
    break;
  case PRI_MIX:    puts("mixture Dirichlet");               break;
  case PRI_STRUCT: puts("structure-based mixtures");        break;
  default: Die("No such strategy");
  }
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (model_style == ERIK_WEIGHT)
    {
      SonnhammerWeights(aseqs, nseqs, ainfo.alen, &wt);
      for (idx = 0; idx < nseqs; idx++)
	{
	  ainfo.sqinfo[idx].weight = wt[idx];
	  ainfo.sqinfo[idx].flags |= SQINFO_WGT;
	}
      free(wt);
    }
  else if (model_style == VORONOI)
    {
      printf("Calculating Voronoi weights. Be patient... this is slow...\n");
      VoronoiWeights(aseqs, nseqs, ainfo.alen, &wt);
      printf("OK, done.\n");
      for (idx = 0; idx < nseqs; idx++)
	{
	  ainfo.sqinfo[idx].weight = wt[idx];
	  ainfo.sqinfo[idx].flags |= SQINFO_WGT;
	}
      free(wt);
    }

  /* Normal or weighted model.
   * Built by maximum likelihood model construction.
   */
  if (model_style == NORMAL ||
      model_style == ERIK_WEIGHT ||
      model_style == VORONOI)
    {
      if (! Maxmodelmaker(aseqs, &ainfo, nseqs, prior, randomseq, mpri, &hmm, &tr))
	Die("Maxmodelmaker failed");
      PriorifyHMM(hmm, prior);
      Renormalize(hmm);
    }
  
  /* Discriminative model, original (1-q) training rule
   * Iterative re-weighting.
   */
  else if (model_style == MINUSQ)
    {
      double *wgt;
      double qscore, weightsum;
      struct hmm_struc *newhmm;

      if (! Maxmodelmaker(aseqs, &ainfo, nseqs, prior, randomseq, mpri, &hmm, &tr))
	Die("Maxmodelmaker failed");
      PriorifyHMM(hmm, prior);
      Renormalize(hmm);

      wgt  = (double *) MallocOrDie (sizeof(double) * nseqs);
      shmm = AllocSearchHMM(hmm->M);

      printf("%10s  %10s  %12s\n", " iteration", "   score  ", "convergence");
      printf("----------  ----------  ------------\n");

      weighted_average = 0.;
      iteration = 0;
      while (1)			/* iterate until we converge */
	{
	  /* Calculate log-odds probabilities that each sequence matches the
	   * current model.
	   */
	  tot_score = 0.0;
	  MakeSearchHMM(hmm, randomseq, shmm);
	  for (idx = 0; idx < nseqs; idx++) 
	    {
	      TraceScore(shmm, aseqs[idx], tr[idx], &score);
	      wgt[idx] = 1.0 / (1.0 + EXP2(score)); /* REQUIRES double precision! */
	      tot_score += score;
	    }
      
	  /* Normalize the weights so that they add up to the number of
	   * training sequences again.
	   */
	  weightsum = 0.;
	  for (idx = 0; idx < nseqs; idx++) 
	    weightsum += wgt[idx];
	  qscore = weightsum;
	  if (weightsum > 0.)
	    for (idx = 0; idx < nseqs; idx++)
	      {
		ainfo.sqinfo[idx].weight = (float) nseqs * (float) (wgt[idx] / weightsum);
		ainfo.sqinfo[idx].flags |= SQINFO_WGT;
	      }
	  else
	    {
	      Warn("no discrimination effect: all scores are very high");
	      for (idx = 0; idx < nseqs; idx++)
		{
		  ainfo.sqinfo[idx].weight = 1.0;
		  ainfo.sqinfo[idx].flags |= SQINFO_WGT;
		}
	    }

	  /* Recount the sequences into a counts-based model, using
	   * new weights
	   */
	  newhmm = AllocHMM(hmm->M);
	  for (idx = 0; idx < nseqs; idx++)
	    TraceCount(newhmm, aseqs[idx], ainfo.sqinfo[idx].weight, tr[idx]);
	  PriorifyHMM(newhmm, prior);
	  Renormalize(newhmm);

	  /* Damp out instability 
	   */
	  HybridizeHMMs(newhmm, hmm, damp_factor);

	  /* Exchange newhmm for old hmm, free old
	   */
	  FreeHMM(hmm);
	  hmm = newhmm;
	  
	  /* Check for convergence. 
	   * This rule checks the current negative log-odds score against
	   * a rolling average of the previous scores. The rolling average
           * damps out slight oscillations which appear as we approach a
           * solution.
	   */
	  qscore = -LOG2(qscore);
	  weighted_average = (9. * weighted_average  + qscore) / 10.;
	  converge_criterion = fabs((weighted_average - qscore) / weighted_average);
	  iteration++;
	  printf("   %4d     %10.2f  %12g\n", 
		 iteration, tot_score / (float) nseqs, converge_criterion);
	  if (converge_criterion < convergence_thresh) break;
	}
      puts("\n");
      free(wgt);
      FreeSearchHMM(shmm);
    }

  /* Discriminative model, using Graeme's maximin training rule
   */
  else if (model_style == MAXIMIN)
    {
      struct hmm_struc *chmm;       /* HMM in counts form */
      
      if (! Maxmodelmaker(aseqs, &ainfo, nseqs, prior, randomseq, mpri, &chmm, &tr))
	Die("Maxmodelmaker failed");
      hmm = HMMDup(chmm);
      PriorifyHMM(hmm, prior);
      Renormalize(hmm);

      printf("%10s  %10s  %10s  %10s  %12s\n", 
	     " iteration", " avg score", " min score", 
	     "  min seq ", "convergence");
      printf("----------  ----------  ----------  ----------  ------------\n");

      shmm = AllocSearchHMM(hmm->M);
      weighted_average = 0.;
      iteration = 0;
      while (1)			/* iterate until we converge */
	{
	  /* Calculate log-odds probabilities that each sequence matches the
	   * current model. Find the worst scoring sequence.
	   */
	  worstscore = FLT_MAX;
	  bestscore  = -FLT_MAX;
	  tot_score  = 0.;
	  MakeSearchHMM(hmm, randomseq, shmm);
	  for (idx = 0; idx < nseqs; idx++) 
	    {
	      TraceScore(shmm, aseqs[idx], tr[idx], &score);
	      if (score < worstscore) { worstscore = score; worstidx = idx; }
	      if (score > bestscore)  bestscore  = score; 
	      tot_score += score;
	    }
	  tot_score /= (double) nseqs;
	  
	  /* Update the counts in the master counts-based model
           * c' = (1-epsilon) c + epsilon * new counts
	   */
				/* reduce counts by 1-epsilon */
	  for (k = 0; k <= chmm->M; k++)
	    {
	      for (x = 0; x < Alphabet_size; x++)
		{ chmm->mat[k].p[x] *= (1.0 - epsilon);
		  chmm->ins[k].p[x] *= (1.0 - epsilon); }
	      for (x = 0; x < 3; x++)
		{
		  chmm->mat[k].t[x] *= (1.0 - epsilon);
		  chmm->ins[k].t[x] *= (1.0 - epsilon);
		  chmm->del[k].t[x] *= (1.0 - epsilon);
		}
	    }
				/* bump counts with low-scoring sequences */
	  for (idx = 0; idx < nseqs; idx++) 
	    {
	      TraceScore(shmm, aseqs[idx], tr[idx], &score);
	      if (score < worstscore + (bestscore - worstscore) * maximin_window)
		TraceCount(chmm, aseqs[idx], epsilon, tr[idx]);
	    }
	  
				/* create new probability-based HMM */
	  FreeHMM(hmm);	  
	  hmm = HMMDup(chmm);
	  PriorifyHMM(hmm, prior);
	  Renormalize(hmm);
	  
	  /* Check for convergence. 
	   * We use a smoothed average of the worst score.
	   */
	  weighted_average = (9. * weighted_average  + worstscore) / 10.;
	  converge_criterion = fabs((weighted_average - worstscore) / weighted_average);
	  iteration++;
	  printf("   %4d     %10.2f  %10.2f  %10s  %f\n", 
		 iteration, tot_score, worstscore, ainfo.sqinfo[worstidx].name, 
		 converge_criterion);
	  if (converge_criterion < convergence_thresh) break;  
	  epsilon *= 0.99;
	}                    
      FreeSearchHMM(shmm);
      puts("\n");
    }

  /* Calculate the final average scores.
   */

  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);
  TraceScore(shmm, aseqs[0], tr[0], &worstscore);
  tot_score = bestscore = worstscore;
  sqsum     = tot_score * tot_score;
  FreeTrace(tr[0]);
  for (idx = 1; idx < nseqs; idx++)
    {
      TraceScore(shmm, aseqs[idx], tr[idx], &score);
      tot_score += score;
      sqsum     += score * score;
      if (score > bestscore)  bestscore  = score;
      if (score < worstscore) worstscore = score;
      FreeTrace(tr[idx]);
    }
  free(tr);

  printf("\nConstructed a hidden Markov model (length %d)\n", hmm->M); 
  printf("Average score:                %10.2f bits\n", 
	 tot_score / (double) nseqs);
  printf("Minimum score:                %10.2f bits\n", worstscore);
  printf("Maximum score:                %10.2f bits\n", bestscore);

				/* Need to protect against small numerical
				   precision problems that feed negative arg to sqrt() */
  if (nseqs > 1)
    {
      stddev = (sqsum - (tot_score * tot_score / (double) nseqs)) / ((double) nseqs - 1.0);
      stddev = (stddev > 0) ? sqrt(stddev) : 0.0;
    }
  else stddev = 0.0;
  printf("Std. deviation:               %10.2f bits\n", stddev);
  printf("Information content:          %10.2f bits\n",
	 HMMInfoContent(hmm, shmm));
  FreeSearchHMM(shmm);

  if ((hmmfp = fopen(hmmfile, "wb")) == NULL) 
    Die("Failed to open %s for writing hmm result\n%s\n", hmmfile, usage); 
  if (save_binary)
    {
      if (! WriteBinaryHMM(hmmfp, hmm))
	Die("failed to save hmm to %s\n", hmmfile);
    }
  else
    if (! WriteHMM(hmmfp, hmm))
      Die("failed to save hmm to %s\n", hmmfile); 
  fclose(hmmfp);
  printf("\nHMM written to file %s\n", hmmfile);
	 
  /* Save the weights?
   */
  if (wtfile != NULL)
    {
      if ((wfp = fopen(wtfile,"w")) == NULL)
	Warn("Failed to open weight save file %s", wtfile);
      else
	{
	  WriteSELEX(wfp, aseqs, nseqs, &ainfo, 50);
	  fclose(wfp);
	  printf("Final weights saved as alignment file %s\n", wtfile);
	}
    }

  FreeHMM(hmm);
  FreePrior(prior);
  FreeAlignment(aseqs, nseqs, &ainfo);
  return 0;
}



