/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* evolve_main.c
 * Build an evolved hmm from a starting alignment.
 * Derived from build_main.c
 *
 * Fri Jun 24 15:57:41 1994
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "states.h"
#include "externs.h" 
#include "squid.h" 
#include "version.h"

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "dfhp:vwA:"

static char usage[]  = "\
Usage: hmm-evolve [-options] <hmmfile output> <alignment file>\n\
   where options are:\n\
   -d        : maximum discrimination (shepherd) model\n\
   -f        : save hmm in flat text format instead of binary\n\
   -h        : print brief help on usage\n\
   -p <scale>: set evolution scale; default 1.0\n\
   -v        : weighted maximum likelihood (Voronoi rule)\n\
   -w        : weighted maximum likelihood (Sonnhammer rule)\n\
   -A <prior>: set architectural prior to <prior> (0-1.0)\n";

static char banner[] = "\
hmm-evolve: hidden Markov model construction from alignment,\n\
            with simulated evolution";

int
main(int argc, char **argv)
{ 
  struct hmm_struc *hmm;        /* calculated hidden Markov model   */ 
  struct hmm_struc *newhmm;     /* tmp or evolved model             */
  struct shmm_s    *shmm;       /* HMM in search form               */
  struct prior_s   *prior;	/* prior HMM model configuration    */
  char **aseqs;                 /* aligned sequences                */
  AINFO  ainfo;			/* info for aseqs                   */
  int    nseqs;			/* number of aseqs                  */ 
  char  *seqfile;               /* sequence file                    */
  int    format;		/* format of sequence file          */
  char  *hmmfile; 		/* OUTPUT: learned hmm              */ 
  FILE  *hmmfp;			/* OUTPUT: pointer to open hmmfile  */
  FILE  *mxfp;                  /* ptr to subst matrix file         */
  int    idx;			/* counter for sequences            */
  double bmx[1][20][20];        /* base substitution matrices       */
  struct trace_s **tr;          /* fake tracebacks of model/aseq alignments */
  double convergence_thresh;	/* convergence threshold for shepherd rules */
  double converge_criterion;
  double weighted_average;
  float  score;
  float  tot_score;
  int    iteration;
  float  worstscore;		/* worst current score */
  float  bestscore;		/* best current score  */
  float  sqsum;			/* for calculating std. dev. of scores */
  double pamscale;		/* evolution scaling factor */
  float *wt;			/* array of sequence weights */
  float  randomseq[MAXABET];	/* random seq model, frequencies */

  float  mpri;
  int    save_binary;
  double damp_factor;		/* damping factor for shepherd rules */
  enum strategy_e { NORMAL, MINUSQ, ERIK_WEIGHT, VORONOI } model_style;
  
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
  pamscale           = 1.0;
  mpri               = 0.85;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'd': model_style  = MINUSQ;       break;
    case 'f': save_binary  = FALSE;        break;
    case 'p': pamscale     = atof(optarg); break;
    case 'v': model_style  = VORONOI;      break;
    case 'w': model_style  = ERIK_WEIGHT;  break;
    case 'A': mpri         = atof(optarg); break;

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
  DefaultRandomModel(randomseq);

  DefaultSimplePrior(&prior);

  /*********************************************** 
   * Create the model 
   ***********************************************/

  printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf(  "Training alignment:                %s\n", seqfile);
  printf(  "Number of sequences:               %d\n", nseqs);
  printf(  "Model output to:                   %s\n", hmmfile);
  printf(  "Model construction strategy:       ");
  switch (model_style) {
  case NORMAL:        puts("Max likelihood");                 break;
  case ERIK_WEIGHT:   puts("Weighted ML (Sonnhammer rule)");  break;
  case VORONOI:       puts("Weighted ML (Voronoi rule)");     break;
  case MINUSQ:        puts("Shepherd (1-q rule)");            break;
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
  
  /* Shepherd model, original (1-q) training rule
   * Iterative re-weighting.
   */
  else if (model_style == MINUSQ)
    {
      double qscore, weightsum;

      if (! Maxmodelmaker(aseqs, &ainfo, nseqs, prior, randomseq, mpri, &hmm, &tr))
	Die("Maxmodelmaker failed");
      PriorifyHMM(hmm, prior);
      Renormalize(hmm);
      shmm = AllocSearchHMM(hmm->M);

      printf("%10s  %10s  %12s\n", " iteration", "  qscore  ", "convergence");
      printf("----------  ----------  ------------\n");

      weighted_average = 0.;
      iteration = 0;
      while (1)			/* iterate until we converge */
	{
	  /* Calculate log-odds probabilities that each sequence matches the
	   * current model.
	   */
	  qscore = 0.0;
	  MakeSearchHMM(hmm, randomseq, shmm);
	  for (idx = 0; idx < nseqs; idx++) 
	    {
	      TraceScore(shmm, aseqs[idx], tr[idx], &score);
	      qscore += 1.0 / (1.0 + EXP2(score));
	      ainfo.sqinfo[idx].weight = 1.0 / (1.0 + EXP2(score));
	      ainfo.sqinfo[idx].flags |= SQINFO_WGT;
	    }
      
	  /* Normalize the weights so that they add up to the number of
	   * training sequences again.
	   */
	  weightsum = 0.;
	  for (idx = 0; idx < nseqs; idx++) 
	    weightsum += ainfo.sqinfo[idx].weight;
	  for (idx = 0; idx < nseqs; idx++)
	    ainfo.sqinfo[idx].weight = nseqs * ainfo.sqinfo[idx].weight / weightsum;
	  
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
		 iteration, qscore, converge_criterion);
	  if (converge_criterion < convergence_thresh) break;
	}
      FreeSearchHMM(shmm);
      puts("\n");
    }

  /* Evolve the HMM
   */
				/* need HMM in counts form */
  newhmm = AllocHMM(hmm->M);
  for (idx = 0; idx < nseqs; idx++)
    TraceCount(newhmm, aseqs[idx], ainfo.sqinfo[idx].weight, tr[idx]);
  FreeHMM(hmm);
  hmm = newhmm;
  PriorifyHMM(hmm, prior);

				/* get the substitution matrices */
  if ((mxfp = fopen("PAM2.smx", "r")) == NULL)
    Die("file open failed for matrix PAM2.smx");
  if (! ParseSubstitutionMatrix(mxfp, bmx[0]))
    Die("parse failed");
  fclose(mxfp);
				/* evolve the hmm */
  EvolveHMM(hmm, bmx, 1, pamscale, &newhmm);
				/* convert to probabilities */
  Renormalize(newhmm);

  /* Calculate the final average scores.
   */

  PriorifyHMM(hmm, prior);
  Renormalize(hmm);
  MakeSearchHMM(hmm, randomseq, shmm);
  TraceScore(shmm, aseqs[0], tr[0], &worstscore);
  tot_score = bestscore = worstscore;
  sqsum     = tot_score * tot_score;
  for (idx = 1; idx < nseqs; idx++)
    {
      TraceScore(shmm, aseqs[idx], tr[idx], &score);
      tot_score += score;
      sqsum     += score * score;
      if (score > bestscore)  bestscore  = score;
      if (score < worstscore) worstscore = score;
    }

  printf("\nUnevolved hidden Markov model (length %d)\n", hmm->M); 
  printf("Average score:                %10.2f bits\n", 
	 tot_score / (double) nseqs);
  printf("Minimum score:                %10.2f bits\n", worstscore);
  printf("Maximum score:                %10.2f bits\n", bestscore);
  printf("Std. deviation:               %10.2f bits\n", 
	 sqrt((sqsum - (tot_score * tot_score / (double) nseqs)) 
	      / ((double) nseqs - 1.0)));
  printf("Information content:          %10.2f bits\n",
	 HMMInfoContent(hmm, shmm));

  /* And now for the evolved model. We can free
   * the traces on this pass
   */
  MakeSearchHMM(newhmm, randomseq, shmm);
  TraceScore(shmm, aseqs[0], tr[0], &worstscore);
  tot_score = bestscore = worstscore;
  sqsum     = tot_score * tot_score;
  for (idx = 1; idx < nseqs; idx++)
    {
      TraceScore(shmm, aseqs[idx], tr[idx], &score);
      tot_score += score;
      sqsum     += score * score;
      if (score > bestscore)  bestscore  = score;
      if (score < worstscore) worstscore = score;
      FreeTrace(tr[idx]);
    }

  printf("\nEvolved hidden Markov model (length %d)\n", hmm->M); 
  printf("Average score:                %10.2f bits\n", 
	 tot_score / (double) nseqs);
  printf("Minimum score:                %10.2f bits\n", worstscore);
  printf("Maximum score:                %10.2f bits\n", bestscore);
  printf("Std. deviation:               %10.2f bits\n", 
	 sqrt((sqsum - (tot_score * tot_score / (double) nseqs)) 
	      / ((double) nseqs - 1.0)));
  printf("Information content:          %10.2f bits\n",
	 HMMInfoContent(newhmm, shmm));
  


  /* Save the HMM
   */
  if ((hmmfp = fopen(hmmfile, "wb")) == NULL) 
    Die("Failed to open %s for writing hmm result\n%s\n", hmmfile, usage); 
  if (save_binary)
    {
      if (! WriteBinaryHMM(hmmfp, newhmm))
	Die("failed to save hmm to %s\n", hmmfile);
    }
  else
    {
      if (! WriteHMM(hmmfp, newhmm))
	Die("failed to save hmm to %s\n", hmmfile); 
    }
  fclose(hmmfp);
  printf("\nEvolved HMM written to file %s\n", hmmfile);
	 
  FreeHMM(hmm);
  FreeHMM(newhmm);
  FreeSearchHMM(shmm);
  FreePrior(prior);
  FreeAlignment(aseqs, nseqs, &ainfo);
  return 0;
}


