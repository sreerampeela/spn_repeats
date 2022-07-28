/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* train_main.c
 * 
 * main() for training an hmm from unaligned sequences.
 * For best results from unaligned sequence, use simulated annealing
 * 
 * Usage: hmmt [-options] <hmmfile out> <seqfile in>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <limits.h>
#include <float.h>
#ifdef NEED_GETOPTH
#include <getopt.h>
#endif


#include "states.h"
#include "externs.h" 
#include "squid.h" 
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char Alphabet[MAXABET];		/* ACGT, for instance       */
int  Alphabet_size;		/* 4 or 20                  */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino    */

static int killnow;		/* flag for SIGINT catching */
static void handle_interrupt(int sig);

#define OPTIONS "a:hi:k:l:o:p:r:s:vA:BP:S:W"

static char usage[]  = "\
Usage: hmmt [-options] <hmmfile output> <seqfile in>\n\
where options are:\n\
   -a <alignfile>    : include (and keep fixed) known alignment in <alignfile>\n\
   -h                : (help) print version and usage info\n\
   -i <hmmfile>      : use HMM in <hmmfile> as starting model\n\
   -k <kT>           : set starting kT for sim annealing run (default 5.0)\n\
   -o <out_afile>    : save final alignment to <out_afile>\n\
   -p <pfile>        : read custom prior from <pfile>\n\
   -r <ramp>         : set multiplier [0 to 1] for sim annealing (default 0.95)\n\
   -s <seed>         : set random number generator seed to <seed>\n\
   -l <length>       : in conjunction with a flat model, set length to l\n\
   -v                : use Viterbi training (default: sim annealing)\n\
 Experimental options:\n\
   -A <prior>        : set architectural prior to <prior> (0-1.0)\n\
   -B                : use full Baum-Welch EM training (default: sim annealing)\n\
   -P <pam>          : use PAM-based ad hoc prior, using matrix in <pam>\n\
   -S <schedulefile> : read simulated annealing schedule from <file>\n\
   -W                : re-weight by the Sonnhammer rule at each iteration\n";

static char banner[] = "hmmt: hidden Markov model training";

static struct hmm_struc *start_from_scratch(char **seqs, int N, int setlength, 
					    struct prior_s *prior);
static struct hmm_struc *start_from_hmm(char *hmmfile);
static int    read_kT(FILE *schedfp, float *ret_kT);

int
main(int argc, char **argv)
{ 
  char   *seqfile;		/* unaligned sequence file                     */
  char  **seqs;			/* training sequences                          */
  char  **seqs1;		/* 1..seqlen sequences, prepped for alignment  */
  SQINFO *sqinfo;		/* array of info for seqs                      */
  int     nseqs;		/* number of unaligned training sequences      */
  char   *fixfile;		/* aligned sequence file (optional)            */
  char  **fix_aseqs;		/* fixed alignment                             */
  AINFO   fix_ainfo;		/* alignment info on fixed seqs                */
  int     fix_nseqs;		/* number of fixed aligned training sequences  */
  struct trace_s  **fixtr;     	/* tracebacks of fixed alignment               */  
  float  *wt;		        /* array of sequence weights */
  struct prior_s   *prior;	/* prior HMM model configuration  */
  struct hmm_struc *hmm;	/* calculated hidden Markov model */ 
  struct hmm_struc *newhmm;	/* re-estimated hmm */ 
  struct shmm_s    *shmm;	/* current model in search form */
  struct vit_s    **mx;		/* Viterbi calculation grid */
  struct trace_s  **tr;		/* tracebacks of unaligned sequences */

  char **aseqs;			/* aligned sequences */

  int    format;		/* format of a sequence file */
  char  *hmmfile;		/* OUTPUT: learned hmm */ 
  FILE  *hmmfp;			/* OUTPUT: pointer to open hmmfile */
  float  score;
  float  tot_score;
  int    idx;
  int    iteration_count; 
  int    maxiter;		/* protection against oscillation */
  AINFO  ainfo;
  float  hmmdistance;		/* measure of parameter change between two models */
  float  randomseq[MAXABET];	/* random sequence model, symbol frequencies */
  enum strategy_e { ANNEALING, BAUMWELCH, VITERBI } strategy;

  double distance_threshold;	/* convergence criterion */
  float  kT;			/* Boltzmann temp factor for simulated annealing */
  float  sa_ramp;               /* multiplier for simulated annealing schedule   */
  char  *inithmmfile;		/* hmm file to read starting model from          */
  int    seed;			/* seed for random number generator              */
  int    setlength;		/* user-preferred length for model               */
  int    do_surgery;		/* use max likelihood model re-estimation        */
  char  *priorfile;		/* file to get customized prior distributions from */
  float  mpri;			/* architectural prior                           */
  char  *schedulefile;          /* customized simulated annealing schedule       */
  FILE  *schedfp;		/* open schedule file                            */
  char  *out_afile;             /* output alignment file                         */
  int    weighted;		/* TRUE to do weight compensation during training*/
  char  *pamfile;		/* name of PAM file to read for PAM-based prior  */
  float  pamwt;			/* weight on PAM-based prior                     */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

				/* set defaults on training parameter configuration */
  distance_threshold  = 0.001;
  inithmmfile         = NULL;
  fixfile             = NULL;
  kT                  = 5.0;	
  sa_ramp             = 0.95;    
  seed                = (int) time ((time_t *) NULL); /* default: "random" seed */
  setlength           = -1;
  do_surgery          = TRUE;
  priorfile           = NULL;
  mpri                = 0.85;
  schedulefile        = NULL;
  strategy            = ANNEALING;
  out_afile           = NULL;
  maxiter             = 500;
  weighted            = FALSE;
  killnow             = FALSE;	/* SIGINT, SIGTERM handling */
  pamfile             = NULL;
  pamwt               = 20.;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'a': fixfile       = optarg;               break;
    case 'i': inithmmfile   = optarg;               break;
    case 'k': kT            = (float) atof(optarg); break;
    case 'l': setlength     = atoi(optarg);         break;
    case 'o': out_afile     = optarg;               break;
    case 'p': priorfile     = optarg;               break;
    case 'r': sa_ramp       = (float) atof(optarg); break;
    case 's': seed          = atoi(optarg);         break;
    case 'v': strategy      = VITERBI;              break;
    case 'A': mpri          = atof(optarg);         break;
    case 'B': strategy      = BAUMWELCH;            break;
    case 'P': pamfile       = optarg;               break;
    case 'S': schedulefile  = optarg;               break;
    case 'W': weighted      = TRUE;                 break;

    case 'h': 
      printf("%s\n version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);
    default:
      Die("Error: unrecognized option %c\n", optc);
    }

  /* Check the options for incompatibilities
   */
  if (strategy == ANNEALING) do_surgery = FALSE;
  if (setlength != -1)       do_surgery = FALSE;      
  if (inithmmfile != NULL && fixfile != NULL)
    Die("Choose only one method for generating start model (-a/-i)");
  if (schedulefile != NULL && strategy != ANNEALING)
    Die("Why are you providing an annealing schedule (-S), if you want -v or -B?"); 

  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);
  
  hmmfile = argv[optind++];
  seqfile = argv[optind]; 
				/* initialize random() */
  sre_srandom(seed);
				/* open SA schedule file */
  if (schedulefile != NULL) {
    if ((schedfp = fopen(schedulefile,"r")) == NULL) 
      Die("Failed to open SA schedule file %s", schedulefile);
    if (! read_kT(schedfp, &kT))
      Die("SA schedule file %s looks empty", schedulefile);
  }

  /* Catch signals gracefully, allowing users to kill runs without 
   * losing models. SIGINT is an interactive interrupt (e.g. ctrl-C).
   * SIGTERM is a 'kill'.
   */
  if (signal(SIGINT,  SIG_IGN) != SIG_IGN) signal(SIGINT,  handle_interrupt);
  if (signal(SIGTERM, SIG_IGN) != SIG_IGN) signal(SIGTERM, handle_interrupt); 

  /*********************************************** 
   * Get unaligned sequence data
   ***********************************************/

  /* The unaligned sequences first
   */
  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }
				/* read the unaligned training seqs from file */
  if (! ReadMultipleRseqs(seqfile, format, &seqs, &sqinfo, &nseqs))
    Die("Failed to read sequence file %s", seqfile);

  if (! DetermineAlphabet(seqs, nseqs))
    Die("Failed to determine alphabet used in seqs in %s", seqfile);
  DefaultRandomModel(randomseq);

  if (priorfile == NULL) DefaultSimplePrior(&prior);
  else                   ReadPrior(priorfile, &prior);

  /* Ad hoc PAM-based prior, extension on a simple Dirichlet prior.
   * Protein only
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

  /*********************************************** 
   * Get fixed alignment data, 
   * and set up starting model from it,
   * if option -a is set.
   ***********************************************/
  
  fix_nseqs = 0;
  if (fixfile != NULL)
    {
      char  **fix_rseqs;		/* raw (unaligned) seqs derived from fix_aseqs */
      char  **all_rseqs;		/* dealigned fixed alignment + unaligned seqs  */
      SQINFO *all_sqinfo;		/* info for all sequences                      */

      if (! SeqfileFormat(fixfile, &format, NULL))
	switch (squid_errno) {
	case SQERR_NOFILE: Die("Alignment file %s could not be opened for reading", fixfile);
	case SQERR_FORMAT: 
	default:           Die("Failed to determine format of alignment file %s", fixfile);
	}
      if (! ReadAlignment(fixfile, format, &fix_aseqs, &fix_nseqs, &fix_ainfo))
	Die("Failed to read fixed alignment file %s", seqfile);
      DealignAseqs(fix_aseqs, fix_nseqs, &fix_rseqs);

      /* Construct the starting model from that alignment
       */
      for (idx = 0; idx < fix_nseqs; idx++)
	s2upper(fix_aseqs[idx]);
      if (! Maxmodelmaker(fix_aseqs, &fix_ainfo, fix_nseqs, prior, randomseq, mpri, &hmm, NULL))
	Die("Maxmodelmaker failed");
      PriorifyHMM(hmm, prior);
      Renormalize(hmm);

      /* Set up all_ stuff, which contains fixed alignment sequences
       * first, then unaligned seqs.
       */ 
      all_rseqs  = (char **)  MallocOrDie ((fix_nseqs+nseqs) * (sizeof(char *)));
      all_sqinfo = (SQINFO *) MallocOrDie ((fix_nseqs+nseqs) * (sizeof(SQINFO)));
      for (idx = 0; idx < fix_nseqs; idx++)
	{
	  all_rseqs[idx] = fix_rseqs[idx];
	  SeqinfoCopy(&(all_sqinfo[idx]), &(fix_ainfo.sqinfo[idx]));
	}
      free(fix_rseqs);
      for (idx = 0; idx < nseqs; idx++)
	{
	  all_rseqs[idx+fix_nseqs] = seqs[idx];
	  SeqinfoCopy(&(all_sqinfo[idx+fix_nseqs]), &(sqinfo[idx]));
	  FreeSequence(NULL, &(sqinfo[idx]));
	}
      free(seqs);
      free(sqinfo);
      seqs   = all_rseqs;
      sqinfo = all_sqinfo;
      nseqs  = fix_nseqs + nseqs;
    }
  
  /* Set up per-sequence weights.
   */
  tr    = (struct trace_s **) MallocOrDie (nseqs * sizeof(struct trace_s *));
  wt    = (float *)           MallocOrDie (nseqs * sizeof(float));
  seqs1 = (char **)           MallocOrDie (nseqs * sizeof(char *));
  for (idx = 0; idx < nseqs; idx++)
    {
      s2upper(seqs[idx]);
      PrepareSequence(seqs[idx], &seqs1[idx], NULL);
      wt[idx] = (sqinfo[idx].flags & SQINFO_WGT) ? sqinfo[idx].weight : 1.0;
    }

  /*********************************************** 
   * Create the initial model, if we haven't already.
   ***********************************************/

  if (inithmmfile != NULL)
    hmm = start_from_hmm(inithmmfile);
  else if (fixfile == NULL)
    hmm = start_from_scratch(seqs, nseqs, setlength, prior);

  if (hmm == NULL)
    Die("Error: failed to create a starting model\n");

  /*********************************************** 
   * Ready to start. Print banner information
   ***********************************************/

  printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);
  printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  printf(  "Training sequences:     %s\n", seqfile);
  if (fixfile != NULL)
    printf("Aligned training seqs:  %s\n", fixfile);
  printf(  "Total no. of sequences: %d ", nseqs);
  if (fixfile != NULL)
    printf("[%d unaligned, %d aligned]\n", nseqs-fix_nseqs, fix_nseqs);                   
  else
    puts("");
  printf(  "Model output to:        %s\n", hmmfile);
  if (out_afile != NULL)
    printf("Alignment output to:    %s\n", out_afile);
  puts("");

  if (fixfile != NULL)
    printf("Initial model:          from alignment in %s (length %d)\n", fixfile, hmm->M);
  else if (inithmmfile != NULL)
    printf("Initial model:          from HMM in file %s (length %d)\n", inithmmfile, hmm->M);
  else
    printf("Initial model:          flat (length %d%s)\n", hmm->M, (setlength != -1) ? " - forced" : "");
    
  if (priorfile != NULL)
    printf("Prior:                  custom, from %s\n", priorfile);
  else
    printf("Prior:                  default (plus-one)\n");

  printf(  "Seed for random():      %d\n", seed);

  switch (strategy) {
  case ANNEALING:
    printf("Strategy:               Simulated annealing\n");
    if (schedulefile != NULL)
      printf("Annealing schedule:     customized schedule %s\n", schedulefile);
    else 
      printf("Annealing schedule:     starting kT = %.2f   ramp = %.2f\n", kT, sa_ramp);
    break;
  case VITERBI:
    printf("Strategy:               Viterbi approximation to Baum-Welch EM\n");
    break;
  case BAUMWELCH:
    printf("Strategy:               Baum-Welch expectation maximization\n");
    break;
  }
  if (weighted)
    printf("Reweighting each iteration, Sonnhammer rule\n");

  printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");


  /*********************************************** 
   * Main training loop 
   ***********************************************/
  printf("%s%s%s%s%s\n",
	 "Iteration   ",
	 "   length   ",
	 "   change   ",
	 "   score    ",
	 (strategy == ANNEALING) ? "     kT     " : "");

  printf("%s%s%s%s%s\n",
	 "---------   ",
	 "   ------   ",
	 "   ------   ",
	 "   -----    ",
	 (strategy == ANNEALING) ? "     --     " : "");
	 
  
  shmm = AllocSearchHMM(hmm->M);
  iteration_count = 1; 
  while (iteration_count <= maxiter && !killnow)
    {
      tot_score = 0.0;
      MakeSearchHMM(hmm, randomseq, shmm);
      
      /* Stage one -- expectation.
       * Align the model to each of the training sequences.
       * Collect tracebacks (Viterbi, SA), or full forward-backward 
       * matrices (Baum-Welch) 
       */
				/* align alignment by Viterbi */
      if (fix_nseqs > 0)
	{
	  ViterbiAlignAlignment(shmm, fix_aseqs, fix_ainfo.alen, fix_nseqs, &fixtr, &score);
	  tot_score += score;
	  for (idx = 0; idx < fix_nseqs; idx++)
	    {
	      tr[idx] = fixtr[idx];
	      DealignTrace(tr[idx], fix_aseqs[idx], fix_ainfo.alen);
	    }
	  free(fixtr);
	}

      if (strategy == VITERBI)
	{
				/* align unaligned sequences */
	  for (idx = fix_nseqs; idx < nseqs; idx++)
	    {
	      if (AlignmentTooBig(sqinfo[idx].len, shmm->M))
		WeeViterbi(shmm, seqs1[idx], sqinfo[idx].len, FALSE, 0., 0., &tr[idx], &score);
	      else
		{
		  if (!ViterbiFill(shmm, seqs1[idx], sqinfo[idx].len, &mx, &score))
		    Die("internal error: ViterbiFill() has failed and aborted\n");
		  if (! ViterbiTrace(mx, shmm, seqs1[idx], sqinfo[idx].len+2, sqinfo[idx].len+1,
				     shmm->M+1, &tr[idx], (int *) NULL, (int *) NULL))
		    Die("error: ViterbiTrace() has failed and aborted\n");
		  FreeViterbiMatrix(mx, sqinfo[idx].len);
		}
	      tot_score += score;
	    }
	}

      else if (strategy == ANNEALING)
	{
	  struct sa_s **sa_mx;
	  struct sa_hmm_s *sahmm;
	  float sascore;

				/* align unaligned seqs */
	  if (( sahmm = CreateSahmm(hmm, kT)) == NULL)
	    Die("internal error: failed to create an SA hmm\n");
	  for (idx = fix_nseqs; idx < nseqs; idx++)
	    {
	      if (! SaFill(sahmm, seqs[idx], &sa_mx))
		Die("internal error: SaFill() has failed and aborted\n");

	      if (! SaTrace(sa_mx, sqinfo[idx].len, sahmm, &tr[idx]))
		{ /* Trace crashes on persistent numerical precision problems. */	
		  Warn("SA trace crash for %s. I'll work around it, don't worry.", sqinfo[idx].name);
		  if (! ViterbiFill(shmm, seqs1[idx], sqinfo[idx].len, &mx, &score))
		    Die("internal error: ViterbiFill() has failed and aborted\n");
		  if (! ViterbiTrace(mx, shmm, seqs1[idx], sqinfo[idx].len+2, sqinfo[idx].len+1,
				     shmm->M+1, &tr[idx], (int *) NULL, (int *) NULL))
		    Die("error: ViterbiTrace() has failed and aborted\n");
		  FreeViterbiMatrix(mx, sqinfo[idx].len);
		}

	      TraceScore(shmm, seqs[idx], tr[idx], &sascore);
	      tot_score += (float) sascore;

	      Free2DArray(sa_mx, sqinfo[idx].len+2);
	    }
	  DestroySahmm(sahmm);
	}

      else if (strategy == BAUMWELCH)
	{
	  struct forback_s **fwd;
	  struct forback_s **bck;
	  float             *scl;
	  
				/* allocate for a counts-based HMM */
	  if ((newhmm = AllocHMM(hmm->M)) == NULL)
	    Die("Failed to allocate for a new HMM\n");
	  
				/* count in the fixed alignment (only an approximation) */
	  for (idx = 0; idx < fix_nseqs; idx++)
	    {
	      TraceCount(newhmm, seqs[idx], wt[idx], tr[idx]);
	      FreeTrace(tr[idx]);
	    }

				/* forward-backward on the unaligned seqs */
	  for (idx = fix_nseqs; idx < nseqs; idx++)
	    {
	      Forward(hmm, seqs[idx], &fwd, &scl);
	      Backward(hmm, seqs[idx], scl, &bck);
	      score = ForwardScore(fwd, sqinfo[idx].len, hmm->M, scl);
	      score -= RandomScore(randomseq, seqs[idx]);
	      ForbackCount(hmm, seqs[idx], sqinfo[idx].len, wt[idx], fwd, bck, scl, newhmm);
	      tot_score += score;
	      
	      free(scl);
	      Free2DArray(fwd, sqinfo[idx].len+2);
	      Free2DArray(bck, sqinfo[idx].len+2);
	    }
	}	  
	  
      /* Stage two -- maximization.
       * Re-estimate a new HMM.
       * Either use tracebacks to create a completely new architecture;
       * or count tracebacks into a new model;
       * or probify an existing counts-based model from Baum-Welch.
       */
      if ((strategy == ANNEALING || strategy == VITERBI) && do_surgery)
	{
	  if (! Traces2Alignment(seqs, sqinfo, nseqs, hmm->M, tr, 0, &aseqs, &ainfo))
	    Die("Traces2Alignment() failed");

	  for (idx = 0; idx < nseqs; idx++)
	    s2upper(aseqs[idx]);

	  if (weighted)
	    {
	      free(wt);
	      SonnhammerWeights(aseqs, nseqs, ainfo.alen, &wt);
	    }

	  if (! Maxmodelmaker(aseqs, &ainfo, nseqs, prior, randomseq, mpri, &newhmm, NULL))
	    Die("Maxmodelmaker() failed");
	  FreeAlignment(aseqs, nseqs, &ainfo);
				/* free old traces */
	  for (idx = 0; idx < nseqs; idx++)
	    FreeTrace(tr[idx]);
	}

      else if (strategy == ANNEALING || strategy == VITERBI)
	{
	  if ((newhmm = AllocHMM(hmm->M)) == NULL)
	    Die("Failed to allocate for a new HMM\n");
	  for (idx = 0; idx < nseqs; idx++)
	    TraceCount(newhmm, seqs[idx], wt[idx], tr[idx]);
	  
	  /* If we're doing structural priors, we still need a new
	   * alignment, so we can re-calculate structural annotation
	   */
	  if (prior->strategy == PRI_STRUCT || weighted)
	    {
	      int    spos, col;
	      float *inp;

	      if (! Traces2Alignment(seqs, sqinfo, nseqs, hmm->M, tr, 0, &aseqs, &ainfo))
		Die("Traces2Alignment() failed");

	      if (weighted)
		{
		  free(wt);
		  SonnhammerWeights(aseqs, nseqs, ainfo.alen, &wt);
		}

	      if (prior->strategy == PRI_STRUCT)
		{
		  AnnotateAlignment(aseqs, nseqs, &ainfo, &inp);
		  memcpy(newhmm->xray, inp, sizeof(float) * NINPUTS);  
		  for (spos = 1, col = 0; col < ainfo.alen; col++)
		    if (!isgap(ainfo.rf[col]))
		      {
			memcpy(newhmm->xray + spos*NINPUTS, inp + col*NINPUTS, sizeof(float)*NINPUTS);
			spos++;
		      }
		  newhmm->flags |= HMM_XRAY;
		  free(inp);
		}
	      FreeAlignment(aseqs, nseqs, &ainfo);
	    }

	  for (idx = 0; idx < nseqs; idx++)
	    FreeTrace(tr[idx]);	
	}

      PriorifyHMM(newhmm, prior);
      Renormalize(newhmm);
#ifdef DEBUG
      VerifyHMM(newhmm);
#endif
      
      if (newhmm->M != hmm->M)
	{
	  FreeSearchHMM(shmm);
	  shmm = AllocSearchHMM(newhmm->M);
	  hmmdistance = FLT_MAX;
	}
      else
	hmmdistance = HMMDistance(newhmm, hmm);

      score = (tot_score / (float) nseqs);
      
      printf(" %4d ::    ", iteration_count); /* what iteration               */
      printf("   %6d   ", newhmm->M);         /* new model length             */
      if (hmm->M == newhmm->M)
	printf("%9.4f   ", hmmdistance);      /* change from old model        */
      else
	printf("     **     ");	              /* (architecture changed)       */
      printf("%8.2f    ", score);             /* average score of old model   */
      if (strategy == ANNEALING) printf("%8.2f    ", kT);  /* annealing only: current temp */
      puts("");

      /* Switch new model to current model
       */
      FreeHMM(hmm); 
      hmm = newhmm; 

      /* Check for convergence.
       */
      if (hmmdistance < distance_threshold)
	{ 
	  if (strategy == ANNEALING)		/* switch to Viterbi before leaving SA */
	    {
	      printf("**** annealed model has converged - switching to Viterbi (kT=0) ****\n");
	      if (setlength == -1) do_surgery = TRUE;	/* switch to max likelihood construction */
	      strategy = VITERBI;
	    }
	  else			/* else, if we converged during Viterbi or Baum-Welch, we're done */
	    break;
	}

      /* Revise annealing parameters
       */
      if (strategy == ANNEALING && schedulefile != NULL) 
	{
	  if (! read_kT(schedfp, &kT)) 
	    {
	      printf("**** end of custom SA schedule  - switching to Viterbi ****\n");
	      if (setlength == -1) do_surgery = TRUE;	
	      strategy = VITERBI;
	      fclose(schedfp);
	    }
	}
      else if (strategy == ANNEALING) 
	{
	  kT *= sa_ramp; 
	  if (kT < 0.10) 
	    {
	      printf("**** kT below threshold - switching to Viterbi ****\n");
	      if (setlength == -1) do_surgery = TRUE;	
	      strategy = VITERBI;
	    }
	}
      iteration_count++; 
    }

  /************************************************** 
   * Save the model 
   **************************************************/

  printf("\nFound a model (length %d) after %d iterations:\n", hmm->M, iteration_count); 

  if ((hmmfp = fopen(hmmfile, "wb")) == NULL) 
    Die("Failed to open %s for writing hmm result\n%s\n", hmmfile, usage); 
  if (! WriteBinaryHMM(hmmfp, hmm))
    Die("Failed to save hmm\n");
  fclose(hmmfp);
  printf("Model written to file %s\n", hmmfile);

  /************************************************** 
   * Make and save a final alignment, if -o is set
   **************************************************/

  if (out_afile != NULL)
    {
      FILE *ofp;

      MakeSearchHMM(hmm, randomseq, shmm);
				/* align alignment */
      if (fix_nseqs > 0)
	{
	  ViterbiAlignAlignment(shmm, fix_aseqs, fix_ainfo.alen, 
				fix_nseqs, &fixtr, &score);
	  for (idx = 0; idx < fix_nseqs; idx++)
	    {
	      tr[idx] = fixtr[idx];
	      DealignTrace(tr[idx], fix_aseqs[idx], fix_ainfo.alen);
	    }
	  free(fixtr);
	}
				/* align unaligned sequences */
      for (idx = fix_nseqs; idx < nseqs; idx++)
	{
	  if (AlignmentTooBig(sqinfo[idx].len, shmm->M))
	    WeeViterbi(shmm, seqs1[idx], sqinfo[idx].len, FALSE, 
		       0., 0., &tr[idx], &score);
	  else
	    {
	      if (! ViterbiFill(shmm, seqs1[idx], sqinfo[idx].len, &mx, &score))
		Die("internal error: ViterbiFill() has failed and aborted\n");
	      if (! ViterbiTrace(mx, shmm, seqs1[idx], sqinfo[idx].len+2, sqinfo[idx].len+1,
				 shmm->M+1, &tr[idx], (int *) NULL, (int *) NULL))
		Die("error: ViterbiTrace() has failed and aborted\n");
	      FreeViterbiMatrix(mx, sqinfo[idx].len);
	    }
	}

      if (! Traces2Alignment(seqs, sqinfo, nseqs, hmm->M, tr, 0, &aseqs, &ainfo))
	Die("Traces2Alignment() failed");
      if ((ofp = fopen(out_afile, "w")) == NULL)
	Die("Failed to open save file %s for final alignment.\n", out_afile);

      if (weighted)
	for (idx = 0; idx < nseqs; idx++)
	  { 
	    ainfo.sqinfo[idx].weight = wt[idx];
	    ainfo.sqinfo[idx].flags |= SQINFO_WGT;
	  }

      WriteSELEX(ofp, aseqs, nseqs, &ainfo, 50);
      printf("Alignment saved in file %s\n", out_afile);
      fclose(ofp);
      FreeAlignment(aseqs, nseqs, &ainfo);
      for (idx = 0; idx < nseqs; idx++)
	FreeTrace(tr[idx]);
    }

  /************************************************** 
   * Garbage collection (for completeness' sake)
   **************************************************/ 

  if (fix_nseqs > 0) FreeAlignment(fix_aseqs, fix_nseqs, &fix_ainfo);
  FreeSearchHMM(shmm);
  FreeHMM(hmm);
  for (idx = 0; idx < nseqs; idx++)
    {
      FreeSequence(seqs[idx], &(sqinfo[idx]));
      free(seqs1[idx]);
    }
  free(sqinfo);
  free(seqs);
  free(seqs1);
  free(wt);
  free(tr);

 /************************************************** 
  * Successful return to invocation environment
  **************************************************/
  return (0); 
}


static struct hmm_struc *
start_from_hmm(char *hmmfile)
{
  struct hmm_struc *hmm;
  FILE *hmmfp;

  if ((hmmfp = fopen(hmmfile, "rb")) == NULL)
    Die("Error: failed to open HMM file %s for reading.\n");
  if ((hmm = ReadHMM(hmmfp)) == NULL)
    Die("Error: failed to parse HMM file %s\n", hmmfile);
  fclose(hmmfp);

  return hmm;
}



static struct hmm_struc *
start_from_scratch(char **seqs, int N, int setlength, struct prior_s *prior)
{
  struct hmm_struc *hmm;        /* RETURN: hmm */
  int M;			/* estimated length of model to build */
  int seqidx;			/* counter for sequences */

				/* make a starting guess for length of model */
  if (setlength < 0)
    {
      M = 0;
      for (seqidx = 0; seqidx < N; seqidx++) 
	M += strlen(seqs[seqidx]);
      M /= N;
    }
  else
    {
      M = setlength;
    }
				/* allocate the hmm */
  hmm = AllocHMM(M);
  if (hmm == NULL) return NULL;

				/* write defaults into the hmm */
  if (!WriteFlatPriorHMM(hmm, prior))
    Die("error: failed to fill hmm with defaults\n");
  Renormalize(hmm);
  
  return hmm;
}

/* Function: read_kT()
 *
 * Purpose:  Read one line from the simulated annealing schedule file.
 *           Ignore comments (lines that start with #) and blank lines.
 *           First field is a float, which is the kT (temperature) value,
 *           typically ranging from 10.0 down to 0.1.
 *            
 *           Returns 0 when the file is empty.
 */
static int
read_kT(FILE *schedfp, float *ret_kT)
{
  char  buffer[256];
  char *sptr;
  float kT;

  do {
    if (fgets(buffer, 256, schedfp) == NULL) return 0; /* out of data */
  } while (*buffer == '#' || (sptr = strtok(buffer, " \t\n")) == NULL);

  kT = atof(sptr);
  if (kT > 10.0) {
    Warn("Trying to melt it? kT reset from %.2f to 10.0", kT);
    kT = 10.0;
  }
  if (kT < 0.1) {
    Warn("kT < 0.1 may cause numerical problems; reset from %.2f to 0.1", kT);
    kT = 0.1;
  }
  *ret_kT = kT;
  return 1;
}
    
/* Function: handle_interrupt()
 * 
 * Purpose:  Handle ctrl-C and kills. Set global flag "killnow",
 *           which will cause this to be the final iteration.
 *           Do not reset the handler, though, so a second signal
 *           will really kill the process.
 */ 
static void 
handle_interrupt(int sig)
{
  if (sig == SIGTERM || sig == SIGINT) 
    {
      killnow = TRUE;
      Warn("Caught a lethal interrupt signal. Will stop after this iteration.");
    }
}


