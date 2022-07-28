/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* search.c
 * 
 * main() for hmms:
 * hmms -- search a database for matches to a hidden Markov model
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#define OPTIONS "hqr:t:FK:N:S"

static char usage[] = "\
Usage: hmms [-options] <hmmfile> <dbfile>\n\
    where available options are:\n\
\n\
    -h            : print short usage and version info, then exit\n\
    -q            : quiet - suppress verbose header info\n\
    -r <rfile>    : read random model from <rfile>\n\
    -t <thresh>   : print only matches scoring above <thresh>\n\
    -F            : do fancy output of model/seq alignment\n\
\n\
   Experimental options for suboptimal alignment output:\n\
    -K <kT>       : set temp factor for -S [1.0]\n\
    -N <num>      : # of suboptimal alignments to sample for -S [1]\n\
    -S            : sample a suboptimal alignment\n";

static char banner[] = "hmms - searching a sequence database for matches to a hidden Markov model";

int
main(int argc, char **argv)
{
  struct hmm_struc *hmm;        /* hmm (probability form)                   */
  struct shmm_s    *shmm;       /* hmm (integer log odds form for searching)*/
  struct trace_s   *tr;		/* traceback from ViterbiTrace()            */
  char        *seqfile;         /* name of database seq file to read        */
  char        *hmmfile;         /* name of HMM file to read                 */
  FILE        *hmmfp;           /* opened HMM file                          */
  SQFILE      *dbfp;            /* opened sequence database file            */
  int          format;		/* format of sequence database              */
  char        *seq;             /* a sequence (0..L-1)                      */
  SQINFO       sqinfo;		/* info for sequence                        */
  char        *seq1;            /* sequence, prepped for alignment (1..L)   */
  int          seqlen;		/* length of sequence                       */
  int          do_fancyoutput;  /* option: fancy model/seq alignment display */
  float        cutoff;
  int          be_quiet;	/* option: suppress verbosity (for piping output) */
  int          do_sample;       /* option: suboptimal alignments                */
  double       kT;		/* temperature factor for suboptimal alignments */
  int          samplenum;	/* number of alignments to sample               */
  struct vit_s **mx;            /* Viterbi DP matrix for an alignment       */
  float        vit_score;       /* score of an alignment                    */
  float        randomseq[MAXABET]; /* random sequence model, freq's         */
  char        *randomfile;

  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /***********************************************
   * Parse the command line
   ***********************************************/
 
  do_fancyoutput   = 0;
  cutoff           = 0.0;
  be_quiet         = FALSE;
  do_sample        = FALSE;
  kT               = 1.0;
  samplenum        = 1;
  randomfile       = NULL;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'q': be_quiet       = TRUE;                 break;
    case 'r': randomfile     = optarg;               break;
    case 't': cutoff         = (float) atof(optarg); break;
    case 'F': do_fancyoutput = TRUE;                 break;
    case 'K': kT             = atof(optarg);         break;
    case 'N': samplenum      = atoi(optarg);         break;
    case 'S': do_sample      = TRUE;                 break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
    }

  if (argc - optind != 2)
    Die("%s\n", usage); 

  hmmfile = argv[optind];
  seqfile = argv[optind+1];
  
  /***********************************************
   * Get the HMM (and sets alphabet)
   ***********************************************/

  if ((hmmfp = fopen(hmmfile, "rb")) == NULL)
    Die("Failed to open HMM file %s", hmmfile);

  if ((hmm = ReadHMM(hmmfp)) == NULL)
    Die("Failed to parse HMM file %s", hmmfile);

  if (fclose(hmmfp) == EOF)
    { fprintf(stderr, "file close failed!?\n"); exit(1); }

  if (randomfile == NULL) DefaultRandomModel(randomseq);
  else                    ReadRandomModel(randomfile, randomseq);

  /***********************************************
   * Open database file
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, "BLASTDB"))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence database %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

  dbfp = SeqfileOpen(seqfile, format, "BLASTDB");
  if (dbfp == NULL)
    Die("Failed to open sequence database file %s\n%s", seqfile, usage);
  
  /***********************************************
   * ready to start; print banner information
   ***********************************************/

  if (! be_quiet)
    {
      printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);

      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      printf(   "HMM file:                %s\n", hmmfile);
      printf(   "Sequence database:       %s\n", seqfile);
      printf(   "Cutoff at score:         %.2f\n", cutoff);
      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /***********************************************
   * Search database, one seq at a time
   ***********************************************/

  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);

  while (ReadSeq(dbfp, format, &seq, &sqinfo))
    {
      s2upper(seq);
      PrepareSequence(seq, &seq1, &seqlen);

      if (AlignmentTooBig(seqlen, shmm->M+1))
	WeeViterbi(shmm, seq1, seqlen, FALSE, 0., 0., &tr, &vit_score);
      else
	{
	  if (! ViterbiFill(shmm, seq1, seqlen, &mx, &vit_score))
	    Die("ViterbiFill failed\n"); 
	  if (! ViterbiTrace(mx, shmm, seq1, seqlen+2, seqlen+1, shmm->M+1, &tr, 
			     (int *) NULL, (int *) NULL))
	    Die("ViterbiTrace() has failed and aborted");
	  FreeViterbiMatrix(mx, seqlen);
	}

      if (vit_score > cutoff)
	{
	  printf("%-6.2f (bits)  Target: %s %s\n", vit_score, sqinfo.name,
		 (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
	  fflush(stdout);

	  if (do_fancyoutput || do_sample)
	    PrintFancyTrace(stdout, shmm, tr, seq, sqinfo.name, 0);
	}
      FreeTrace(tr);


      /* Provide additional sampled suboptimal alignments, if asked
       */
      if (do_sample && vit_score > cutoff)
	{
	  struct sa_s **sa_mx;
	  struct sa_hmm_s *sahmm;
	  float  score;
	  int i;
	  
	  sre_srandom((int) time((time_t *) NULL));

	  if (( sahmm = CreateSahmm(hmm, kT)) == NULL)
	    Die("internal error: failed to create an SA hmm\n");
	  	  
	  if (! SaFill(sahmm, seq, &sa_mx))
	    Die("internal error: SaFill() has failed and aborted\n");

	  for (i = 1; i <= samplenum; i++)
	    {
	      if (! SaTrace(sa_mx, seqlen, sahmm, &tr))
		Die("error: SaTrace() has failed and aborted\n");

	      TraceScore(shmm, seq, tr, &score);
	      printf("Suboptimal #%d: %-6.2f (bits)  Target: %s\n", i, score, sqinfo.name);
	      PrintFancyTrace(stdout, shmm, tr, seq, sqinfo.name, 0);
	    }

	  FreeTrace(tr);
	  Free2DArray(sa_mx, seqlen+2);
	  DestroySahmm(sahmm);
	}
      FreeSequence(seq, &sqinfo);
      free(seq1);
    }

  /**************************************************
   * Successful return to invocation environment
   **************************************************/
  SeqfileClose(dbfp);
  FreeHMM(hmm);
  FreeSearchHMM(shmm);
  return (0);
}



