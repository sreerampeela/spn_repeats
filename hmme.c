/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmme.c
 * Emit sequences from a model.
 * 
 * Given an hmm, emit a number of sequences.
 *             
 * Tue Nov 24 10:16:29 1992
 * radical alterations Wed Jul 28 12:12:57 1993: hmme
 *   used to be "hmm evaluation" and did some pretty
 *   useless stuff.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "squid.h"
#include "states.h"
#include "externs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#define OPTIONS "bfhn:o:qr:s:B"

static char usage[] = "\
Usage: hmme [-options] <hmmfile>\n\
   where available options are:\n\
   -b            : emit single best (most probable) sequence\n\
   -f            : write in FASTA format, not SELEX\n\
   -h            : print short help and version info\n\
   -n <number>   : emit <number> sequences\n\
   -o <outfile>  : save sequences to outfile\n\
   -q            : quiet - suppress verbose banner\n\
   -r <rfile>    : read random model from <rfile>\n\
   -s <seed>     : set seed for random()\n\
\n\
   Experimental options:\n\
   -B            : like -b, but emit highest scoring sequence\n";

static char banner[] = "hmme -- sequence emission from a hidden Markov model";

int
main(int argc, char **argv)
{
  char   *hmmfile;              /* HMM file to open                */
  FILE   *hmmfp;                /* opened hmm file pointer         */
  struct hmm_struc *hmm;        /* the hidden Markov model         */
  int     idx;			/* counter over sequences          */
  char  **rseqs;		/* emitted sequences               */
  SQINFO *sqinfo;		/* array of info for rseqs         */
  char  **aseqs;                /* alignment of emitted sequences  */
  AINFO   ainfo;		/* info for aseqs                  */
  struct trace_s **tr;          /* tracebacks from emitted seqs    */
  float   score;                /* score of emitted sequence       */
  float   randomseq[MAXABET];	/* random sequence model, freqs    */
  char   *randomfile;

  int          optc;
  extern int   optind;
  extern char *optarg;
  
  int    emitnum;		/* option: number of sequences to emit */
  int    seed; 
  char  *outfile;
  FILE  *ofp;
  int    do_best;
  int    do_rd_best;
  int    be_quiet;
  int    do_fasta;

  /***********************************************
   * Parse the command line
   ***********************************************/

  emitnum    = 10;
  seed       = (int) time ((time_t *) NULL);  /* default: "random" seed */
  outfile    = NULL;
  do_best    = FALSE;
  do_rd_best = FALSE;
  be_quiet   = FALSE;
  do_fasta   = FALSE;
  randomfile = NULL;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'b': do_best    = TRUE;         break;
    case 'f': do_fasta   = TRUE;         break;
    case 'n': emitnum    = atoi(optarg); break;
    case 'o': outfile    = optarg;       break;
    case 'q': be_quiet   = TRUE;         break;
    case 'r': randomfile = optarg;       break;
    case 's': seed       = atoi(optarg); break;
    case 'B': do_rd_best = TRUE;         break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      fprintf(stderr, "Unrecognized option: -%c\n%s\n", optc, usage);
      exit(1);
    }

  if (argc - optind != 1)
    { fprintf(stderr, "%s\n", usage); exit(1);  }

  hmmfile = argv[optind];

  sre_srandom(seed);

  if (outfile == NULL) ofp = stdout;
  else if ((ofp = fopen(outfile, "w")) == NULL)
    Die ("Failed to open output file %s", outfile);

  /***********************************************
   * Input of HMM (and sets alphabet)
   ***********************************************/

  if ((hmmfp = fopen(hmmfile, "rb")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);

  if ((hmm = ReadHMM(hmmfp)) == NULL)
    Die("failed to parse HMM file %s", hmmfile);

  if (fclose(hmmfp) == EOF)
    Die("file close failed!?");

  if (randomfile == NULL) DefaultRandomModel(randomseq);
  else                    ReadRandomModel(randomfile, randomseq);

  /***********************************************
   * Allocations
   ***********************************************/

  if (do_best || do_rd_best) emitnum = 1;
  if ((rseqs      = (char **) malloc (emitnum * sizeof(char *))) == NULL ||
      (sqinfo     = (SQINFO *) malloc (emitnum * sizeof(SQINFO))) == NULL ||
      (tr = (struct trace_s **) malloc (emitnum * sizeof(struct trace_s *))) == NULL)
    Die("malloc failed");

  /***********************************************
   * Emit and align the sequences, then print 'em
   ***********************************************/

  if (! be_quiet)
    {
      printf("%s\n", banner);
      printf("     version %s, %s\n", RELEASE, RELEASEDATE);
      puts("");
    }

  /* There are two versions of the best sequence code.
   * Mine produces the best *state* sequence, then emits the most likely
   * symbol sequence from it. Richard Durbin's version emits the 
   * most likely sequence relative to a random model.
   */
  if (do_best || do_rd_best)
    {
      struct shmm_s *shmm;

      if (do_rd_best)
	RD_EmitBestSequence(hmm, randomseq, &(rseqs[0]), &(tr[0]));
      else
	EmitBestSequence(hmm, &(rseqs[0]), &(tr[0]));

      shmm = AllocSearchHMM(hmm->M);
      MakeSearchHMM(hmm, randomseq, shmm);
      TraceScore(shmm, rseqs[0], tr[0], &score);

      if (do_rd_best)
	{
	  strcpy(sqinfo[0].name, "best-score");
	  sprintf(sqinfo[0].desc, "Best scoring sequence (score %.2f bits)", score);
	  sqinfo[0].flags = SQINFO_NAME | SQINFO_DESC;
	}
      else
	{
	  strcpy(sqinfo[0].name, "most-likely");
	  sprintf(sqinfo[0].desc, "Most probable sequence (score %.2f bits)", score);
	  sqinfo[0].flags = SQINFO_NAME | SQINFO_DESC;
	}

      WriteSeq(ofp, kPearson, rseqs[0], &(sqinfo[0]));

      FreeSequence(rseqs[0], &(sqinfo[0]));
      FreeTrace(tr[0]);
      FreeSearchHMM(shmm);
    }
  else
    {
      for (idx = 0; idx < emitnum; idx++)
	if (! EmitSequence(hmm, &rseqs[idx], &tr[idx]))
	  Die("failed to emit sequence number %d", idx);

				/* give them names */
      for (idx = 0; idx < emitnum; idx++)
	{
	  sqinfo[idx].flags = SQINFO_NAME;
	  sprintf(sqinfo[idx].name, "seq%d", idx);
	}

      if (do_fasta)
	{
	  for (idx = 0; idx < emitnum; idx++)
	    WriteSeq(ofp, kPearson, rseqs[idx], &(sqinfo[idx]));
	}
      else /* SELEX alignment */
	{
	  if (! Traces2Alignment(rseqs, sqinfo, emitnum, hmm->M, tr, FALSE, &aseqs, &ainfo))
	    Die("Traces2Alignment failed");
	  if (! WriteSELEX(ofp, aseqs, emitnum, &ainfo, 50))
	    Die("WriteSELEX failed");
	  FreeAlignment(aseqs, emitnum, &ainfo);
	}
      
      for (idx = 0; idx < emitnum; idx++)
	{
	  FreeSequence(rseqs[idx], &(sqinfo[idx]));
	  FreeTrace(tr[idx]);
	}
    }

  /***********************************************
   * Exit
   ***********************************************/

  if (outfile != NULL) fclose(ofp);
  free(tr);
  free(rseqs);
  free(sqinfo);
  FreeHMM(hmm);
  return 0;
}


