/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmma.c
 * modified Wed Jul 28 12:35:11 1993: alignment routine broken out to align.c
 * 
 * Main driver for hmma: multiple sequence alignment using HMM's
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#define OPTIONS "hmo:qr:s:RX"
static char usage[] = "\
Usage: hmma [-options] <hmmfile> <seqfile>,\n\
    where options are:\n\
\n\
    -h        : print short help and usage info\n\
    -m        : only print symbols aligned to MATCH states\n\
    -o <file> : save alignment in <file> in SELEX format\n\
    -q        : quiet - suppress verbose banner\n\
    -r <rfile>: read random model from <rfile>\n\
    -s <file> : save individual scores in <file>\n\
\n\
    and additional experimental options are:\n\
    -R        : ragged alignment, allowing Smith/Waterman fragments\n";

static char banner[] = "hmma  -- multiple sequence alignment using a hidden Markov model";

int
main(int argc, char **argv)
{
  struct hmm_struc *hmm;        /* hmm model to align with               */
  struct shmm_s    *shmm;       /* HMM in search form                    */
  struct trace_s  **tr;         /* tracebacks for each sequence          */
  struct vit_s    **mx;         /* storage for viterbi matrix            */
  char        *seqfile;         /* name of SELEX-format seq file to read */
  int          format;		/* format of sequence file               */
  char        *hmmfile;         /* name of HMM file to read              */
  FILE        *hmmfp;           /* opened HMM file                       */
  char       **rseqs;           /* de-aligned, raw sequence set          */
  SQINFO      *sqinfo;          /* info assoc. with sequences (0..N-1)   */
  char       **aseqs;           /* calculated new alignment              */
  AINFO        ainfo;		/* alignment info                        */
  char        *seq1;            /* prepared seq, with 1..L offset        */
  int          N;		/* number of sequences                   */
  int          idx;		/* counter over sequences                */
  float        score;		/* score of model vs. 1 seq              */
  float        tot_score;	/* score of model vs. all seqs           */
  int          L;		/* length of a sequence                  */
  FILE        *sfp;		/* ptr to open scorefile                 */
  int          swi, swj;	/* start, end point of a Smith/Waterman match 1..L */
  int          kstart, kend;	/* start, end of S/W match on the model 1..M  */
  float        P1,P2,P3;	/* controlling parameters for Smith/Waterman */
  float        randomseq[MAXABET]; /* random seq model, frequencies      */

  char        *randomfile;
  char        *outfile;         /* output alignment file */
  FILE        *outfp;           /* opened output file    */  
  char        *scorefile;	/* output score file     */
  int          matchonly;	/* TRUE to print only match-aligned symbols */
  int          be_quiet;	/* TRUE to suppress verbose banner */
  int          ragged;		/* TRUE to do fragment alignments by Smith/Waterman */

  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /***********************************************
   * Parse the command line
   ***********************************************/
 
  randomfile= NULL;
  outfile   = NULL;
  scorefile = NULL;
  matchonly = FALSE;
  be_quiet  = FALSE;
  ragged    = FALSE;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'm': matchonly = ! matchonly; break;
    case 'o': outfile   = optarg;      break;
    case 'q': be_quiet  = TRUE;        break;
    case 'r': randomfile= optarg;      break;
    case 's': scorefile = optarg;      break;

    case 'R': ragged    = TRUE;        break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
    }

  if (argc - optind != 2)
    { fprintf(stderr, "%s\n", usage); exit(1);  }

  hmmfile = strdup(argv[optind]);
  seqfile = strdup(argv[optind+1]);
  
  /***********************************************
   * Get the HMM (also sets alphabet)
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
   * Get sequence data
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

				/* read the seqs from file, ignore weights if any */
  if (! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &N))
    Die("Failed to read sequence file %s", seqfile);

  for (idx = 0; idx < N; idx++)
    s2upper(rseqs[idx]);

  /***********************************************
   * Allocations and output file opening
   ***********************************************/
  
  if ((tr = (struct trace_s **) malloc (sizeof(struct trace_s) * N)) == NULL)
    Die("malloc failed");

  if (scorefile != NULL)
    if ((sfp = fopen(scorefile, "w")) == NULL)
      Die("Failed to open score file %s", scorefile);


  /***********************************************
   * Print banner information
   ***********************************************/

  if (! be_quiet)
    {
      printf("%s\n    version %s, %s\n", banner, RELEASE, RELEASEDATE);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      printf("Model:    \t%s\n", hmmfile);
      printf("Sequences:\t%s   (%d total sequences)\n", seqfile, N);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      puts("");
    }

  /***********************************************
   * Gather tracebacks, then construct alignment
   ***********************************************/

  if (ragged && scorefile != NULL)
    {
      fprintf(sfp, "Score   seq-f seq-t hmm-f hmm-t Name and description\n");
      fprintf(sfp, "-----   ----- ----- ----- ----- --------------------\n"); 
    }

  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);
  P2 = 0.5;
  P3 = 1.0 / (float) (hmm->M-1);

  tot_score = 0.0;
  for (idx = 0; idx < N; idx++)
    {
      PrepareSequence(rseqs[idx], &seq1, &L);
      P1 = (float) L / (float) (L+1);

      if (ragged)
	{
	  if (AlignmentTooBig(L, shmm->M))
	    {
				/* the TRUE flag toggles Smith/Waterman behavior */
	      WeeViterbi(shmm, seq1, L, TRUE, P2, P3, &tr[idx], &score);
	      swi    = 1; 
	      swj    = L;
	      kstart = tr[idx]->nodeidx[1];
	      kend   = tr[idx]->nodeidx[tr[idx]->tlen-2];
	    }
	  else
	    {
				/* the TRUE flag forces all the sequence to align */
	      if (!SWViterbi(shmm, seq1, L, P1, P2, P3, TRUE,
			     &swi, &swj, &kstart, &kend, &score, &tr[idx]))
		Die("SWViterbi() failed");
	    }
	  tot_score += score;
	  
	  if (scorefile != NULL) 
	    fprintf(sfp, "%-6.2f  %5d %5d %5d %5d %s %s\n",
		    score, swi, swj, kstart, kend, 
		    sqinfo[idx].name,
		    (sqinfo[idx].flags & SQINFO_DESC) ? sqinfo[idx].desc : "");
	}
      else
	{
	  if (AlignmentTooBig(L, shmm->M))
	    {
				/* FALSE flag toggles Needleman/Wunsch behavior */
	      WeeViterbi(shmm, seq1, L, FALSE, 0., 0., &tr[idx], &score);
	    }
	  else
	    {
	      if (! ViterbiFill(shmm, seq1, L, &mx, &score))
		Die("ViterbiFill failed");
	      if (! ViterbiTrace(mx, shmm, seq1, L+2, L+1, shmm->M+1, &tr[idx], NULL, NULL))
		Die("ViterbiTrace() has failed and aborted");
	      FreeViterbiMatrix(mx, L);
	    }

	  tot_score += score;
      
	  if (scorefile != NULL) 
	    fprintf(sfp, "%-8.3f   : %s\n", score, sqinfo[idx].name);
	}

      free(seq1);
    }
  FreeSearchHMM(shmm);

  if (! Traces2Alignment(rseqs, sqinfo, N, hmm->M, tr, matchonly, &aseqs, &ainfo))
    Die("Traces2Alignment failed");
  
  /***********************************************
   * Print the alignment
   ***********************************************/

  if (outfile != NULL && (outfp = fopen(outfile, "w")) != NULL)
    {
      WriteSELEX(outfp, aseqs, N, &ainfo, 50);
      printf("Alignment saved in file %s\n", outfile);
      fclose(outfp);
    }
  else
    WriteSELEX(stdout, aseqs, N, &ainfo, 50);

  /**************************************************
   * Garbage collection (for completeness' sake)
   **************************************************/
  
  FreeAlignment(aseqs, N, &ainfo);
  for (idx = 0; idx < N; idx++)
    {
      FreeTrace(tr[idx]);
      FreeSequence(rseqs[idx], &(sqinfo[idx]));
    }
  free(sqinfo);
  free(tr);
  if (scorefile != NULL) fclose(sfp);

  /**************************************************
   * Successful return to invocation environment
   **************************************************/

  if (! be_quiet)
    printf("Alignment of %d sequences: score %.2f\n\n", N, tot_score / (float) N);
  return (0);
}

