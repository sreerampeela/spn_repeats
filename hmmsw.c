/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmsw.c - Smith-Waterman alignment of HMM's to sequences
 *           finds single best match of a submodel to a subsequence
 *
 * SRE, Mon Jan 25 13:43:45 1993
 * v1.1: SRE, Thu Jul 15 18:27:32 1993
 *
 * modified from search.c
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "states.h"

#include "squid.h"
#include "externs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

/* For linear-memory search and linear-memory alignment reconstruction
 */
static int linear_match(struct shmm_s *shmm, char *seq, int L, int i, int j, float score);
static SQINFO sqinfo;
static int    is_revcomp;
static int    fancy_output;


#define OPTIONS "chqr:t:F"

static char usage[] = "\
Usage: hmmsw [-options] <hmmfile> <dbfile>\n\
    where available options are:\n\
    -c          : search complementary strand too (DNA only)\n\
    -h          : print short usage and version info\n\
    -q          : quiet - suppress verbose header info\n\
    -r <rfile>  : read random model from <rfile>\n\
    -t <thresh> : only report matches above a score of <cutoff>\n\
    -F          : fancy BLAST-style alignment output of matches\n";

static char banner[] = "\
hmmsw - Smith-Waterman local searching of a sequence database\n\
        for best matches to a hidden Markov model";

int
main(int argc, char **argv)
{
  struct hmm_struc *hmm;        /* hmm                                      */
  struct shmm_s    *shmm;       /* search model, integer log odds form      */
  char        *seqfile;         /* name of database seq file to read        */
  char        *hmmfile;         /* name of HMM file to read                 */
  FILE        *hmmfp;           /* opened HMM file                          */
  SQFILE      *dbfp;            /* opened sequence database file            */
  int          format;		/* format of sequence database              */
  int          seqlen;          /* length of sequence                       */
  char        *seq;             /* sequence to search against               */
  char        *seq1;            /* prepared sequence 1..L                   */
  char        *rev;             /* reverse complement of seq                */
  float        cutoff;		/* cutoff threshold for reporting matches   */
  int          best_i;		/* start of match on sequence, 1..L         */
  int          best_j;		/* end of match on sequence, 1..L           */
  int          best_kstart;     /* start of match on model, 1..M            */
  int          best_kend;	/* end of match on model, 1..M              */
  float        best_score;	/* score of match, bits                     */
  struct trace_s *best_tr;      /* traceback of match                       */
  float        randomseq[MAXABET]; /* random sequence model, freq's         */
  char        *randomfile;

  float        P1;		/* controls expected spacing between hits   */
  float        P2;		/* controls starting at beginning state     */
  float        P3;              /* controls ending from internal states     */
  
  int          complement_too;  /* option: do complementary strand too      */
  int          be_quiet;	/* option: suppress verbosity (for piping)  */

  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /***********************************************
   * Parse the command line
   ***********************************************/
 
  complement_too = FALSE;
  be_quiet       = FALSE;
  cutoff         = 0.0;
  fancy_output   = FALSE;
  randomfile     = NULL;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {
      
    case 'c': complement_too = TRUE;                 break;
    case 'q': be_quiet       = TRUE;                 break;
    case 'r': randomfile     = optarg;               break;
    case 't': cutoff         = (float) atof(optarg); break;
    case 'F': fancy_output   = TRUE;                 break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
    }

  if (argc - optind != 2) Die("%s\n", usage);
  hmmfile = argv[optind];
  seqfile = argv[optind+1];
  
  /***********************************************
   * Get the HMM (and sets alphabet)
   ***********************************************/

  if ((hmmfp = fopen(hmmfile, "rb")) == NULL)
    Die("HMM file %s open failed", hmmfile);

  if ((hmm = ReadHMM(hmmfp)) == NULL)
    Die("failed to parse HMM file %s", hmmfile);

  if (fclose(hmmfp) == EOF)
    { fprintf(stderr, "file close failed!?\n"); exit(1); }

  if (randomfile == NULL) DefaultRandomModel(randomseq);
  else                    ReadRandomModel(randomfile, randomseq);

  /***********************************************
   * Open database file
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, "BLASTDB"))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

  dbfp = SeqfileOpen(seqfile, format, "BLASTDB");
  if (dbfp == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
  
  /***********************************************
   * ready to start; print banner information
   ***********************************************/

  if (! be_quiet)
    {
      printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);

      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      printf(   "HMM file:                 %s\n", hmmfile);
      printf(   "Sequence database:        %s\n", seqfile);
      printf(   "Search strategy:          Smith-Waterman\n");
      printf(   "Cutoff at score:          %.2f\n", cutoff);
      if (complement_too)
	printf(   "Search complement too:    yes\n");
      if (fancy_output)
	printf(   "Alignment output:         yes\n");
      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
      printf("Score   seq-f seq-t hmm-f hmm-t Name and description\n");
      printf("-----   ----- ----- ----- ----- --------------------\n");  
    }

  /***********************************************
   * Search database, one seq at a time
   ***********************************************/

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif

  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);
  P2 = 0.5;
  P3 = 1.0 / (float) (hmm->M - 1);

  while (ReadSeq(dbfp, format, &seq, &sqinfo))
    {
      s2upper(seq);
      PrepareSequence(seq, &seq1, &seqlen);
      P1 = (float) seqlen / (float) (seqlen + 1);

      if (AlignmentTooBig(seqlen, shmm->M))
				/* TRUE flag sets single-hit behavior */
	FragViterbi(shmm, seq1, seqlen, TRUE, P1,P2,P3, cutoff, linear_match);
      else
	{
	  if (! SWViterbi(shmm, seq1, seqlen, P1,P2,P3, FALSE, &best_i, &best_j, 
		      &best_kstart, &best_kend, &best_score, &best_tr))
	    Die("SWViterbi() failed");
	  
	  if (best_score > cutoff)
	    {
	      printf("%-6.2f  %5d %5d %5d %5d %s %s\n",
		     best_score, best_i, best_j, best_kstart, best_kend, 
		     sqinfo.name,
		     (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
	      fflush(stdout);

	      if (fancy_output)		/* print alignment */
		PrintFancyTrace(stdout, shmm, best_tr, seq, sqinfo.name, best_i-1);
	    }

	  FreeTrace(best_tr);
	}
      free(seq1);

      if (complement_too)
	{
	  if ((rev = (char *) malloc ((seqlen+1) * sizeof(char))) == NULL)
	    Die("malloc failed");
	  revcomp(rev, seq);
	  PrepareSequence(rev, &seq1, &seqlen);
	  
	  if (AlignmentTooBig(seqlen, shmm->M))
	    FragViterbi(shmm, seq1, seqlen, TRUE, P1,P2,P3, cutoff, linear_match);
	  else
	    {
	      if (! SWViterbi(shmm, seq1, seqlen, P1,P2,P3, FALSE, &best_i, &best_j, 
			      &best_kstart, &best_kend, &best_score, &best_tr))
		Die("SWViterbi() failed");
	  
	      if (best_score > cutoff)
		{
		  printf("%-6.2f  %5d %5d %5d %5d %s %s\n",
			 best_score, seqlen - best_i + 1, seqlen - best_j + 1,
			 best_kstart, best_kend, sqinfo.name,
			 (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
		  fflush(stdout);

		  if (fancy_output)		/* print alignment */
		    PrintFancyTrace(stdout, shmm, best_tr, rev, sqinfo.name, best_i-1);
		}
	      FreeTrace(best_tr);
	    }
	  free(rev);
	  free(seq1);
	}
      FreeSequence(seq, &sqinfo);
    }

  /**************************************************
   * Successful return to invocation environment
   **************************************************/
  FreeHMM(hmm);
  FreeSearchHMM(shmm);
  SeqfileClose(dbfp);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "No memory leaks, sir.\n");
#endif
    

  return (0);
}


/* Function: linear_match()
 * 
 * Purpose:  Score reporting function passed to fragviterbi.c, when
 *           doing linear-memory best hit Smith/Waterman
 *           Prototyped in fragviterbi.c: do not change arguments
 *           Modified from fragmatch(), in hmmfs.c
 *           
 * Args:     shmm  - model, integer log odds form
 *           seq   - prepared sequence 1..L
 *           L     - length of sequence
 *           i     - start point of match (1..L)
 *           j     - end point of match (1..L)
 *           score - score of the match, in bits
 *                   
 * Return:   1 on success, 0 on failure                  
 */
static int 
linear_match(struct shmm_s *shmm, char *seq, int L, int i, int j, float score)
{
  struct trace_s *tr;           /* traceback of match             */
  char *ssq;                    /* extracted subseq               */
  char *ssq1;                   /* subseq, prepared for alignment */
  int   matchlen;               /* length of match                */
  int   showfrom, showto;	/* 1..L coords to be printed      */
  int   kstart, kend;		/* start, end on model, 1..M      */
  int   best_i, best_j;         /* best start, end on subseq (should be 1,matchlen) */
  float best_score;            /* score of SW match to subseq    */
  float P1,P2,P3;		/* controlling parameters for SW */

  if (is_revcomp)
    { showfrom = L - i + 1;
      showto   = L - j + 1;
    }
  else
    { showfrom = i;
      showto   = j;
    }
  
  /* Extract the subsequence that matched, and realign using
   * SWViterbi to get the best match. 
   */
  matchlen = j - i + 1;
  if ((ssq = (char *) malloc (sizeof(char) * (matchlen+1))) == NULL)
    Die("malloc failed");
  strncpy(ssq, seq + i, matchlen);
  ssq[matchlen] = '\0';
  PrepareSequence(ssq, &ssq1, &matchlen);
  P1 = (float) L / (float) (L+1); 
  P2 = 0.5;
  P3 = 1.0 / (float) (shmm->M-1);

  if (AlignmentTooBig(L, shmm->M))
    {
      WeeViterbi(shmm, ssq1, matchlen, TRUE, P2, P3, &tr, &best_score);
      best_i = 1; 
      best_j = matchlen;
      kstart = tr->nodeidx[1];
      kend   = tr->nodeidx[tr->tlen-2];
    }
  else
				/* TRUE forces whole sequence to align */
    if (! SWViterbi(shmm, ssq1, matchlen, P1,P2,P3, TRUE, &best_i, &best_j, &kstart, &kend,
		    &best_score, &tr))
      Die("SWViterbi() failed");

  if (best_i != 1 || best_j != matchlen)
    Warn("So much for that idea. new i,j = %d,%d\n", best_i, best_j);

  /* Print the hit
   */
  printf("%-6.2f  %5d %5d %5d %5d %s %s\n",
	 score, showfrom, showto, kstart, kend, sqinfo.name,
	 (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");


  if (fancy_output)		/* print alignment */
    PrintFancyTrace(stdout, shmm, tr, seq+1, sqinfo.name, i-1);

  free(ssq);
  free(ssq1);
  FreeTrace(tr);
  return 1;
}

