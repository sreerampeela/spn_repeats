/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmfs.c - multiple-hit Smith-Waterman alignment of HMM's to sequences
 * SRE, Fri Jul 22 13:25:22 1994
 *
 * modified from hmmsw.c
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

static int fragmatch(struct shmm_s *shmm, char *seq, int L, int i, int j, float score);

static SQINFO sqinfo;
static int    is_revcomp;
static int    fancy_output;
static int    small_memory;

#define OPTIONS "chqr:t:FS"

static char usage[] = "\
Usage: hmmfs [-options] <hmmfile> <dbfile>\n\
   where available options are:\n\
    -c          : search complementary strand too (DNA only)\n\
    -h          : print short usage and version info\n\
    -q          : quiet - suppress verbose header info\n\
    -r <rfile>  : read random model from <rfile>\n\
    -t <thresh> : only report matches above a score of <cutoff>\n\
    -F          : fancy BLAST-style alignment output of matches\n";

static char banner[] = "\
hmmfs - multiple-hit Smith-Waterman local searching of a sequence database\n\
        for matches to a hidden Markov model";

int
main(int argc, char **argv)
{
  struct hmm_struc *hmm;        /* hmm model to search with                 */
  struct shmm_s    *shmm;       /* integer log-odds form model              */
  char        *seqfile;         /* name of database seq file to read        */
  int          format;		/* format of seqfile (i.e. kPearson)        */
  char        *hmmfile;         /* name of HMM file to read                 */
  FILE        *hmmfp;           /* opened HMM file                          */
  SQFILE      *dbfp;            /* opened sequence database file            */
  char        *seq;             /* sequence to search against               */
  char        *seq1;            /* prepared sequence for alignment          */
  char        *rev;             /* reverse complement of seq                */
  int          seqlen;
  float        randomseq[MAXABET]; /* random sequence model, freq's         */
  char        *randomfile;

  float        P1;		/* controls expected spacing between hits   */
  float        P2;		/* controls starting at beginning state     */
  float        P3;              /* controls ending from internal states     */

  float        cutoff;		/* cutoff threshold for reporting matches   */
  int          complement_too;  /* option: do complementary strand too      */
  int          be_quiet;	/* option: suppress verbosity (for piping)  */

  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_size(&histid1);
  printf("memory debug ON\n");
#endif

  /***********************************************
   * Parse the command line
   ***********************************************/
 
  complement_too = FALSE;
  be_quiet       = FALSE;
  cutoff         = 0.0;
  fancy_output   = FALSE;
  small_memory   = FALSE;
  randomfile     = NULL;

  P1 = 1000. / 1001.;
  P2 = 0.5;			/* some tendency to match at start */
  P3 = 1.0;			/* later divided by M - 1 */

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {
      
    case 'c': complement_too = TRUE;                 break;
    case 'q': be_quiet       = TRUE;                 break;
    case 'r': randomfile     = optarg;               break;
    case 't': cutoff         = (float) atof(optarg); break;
    case 'F': fancy_output   = TRUE;                 break;
    case 'S': small_memory   = TRUE;                 break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
    }

  if (argc - optind != 2) Die("%s\n", usage);
  if (small_memory && fancy_output)
    Die("Small memory (-S) and fancy output (-F) options are incompatible.");

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
    Die("file close failed!?\n");

  if (randomfile == NULL) DefaultRandomModel(randomseq);
  else                    ReadRandomModel(randomfile, randomseq);

  P3 /= (float) (hmm->M-1);

  /***********************************************
   * Open database file
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, "BLASTDB"))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Sequence file %s could not be opened for reading", seqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of sequence file %s", seqfile);
    }

  if ((dbfp = SeqfileOpen(seqfile, format, "BLASTDB")) == NULL)
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
      printf(   "Search strategy:          multiple-hit Smith-Waterman\n");
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


  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);

#ifdef EXTREME_DEBUG
  PrintSearchHMM(stdout, shmm);
#endif /* EXTREME_DEBUG */

  while (ReadSeq(dbfp, format, &seq, &sqinfo))
    {
      is_revcomp = FALSE;
      s2upper(seq);
      PrepareSequence(seq, &seq1, &seqlen);
				/* FALSE flag sets multiple-hit behavior */
      if (! FragViterbi(shmm, seq1, seqlen, FALSE, P1,P2,P3, cutoff, fragmatch))
	Die("FragViterbi failed");
      free(seq1);

      /* Do the reverse complement, if asked
       */
      if (complement_too)
	{
	  if ((rev = (char *) malloc ((seqlen+1) * sizeof(char))) == NULL)
	    Die("malloc failed");
	  revcomp(rev, seq);

	  is_revcomp = TRUE;
	  PrepareSequence(rev, &seq1, &seqlen);
	  if (! FragViterbi(shmm, seq1, seqlen, FALSE, P1,P2,P3, cutoff, fragmatch))
	    Die("FragViterbi failed");
	  free(seq1);
	  free(rev);
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
    {
      printf("Bad memory, boss.\n");
      malloc_list(2, histid1, histid2);
    }
  else
    fprintf(stderr, "No memory leaks, sir.\n");
#endif
    
  return (0);
}



/* Function: fragmatch()
 * 
 * Purpose:  Score reporting function passed to fragviterbi.c
 *           Prototyped in fragviterbi.c: do not change arguments
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
fragmatch(struct shmm_s *shmm, char *seq, int L, int i, int j, float score)
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
    if (! SWViterbi(shmm, ssq1, matchlen, P1,P2,P3, TRUE, 
		    &best_i, &best_j, &kstart, &kend,
		    &best_score, &tr))
      Die("SWViterbi() failed");

  if (best_i != 1 || best_j != matchlen)
    Warn("So much for that idea. new i,j = %d,%d\n", best_i, best_j);

  /* Print the hit
   */
  printf("%-6.2f  %5d %5d %5d %5d %s %s\n",
	 score, showfrom, showto, kstart, kend, sqinfo.name,
	 (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
  fflush(stdout);

  if (fancy_output)		/* print alignment */
    PrintFancyTrace(stdout, shmm, tr, seq+1, sqinfo.name, i-1);

  free(ssq);
  free(ssq1);
  FreeTrace(tr);
  return 1;
}
