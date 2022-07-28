/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmls.c
 * Tue Aug 17 11:15:10 1993
 * modified from search.c
 *
 * main() for hmmls -- search long sequences for local matches to a hidden Markov model
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#define OPTIONS "bchFqr:t:w:"

static char usage[] = "\
Usage: hmmls [-options] <hmmfile> <dbfile>\n\
    where available options are:\n\
    -b          : report only single best match per sequence\n\
    -c          : search complementary strand too (DNA only)\n\
    -h          : print short help and version info, then exit\n\
    -q          : quiet - suppress verbose banner (for piping output)\n\
    -r <rfile>  : read random model from <rfile>\n\
    -t <thresh> : print only matches scoring above <thresh> (default 0)\n\
    -w <window> : set max length of match to <window> (default 1000)\n\
    -F          : fancy BLAST-style alignment output of matches\n";

static char banner[] = "hmmls - search long sequences for local matches to a hidden Markov model";

static SQINFO sqinfo;
static int   ext_revcomp;	/* True if we're doing the reverse strand */
static char *ext_seq;           /* make sequence available to print_hit() */
static char *ext_seq1;		/* sequence, prepped for alignment, 1..L  */
static int   ext_seqlen;	/* length of sequence */
static struct shmm_s *ext_shmm; /* scoring HMM, available to print_hit() */
static int   ext_fancyoutput;   /* True if we're doing fancy alignment output */

static int ext_beststart;
static int ext_bestend;
static int ext_bestscore;



static int dbmatch(struct vit_s **mx, int window, int M, int from_j, int from_k);
static int bestmatch(struct vit_s **mx, int window, int M, int from_j, int from_k);
static int print_hit(int from, int to, double score);

int
main(int argc, char **argv)
{
  struct hmm_struc *hmm;        /* hmm model to search with           */
  char        *seqfile;         /* name of database seq file to read  */
  char        *hmmfile;         /* name of HMM file to read           */
  FILE        *hmmfp;           /* opened HMM file                    */
  SQFILE      *dbfp;            /* opened sequence database file      */
  int          format;		/* format of sequence database        */
  char        *seq;             /* sequence to search against         */
  char        *rev;             /* reverse complement of seq          */
  int          window;	        /* OPTION: search matrix dimension    */
  float        cutoff;		/* OPTION: cutoff for score reporting */
  float        randomseq[MAXABET]; /* random sequence model, freq's   */
  char        *randomfile;

  int          do_complement;
  int          do_best;
  int          be_quiet;
  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */
#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif
  

  /***********************************************
   * Parse the command line
   ***********************************************/
 
  cutoff           = 0.0;
  window           = 1000;
  do_complement    = FALSE;
  ext_fancyoutput  = FALSE;
  do_best          = FALSE;
  be_quiet         = FALSE;
  randomfile       = NULL;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'b': do_best        = TRUE;         break;
    case 'c': do_complement  = TRUE;         break;
    case 'F': ext_fancyoutput= TRUE;         break;
    case 'q': be_quiet       = TRUE;         break;
    case 'r': randomfile     = optarg;       break;
    case 't': cutoff = (float) atof(optarg); break;
    case 'w': window = atoi(optarg);         break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
      exit(1);
    }

  if (argc - optind != 2)
    { fprintf(stderr, "%s\n", usage); exit(1);  }

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
    case SQERR_NOFILE: Die("Alignment file %s could not be opened for reading", seqfile);
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
      printf("%s\n     version %s, %s\n\n", banner, RELEASE, RELEASEDATE);

      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      printf("HMM file:                   %s\n", hmmfile);
      printf("Sequence database:          %s\n", seqfile);
      printf("Report scores above:        %.2f\n", cutoff);
      printf("Scan window size:           %d\n", window);
      printf("Do complementary strand:    %s\n", do_complement ? "yes" : "no");
      printf("Fancy alignment output:     %s\n", ext_fancyoutput ? "yes" : "no");
      if (do_best)
	printf("[Printing single best match per sequence]\n");
      else
	printf("[Printing multiple non-overlapping hits per sequence]\n");
      printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /***********************************************
   * Search database, one seq at a time
   ***********************************************/

  ext_shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, ext_shmm);
  while (ReadSeq(dbfp, format, &seq, &sqinfo))
    {
#ifdef MEMDEBUG
      orig_size = malloc_size(&histid1);
#endif
      s2upper(seq);
      PrepareSequence(seq, &ext_seq1, &ext_seqlen);
				/* pass dbmatch the name of the seq through
				   a static global. */
      ext_revcomp       = False;
      ext_seq           = seq;
      if (do_best) 
	{
	  ext_bestscore = -99999999;
	  if (DBViterbi(ext_shmm, ext_seq1, ext_seqlen, window, cutoff, bestmatch) == 0)
	    Die("viterbi fill step failed and returned\n");
	  if ((double) (ext_bestscore / INTSCALE) > cutoff)
	    print_hit(ext_beststart, ext_bestend, 
		      (double) (ext_bestscore / INTSCALE));
	}
      else
	{
	  if (DBViterbi(ext_shmm, ext_seq1, ext_seqlen, window, cutoff, dbmatch) == 0)
	    Die("viterbi fill step failed and returned\n");

	  /* If we're filtering overlaps, clean up now.
	   */
	  ReportScanHit(-1, -1, 0.0, print_hit);
	}
      free(ext_seq1);

      /* Do the reverse complement, if asked
       */
      if (do_complement)
	{
	  ext_revcomp = True;
	  if ((rev = (char *) malloc ((ext_seqlen+1) * sizeof(char))) == NULL)
	    Die("malloc failed");
	  revcomp(rev, seq);
	  ext_seq = rev;
	  PrepareSequence(rev, &ext_seq1, &ext_seqlen);

	  if (do_best)
	    {
	      ext_bestscore = -99999999;
	      if (DBViterbi(ext_shmm, ext_seq1, ext_seqlen, window, cutoff, bestmatch) == 0)
		Die("viterbi fill step failed and returned\n");
	      if ((double) (ext_bestscore / INTSCALE) > cutoff)
		print_hit(ext_beststart, ext_bestend, 
			  (double) (ext_bestscore / INTSCALE));
	    }
	  else
	    {
	      if (DBViterbi(ext_shmm, ext_seq1, ext_seqlen, window, cutoff, dbmatch) == 0)
		Die("viterbi fill step failed and returned\n");
	      ReportScanHit(-1, -1, 0.0, print_hit);
	    }
	  free(ext_seq1);
	  free(rev);
	}
      FreeSequence(seq, &sqinfo);

#ifdef MEMDEBUG
      if ((current_size = malloc_size(&histid2)) > orig_size)
	{
	  printf("Bad memory: current size %ld, original size %ld\n", 
		 current_size, orig_size);
	  malloc_list(2, histid1, histid2);
	}
#endif
    }

  /**************************************************
   * Successful return to invocation environment
   **************************************************/
  
  SeqfileClose(dbfp);
  return (0);
}



/* Function: dbmatch()
 * 
 * Purpose:  Score reporting function passed to dbviterbi.c. 
 * 
 * Args:     (Prototyped by dbviterbi.c. Do not alter.)
 *           mx     - scoring matrix filled by dbviterbi
 *           window - length of scrolled mx window, in rows
 *           M      - width of mx window, in columns (size of model)
 *           from_j - endpoint of match on sequence
 *           from_k - endpoint of match on model
 *                    
 * Return:   1 on success, 0 on failure.
 */         
static int
dbmatch(struct vit_s **mx, int window, int M, int from_j, int from_k)
{
  int   istart;

  /* A quick check. If our current score is worse than the score
   * of the previous nucleotide, don't bother. Only do this for
   * full model matches, where from_k = hmm->M+1
   */
  if (from_k == M+1 &&
      mx[from_j%window][from_k].score_m < mx[(from_j-1)%window][from_k].score_m)
    return 1;

  /* Traceback to get start point. Because END and BEGIN states 
   * are treated as MATCHes, istart and from_j are both 1 more than 
   * they should be.
   * If the traceback fails (returns 0) this means it was too long
   * for the window. We *assume* that this is because we're looking
   * at a long tail (the Alu bug) and we quietly return.
   */
  if (! ViterbiTrace(mx, ext_shmm, ext_seq1, window, from_j, from_k, 
		     (struct trace_s **) NULL, &istart, (int *) NULL))
    return 1;

  /* If we're filtering out overlapping scores, report the hit to
   * the filtering package (see scorestack.c); else print it now.
   */
  ReportScanHit(istart+1, from_j-1, (double) (mx[from_j%window][from_k].score_m / INTSCALE), 
		print_hit);
  return 1;
}


/* Function: bestmatch()
 * 
 * Purpose:  Score reporting function passed to dbviterbi.c. 
 * 
 * Args:     (Prototyped by dbviterbi.c. Do not alter.)
 *           mx     - scoring matrix filled by dbviterbi
 *           window - length of scrolled mx window, in rows
 *           M      - width of mx window, in columns (size of model)
 *           from_j - endpoint of match on sequence
 *           from_k - endpoint of match on model
 *                    
 * Return:   1 on success, 0 on failure.
 */
/*ARGSUSED*/
static int
bestmatch(struct vit_s **mx, int window, int M, int from_j, int from_k)
{
  int   istart;

  if (mx[from_j%window][from_k].score_m > ext_bestscore)
    {
      if (! ViterbiTrace(mx, ext_shmm, ext_seq1, window, from_j, from_k, 
			 (struct trace_s **) NULL, &istart, (int *) NULL))
	{			/* trace probably too long. Return quietly. */
	  ext_beststart = -1;
	  ext_bestend   = from_j - 1;
	  return 1;
	}
      ext_beststart = istart + 1;
      ext_bestend   = from_j - 1;
      ext_bestscore = mx[from_j%window][from_k].score_m;
    }
  return 1;
}


/* Function: print_hit()
 * 
 * Purpose:  Record a hit. We save the start and end points, the score,
 *           and the name of the sequence. 
 *           
 * Args:     (Prototyped by scorestack.c module. Do not alter)
 *           from  - start point of match
 *           to    - end point of match
 *           score - score of match
 *           
 * Return:   1 on success, 0 on failure.
 */
static int
print_hit(int from, int to, double score)
{
  int matchlen;			/* length of ext_seq that matches */
  char *ssq;			/* subsequence that matches       */
  char *ssq1;			/* prepared sequence              */
  struct vit_s **mx;
  struct trace_s        *tr;	/* traceback */
  int   showfrom, showto;	/* top strand coords */
  float  tmp;
  
				/* convert coords to top strand */
  if (ext_revcomp)
    {
      showfrom = ext_seqlen - from + 1;
      showto   = ext_seqlen - to   + 1;
    }
  else 
    {
      showfrom = from;
      showto   = to;
    }

  printf("%-6.2f (bits) f:%5d t:%5d Target: %s %s\n",
	 score, showfrom, showto, sqinfo.name,
	 (sqinfo.flags & SQINFO_DESC) ? sqinfo.desc : "");
  fflush(stdout);

  if (ext_fancyoutput)
    {
				/* extract subseq. to, from are 1..seqlen */
      matchlen = to - from + 1;
      if ((ssq = (char *) malloc (sizeof(char) * (matchlen + 1))) == NULL)
	Die("malloc failed");
      strncpy(ssq, ext_seq+from-1, matchlen);
      ssq[matchlen] = '\0';
      PrepareSequence(ssq, &ssq1, &matchlen);
				/* align to model and get traceback */
      if (! ViterbiFill(ext_shmm, ssq1, matchlen, &mx, &tmp))
	Die("ViterbiFill() failed\n"); 

      if (! ViterbiTrace(mx, ext_shmm, ssq1, matchlen+2, matchlen+1, ext_shmm->M+1, 
			 &tr, NULL, NULL))
	Die("ViterbiTrace() failed");

				/* print alignment */
      PrintFancyTrace(stdout, ext_shmm, tr, ext_seq, sqinfo.name, from-1);
      FreeTrace(tr);
      FreeViterbiMatrix(mx, matchlen);
      free(ssq);
      free(ssq1);
    }
  return 1;
}
