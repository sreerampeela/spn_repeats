/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* align_homologues.c
 * Adapted from hmma, Sat Jan 21 15:22:48 1995
 * 
 * This is not really an HMM program. It is a quick hack for a
 * particular task, using warped HMM code.
 *
 * The task is to align a set of homologues to a given alignment,
 * by pairwise alignment to the nearest relative in the alignment.
 * The original alignment itself is kept fixed. The final alignment
 * contains both the original aligned sequences and the homologues.
 * The homologues are put underneath their best relative in
 * the original alignment.
 *
 * The scoring parameters used for the alignment are PAM scores,
 * with heuristic gap open and gap extend costs -- not HMM
 * probabilities!
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "states.h"
#include "externs.h"
#include "squid.h"
#include "version.h"

char  Alphabet[MAXABET];		/* ACGT, for instance     */
int   Alphabet_size;		/* 4 or 20                */
int   Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* pam.h (from SWAT)
 * Default PAM matrix for searching
 *
 * This is the BLOSUM45 matrix of Henikoff & Henikoff,
 * with rows and columns rearranged into alphabetical order.
 * #  Matrix made by matblas from blosum45.iij
 * #  * column uses minimum score
 * #  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
 * #  Blocks Database = /data/blocks_5.0/blocks.dat
 * #  Cluster Percentage: >= 45
 * #  Entropy =   0.3795, Expected =  -0.2789
 */
static int default_pam[27][27] = {
{  5,-1,-1,-2,-1,-2, 0,-2,-1,0,-1,-1,-1,-1,0,-1,-1,-2, 1, 0,0, 0,-2, 0,-2,-1,-5 },
{ -1, 4,-2, 5, 1,-3,-1, 0,-3,0, 0,-3,-2, 4,0,-2, 0,-1, 0, 0,0,-3,-4,-1,-2, 2,-5 },
{ -1,-2,12,-3,-3,-2,-3,-3,-3,0,-3,-2,-2,-2,0,-4,-3,-3,-1,-1,0,-1,-5,-2,-3,-3,-5 },
{ -2, 5,-3, 7, 2,-4,-1, 0,-4,0, 0,-3,-3, 2,0,-1, 0,-1, 0,-1,0,-3,-4,-1,-2, 1,-5 },
{ -1, 1,-3, 2, 6,-3,-2, 0,-3,0, 1,-2,-2, 0,0, 0, 2, 0, 0,-1,0,-3,-3,-1,-2, 4,-5 },
{ -2,-3,-2,-4,-3, 8,-3,-2, 0,0,-3, 1, 0,-2,0,-3,-4,-2,-2,-1,0, 0, 1,-1, 3,-3,-5 },
{  0,-1,-3,-1,-2,-3, 7,-2,-4,0,-2,-3,-2, 0,0,-2,-2,-2, 0,-2,0,-3,-2,-1,-3,-2,-5 },
{ -2, 0,-3, 0, 0,-2,-2,10,-3,0,-1,-2, 0, 1,0,-2, 1, 0,-1,-2,0,-3,-3,-1, 2, 0,-5 },
{ -1,-3,-3,-4,-3, 0,-4,-3, 5,0,-3, 2, 2,-2,0,-2,-2,-3,-2,-1,0, 3,-2,-1, 0,-3,-5 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0 },
{ -1, 0,-3, 0, 1,-3,-2,-1,-3,0, 5,-3,-1, 0,0,-1, 1, 3,-1,-1,0,-2,-2,-1,-1, 1,-5 },
{ -1,-3,-2,-3,-2, 1,-3,-2, 2,0,-3, 5, 2,-3,0,-3,-2,-2,-3,-1,0, 1,-2,-1, 0,-2,-5 },
{ -1,-2,-2,-3,-2, 0,-2, 0, 2,0,-1, 2, 6,-2,0,-2, 0,-1,-2,-1,0, 1,-2,-1, 0,-1,-5 },
{ -1, 4,-2, 2, 0,-2, 0, 1,-2,0, 0,-3,-2, 6,0,-2, 0, 0, 1, 0,0,-3,-4,-1,-2, 0,-5 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0 },
{ -1,-2,-4,-1, 0,-3,-2,-2,-2,0,-1,-3,-2,-2,0, 9,-1,-2,-1,-1,0,-3,-3,-1,-3,-1,-5 },
{ -1, 0,-3, 0, 2,-4,-2, 1,-2,0, 1,-2, 0, 0,0,-1, 6, 1, 0,-1,0,-3,-2,-1,-1, 4,-5 },
{ -2,-1,-3,-1, 0,-2,-2, 0,-3,0, 3,-2,-1, 0,0,-2, 1, 7,-1,-1,0,-2,-2,-1,-1, 0,-5 },
{  1, 0,-1, 0, 0,-2, 0,-1,-2,0,-1,-3,-2, 1,0,-1, 0,-1, 4, 2,0,-1,-4, 0,-2, 0,-5 },
{  0, 0,-1,-1,-1,-1,-2,-2,-1,0,-1,-1,-1, 0,0,-1,-1,-1, 2, 5,0, 0,-3, 0,-1,-1,-5 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0 },
{  0,-3,-1,-3,-3, 0,-3,-3, 3,0,-2, 1, 1,-3,0,-3,-3,-2,-1, 0,0, 5,-3,-1,-1,-3,-5 },
{ -2,-4,-5,-4,-3, 1,-2,-3,-2,0,-2,-2,-2,-4,0,-3,-2,-2,-4,-3,0,-3,15,-2, 3,-2,-5 },
{  0,-1,-2,-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,0,-1,-1,-1, 0, 0,0,-1,-2,-1,-1,-1,-5 },
{ -2,-2,-3,-2,-2, 3,-3, 2, 0,0,-1, 0, 0,-2,0,-3,-1,-1,-2,-1,0,-1, 3,-1, 8,-2,-5 },
{ -1, 2,-3, 1, 4,-3,-2, 0,-3,0, 1,-2,-1, 0,0,-1, 4, 0, 0,-1,0,-3,-2,-1,-2, 4,-5 },
{ -5,-5,-5,-5,-5,-5,-5,-5,-5,0,-5,-5,-5,-5,0,-5,-5,-5,-5,-5,0,-5,-5,-5,-5,-5, 1 },
};

#define OPTIONS "ho:ql:u:"
static char usage[] = "\
Usage: align_homologues [-options] <alignment> <seqfile>,\n\
    where options are:\n\
\n\
    -h        : print short help and usage info\n\
    -o <file> : save alignment in <file> in SELEX format\n\
    -q        : quiet - suppress verbose banner\n\
    -l <id%>  : (lower limit) exclude homologues below this overall % identity [0-100]\n\
    -u <id%>  : (upper limit) exclude homologues above this overall % identity [0-100]\n";

static char banner[] = "align_homologues -- multiple alignment of homologues to a trusted alignment";

static void           defaultPAM(int ***ret_pam);
static struct shmm_s *single_query_hmm(char *rseq, int len, int **pam, int gop, int gex);
static float          pairwise_trace_identity(char *orig, struct trace_s *otr, 
                                              char *homo, struct trace_s *htr);
struct trace_s *      canonize_trace(char *aseq, int alen, int rlen, struct trace_s *tr);

int
main(int argc, char **argv)
{
  char            *aseqfile;    /* name of aligned seq file                   */
  int              aseq_fmt;    /* format of aseqfile                         */
  char           **orig_aseq;	/* original alignment                         */
  char           **orig_rseq;	/* dealigned sequences                        */
  AINFO            orig_ainfo;  /* original alignment info                    */
  int              orig_nseq;   /* # of aligned sequences                     */
  struct shmm_s  **shmm;        /* array of search HMMs for the aligned seqs  */
  struct trace_s **orig_tr;     /* tracebacks for each aligned seq            */
  int              idx;		/* counter over sequences                     */
  int              aidx;	/* counter over aligned sequences             */
  float            lower_id;    /* identity below which to exclude [0-1]      */ 
  float            upper_id;    /* identity above which to exclude [0-1]      */

  char            *homofile;	/* name of unaligned homologue sequence file  */
  int              homofmt;	/* format of unaligned homologue seq file     */
  char           **homo_rseq;	/* unaligned homologues                       */
  SQINFO          *homo_sqinfo; /* info for unaligned homologues              */
  int              homo_nseq;	/* number of unaligned homologues             */

  float           *homo_id;	/* fractional identity of aseqs vs. homologue */
  struct trace_s **homo_tr;	/* homologue tracebacks against best aseq     */
  int             *homo_aseq;	/* index of aseq aligned to; -1 for exclusion */

  struct trace_s **final_tr;	/* arrayed aseq + homologue traces            */
  SQINFO          *final_info;	/* arrayed aseq + homologues sqinfo           */
  char           **final_rseq;  /* raw sequences                              */
  char           **final_aseq;  /* final alignment                            */
  AINFO            final_ainfo; /* final alignment optional info              */
  int              final_nseq;  /* number of seqs included in final alignment */

  int              gop, gex;    /* gap open, gap extend penalties             */
  char             pammatrix[32];/* name of PAM matrix file                   */
  int            **pam;		/* PAM matrix, 27x27, alphabetic + *          */
  
  char            *outfile;     /* output alignment file                      */
  FILE            *outfp;       /* opened output file                         */  
  int              be_quiet;	/* TRUE to suppress verbose banner            */

  int              optchar;	/* option character                           */
  extern char     *optarg;      /* for getopt()                               */
  extern int       optind;	/* for getopt()                               */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_size(&histid1);
#endif


  /***********************************************
   * Parse the command line
   ***********************************************/
 
  outfile    = NULL;
  be_quiet   = FALSE;
  lower_id   = 0.0;
  upper_id   = 100.0;
  gop        = -14;
  gex        = -1;
  strcpy(pammatrix, "BLOSUM45");

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'o': outfile  = optarg;               break;
    case 'q': be_quiet = TRUE;                 break;
    case 'l': lower_id = atof(optarg) / 100.0; break;
    case 'u': upper_id = atof(optarg) / 100.0; break;

    case 'h': 
      printf("%s\n   version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage); 
      exit(0);
    default:
      Die("Unrecognized option: -%c\n%s\n", optchar, usage);
    }

  if (argc - optind != 2) 
    Die ("Incorrect number of command line arguments (need two)\n%s\n", usage); 

  aseqfile = argv[optind++];
  homofile = argv[optind];
  
  /***********************************************
   * Print banner information
   ***********************************************/

  if (! be_quiet)
    {
      printf("%s\n    version %s, %s\n", banner, RELEASE, RELEASEDATE);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      printf("Trusted alignment:    %s\n", aseqfile);
      printf("Homologues:           %s\n", homofile);
      printf("PAM matrix:           %s\n", pammatrix);
      printf("Gap open/extend:      %d/%d\n", gop, gex);
      if (lower_id > 0.0)
	printf("Exclude distant homologues below %3.0f%% identity\n", lower_id * 100.0);
      if (upper_id < 1.0)
	printf("Exclude close homologues above   %3.0f%% identity\n", upper_id * 100.0);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      puts("");
    }

  /***********************************************
   * Get the aligned sequences: assure upper case and determine alphabet
   ***********************************************/

  if (! SeqfileFormat(aseqfile, &aseq_fmt, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Alignment file %s could not be opened for reading", aseqfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of alignment file %s", aseqfile);
    }

  if (! ReadAlignment(aseqfile, aseq_fmt, &orig_aseq, &orig_nseq, &orig_ainfo))
    Die("Failed to parse aligned sequence file %s", aseqfile);

  for (idx = 0; idx < orig_nseq; idx++)
    s2upper(orig_aseq[idx]);

  if (! DetermineAlphabet(orig_aseq, orig_nseq))
    Die("Failed to determine alphabet used in seqs in %s", aseqfile);

  /***********************************************
   * Get unaligned homologues, assure upper case, assume same alphabet
   ***********************************************/

  if (! SeqfileFormat(homofile, &homofmt, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Homologue sequence file %s could not be opened for reading", homofile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of homologue sequence file %s", homofile);
    }

  if (! ReadMultipleRseqs(homofile, homofmt, &homo_rseq, &homo_sqinfo, &homo_nseq))
    Die("Failed to read sequence file %s", homofile);

  for (idx = 0; idx < homo_nseq; idx++)
    s2upper(homo_rseq[idx]);

  /***********************************************
   * Allocations and setup...
   ***********************************************/

  homo_id   = (float *)           MallocOrDie (sizeof(struct trace_s *) * homo_nseq);
  homo_tr   = (struct trace_s **) MallocOrDie (sizeof(struct trace_s *) * homo_nseq);
  homo_aseq = (int *)             MallocOrDie (sizeof(struct trace_s *) * homo_nseq);

				/* Get PAM matrix to score with */
  defaultPAM(&pam);

				/* Build search models from alignment */
  shmm = (struct shmm_s **) MallocOrDie (orig_nseq * sizeof(struct shmm_s *));
  DealignAseqs(orig_aseq, orig_nseq, &orig_rseq);
  for (idx = 0; idx < orig_nseq; idx++)
    shmm[idx] = single_query_hmm(orig_rseq[idx], strlen(orig_rseq[idx]), pam, gop, gex);

				/* build dummy traces for dealigned orig seqs*/
  orig_tr   = (struct trace_s **) MallocOrDie (sizeof(struct trace_s *) * orig_nseq);
  for (idx = 0; idx < orig_nseq; idx++)
    {
      int k;

      AllocTrace(shmm[idx]->M+2, &(orig_tr[idx]));
      orig_tr[idx]->tlen = shmm[idx]->M+2;
      for (k = 0; k <= shmm[idx]->M+1; k++)
	{
	  orig_tr[idx]->nodeidx[k]   = k;
	  orig_tr[idx]->statetype[k] = MATCH;
	  orig_tr[idx]->rpos[k]      = k-1;
	}
      orig_tr[idx]->rpos[shmm[idx]->M+1] = -1;
    }

  /***********************************************
   * Gather tracebacks of homologues against best matching aligned seq
   ***********************************************/

  for (idx = 0; idx < homo_nseq; idx++)
    {
      char           *seq1;
      int             L;
      struct trace_s *test_tr;
      float           ident;
      float           score;

      PrepareSequence(homo_rseq[idx], &seq1, &L);
      
      homo_id[idx]   = -0.01;	/* id must be positive, so one must be > than this  */
      homo_tr[idx]   = NULL;
      homo_aseq[idx] = -1;
      for (aidx = 0; aidx < orig_nseq; aidx++)
	{
	  /* Align homologue to each aligned seq.
	   * Save traces, scores.
	   */
	  if (AlignmentTooBig(L, shmm[aidx]->M))
	    {				/* FALSE flag toggles Needleman/Wunsch behavior */
	      WeeViterbi(shmm[aidx], seq1, L, FALSE, 0., 0., &test_tr, &score);
	    }
	  else
	    {
	      struct vit_s **mx;

	      if (! ViterbiFill(shmm[aidx], seq1, L, &mx, &score))
		Die("ViterbiFill failed");
	      if (! ViterbiTrace(mx, shmm[aidx], seq1, L+2, L+1, shmm[aidx]->M+1, &test_tr, NULL, NULL))
		Die("ViterbiTrace() has failed and aborted");
	      FreeViterbiMatrix(mx, L);
	    }

	  /* Calculate percent identity
	   */
	  ident = pairwise_trace_identity(orig_rseq[aidx], orig_tr[aidx], 
					  homo_rseq[idx], test_tr);
	  if (ident > homo_id[idx])
	    {
	      if (homo_tr[idx] != NULL) FreeTrace(homo_tr[idx]);
	      homo_id[idx]   = ident;
	      homo_tr[idx]   = test_tr;
	      homo_aseq[idx] = aidx;
	    }
	  else
	    FreeTrace(test_tr);
	}
      free(seq1);

      printf("%15.15s -- %3.0f%% identity to %15.15s",
	     homo_sqinfo[idx].name, homo_id[idx] * 100.0,
	     orig_ainfo.sqinfo[homo_aseq[idx]].name);
      if      (homo_id[idx] < lower_id) printf("  (rejected, too distant)\n");
      else if (homo_id[idx] > upper_id) printf("  (rejected, too close)\n");
      else printf("\n");
    }
  puts("");

  /***********************************************
   * Add deletes to tracebacks -- "canonize" them
   * make all traces relative to the alignment columns, rather than 
   * raw sequence models.
   * Arrange tracebacks for a single call to Traces2Alignment,
   * making the final multiple alignment. 
   ***********************************************/

  final_tr   = (struct trace_s **) MallocOrDie ((orig_nseq+homo_nseq) * sizeof(struct trace_s *));
  final_info = (SQINFO *) MallocOrDie ((orig_nseq+homo_nseq) * sizeof(SQINFO));
  final_rseq = (char **) MallocOrDie((orig_nseq+homo_nseq) * sizeof(char *));

  final_nseq = 0;
  for (aidx = 0; aidx < orig_nseq; aidx++)
    {
				/* original aligned seq definitely goes in final alignment */
      final_tr[final_nseq] = canonize_trace(orig_aseq[aidx], orig_ainfo.alen, 
					    strlen(orig_rseq[aidx]), orig_tr[aidx]);

      final_rseq[final_nseq] = orig_rseq[aidx];
      SeqinfoCopy(&(final_info[final_nseq]), &(orig_ainfo.sqinfo[aidx]));
      final_nseq++;
      
				/* pull out homologues that align best to this aseq */
      for (idx = 0; idx < homo_nseq; idx++)
	{
	  if (homo_aseq[idx] == aidx && homo_id[idx] >= lower_id && homo_id[idx] <= upper_id)
	    {
	      final_tr[final_nseq] = canonize_trace(orig_aseq[aidx], orig_ainfo.alen,
						    strlen(orig_rseq[aidx]), homo_tr[idx]);
	      final_rseq[final_nseq] = homo_rseq[idx];
	      SeqinfoCopy(&(final_info[final_nseq]), &(homo_sqinfo[idx]));
	      final_nseq++;
	    }
	}
#ifdef MEMDEBUG
      malloc_chain_check(1);
#endif
    }

  if (! Traces2Alignment(final_rseq, final_info, final_nseq, orig_ainfo.alen, final_tr, FALSE, 
			 &final_aseq, &final_ainfo))
    Die("Traces2Alignment failed");
  
  /***********************************************
   * Print the alignment
   ***********************************************/

  if (outfile != NULL && (outfp = fopen(outfile, "w")) != NULL)
    {
      WriteSELEX(outfp, final_aseq, final_nseq, &final_ainfo, 50);
      printf("Alignment saved in file %s\n", outfile);
      fclose(outfp);
    }
  else
    WriteSELEX(stdout, final_aseq, final_nseq, &final_ainfo, 50);

  /**************************************************
   * Garbage collection (for completeness' sake)
   **************************************************/
  
  FreeAlignment(orig_aseq, orig_nseq, &orig_ainfo);
  FreeAlignment(final_aseq, final_nseq, &final_ainfo);

  for (idx = 0; idx < homo_nseq; idx++)
    {
      FreeSequence(homo_rseq[idx], &(homo_sqinfo[idx]));
      FreeTrace(homo_tr[idx]);
    }
  free(homo_rseq);
  free(homo_sqinfo);
  free(homo_id);
  free(homo_tr);
  free(homo_aseq);

  for (aidx = 0; aidx < orig_nseq; aidx++)
    {
      FreeSearchHMM(shmm[aidx]);
      free(orig_rseq[aidx]);
      FreeTrace(orig_tr[aidx]);
    }
  free(shmm);
  free(orig_rseq);
  free(orig_tr);

  for (idx = 0; idx < final_nseq; idx++)
    {
      if (final_info[idx].flags & SQINFO_SS)   free(final_info[idx].ss);
      if (final_info[idx].flags & SQINFO_FREE) free(final_info[idx].free);
      FreeTrace(final_tr[idx]);
    }
  free(final_info);
  free(final_tr);
  free(final_rseq);

  Free2DArray(pam, 27);

  
  /**************************************************
   * Successful return to invocation environment
   **************************************************/
  
#ifdef MEMDEBUG
  malloc_chain_check(1);
  current_size = malloc_size(&histid2);
  
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    printf("No memory leaks, sir. (Don't look so surprised.)\n");
#endif
  return (0);
}


/* Function: defaultPAM()
 * 
 * Purpose:  Borrowed from swat - load the default PAM matrix from default_pam.
 *           Return pointer to a 2D 27x27 alphabetic order matrix.
 *           Default matrix happens to be BLOSUM45.
 */
static void
defaultPAM(int ***ret_pam)
{
  int **pam;
  int   i,j;

  if ((pam = (int **) malloc (sizeof(int *) * 27)) == NULL)
    Die("malloc failed");
  for (i = 0; i < 27; i++)
    if ((pam[i] = (int *) malloc (sizeof(int) * 27)) == NULL)
      Die("malloc failed");

  for (i = 0; i < 27; i++)
    for (j = 0; j < 27; j++)
      pam[i][j] = default_pam[i][j];

  *ret_pam = pam;
}



/* Function: single_query_hmm()
 * 
 * Purpose:  Create an HMM from a single query sequence.
 *           Use PAM scores and ad hoc gap penalties.
 *           This is for using HMM code to do classical-style
 *           alignment.
 *           
 * Args:     rseq:   query sequence. 0..len-1. All upper case.
 *           len:    length of rseq
 *           pam:    27x27 symmetric integer array, A..Z + *.
 *           gop:    gap open penalty. Must be negative.
 *           gex:    gap extend penalty. Must be negative.
 *           
 * Return:   ptr to allocated search HMM (integer log odds form).
 *           Caller must free.
 */
static struct shmm_s *
single_query_hmm(char *rseq, int len, int **pam, int gop, int gex)
{
  struct shmm_s *ersatz;
  int k;
  int sym, x;

  ersatz = AllocSearchHMM(len);
  for (k = 0; k <= len; k++)
    {
      if (k > 0)		/* Match symbol */
	{
				/* what's the symbol here? */
	  if      (isupper((int)rseq[k-1]))  sym = rseq[k-1] - 'A'; 
	  else { 
	    Warn("Bogus symbol %c in query sequence\n", rseq[k-1]);
	    sym = 'X'-'A';	/* treat as an X */
	  }
				/* Put PAM row for that symbol in */
	  for (x = 0; x < 26; x++)
	    ersatz->m_emit[x][k] = pam[sym][x];
	}  
      
				/* Inserts -- hardwired to zero score */
      for (x = 0; x < 26; x++)
	ersatz->i_emit[x][k] = 0;

				/* gap open, gap extends */
      ersatz->t[k*9 + Tdd] = gex;
      ersatz->t[k*9 + Tdi] = gop;
      ersatz->t[k*9 + Tdm] = 0;
      ersatz->t[k*9 + Tid] = gop;
      ersatz->t[k*9 + Tii] = gex;
      ersatz->t[k*9 + Tim] = 0;
      ersatz->t[k*9 + Tmd] = gop;
      ersatz->t[k*9 + Tmi] = gop;
      ersatz->t[k*9 + Tmm] = 0;
    }
  ersatz->flags = 0;
  return ersatz;
}

/* Function: pairwise_trace_identity()
 * 
 * Purpose:  Given tracebacks, calculate pairwise identity between
 *           two aligned sequences.
 */
static float
pairwise_trace_identity(char *orig, struct trace_s *otr, char *homo, struct trace_s *htr)
{
  char   *rseq[2];
  SQINFO  sqinfo[2];
  struct trace_s *tr[2];
  char  **aseq;
  AINFO   ainfo;
  float   ident;

  rseq[0] = orig;
  rseq[1] = homo;
  sqinfo[0].flags = sqinfo[1].flags = 0; /* sqinfo can be blank */
  tr[0]   = otr;
  tr[1]   = htr;
  Traces2Alignment(rseq, sqinfo, 2, strlen(orig), tr, FALSE, &aseq, &ainfo);
  
  s2upper(aseq[0]);
  s2upper(aseq[1]);
  PairwiseIdentity(aseq[0], aseq[1], FALSE, NULL, NULL, &ident);
  
  FreeAlignment(aseq, 2, &ainfo);
  return ident;
}
  
/* Function: canonize_trace()
 * 
 * Purpose:  Make a trace against a raw sequence model be relative to
 *           an aligned sequence model. 
 *           
 *           Searches are done against several different models built
 *           from dealigned sequences. To reconstruct the final alignment,
 *           all these traces are made consistent with a single model.
 *           The single "canonical" model is one in which all the
 *           columns of the original alignment are defined as MATCH
 *           columns.
 *           
 * Args:     aseq - aligned sequence the search model was built from
 *           alen - length of aseq (inclusive of gaps)
 *           rlen - length of rseq derived from aseq
 *           tr   - trace for rseq aligned to aseq-derived search model
 *           
 * Returns:  Returns a "canonized" traceback relative to the original
 *           alignment.
 */ 
struct trace_s *
canonize_trace(char *aseq, int alen, int rlen, struct trace_s *tr)
{
  struct trace_s *atr;          /* alignment-relative trace               */
  int tlen;			/* length of new alignment-relative trace */
  int otpos, ntpos;		/* position in old, new trace             */
  int apos;			/* position in aligned seq                */

				/* alen-rlen = DELETE states to add */
  tlen = tr->tlen + (alen - rlen);
  AllocTrace(tlen, &atr);
  
  otpos = ntpos = 0;
  apos = 0;
  while (otpos < tr->tlen)
    {
      atr->nodeidx[ntpos]  = apos;
      atr->statetype[ntpos]= tr->statetype[otpos];
      atr->rpos[ntpos]     = tr->rpos[otpos];

      if (tr->statetype[otpos] != INSERT && otpos != tr->tlen-1) apos++;
      otpos++; 
      ntpos++;

      while (isgap(aseq[apos-1]))
	{
				/* fold inserts into delete spaces as fake matches*/
	  if (tr->statetype[otpos] == INSERT)
	    {
	      atr->nodeidx[ntpos]  = apos;
	      atr->statetype[ntpos]= MATCH;
	      atr->rpos[ntpos]     = tr->rpos[otpos];
	      apos++;
	      ntpos++; otpos++;
	    }
	  else			/* add a delete to the trace */
	    {
	      atr->nodeidx[ntpos]   = apos;
	      atr->statetype[ntpos] = DELETE;
	      atr->rpos[ntpos]      = -1;
	      apos++;
	      ntpos++;
	    }
	}
    }
  if (ntpos > tlen) Die("Final canonized trace longer than estimate, %d > %d", ntpos, tlen);
  atr->tlen = ntpos;
#ifdef MEMDEBUG
  malloc_chain_check(1);
#endif
  return atr;
}
