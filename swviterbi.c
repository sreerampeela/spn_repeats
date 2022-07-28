/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* 
 * swviterbi.c -- Smith-Waterman alignment of Haussler HMM's to sequences;
 * modified from dbviterbi.c
 * 
 * Search for local alignments (part of sequence against part
 * of HMM). Reports only the single best match. For multiple
 * matches, see fragviterbi.c.
 * 
 * 
 * The prototype for this function is:
 *     int (*gotone_f)(mx, window, hmm, i, k, ret_thresh)
 *         struct swvit_s **mx;
 *         int              window;
 *         struct hmm_struc *hmm;
 *         int              i        (start trace from row i)
 *         int              k        (start trace from column k)
 *         float           *ret_thresh; 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <math.h>

#include "states.h"
#include "squid.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: SWViterbi()
 * 
 * Perform the Viterbi dynamic programming calculation of aligning and 
 * scoring an HMM against a sequence. Restrict the current matrix to
 * a window, which is "scrolling" across the target sequence. 
 * 
 * Upon finding a significant match, trace it back. Keep the best trace.
 * When done, report the best trace, best score, etc. back to the caller.
 * 
 * Args: shmm       - model, integer log odds form
 *       seq        - prepared sequence, 1..L
 *       L          - length of sequence
 *       P1         - controls spacing between hits -- 1/P1 = expected spacing
 *       P2         - controls how often start is preferred -- P2=1/M = flat expectation
 *       P3         - controls how easy it is to leave the model
 *       fullseq    - TRUE to force whole sequence to align
 *       ret_i      - RETURN: start point on seq, 1..L
 *       ret_j      - RETURN: end point on seq, 1..L
 *       ret_kstart - RETURN: start point on model, 1..M
 *       ret_kend   - RETURN: end point on model, 1..M
 *       ret_score  - RETURN: score of best match
 *       ret_tr     - RETURN: traceback of best match
 *
 * Return 1 on success, 0 on failure.
 * ret_tr is alloc'ed here, must be free'd by caller
 */
int
SWViterbi(struct shmm_s *shmm, char *seq, int L, 
	  float P1, float P2, float P3, int fullseq,
	  int *ret_i, int *ret_j, int *ret_kstart, int *ret_kend, 
	  float *ret_score, struct trace_s **ret_tr)
{
  struct vit_s   **mx;          /* the viterbi calculation grid                 */
  int              score;	/* tmp variable for score                       */
  int              i;		/* counter for sequence position: 0,1..L        */
  int              k;		/* counter for model position: 0,1..M           */
  struct vit_s    *thisrow;     /* (optimization) ptr to current row            */
  struct vit_s    *nextrow;     /* (optimization) ptr to next row               */
  struct vit_s    *init_row;    /* (optimization) generic initialized row       */
  int   init_row_size;          /* size of init_row in bytes                    */
  int   max;                    /* keeps track of maximum score in row          */
  int   k_of_max;               /* keeps track of where max is                  */
  int   next_k;                 /* (optimization) k+1                           */
  int   i_symidx;	        /* index of symbol in alphabet (optimization)   */
  int  *m_emit;			/* ptr into match scores for current symbol, 0..m  */
  int  *i_emit;                 /* ptr into insert scores for current symbol, 0..m */
  int  *tptr;                   /* ptr into transition scores, (m+1)*9          */
  int  *startcost;		/* cost to start a S/W match at next row        */
  int  *endcost;		/* cost to end a S/W match at this row          */
  int   besti;			/* best start point                             */
  int   bestj;                  /* best end point                               */
  int   best_kstart;            /* best model start point                       */
  int   best_kend;              /* best model end point                         */
  int   best_score;             /* best score, scaled integer form              */
  struct trace_s *best_tr;      /* traceback of best hit                        */

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* allocate a scrolling window for the calculation matrix,
				   which is 0..window-1 rows by 0..M+1 cols */
  if (( mx = (struct vit_s **) malloc (sizeof(struct vit_s) * (L+2)) ) == NULL) 
    Die("malloc failed"); 
  for (i = 0; i < L+2 ; i++)
    if ((mx[i] = (struct vit_s *) calloc (shmm->M+2, sizeof(struct vit_s))) == NULL) 
      Die("malloc failed"); 
  if ((startcost = (int *) malloc (sizeof(int) * (shmm->M+2))) == NULL ||
      (endcost   = (int *) malloc (sizeof(int) * (shmm->M+2))) == NULL)
    Die("malloc failed");
  
  /********************************************
   * Initialization
   ********************************************/

  /* Create position-specific start and end costs
   */
  for (k = 2; k <= shmm->M; k++)
    startcost[k] = (int) (INTSCALE * (LOG2((1.0-P1)/2.0) + 
				      LOG2((1.0-P2)/(double)(shmm->M-1))));
  startcost[0]         = -99999999;
  startcost[1]         = (int) (INTSCALE * (LOG2((1.0-P1)/2.0) + LOG2(P2)));
  startcost[shmm->M+1] = -99999999;
  
  for (k = 1; k < shmm->M; k++)
    endcost[k] = (int) (INTSCALE * LOG2(P3));
  endcost[0]       = -99999999;
  endcost[shmm->M] = 0;
  
				/* create init_row, an initialized row of scores
				   that can be slapped in place by a memcpy() */
  init_row_size = (shmm->M + 2) * sizeof(struct vit_s);
  if ((init_row = (struct vit_s *) malloc (init_row_size)) == NULL)
    Die("malloc failed");
  for (k = 0; k <= shmm->M + 1; k++)
    {
      init_row[k].score_m = -99999999;
      init_row[k].score_d = -99999999;
      init_row[k].score_i = -99999999;
    }
				/* set up the top rows */
  memcpy(mx[0], init_row, init_row_size);

  /********************************************
   * Recursion: fill in the mx matrix
   ********************************************/

  best_score = -99999999;
  for (i = 0; i <= L; i++)
    {
      thisrow = mx[i];
      nextrow = mx[i + 1];
				/* initialize the values in the next row. */
				/* note the Smith/Waterman statistical corrections,
				   startcost and endcost; see notes */
      memcpy(nextrow, init_row, init_row_size);
      if (!fullseq || i == 0)
	for (k = 1; k <= shmm->M; k++)
	  nextrow[k].score_m = startcost[k];

				/* set up HMM scoring ptrs for this row */
      i_symidx = seq[i] - 'A';
      m_emit = shmm->m_emit[i_symidx] + 1;
      i_emit = shmm->i_emit[i_symidx] + 1;
      tptr   = shmm->t + 9;

      max = -99999999;
      for (k = 1; k <= shmm->M; k++)
	{ /* begin inner loop -- this is where all the time is spent */
	  next_k = k + 1;

				/* add in emission scores to the current cell. */
	  thisrow[k].score_m += *m_emit; m_emit++;
	  thisrow[k].score_i += *i_emit; i_emit++;

				/* deal with transitions out of delete state */
	  if (fullseq || thisrow[k].score_d > startcost[k])
	    {
	      if ((score = thisrow[k].score_d + *tptr) > thisrow[next_k].score_d)
		thisrow[next_k].score_d = score;
	      if ((score = thisrow[k].score_d + *(tptr+1)) > nextrow[k].score_i)
		nextrow[k].score_i = score;
	      if ((score = thisrow[k].score_d + *(tptr+2)) > nextrow[next_k].score_m)
		nextrow[next_k].score_m = score;
	    }
	  tptr += 3;

				/* deal with transitions out of insert state */
	  if (fullseq || thisrow[k].score_i > startcost[k])
	    {
	      if ((score = thisrow[k].score_i + *tptr) > thisrow[next_k].score_d)
		thisrow[next_k].score_d = score;
	      if ((score = thisrow[k].score_i + *(tptr+1)) > nextrow[k].score_i)
		nextrow[k].score_i = score;
	      if ((score = thisrow[k].score_i + *(tptr+2)) > nextrow[next_k].score_m)
		nextrow[next_k].score_m = score;
	    }
	  tptr+= 3;

				/* deal with transitions out of match state */
	  if (fullseq || thisrow[k].score_m > startcost[k])
	    {
	      if ((score = thisrow[k].score_m + *tptr) > thisrow[next_k].score_d)
		thisrow[next_k].score_d = score;
	      if ((score = thisrow[k].score_m + *(tptr+1)) > nextrow[k].score_i)
		nextrow[k].score_i = score;
	      if ((score = thisrow[k].score_m + *(tptr+2)) > nextrow[next_k].score_m)
		nextrow[next_k].score_m = score;
	    }
	  tptr += 3;

				/* look at our current score and decide if it's a maximum */
	  if (!fullseq || i == L)
	    if (thisrow[k].score_m + endcost[k] > max)
	      {
		max      = thisrow[k].score_m + endcost[k];
		k_of_max = k;
	      }
	}
				/* check the best score in this row */
      if (max > best_score)
	{ 
	  bestj      = i;
	  best_kend  = k_of_max;
	  best_score = max;
	}
    }

  /* Traceback of the best hit
   */
  if (! ViterbiTrace(mx, shmm, seq, L+2, bestj, best_kend, 
		     &best_tr, &besti, &best_kstart))
    Die("Traceback failed in SWViterbi()\n");


  /********************************************
   * Garbage collection
   ********************************************/
  
  Free2DArray(mx, (L+2));
  free(init_row);
  free(startcost);
  free(endcost);

  /********************************************
   * Return
   ********************************************/
  *ret_i      = besti;
  *ret_j      = bestj;
  *ret_kstart = best_kstart;
  *ret_kend   = best_kend;
  *ret_score  = (float) best_score / INTSCALE;
  *ret_tr     = best_tr;
  return 1;
}

