/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* fragviterbi.c -- Smith-Waterman alignment of Haussler HMM's to sequences;
 *                  allowing multiple non-overlapping hits per sequence.
 * modified from swviterbi.c
 * 
 * Uses a rolling pointer trick to minimize memory costs; only two
 * rows need to be kept. Keeps traceback pointers so it can reconstruct
 * a series of matches.
 *
 * Can do either multiple hit, or single best hit Smith/Waterman. Normally,
 * one would use SWViterbi() for single best hit; the advantage of FragViterbi()
 * is that it can recover the bounds of the best Smith/Waterman hit on the
 * sequence in linear memory. It is not as convenient to use, however,
 * since the match is reported back via a function call, (*gotone_f)().
 * 
 * The matches are reordered and then passed one at a time to
 * (*gotone_f)() to report the hit. The sequence and the coords
 * are both 1..L.
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

/* Function: FragViterbi()
 * 
 * Perform the Viterbi dynamic programming calculation of aligning and 
 * scoring an HMM against a sequence. Restrict the current matrix to
 * a window, which is "scrolling" across the target sequence.
 * 
 * Args: shmm       - model, integer log odds form
 *       seq        - sequence 1..L, upper case, no ambiguous code
 *       L          - length of sequence
 *       singlehit  - TRUE to do single, rather than multiple S/W hits
 *       P1         - controls spacing between hits -- 1/P1 = expected spacing
 *       P2         - controls how often start is preferred -- P2=1/M = flat expectation
 *       P3         - controls how easy it is to leave the model
 *       thresh     - scores above this are reported
 *       gotone_f() - caller-supplied score reporting function 
 * 
 * Return 1 on success, 0 on failure.
 */
int
FragViterbi(struct shmm_s *shmm, char *seq, int L, int singlehit,
	    float P1, float P2, float P3, float thresh,
	    int (*gotone_f)(struct shmm_s *, char *, int,int,int,float))
{
  struct fvit_s   *mx[2];       /* the viterbi calculation grid, 2 rows     */
  int              score;	/* tmp variable for score                   */
  int              i;		/* counter for sequence position: 0,1..L    */
  int              k;		/* counter for model position: 0,1..M       */
  struct fvit_s   *thisrow;     /* (optimization) ptr to current row        */
  struct fvit_s   *nextrow;     /* (optimization) ptr to next row           */
  struct fvit_s   *init_row;    /* (optimization) generic initialized row   */
  int   init_row_size;          /* size of init_row in bytes                */
  int   next_k;                 /* (optimization) k+1                       */
  int  *m_emit;			/* ptr to match scores, 0..m                */
  int  *i_emit;                 /* ptr into insert scores, 0..m             */
  int  *tptr;                   /* ptr into transition scores, (m+1)*9      */
  int   ithresh;                /* integer version of thresh                */
  int  *startcost;		/* costs to start a S/W match at next row   */
  int  *endcost;                /* costs to terminate a S/W match           */
  int  *tracei;                 /* traceback ptr kept in end state          */
  int  *tracej;                 /* traceback ptr kept in end state          */
  int  *escores;                /* scores at END state, saved               */
  int   curri, currj, currend;	/* used during traceback                    */
  struct intstack_s *stack;     /* used for reversing traceback to 5'->3'   */
  int   bestsc, besti, bestj;   /* for best single-hit score                */

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* integer version of the threshold */
  ithresh = (int) (thresh * INTSCALE);

				/* allocate two rows for the calculation matrix, 
				   each 0..M+1 columns wide */
  for (i = 0; i < 2 ; i++)
    if ((mx[i]=(struct fvit_s *)calloc(shmm->M+2,sizeof(struct fvit_s)))==NULL)
      Die("malloc failed"); 

  if (! singlehit)
    {
      if ((tracei     = (int *) malloc (sizeof(int) * (L+2)))      == NULL ||
	  (tracej     = (int *) malloc (sizeof(int) * (L+2)))      == NULL ||
	  (escores    = (int *) malloc (sizeof(int) * (L+2)))      == NULL)
	Die("malloc failed");
    }
  
  if ((startcost  = (int *) malloc (sizeof(int) * (shmm->M+2)))== NULL ||
      (endcost    = (int *) malloc (sizeof(int) * (shmm->M+1)))== NULL)
    Die("malloc failed");
  
  /********************************************
   * Initialization
   ********************************************/
      
  /* Costs to start at various positions in the model are biased
   * towards a start at match 1 by parameter P2.
   */
  for (k = 2; k <= shmm->M; k++)
    startcost[k]       = (int) (INTSCALE * (LOG2((1.0-P1)/2.0) + 
					    LOG2((1.0-P2)/(double)(shmm->M-1))));
  startcost[0]         = -99999999;
  startcost[1]         = (int) (INTSCALE * (LOG2((1.0-P1)/2.0) + LOG2(P2)));
  startcost[shmm->M+1] = -99999999;

  /* Costs to end at various positions in the model are set by
   * the parameter P3; the final match, however, has probability 1.0
   */
  for (k = 1; k < shmm->M; k++)
    endcost[k] = (int) (INTSCALE * LOG2(P3));
  endcost[0]       = -99999999;
  endcost[shmm->M] = 0;

  /* create init_row, an initialized row of 
   * scores that can be slapped in place by
   * a memcpy() 
   */
  init_row_size = (shmm->M + 2) * sizeof(struct fvit_s);
  if ((init_row = (struct fvit_s *) malloc (init_row_size)) == NULL)
    Die("malloc failed");
  for (k = 0; k <= shmm->M+1; k++)
    {
      init_row[k].score_m = startcost[k];
      init_row[k].score_d = -99999999;
      init_row[k].score_i = -99999999;
      init_row[k].tback_m = -1;
      init_row[k].tback_d = -1;
      init_row[k].tback_i = -1;
    }
  if (! singlehit)
    {
      escores[0]  = -99999999; /* can't match to null sequence */
      tracei[0]   = -1;
      tracej[0]   = -1;
    }
  else
    bestsc = -99999999;

				/* set up the top row */
  memcpy(mx[0], init_row, init_row_size);
  for (k = 0; k <= shmm->M+1; k++)
    mx[0][k].score_m = -99999999;
  mx[0][shmm->M+1].score_m = 0;
  
  /********************************************
   * Recursion: fill in the mx matrix
   ********************************************/

  for (i = 0; i <= L; i++)
    {
			/* calculate row indices in the scrolling window */
      thisrow = mx[i % 2];
      nextrow = mx[(i + 1) % 2];

			/* init the next row. */
      memcpy(nextrow, init_row, init_row_size);
      for (k = 0; k <= shmm->M; k++)
	nextrow[k].tback_m = i+1;
      
      /* Check for connections from END->MATCH in this row; these
       * are restarts, paths that include previous matches, and are an
       * alternative to starting a completely fresh path. 
       * 
       */
      if (! singlehit)
	{
				/* >0 test *assumes* start costs more than
				 * best match emission gain  */
	  if (thisrow[shmm->M+1].score_m > 0) 
	    for (k = 1; k <= shmm->M; k++)
	      {
		score = thisrow[shmm->M+1].score_m + startcost[k];
		if (score > thisrow[k].score_m)
		  {
		    thisrow[k].score_m = score;
		    thisrow[k].tback_m = i;
		  }
	      }
	}

			/* set up HMM scoring ptrs for this row */
      m_emit = shmm->m_emit[seq[i]-'A'];
      i_emit = shmm->i_emit[seq[i]-'A'];
      tptr   = shmm->t;

			/* begin inner loop -- this is where all the time is spent */
      for (k = 0; k < shmm->M; k++)
	{ 
	  next_k = k + 1;

			/* add in emission scores to the current cell. */
	  thisrow[k].score_m += *m_emit; m_emit++;
	  thisrow[k].score_i += *i_emit; i_emit++;

			/* deal with transitions out of delete state */
	  if ((score = thisrow[k].score_d + *tptr) > thisrow[next_k].score_d)
	    {
	      thisrow[next_k].score_d = score; 
	      thisrow[next_k].tback_d = thisrow[k].tback_d;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_d + *tptr) > nextrow[k].score_i)
	    {
	      nextrow[k].score_i = score;
	      nextrow[k].tback_i = thisrow[k].tback_d;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_d + *tptr) > nextrow[next_k].score_m)
	    {
	      nextrow[next_k].score_m = score;
	      nextrow[next_k].tback_m = thisrow[k].tback_d;
	    }
	  tptr++;

			/* deal with transitions out of insert state */
	  if ((score = thisrow[k].score_i + *tptr) > thisrow[next_k].score_d)
	    {
	      thisrow[next_k].score_d = score;
	      thisrow[next_k].tback_d = thisrow[k].tback_i;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[k].score_i)
	    {
	      nextrow[k].score_i = score;
	      nextrow[k].tback_i = thisrow[k].tback_i;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[next_k].score_m)
	    {
	      nextrow[next_k].score_m = score;
	      nextrow[next_k].tback_m = thisrow[k].tback_i;
	    }
	  tptr++;

			/* deal with transitions out of match state */
	  if ((score = thisrow[k].score_m + *tptr) > thisrow[next_k].score_d)
	    {
	      thisrow[next_k].score_d = score;
	      thisrow[next_k].tback_d = thisrow[k].tback_m;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[k].score_i)
	    {
	      nextrow[k].score_i = score;
	      nextrow[k].tback_i = thisrow[k].tback_m;
	    }
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[next_k].score_m)
	    {
	      nextrow[next_k].score_m = score;
	      nextrow[next_k].tback_m = thisrow[k].tback_m;
	    }
	  tptr++;
	}

      /* Termination: k==M case.
       * insert M is irrelevant (terminal inserts are useless)
       * delete M is irrelevant (we'd stop the match at a previous match state instead)
       * the only thing we need to worry about is transitions into the
       * end state (MAT M+1) from all the match states 1..M.
       */
				/* add in emission scores to current cell M */
      thisrow[shmm->M].score_m += *m_emit; m_emit++;

      /* Multiple-hit ends: keep traceback info in tracei, tracej, escores
       */
      if (! singlehit)
	{ 
				/* transitions out of all match states to next row's end */
	  for (k = 1; k <= shmm->M; k++)
	    {
	      if ((score = thisrow[k].score_m + endcost[k]) > nextrow[shmm->M+1].score_m)
		{
		  nextrow[shmm->M+1].score_m = score;
		  tracei[i+1]                = thisrow[k].tback_m;
		  tracej[i+1]                = i;
		  escores[i+1]               = score;
		}
	    }

				/* alternatively, previous end may be better */
	  if (thisrow[shmm->M+1].score_m > nextrow[shmm->M+1].score_m)
	    {
	      nextrow[shmm->M+1].score_m = thisrow[shmm->M+1].score_m;
	      tracei[i+1]                = tracei[i];
	      tracej[i+1]                = tracej[i];
	      escores[i+1]               = thisrow[shmm->M+1].score_m;
	    }
	}

      /* else, single-hit ends.
       */
      else
	{
	  for (k = 1; k <= shmm->M; k++)
	    if ((score = thisrow[k].score_m + endcost[k]) > bestsc)
	      {
		bestsc = nextrow[shmm->M+1].score_m = score;
		besti  = thisrow[k].tback_m;
		bestj  = i;
	      }
	}
	  
#ifdef EXTREME_DEBUG
      /* Debugging: print this row
       */
      printf("%d MAT: ",i);
      for (k = 0; k <= shmm->M+1; k++)
	printf("%10d ", thisrow[k].score_m);
      printf("%10d", escores[i]);
      puts("");

      printf("%c INS: ", seq[i]);
      for (k = 0; k <= shmm->M+1; k++)
	printf("%10d ", thisrow[k].score_i);
      printf("%10d tracei", tracei[i]);
      puts("");
      
      printf("  DEL: ");
      for (k = 0; k <= shmm->M+1; k++)
	printf("%10d ", thisrow[k].score_d);
      printf("%10d tracej", tracej[i]);
      puts("");
      puts("");
#endif /* DEBUG */

    }

#ifdef EXTREME_DEBUG
  /* Debugging: print next row
   */
  printf("%d MAT: ",i);
  for (k = 0; k <= shmm->M+1; k++)
    printf("%10d ", nextrow[k].score_m);
  printf("%10d", escores[i]);
  puts("");

  printf("%c INS: ", seq[i]);
  for (k = 0; k <= shmm->M+1; k++)
    printf("%10d ", nextrow[k].score_i);
  printf("%10d tracei", tracei[i]);
  puts("");
      
  printf("  DEL: ");
  for (k = 0; k <= shmm->M+1; k++)
    printf("%10d ", nextrow[k].score_d);
  printf("%10d tracej", tracej[i]);
  puts("");
  puts("");
#endif /* DEBUG */

 /********************************************
  * Traceback - collect the matches.
  ********************************************/
 
  if (singlehit)
    {
      if (bestsc >= ithresh)
	if (!(*gotone_f)(shmm, seq, L, besti, bestj, (float) bestsc / INTSCALE))
	  Warn("Caller ignored a match at %d,%d, score %.1f\n", 
	       besti, bestj, (float) bestsc/INTSCALE);
    }
  else
    {
      /* We push j,i,score onto an integer stack in 3'->5' backwards order;
       * then pop off score,i,j in 5'->3' order and report them to caller
       */
      stack   = InitIntStack();
      curri   = tracei[L+1];
      currj   = tracej[L+1];
      currend = L+1;
      while (curri != -1 && currj != -1)
	{
	  score = escores[currend] - escores[curri];
	  PushIntStack(stack, score);
	  PushIntStack(stack, curri);
	  PushIntStack(stack, currj);

	  currend = curri;
	  curri   = tracei[currend];
	  currj   = tracej[currend];
	}

      while (PopIntStack(stack, &currj))
	{
	  PopIntStack(stack, &curri);
	  PopIntStack(stack, &score);
      
	  if (score >= ithresh)
	    if (!(*gotone_f)(shmm, seq, L, curri, currj, (float) score / INTSCALE))
	      Warn("Caller ignored a match at %d,%d, score %.1f\n", 
		   curri, currj, (float) score/INTSCALE);
	}
      FreeIntStack(stack);
    }
  /********************************************
   * Garbage collection
   ********************************************/
  
  free(mx[0]);
  free(mx[1]);
  free(init_row);
  free(startcost);
  free(endcost);
  if (!singlehit)
    {
      free(tracei);
      free(tracej);
      free(escores);
    }


  /********************************************
   * Return
   ********************************************/
  return 1;
}

