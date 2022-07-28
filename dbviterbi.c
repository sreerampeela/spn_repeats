/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* dbviterbi.c
 * modified from viterbi.c
 * 
 * Search a long sequence for complete matches to an HMM.
 * Uses a rolling row trick to avoid enormous memory costs.
 * Action taken on a match depends on a passed function.
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "states.h"
#include "externs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: DBViterbi()
 * 
 * Perform the Viterbi dynamic programming calculation of aligning and 
 * scoring an HMM against a sequence. Restrict the current matrix to
 * a window, which is "scrolling" across the target sequence. Upon
 * finding a significant match, as defined by the caller, pass the current
 * matrix and the coordinates of the high-scoring cell to a caller-defined
 * function.
 * 
 * Return 1 on success, 0 on failure.
 * 
 */
int
DBViterbi(struct shmm_s *shmm,     /* HMM in search form                */
	  char          *seq,      /* sequence, 1..L                    */
	  int            L,        /* length of sequence                */
	  int            window,   /* size of scrolling window, rows    */
	  float          thresh,   /* scores above this are significant */
	  int          (*gotone_f)(struct vit_s **,int,int,int,int))
{
  struct vit_s **mx;		/* the viterbi calculation grid */
  struct vit_s  *thisrow;       /* ptr into current calculation row */
  struct vit_s  *nextrow;       /* ptr into next calculation row    */
  int   score;	                /* tmp variable for scores */
  int   i;			/* counter for sequence position: 0,1..l */
  int   k;			/* counter for model position: 0,1..m    */
  int   i_symidx;		/* index of symbol in alphabet (optimization) */
  int  *m_emit;			/* ptr into match scores for current symbol, 0..m  */
  int  *i_emit;                 /* ptr into insert scores for current symbol, 0..m */
  int  *tptr;                   /* ptr into transition scores, (m+1)*9 */
  int   ithresh;                /* integer version of thresh */

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* integer version of the threshold */
  ithresh = (int) (thresh * INTSCALE);

				/* if window is bigger than the sequence, just
				   do the whole sequence at once*/
  if (L+2 < window) window = L+2;
				/* allocate a scrolling window for the calculation matrix,
				   which is 0..window-1 rows by 0..m+1 cols */
  if (( mx = (struct vit_s **) malloc (sizeof(struct vit_s) * window) ) == NULL) 
    Die("memory failure allocating viterbi matrix\n");
  for (i = 0; i < window ; i++)
    if ((mx[i] = (struct vit_s *) malloc (sizeof(struct vit_s) * (shmm->M+2))) == NULL) 
      Die("memory failure allocating viterbi matrix, row %d\n", i);
  
  /********************************************
   * Initialization
   ********************************************/
      
				/* set up the 0,0 cell */
  mx[0][0].score_m = 0;		/* 0,0 cell, by definition, was produced by a match state */
  mx[0][0].score_i = -99999999;
  mx[0][0].score_d = -99999999;

				/* initialize the top row */
  for (k = 1; k <= shmm->M+1; k++)
    {
      mx[0][k].score_m = -99999999;
      mx[0][k].score_i = -99999999;
    }

  /********************************************
   * recursion: fill in the mx matrix
   ********************************************/

  for (i = 0; i <= L; i++)
    {
				/* calculate row indices in the scrolling window */
      thisrow = mx[i%window];
      nextrow = mx[(i+1)%window];
				/* initialize in next row */
      nextrow[0].score_m =  0;
      nextrow[0].score_d = -99999999;

				/* setup HMM scoring ptrs for this row */
      i_symidx = seq[i] - 'A';
      m_emit = shmm->m_emit[i_symidx];
      i_emit = shmm->i_emit[i_symidx];
      tptr   = shmm->t;

      for (k = 0; k <= shmm->M; k++)
	{
				/* add in emission scores to the current cell. */
	  if (k > 0) thisrow[k].score_m += *m_emit; 
	  m_emit++;
	  thisrow[k].score_i += *i_emit; 
	  i_emit++;
				/* deal with transitions out of delete state */
				/* note: tptr assumes order dd,di,dm,id,ii,im,md,mi,mm */
				/* we initialize with these */
	  thisrow[k+1].score_d = thisrow[k].score_d + *tptr; tptr++;
	  nextrow[k].score_i   = thisrow[k].score_d + *tptr; tptr++;
	  nextrow[k+1].score_m = thisrow[k].score_d + *tptr; tptr++;

				/* deal with transitions out of insert state */
	  if ((score = thisrow[k].score_i + *tptr) > thisrow[k+1].score_d)
	    thisrow[k+1].score_d = score;
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[k].score_i)
	    nextrow[k].score_i = score;
	  tptr++;
	  if ((score = thisrow[k].score_i + *tptr) > nextrow[k+1].score_m)
	    nextrow[k+1].score_m = score;
	  tptr++;
	  
			/* deal with transitions out of match state */
	  if ((score = thisrow[k].score_m + *tptr) > thisrow[k+1].score_d)
	    thisrow[k+1].score_d = score;
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[k].score_i)
	    nextrow[k].score_i = score;
	  tptr++;
	  if ((score = thisrow[k].score_m + *tptr) > nextrow[k+1].score_m)
	    nextrow[k+1].score_m = score;
	  tptr++;
	}

      /* Now check to see if best alignment ending at the current
       * row is significant; if it is, report it to caller 
       */
      if (thisrow[shmm->M+1].score_m  > ithresh)
	if (!(*gotone_f)(mx, window, shmm->M, i, shmm->M+1))
	  fprintf(stderr, "caller ignored report of a match!\n");
    }

				/* check final row. */
  thisrow = mx[(L+1)%window];
  if (thisrow[shmm->M+1].score_m  > ithresh)
    if (!(*gotone_f)(mx, window, shmm->M, i, shmm->M+1))
      fprintf(stderr, "caller ignored report of a match!\n");

  /********************************************
   * Garbage collection and return
   ********************************************/
#ifdef EXTREME_DEBUG
  PrintViterbiMatrix(mx, seq, L, shmm->M);
#endif
  
  Free2DArray(mx, window);
  return 1;
}
