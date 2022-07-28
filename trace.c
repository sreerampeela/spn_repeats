/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* trace.c
 * Tue Mar  8 16:29:22 1994
 *
 * Stuff to do with tracebacks.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "squid.h"
#include "states.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: AllocTrace(), ReallocTrace(), FreeTrace()
 * 
 * Purpose:  allocation and freeing of traceback structures
 */
void
AllocTrace(int tlen, struct trace_s **ret_tr)
{
  struct trace_s *tr;
  
  if ((tr = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    Die("malloc failed");
  if ((tr->nodeidx   = (int *)  malloc (sizeof(int)  * tlen)) == NULL ||
      (tr->statetype = (char *) malloc (sizeof(char) * tlen)) == NULL ||
      (tr->rpos      = (int *)  malloc (sizeof(int)  * tlen)) == NULL)
    Die("malloc failed");
  *ret_tr = tr;
}
void
ReallocTrace(struct trace_s *tr, int tlen)
{
  if ((tr->nodeidx   = (int *)  realloc (tr->nodeidx,   tlen * sizeof(int)))  == NULL ||
      (tr->statetype = (char *) realloc (tr->statetype, tlen * sizeof(char))) == NULL ||
      (tr->rpos      = (int *)  realloc (tr->rpos,      tlen * sizeof(int)))  == NULL)
    Die("malloc failed");
}
void 
FreeTrace(struct trace_s *tr)
{
  free(tr->nodeidx);
  free(tr->statetype);
  free(tr->rpos);
  free(tr);
}



/* Function: ReverseTrace()
 * 
 * Purpose:  Tracebacks are more easily constructed backwards,
 *           using overallocated trace_s structures. Here we
 *           reverse the arrays of a traceback and give some
 *           extra memory back to malloc. 
 *           
 * Arguments: tr   - the traceback to reverse
 *            tlen - actual length of the traceback.
 *                   This means that END is at 0 and BEGIN is at
 *                   tlen-1.              
 */
void
ReverseTrace(struct trace_s *tr, int tlen)
{
  int  *new_nodeidx;
  char *new_statetype;
  int  *new_rpos;
  int   opos;
  int   npos;

				/* allocate for new reversed arrays */
  if ((new_nodeidx   = (int *)  malloc (sizeof(int)  * tlen)) == NULL ||
      (new_statetype = (char *) malloc (sizeof(char) * tlen)) == NULL ||
      (new_rpos      = (int *)  malloc (sizeof(int)  * tlen)) == NULL)
    Die("malloc failed");
				/* reverse the arrays */
  for (npos = 0, opos = tlen-1; npos < tlen; npos++, opos--)
    {
      new_nodeidx[npos]   = tr->nodeidx[opos];
      new_statetype[npos] = tr->statetype[opos];
      new_rpos[npos]      = tr->rpos[opos];
    }
				/* free old, switch in the new */
  free(tr->nodeidx);   tr->nodeidx   = new_nodeidx;
  free(tr->statetype); tr->statetype = new_statetype;
  free(tr->rpos);      tr->rpos      = new_rpos;
  tr->tlen = tlen;
}

/* Function: PrintTrace()
 * 
 * Purpose:  Debugging. Print out a traceback.
 */
void
PrintTrace(struct trace_s *tr)
{
  int j;

  printf("=== Traceback === (length %d)\n", tr->tlen);
  printf("Node indexes:  ");
  for (j = 0; j < tr->tlen; j++) printf("%2d ", tr->nodeidx[j]);
  printf("\nState types:   ");
  for (j = 0; j < tr->tlen; j++) printf("%2d ", tr->statetype[j]);
  printf("\nSeq positions: ");
  for (j = 0; j < tr->tlen; j++) printf("%2d ", tr->rpos[j]);
  puts("");
}


/* Function: TraceCount()
 * 
 * Purpose:  Count a traceback into a count-based HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 * Args:     hmm   - counts-based HMM
 *           seq   - the sequence that the traceback aligns to the HMM (0..L-1)
 *           wt    - weight on the sequence
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
void
TraceCount(struct hmm_struc *hmm, char *seq, float wt, struct trace_s *tr)
{
  int spos;			/* position in tr */
  int rpos;                     /* symbol position in seq */
  
  for (spos = 0; spos < tr->tlen; spos++)
    {
      rpos = tr->rpos[spos];

      /* Emission counts
       */
      if (tr->nodeidx[spos] >  0 && 
	  tr->nodeidx[spos] <= hmm->M &&
	  tr->statetype[spos] == MATCH)
	CountSymbol(seq[rpos], wt, hmm->mat[tr->nodeidx[spos]].p);
      else if (tr->statetype[spos] == INSERT)
	CountSymbol(seq[rpos], wt, hmm->ins[tr->nodeidx[spos]].p);

      /* State transition counts
       */
      if (spos < tr->tlen-1) {
	switch (tr->statetype[spos]) {
	case MATCH:
	  switch (tr->statetype[spos+1]) 
	    {
	    case MATCH:  hmm->mat[tr->nodeidx[spos]].t[MATCH]  += wt; break;
	    case INSERT: hmm->mat[tr->nodeidx[spos]].t[INSERT] += wt; break;
	    case DELETE: hmm->mat[tr->nodeidx[spos]].t[DELETE] += wt; break;
	    default: Die("unrecognized statetype %d", tr->statetype[spos+1]);
	    }
	  break;

	case INSERT:
	  switch (tr->statetype[spos+1]) 
	    {
	    case MATCH:  hmm->ins[tr->nodeidx[spos]].t[MATCH]  += wt; break;
	    case INSERT: hmm->ins[tr->nodeidx[spos]].t[INSERT] += wt; break;
	    case DELETE: hmm->ins[tr->nodeidx[spos]].t[DELETE] += wt; break;
	    default: Die("unrecognized statetype %d", tr->statetype[spos+1]);
	    }
	  break;

	case DELETE:
	  switch (tr->statetype[spos+1]) 
	    {
	    case MATCH:  hmm->del[tr->nodeidx[spos]].t[MATCH]  += wt; break;
	    case INSERT: hmm->del[tr->nodeidx[spos]].t[INSERT] += wt; break;
	    case DELETE: hmm->del[tr->nodeidx[spos]].t[DELETE] += wt; break;
	    default: Die("unrecognized statetype %d", tr->statetype[spos+1]);
	    }
	  break;
	  
	default: Die("Unrecognized statetype %d", tr->statetype[spos]);
	}
      }
    }
}


/* Function: TraceScore()
 * 
 * Purpose:  Calculate a score from a traceback. Used for the emit.c
 *           functions, and eventually maxmodelmaker.c.
 *           
 * Args:     shmm      - search form HMM structure
 *           seq       - sequence 0..len-1
 *           tr        - traceback that aligns seq to hmm
 *           ret_score - RETURN: the score of that alignment
 *           
 * Return:   1 on success, 0 on failure.          
 */
int
TraceScore(struct shmm_s  *shmm, 
	   char           *seq, 
	   struct trace_s *tr,
	   float          *ret_score)
{
  int    pos;			/* position in seq */
  int    sym;
  int    spos;			/* position in state sequence */
  int    score;
  int    k;

  score = 0;
  for (spos = 0; spos < tr->tlen; spos++)
    {
      pos = tr->rpos[spos];

      /* Transition cost
       */
      if (spos > 0) {
	k = tr->nodeidx[spos-1];
	switch (tr->statetype[spos-1]) {
	case MATCH:
	  switch (tr->statetype[spos]) {
	  case MATCH:  score += shmm->t[k*9 + Tmm]; break;
	  case DELETE: score += shmm->t[k*9 + Tmd]; break;
	  case INSERT: score += shmm->t[k*9 + Tmi]; break;
	  default: Die("unrecognized statetype %d\n", spos);
	  }
	  break;
	case DELETE:
	  switch (tr->statetype[spos]) {
	  case MATCH:  score += shmm->t[k*9 + Tdm]; break;
	  case DELETE: score += shmm->t[k*9 + Tdd]; break;
	  case INSERT: score += shmm->t[k*9 + Tdi]; break;
	  default: Die("unrecognized statetype %d\n", spos);
	  }
	  break;
	case INSERT:
	  switch (tr->statetype[spos]) {
	  case MATCH:  score += shmm->t[k*9 + Tim]; break;
	  case DELETE: score += shmm->t[k*9 + Tid]; break;
	  case INSERT: score += shmm->t[k*9 + Tii]; break;
	  default: Die("unrecognized statetype %d\n", spos);
	  }
	  break;
	default: Die("unrecognized statetype %d\n", spos);
	}
      }

      /* Emission cost
       */
      k = tr->nodeidx[spos];
      if (tr->statetype[spos] == MATCH && 
	  tr->nodeidx[spos] <= shmm->M &&
	  tr->nodeidx[spos] > 0)
	{
	  sym = isupper((int) seq[pos]) ? seq[pos] - 'A' : seq[pos] - 'a';
	  score += shmm->m_emit[sym][k];
	  pos++;
	}
      if (tr->statetype[spos] == INSERT)
	{
	  sym = isupper((int) seq[pos]) ? seq[pos] - 'A' : seq[pos] - 'a';
	  score += shmm->i_emit[sym][k];
	  pos++;
	}
    }	  
  *ret_score = (float) (score / INTSCALE);
  return 1;
}


/* Function: DealignTrace()
 * 
 * Purpose:  Take a traceback relative to an aligned sequence
 *           (such as the fake tracebacks produced by Maxmodelmaker())
 *           and make it relative to the raw sequence.
 *
 *           Robust against S/W (local alignment) traces.
 *           
 * Args:     tr    - traceback to dealign
 *           aseq  - aligned sequence corresponding to traceback
 *           alen  - length of aseq
 */
void
DealignTrace(struct trace_s *tr, char *aseq, int alen)
{
  int *rmap;			/* position in raw sequence, 0..alen-1 */
  int  rpos, apos;
  int  tpos;			/* position in traceback */

  /* Construct a mapping of aligned seq to raw sequence coords
   */
  rmap = (int *) MallocOrDie (sizeof(int) * alen);
  for (apos = rpos = 0; apos < alen; apos++)
    if (!isgap(aseq[apos])) 
      rmap[apos] = rpos++;
    else
      rmap[apos] = -1;
  
  /* Dealign the trace
   */
  for (tpos = 0; tpos < tr->tlen; tpos++)
    if (tr->rpos[tpos] != -1) 
      tr->rpos[tpos] = rmap[tr->rpos[tpos]];
  
  free(rmap);
}
