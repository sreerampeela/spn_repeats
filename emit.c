/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* emit.c
 * emission of raw sequence and tracebacks from a model
 * also, information content calculations
 *
 * SRE, Wed Jul 28 12:53:49 1993
 */

#include <stdio.h>
#include <stdlib.h>


#include "states.h"
#include "externs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: EmitSequence()
 * 
 * Purpose:  Given a model, emit a sequence from it. Also keeps traceback
 *           info if the caller wants it.
 *           
 * Args:     hmm          - the model 
 *           ret_seq      - RETURN: emitted sequence
 *           ret_tr       - RETURN: "traceback" of seq/model alignment
 *                            pass NULL if you don't want a trace
 *
 * Return:   1 on success, 0 on failure.
 *           Caller is reponsible for freeing ret_seq and ret_tr.
 */
int
EmitSequence(struct hmm_struc *hmm, char **ret_seq, struct trace_s **ret_tr)
{
  char *seq;
  struct trace_s *tr;
  int curridx;			/* current state index    */
  int currtype;			/* current state type     */
  int seqlen;			/* sequence length so far */
  int alloc_pathlen;		/* pathlen alloced so far */
  int alloc_seqlen;		/* seqlen alloced do far  */
  int   x;			/* counter for symbols     */
  char  sym;			/* emitted symbol          */

  /* Initial allocations. We use a blocksize of 64, and
   * realloc as the trace and the sequence grow.
   */
  AllocTrace(64, &tr);
  if ((seq       = (char *) malloc (64 * sizeof(char))) == NULL)
    Die("malloc failed");
  alloc_seqlen  = 64;
  alloc_pathlen = 64;

  curridx          = 0;
  currtype         = BEGIN;
  seqlen           = 0;
  tr->nodeidx[0]   = curridx;
  tr->statetype[0] = currtype;
  tr->rpos[0]      = -1;
  tr->tlen         = 1;
  while (curridx != hmm->M+1)
    {
      /* Transit to a new state, according to state transition
       * probabilities. Too wordy, but we're locked into this
       * unwieldy hmm structure declaration.
       */
      if (currtype == MATCH)
	{
	  currtype = FChoose(hmm->mat[curridx].t, 3);
	  if (currtype != INSERT) curridx++;
	}
      else if (currtype == INSERT)
	{
	  currtype = FChoose(hmm->ins[curridx].t, 3);
	  if (currtype != INSERT) curridx++;
	}
      else if (currtype == DELETE)
	{
	  currtype = FChoose(hmm->del[curridx].t, 3);
	  if (currtype != INSERT) curridx++;
	}

      /* Record the new state in the path. Make sure to
       * expand allocations if necessary.
       */
      tr->nodeidx[tr->tlen]   = curridx;
      tr->statetype[tr->tlen] = currtype;
      tr->tlen++;

      if (tr->tlen == alloc_pathlen)
	{
	  alloc_pathlen += 64;
	  ReallocTrace(tr, alloc_pathlen);
	}

      /* If new currtype is an emitting state, choose a symbol
       * according to the appropriate probability distribution
       */
      if (currtype == MATCH && curridx <= hmm->M)
	{
	  x = FChoose(hmm->mat[curridx].p, Alphabet_size);
	  if (x == Alphabet_size) 
	    Die("failed to choose symbol");
	  sym = Alphabet[x];
	}
      else if (currtype == INSERT)
	{
	  x = FChoose(hmm->ins[curridx].p, Alphabet_size);
	  if (x == Alphabet_size) 
	    Die("failed to choose symbol");
	  sym = sre_tolower((int) Alphabet[x]);
	}


      /* If currtype is an emitting state, insert the new symbol
       * we just chose into the sequence.
       */
      if ((currtype == MATCH || currtype == INSERT) && curridx <= hmm->M)
	{
	  if (seqlen+1 >= alloc_seqlen)
	    {
	      if ((seq = (char *) realloc (seq, (alloc_seqlen+64) * sizeof(char))) == NULL)
		Die("realloc failed");
	      alloc_seqlen += 64;
	    }
	  seq[seqlen]           = sym;
	  tr->rpos[tr->tlen-1]  = seqlen;
	  seqlen++;
	}
      else
	tr->rpos[tr->tlen-1] = -1;
    }
  
  /* Clean up: null terminate sequence, set passed ptrs for return
   */
  seq[seqlen] = '\0';

  *ret_seq = seq;
  if (ret_tr == NULL) FreeTrace(tr); else *ret_tr = tr;
  return 1;
}



/* Function: EmitBestSequence()
 * 
 * Purpose:  Given a model, emit the best sequence from it. Also keeps 
 *           traceback info if the caller wants it.
 * 
 *           This is the best sequence by Viterbi probability, 
 *           max P(seq | model). See RD_EmitBestSequence() for a 
 *           function that gives the best *scoring* sequence,
 *           max (P(seq | model) / P(seq | random)).
 *           
 * Args:     hmm          - the model 
 *           ret_seq      - RETURN: emitted sequence
 *           ret_tr       - RETURN: "traceback" of the model/seq alignment.
 *                            Pass NULL if this info isn't wanted.
 *
 * Return:   1 on success, 0 on failure.
 */
int
EmitBestSequence(struct hmm_struc *hmm,
		 char            **ret_seq,
		 struct trace_s  **ret_tr)
{
  struct trace_s *tr;
  struct vit_s   *mx;           /* 1D viterbi matrix      */
  char    *rev_m;               /* traceback pointers for matches */         
  char    *rev_i;               /* traceback pointers for inserts */
  char    *rev_d;               /* traceback pointers for deletes */
  char    *seq;
  int      x;			/* counter for symbols     */
  int      k;			/* counter for nodes       */
  double   sc;			/* temporary score         */
  int      bestsym;
  int      rpos;		/* position in a sequence  */
  int      spos;		/* position in a traceback */
  int      tlen;		/* length of the traceback */

  /* Initial allocations. 
   */
  if ((mx = (struct vit_s *) malloc (sizeof(struct vit_s) * (hmm->M+2))) == NULL ||
      (rev_m = (char *) malloc (sizeof(char) * (hmm->M+2))) == NULL ||
      (rev_d = (char *) malloc (sizeof(char) * (hmm->M+2))) == NULL ||
      (rev_i = (char *) malloc (sizeof(char) * (hmm->M+2))) == NULL)
    Die("malloc failed");

  /* Fill stage.
   */
  mx[0].score_m    = 0;
  mx[0].score_d    = -99999999;
  rev_m[0] = rev_d[0] = -1;
  for (k = 0; k <= hmm->M; k++)
    {
      /* Transits out of match state (init with these)
       */
      mx[k+1].score_m = mx[k].score_m + (int) (INTSCALE * LOG2(hmm->mat[k].t[MATCH]));
      mx[k].score_i   = mx[k].score_m + (int) (INTSCALE * LOG2(hmm->mat[k].t[INSERT]));
      mx[k+1].score_d = mx[k].score_m + (int) (INTSCALE * LOG2(hmm->mat[k].t[DELETE]));
      rev_m[k+1] = rev_i[k] = rev_d[k+1] = MATCH;
      
      /* Check transits out of delete state
       */
      if ((sc = mx[k].score_d + (int) (INTSCALE * LOG2(hmm->del[k].t[MATCH]))) > mx[k+1].score_m)
	{ mx[k+1].score_m = sc; rev_m[k+1] = DELETE; }
      if ((sc = mx[k].score_d + (int) (INTSCALE * LOG2(hmm->del[k].t[INSERT]))) > mx[k].score_i)
	{ mx[k].score_i = sc; rev_i[k] = DELETE; }  
      if ((sc = mx[k].score_d + (int) (INTSCALE * LOG2(hmm->del[k].t[DELETE]))) > mx[k+1].score_d)
	{ mx[k+1].score_d = sc; rev_d[k+1] = DELETE; }

      /* Check transits out of insert state
       * (looping on insert state is never good)
       */
      if ((sc = mx[k].score_i + (int) (INTSCALE * LOG2(hmm->ins[k].t[MATCH]))) > mx[k+1].score_m)
	{ mx[k+1].score_m = sc; rev_m[k+1] = INSERT; }
      if ((sc = mx[k].score_i + (int) (INTSCALE * LOG2(hmm->ins[k].t[DELETE]))) > mx[k+1].score_d)
	{ mx[k+1].score_d = sc; rev_d[k+1] = INSERT; }
    }
      
  /* Traceback
   * We know the traceback must be <= 2M+1, plus 2 for BEGIN and END.
   */
  AllocTrace(hmm->M*2+3, &tr);
  k    = hmm->M+1;
  tr->statetype[0] = MATCH; /* END */
  tr->nodeidx[0]   = k;
  spos = 1;
  while (k > 0 || tr->statetype[spos-1] == INSERT)
    {
      switch (tr->statetype[spos-1]) {
      case MATCH:
	tr->statetype[spos] = rev_m[k];
	tr->nodeidx[spos]   = --k;
	break;
      case INSERT:
	tr->statetype[spos] = rev_i[k];
	tr->nodeidx[spos]   = k;
	break;
      case DELETE:
	tr->statetype[spos] = rev_d[k];
	tr->nodeidx[spos]   = --k;
	break;
      }
      spos++;
    }
  tlen = spos;
  ReverseTrace(tr, tlen);	/* flip over to 0..tlen-1 order */


  /* Now, from that trace, emit a sequence.
   * We know the sequence is <= tlen-2 characters long.
   */
  if ((seq = (char *) malloc ((tlen-1) * sizeof(char))) == NULL)
    Die("malloc failed");
  for (spos = 1, rpos = 0; spos < tlen-1; spos++)
    {
      switch (tr->statetype[spos]) {
      case MATCH:
	bestsym = 0;
	for (x = 1; x < Alphabet_size; x++)
	  if (hmm->mat[tr->nodeidx[spos]].p[x] > hmm->mat[tr->nodeidx[spos]].p[bestsym])
	    bestsym = x;
	seq[rpos]      = Alphabet[bestsym];
	tr->rpos[spos] = rpos;
	rpos++;
	break;
	
      case INSERT:
	bestsym = 0;
	for (x = 1; x < Alphabet_size; x++)
	  if (hmm->ins[tr->nodeidx[spos]].p[x] > hmm->ins[tr->nodeidx[spos]].p[bestsym])
	    bestsym = x;
	seq[rpos]      = sre_tolower((int) Alphabet[bestsym]);
	tr->rpos[spos] = rpos;
	rpos++;
	break;

      case DELETE:
	tr->rpos[spos] = -1;
	break;
      }
    }

  /* Clean up: null terminate sequence, set passed ptrs for return
   */
  tr->rpos[0]      = -1;
  tr->rpos[tlen-1] = rpos;
  seq[rpos]        = '\0';
  *ret_seq         = seq;
  if (ret_tr == NULL) FreeTrace(tr); else *ret_tr = tr;
  free(mx); free(rev_m); free(rev_i); free(rev_d);
  return 1;
}


/* Function: RD_EmitBestSequence()
 * 
 * Purpose:  Given a model, emit the sequence with highest Viterbi match score.
 *            Also keeps traceback info if the caller wants it.
 *           (Richard and I coincidentally wrote versions of EmitBestSequence()
 *            on the same day. This is Richard's version. The difference
 *            is that Richard's calculates the best scoring sequence, using
 *            log odds, and I calculate the highest probability sequence)
 *           
 * Args:     hmm          - the model 
 *           randomseq    - random sequence model
 *           ret_seq      - RETURN: emitted sequence
 *           ret_tr       - RETURN: "traceback" of model/seq alignment
 *                            pass NULL if this info isn't wanted
 *
 * Return:   1 on success, 0 on failure.
 *           Caller is responsible for freeing ret_seq and ret_tr.
 */
int
RD_EmitBestSequence(struct hmm_struc *hmm,
		    float            *randomseq,
		    char            **ret_seq,
		    struct trace_s  **ret_tr)
{
  struct shmm_s    *shmm;       /* model, log odds scoring form */
  struct vit_s     *rx ;	/* row of viterbi cells - one per state */
  int   i, k, x ;		/* indices */
  int   score, best;		/* tmp variables */
  int   next_k, xbest ;		/* what they say they are*/
  char  *back_d, *back_i, *back_m ; /* back pointers */
  struct trace_s *tr;           /* traceback, initially written backwards */
  char *seq;			/* sequence to be returned */
  int   N, L;			/* lengths of optimal state path, sequence */
  int  *tptr;

  shmm = AllocSearchHMM(hmm->M);
  MakeSearchHMM(hmm, randomseq, shmm);

  if ((rx = (struct vit_s *) malloc (sizeof(struct vit_s) * (shmm->M + 2) )) == NULL) 
    Die ("memory failure allocating viterbi row") ;

  if (!(back_d = (char*) malloc (sizeof(char) * (shmm->M+2))))
    Die ("memory failure allocating back pointers") ;
  if (!(back_i = (char*) malloc (sizeof(char) * (shmm->M+2))))
    Die ("memory failure allocating back pointers") ;
  if (!(back_m = (char*) malloc (sizeof(char) * (shmm->M+2))))
    Die ("memory failure allocating back pointers") ;

	/* first run forward in a bastardisation of ViterbiFill,
	   but in this case keep back pointers, for simplicity
	*/
				/* initialize the first cell 0,0 */
  rx[0].score_m = 0;
  rx[0].score_d = -99999999;
  rx[0].score_i = -99999999;

  tptr = shmm->t;
  for (k = 0; k <= shmm->M; k++, tptr += 9)
    { 
      next_k = k+1 ;
				/* match emission first */
      if (k > 0)
	{ best = -99999999;
	  for (x = 0; x < 26; x++)
	    if (shmm->m_emit[x][k] > best)
	      best = shmm->m_emit[x][k];
	  rx[k].score_m += best ;
	}

      /* The exact order of the next nine score calculations is
       * totally dependent on the arrangement of state transitions
       * in the search model (tptr)
       */
				/* deal with transitions out of delete state */
				/* we initialize with these */
      rx[next_k].score_d = rx[k].score_d + *(tptr + Tdd); back_d[next_k] = DELETE;
      rx[k].score_i      = rx[k].score_d + *(tptr + Tdi); back_i[k]      = DELETE; 
      rx[next_k].score_m = rx[k].score_d + *(tptr + Tdm); back_m[next_k] = DELETE; 

				/* deal with transitions out of match state */
				/* must do match before insert in this case */
      if ((score = rx[k].score_m + *(tptr + Tmd)) > rx[next_k].score_d)
	{ rx[next_k].score_d = score; back_d[next_k] = MATCH ; }
      if ((score = rx[k].score_m + *(tptr + Tmi)) > rx[k].score_i)
	{ rx[k].score_i = score; back_i[k] = MATCH ; }
      if ((score = rx[k].score_m + *(tptr + Tmm)) > rx[next_k].score_m)
	{ rx[next_k].score_m = score; back_m[next_k] = MATCH ; }

				/* now do insert emission */
      best = -99999999;
      for (x = 0; x < 26; x++)
	if (shmm->i_emit[x][k] > best)
	  best = shmm->i_emit[x][k];
      rx[k].score_i += best;

				/* deal with transitions out of insert state */
      if ((score = rx[k].score_i + *(tptr + Tid)) > rx[next_k].score_d)
	{ rx[next_k].score_d = score; back_d[next_k] = INSERT ; }
/* Avoid multiple inserts - if they have positive score, we loop forever
      if ((score = rx[k].score_i + *tptr) > rx[k].score_i)
	rx[k].score_i = score;
*/
      if ((score = rx[k].score_i + *(tptr + Tim)) > rx[next_k].score_m)
	{ rx[next_k].score_m = score; back_m[next_k] = INSERT ; }
    }

	/* now trace back using a our pointers, with help from ViterbiTrace code */
				/* we know N <= 2M+1, so alloc accordingly */
				/* leaving space for dummy BEGIN and END  */
  AllocTrace(2*hmm->M + 3, &tr);
				/* start in dummy state at END */
  tr->nodeidx[0]   = hmm->M+1;
  tr->statetype[0] = MATCH;
  k = hmm->M+1;
  N = 1;

  /**************************************************
   * Traceback
   **************************************************/

  while (k > 0 || tr->statetype[N-1] == INSERT)
    {				/* if we look back one position in tmp_statetype,
				   we find out what substate we're supposed to use */
      switch (tr->statetype[N-1]) {

      case MATCH:
	tr->statetype[N] = back_m[k] ;
	tr->nodeidx[N]   = --k ;
	break ;

      case DELETE:
	tr->statetype[N] = back_d[k] ;
	tr->nodeidx[N]   = --k ;
	break ;

      case INSERT:
	tr->statetype[N] = back_i[k] ;
	tr->nodeidx[N]   = k ;
      }
      N++ ;
    }
  
  /**************************************************
   * Reverse the state path, and find the emitted sequence
   **************************************************/

  ReverseTrace(tr, N);
  L = 0 ;
  if (!(seq = (char*) malloc (N)))
    Die ("Couldn't allocate memory for the sequence") ;
  for (i = 0; i < N ; i++)
    { 
      k = tr->nodeidx[i];
      if (i == 0)        { tr->rpos[i] = -1; continue; } /* BEGIN */
      if (k == hmm->M+1) { tr->rpos[i] = L;  continue; } /* END   */
      
      if (tr->statetype[i] == MATCH)
	{ best = shmm->m_emit[0][k] ; xbest = 0 ;
	  for (x = 1; x < 26; x++)
	    if (shmm->m_emit[x][k] > best)
	      { best =  shmm->m_emit[x][k]; xbest = x ; }
	  tr->rpos[i] = L;
	  seq[L++] = (char) 'A' + xbest;

	}
      else if (tr->statetype[i] == INSERT)
	{ best = shmm->i_emit[0][k] ; xbest = 0 ;
	  for (x = 1; x < 26; x++)
	    if (shmm->i_emit[x][k] > best)
	      { best = shmm->i_emit[x][k]; xbest = x ; }
	  tr->rpos[i] = L;
	  seq[L++] = sre_tolower((char) 'A' + xbest) ;
	}
    }
  seq[L] = '\0';

  /**************************************************
   * Garbage collection and return
   **************************************************/

  free (rx) ;
  free (back_i) ;
  free (back_d) ;
  free (back_m) ;

  *ret_seq = seq;
  if (! ret_tr) FreeTrace(tr); else *ret_tr = tr;
  FreeSearchHMM(shmm);
  return 1;
}



/* Function: HMMInfoContent()
 * 
 * Purpose:  Calculate the information content of an HMM,
 *           by a sampling numerical integration.
 *           
 * Args:     hmm   - model, in probability form           
 *           shmm  - model, in search (integer log-odds) form
 *           
 * Return:   information content of model, in bits.          
 */
double
HMMInfoContent(struct hmm_struc *hmm, struct shmm_s *shmm)
{
  char             *seq;            /* an emitted (sampled) sequence */
  struct trace_s   *tr;             /* trace of an emitted sequence  */
  float             score;
  float             info;
  int               i;
  int               maxi = 500; /* number of sampled sequences */
  
  /* Do a sampling integration of the model info content.
   * This is a gross approximation, because we do not re-align.
   */
  info = 0.0;
  for (i = 0; i < maxi; i++)
    {
      EmitSequence(hmm, &seq, &tr);
      TraceScore(shmm, seq, tr, &score);
      info    += score;
      free(seq);
      FreeTrace(tr);
    }
  return info / (double) maxi;
}



