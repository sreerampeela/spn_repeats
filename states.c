/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* states.c
 * 
 * alloc, free, and initialization of state structures
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "squid.h"
#include "states.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


struct hmm_struc *
AllocHMM(int M)               		/* length of model to make */
{
  struct hmm_struc *hmm;        /* RETURN: blank HMM */
  int idx, k, ts;		/* counters */

				/* alloc for a new_hmm */
  if ((hmm = (struct hmm_struc *) malloc (sizeof(struct hmm_struc))) == NULL)
    return NULL;
  hmm->M = M;
				/* malloc and initialize optional info */
  if ((hmm->ref  = (char *) malloc ((M+2) * sizeof(char))) == NULL ||
      (hmm->cs   = (char *) malloc ((M+2) * sizeof(char))) == NULL ||
      (hmm->xray = (float *) malloc ((M+2) * sizeof(float) * NINPUTS)) == NULL)
    return NULL;
  hmm->flags = 0;

				/* malloc states; allow room for 2 dummies,
				   BEGIN and END */
  hmm->ins = (struct basic_state *) malloc (sizeof(struct basic_state) * (hmm->M + 2));
  hmm->del = (struct basic_state *) malloc (sizeof(struct basic_state) * (hmm->M + 2));
  hmm->mat = (struct basic_state *) malloc (sizeof(struct basic_state) * (hmm->M + 2));
  if (hmm->ins == NULL || hmm->del == NULL || hmm->mat == NULL) return NULL;

				/* initialize the counts to zero */
  for (k = 0; k <= hmm->M+1; k++)
    {
      for (ts = 0; ts < 3; ts++)
	{
	  hmm->mat[k].t[ts] = 0.0;
	  hmm->ins[k].t[ts] = 0.0;
	  hmm->del[k].t[ts] = 0.0;
	}
      for (idx = 0; idx < Alphabet_size; idx++)
	{
	  hmm->mat[k].p[idx]   = 0.0;
	  hmm->ins[k].p[idx]   = 0.0;
	  hmm->del[k].p[idx]   = 0.0;
	}
    }
  return hmm;
}

/* Function: WriteFlatPriorHMM()
 * 
 * Purpose:  Fill an HMM with expected probabilities according
 *           to a given prior. Used to construct "flat" initial
 *           models for hmmt.
 */
int
WriteFlatPriorHMM(struct hmm_struc *hmm, struct prior_s *prior)
{
  int    k;			/* counter across model                */
  int    q;			/* counter over mixtures               */
  int    x;			/* counter over symbols or transitions */
  float  malpha;		/* alpha for mixture                   */
  float  ialpha;		/* alpha for insert mixture            */
  float  dalpha;		/* alpha for delete mixture            */

  for (k = 0; k <= hmm->M; k++)
    {
				/* xray info for structure prior */
      if (prior->strategy == PRI_STRUCT)
	{
	  hmm->xray[k*NINPUTS + XRAY_bias] = 1.0;
	  hmm->xray[k*NINPUTS + XRAY_E]    = 0.0;
	  hmm->xray[k*NINPUTS + XRAY_H]    = 0.0;
	  hmm->xray[k*NINPUTS + XRAY_SA]   = 0.0;
	}
				/* match symbol emissions */
      for (x = 0; x < Alphabet_size; x++) 
	hmm->mat[k].p[x] = 0.0;
      if (k > 0)
	for (q = 0; q < prior->mnum; q++)
	  {
	    if (prior->strategy == PRI_STRUCT) 
	      prior->mq[q] = 1.0 / prior->mnum;
	    malpha = 0.0;
	    for (x = 0; x < Alphabet_size; x++)
	      malpha += prior->mat[q][x];
	    for (x = 0; x < Alphabet_size; x++)
	      hmm->mat[k].p[x] += prior->mq[q] * prior->mat[q][x] / malpha;
	  }
				/* insert emissions */
      for (x = 0; x < Alphabet_size; x++) 
	hmm->ins[k].p[x] = 0.0;
      for (q = 0; q < prior->inum; q++)
	{
	  if (prior->strategy == PRI_STRUCT) 
	    prior->iq[q] = 1.0 / prior->inum;
	  ialpha = 0.0;
	  for (x = 0; x < Alphabet_size; x++)
	    ialpha += prior->ins[q][x];	      
	  for (x = 0; x < Alphabet_size; x++)
	    hmm->ins[k].p[x] += prior->iq[q] * prior->ins[q][x] / ialpha;
	}

				/* state transitions  */
      for (x = 0; x < 3; x++)
	hmm->mat[k].t[x] = hmm->ins[k].t[x] = hmm->del[k].t[x] = 0.0;
      for (q = 0; q < prior->tnum; q++)
	{
	  if (prior->strategy == PRI_STRUCT) 
	    prior->tq[q] = 1.0 / prior->tnum;
	  malpha = ialpha = dalpha = 0.0;
	  for (x = 0; x < 3; x++)
	    {
	      malpha += prior->tm[q][x];
	      ialpha += prior->ti[q][x];
	      dalpha += prior->td[q][x];
	    }
	  for (x = 0; x < 3; x++)
	    {
	      hmm->mat[k].t[x] += prior->tq[q] * prior->tm[q][x] / malpha;
	      hmm->ins[k].t[x] += prior->tq[q] * prior->ti[q][x] / ialpha;
	      if (k > 0) hmm->del[k].t[x] += prior->tq[q] * prior->td[q][x] / dalpha;
	    }
	}
    }
				/* the final state never transits to d+1 */
  hmm->mat[hmm->M].t[DELETE] = 0.0;
  hmm->ins[hmm->M].t[DELETE] = 0.0;
  hmm->del[hmm->M].t[DELETE] = 0.0;
  Renormalize(hmm);
  return 1;
}


/* Function: HMMDup()
 * 
 * Purpose:  Create a duplicate copy of an HMM.
 * 
 * Return:   Pointer to the duplicate. 
 *           Caller is responsible for free'ing the duplicate.
 */
struct hmm_struc *
HMMDup(struct hmm_struc *hmm)
{
  struct hmm_struc *newhmm;
  int k, x, ts;

  if ((newhmm = AllocHMM(hmm->M)) == NULL)
    Die("AllocHMM() failed");

  newhmm->flags = hmm->flags;
  if (hmm->flags & HMM_REF)  
    strcpy(newhmm->ref, hmm->ref);
  if (hmm->flags & HMM_CS)   
    strcpy(newhmm->cs, hmm->cs);
  if (hmm->flags & HMM_XRAY) 
    memcpy(newhmm->xray, hmm->xray, NINPUTS * (hmm->M+2) * sizeof(float));

  for (k = 0; k <= hmm->M+1; k++)
    {
			/* copy transition T's */
      for (ts = 0; ts < 3; ts++)
	{
	  newhmm->mat[k].t[ts] = hmm->mat[k].t[ts];
	  newhmm->ins[k].t[ts] = hmm->ins[k].t[ts];
	  newhmm->del[k].t[ts] = hmm->del[k].t[ts];
	}
				/* copy symbol P table and regularizers */
      for (x = 0; x < Alphabet_size; x++) 
	{
	  newhmm->mat[k].p[x]   = hmm->mat[k].p[x];
	  newhmm->ins[k].p[x]   = hmm->ins[k].p[x];
	}
    }
  return newhmm;
}


int
FreeHMM(struct hmm_struc *hmm)
{
  if (hmm == NULL) return 0;
  free(hmm->ref);
  free(hmm->cs);
  free(hmm->xray);
  if (hmm->mat != NULL)  free (hmm->mat);
  if (hmm->ins != NULL)  free (hmm->ins);
  if (hmm->del != NULL)  free (hmm->del);
  free(hmm);
  return 1;
}


/* Function: CountSymbol()
 * 
 * Given an observed symbol, and a number of counts to
 * distribute (typically just 1.0), bump the appropriate counter(s).
 * 
 * This is completely trivial only so long as the symbols
 * always come from the expected alphabet; since we also
 * have to deal with degenerate symbols for both nucleic
 * acid and protein languages, we make a function to deal
 * with this.
 *
 * Returns 1 on success and bumps the necessary counters.
 * Returns 0 on failure and bumps each counter evenly, as
 * if it saw a completely ambiguous symbol; this lets
 * the caller silently accept garbage symbols, if it cares to.
 */
int
CountSymbol(char   sym,		/* observed symbol                        */
	    double wt,          /* number of counts to distribute (1.0)   */
	    float *counters)    /* array of 4 or 20 counters to increment */
{
  char *alphptr;                /* pointer into symbol in hmm->alphabet         */
  int   status;			/* RETURN: status; did we recognize the symbol? */
  int   i;

				/* trivial case: symbol is in alphabet */
  if ((alphptr = strchr(Alphabet, sym)) != NULL)
    {
      counters[alphptr - Alphabet] += wt;
      return 1;
    }
				/* non trivial case: symbol not in alphabet;
				   either degenerate symbol, or it's garbage */
  status = 1;
  if (Alphabet_type == kAmino)
    {
      switch (sym) {
      case 'B': 
	counters[SYMIDX('N')] += wt * 0.5;
	counters[SYMIDX('D')] += wt * 0.5;
	break;
      case 'Z':
	counters[SYMIDX('Q')] += wt * 0.5;
	counters[SYMIDX('E')] += wt * 0.5;
	break;
      default:
	Warn("unrecognized character %c (%d) in sequence\n", sym, (int) sym);
	status = 0;
				/* break thru to case 'X' */
      case 'X':
	for (i = 0; i < Alphabet_size; i++)
	  counters[i] += wt / (float) Alphabet_size;
	break;
      }
    }

  else if (Alphabet_type == kDNA || Alphabet_type == kRNA)
    {
				/* Deal with IUPAC code degeneracies. 
				   WARNING: Expects that the alphabet
				   is "ACGT" or "ACGU"; any other order 
				   will break this code! */
      switch (sym) {
      case 'B': counters[1] += wt/3.0; counters[2] += wt/3.0; counters[3] += wt/3.0; break;
      case 'D': counters[0] += wt/3.0; counters[2] += wt/3.0; counters[3] += wt/3.0; break;
      case 'H': counters[0] += wt/3.0; counters[1] += wt/3.0; counters[3] += wt/3.0; break;
      case 'K': counters[2] += wt/2.0; counters[3] += wt/2.0;                        break;
      case 'M': counters[0] += wt/2.0; counters[1] += wt/2.0;                        break;
      case 'R': counters[0] += wt/2.0; counters[2] += wt/2.0;                        break;
      case 'S': counters[1] += wt/2.0; counters[2] += wt/2.0;                        break;
      case 'T': counters[3] += wt;                                                   break;
      case 'U': counters[3] += wt;                                                   break;
      case 'V': counters[0] += wt/3.0; counters[1] += wt/3.0; counters[2] += wt/3.0; break;
      case 'W': counters[0] += wt/2.0; counters[3] += wt/2.0;                        break;
      case 'Y': counters[1] += wt/2.0; counters[3] += wt/2.0;                        break;
      default:
	Warn("unrecognized character %c (%d) in sequence\n", sym, (int) sym);
	status = 0;
				/* break thru to case 'N' */
      case 'N':
	for (i = 0; i < Alphabet_size; i++)
	  counters[i] += wt / (float) Alphabet_size;
	break;
      }
    }

  else
    {
      status = 0;
      Warn("unrecognized character %c (%d) in sequence\n", sym, (int) sym);
      for (i = 0; i < Alphabet_size; i++)
	counters[i] += wt / (float) Alphabet_size;
    }

  return status;
}





/* Function: HMMDistance()
 * 
 * Purpose:  Test two models for how different they are, using
 *           a simple squared difference measure on all homologous
 *           parameters. They must have the same architecture:
 *           i.e. check that newhmm->M == oldhmm->M before calling.
 *           
 * Args:     newhmm  - new HMM, probability form
 *           oldhmm  - old HMM, probability form
 *           
 * Return:   distance.
 */
float
HMMDistance(struct hmm_struc *newhmm, struct hmm_struc *oldhmm)
{
  int    k,x, ts;
  float  distance = 0.0;

  for (k = 0; k <= newhmm->M; k++)
    {
				/* state transition distances */
      if (k > 0)
	{			
	  for (ts = 0; ts < 3; ts++)
	    distance += SQR( 100. * (newhmm->del[k].t[ts] - oldhmm->del[k].t[ts]));
	}
      for (ts = 0; ts < 3; ts++)
	distance += SQR( 100. * (newhmm->mat[k].t[ts] - oldhmm->mat[k].t[ts]));
      for (ts = 0; ts < 3; ts++)
	distance += SQR( 100. * (newhmm->ins[k].t[ts] - oldhmm->ins[k].t[ts]));
      
				/* symbol emission distances */
      if (k > 0)
	for (x = 0; x < Alphabet_size; x++)
	  distance += SQR( 100. * (newhmm->mat[k].p[x] - oldhmm->mat[k].p[x]));
      for (x = 0; x < Alphabet_size; x++)
	  distance += SQR( 100. * (newhmm->ins[k].p[x] - oldhmm->ins[k].p[x]));
    }
  distance = sqrt(distance) / newhmm->M;
  return distance;
}



/* Function: VerifyHMM()
 * 
 * Purpose:  Make sure the probabilities in an HMM sum to 1
 *           where they should.
 * 
 * Return:   1 on success, 0 on failure.
 *           Prints warnings if the HMM is not OK.          
 */
void
VerifyHMM(struct hmm_struc *hmm)
{
  int    k;
  double tolerance = 0.00001;
  int    bad = 0;

  for (k = 0; k <= hmm->M; k++)
    {
      if ((k > 0 && fabs(FSum(hmm->del[k].t, 3) - 1.0) > tolerance) ||
	  (fabs(FSum(hmm->mat[k].t, 3) - 1.0) > tolerance)          ||
          (fabs(FSum(hmm->ins[k].t, 3) - 1.0) > tolerance)          ||
	  (k > 0 && fabs(FSum(hmm->mat[k].p, Alphabet_size) - 1.0) > tolerance) ||
	  (fabs(FSum(hmm->ins[k].p, Alphabet_size) - 1.0) > tolerance))
	bad++;
    }

  if (bad) Die("That HMM is no good, pal. %d distributions don't sum to 1.0", bad);
}




/* Function: Renormalize()
 * 
 * Normalize all P distributions so they sum to 1.
 * P distributions that are all 0, or contain negative
 * probabilities, are left untouched.
 * 
 * Returns 1 on success, or 0 on failure.
 */
void
Renormalize(struct hmm_struc *hmm)
{
  int    k;			/* counter for states                  */

  for (k = 0; k <= hmm->M ; k++)
    {
				/* match state transition frequencies */
      FNorm(hmm->mat[k].t, 3);
      FNorm(hmm->ins[k].t, 3);
      if (k > 0) FNorm(hmm->del[k].t, 3);

      if (k > 0) FNorm(hmm->mat[k].p, Alphabet_size);
      FNorm(hmm->ins[k].p, Alphabet_size);
    }
}


/* Function: HybridizeHMMs()
 * 
 * Purpose:  Stabilization of the (1-q) shepherd training rule.
 *           Graeme Mitchison has found that training can
 *           be stabilized by slowing down the rate of change of parameters.
 *           We do this with a term, damp_factor, which varies from 0 to 1.
 *           New probabilities are a combination of the old probability
 *           and the reestimated probability: 
 *             damp_factor * old + (1- damp_factor) * new
 *           The newhmm and oldhmm must have identical architectures.
 *             
 * Args:     newhmm      - HMM to be damped (probability form)
 *           oldhmm      - old HMM          (probability form)
 *           damp_factor - how much weight on old hmm parameters (0.95 works)
 *          
 * Return:   (void)
 *           newhmm probabilities are changed according to the damping rule.
 */        
void
HybridizeHMMs(struct hmm_struc *newhmm, struct hmm_struc *oldhmm, double damp_factor)
{
  int k;			/* counter for model nodes */
  int x;			/* counter for alphabet symbols */
  int ts;

  for (k = 0; k <= newhmm->M; k++)
    {
      for (ts = 0; ts < 3; ts++)
	{
	  newhmm->mat[k].t[ts] = damp_factor * oldhmm->mat[k].t[ts] +
	                         (1-damp_factor) * newhmm->mat[k].t[ts];
	  newhmm->ins[k].t[ts] = damp_factor * oldhmm->ins[k].t[ts] +
	                         (1-damp_factor) * newhmm->ins[k].t[ts];
	  newhmm->del[k].t[ts] = damp_factor * oldhmm->del[k].t[ts] +
	                         (1-damp_factor) * newhmm->del[k].t[ts];
	}
      for (x = 0; x < Alphabet_size; x++)
	{
	  newhmm->mat[k].p[x] = damp_factor * oldhmm->mat[k].p[x] +
	                     (1-damp_factor) * newhmm->mat[k].p[x];
	  newhmm->ins[k].p[x] = damp_factor * oldhmm->ins[k].p[x] +
	                     (1-damp_factor) * newhmm->ins[k].p[x];
	}
    }
  Renormalize(newhmm);		/* be safe */
}
