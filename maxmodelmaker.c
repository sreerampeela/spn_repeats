/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* maxmodelmaker.c
 * SRE, Mon Aug 16 09:57:38 1993
 * 
 * Maximum likelihood HMM construction, given a multiple
 * sequence alignment.
 * 
 * Assumes that insert emission scores are always zero.
 *
 * Supercedes the ad hoc model construction procedure
 * of modelmaker.c and the ad hoc model surgery procedure
 * of surgery.c
 * 
 **********************************************************
 *
 * Given a multiple sequence alignment of N columns, we
 * construct a maximally likely HMM, using a recursive
 * calculation. We can't do this by straightforward dynamic
 * programming because when we assign a new column to an
 * INS state, we change transition and emission probabilities
 * the INS state if we've assigned one or more previous
 * columns to it.
 * 
 * The recursion is as follows. Consider a column j. We
 * can account for column j by 1) MAT and DEL states or
 * 2) INS states. We will calculate and save the score
 * for emitting it with MAT states. This is:
 *
 *  S(j) =  max     { S(i) 
 *           0<=i<j   + log P(seqs i+1..j-1 | INS)  
 *                    + log P(seqs j | MAT)
 *                    + log T(seqs i..j | t)
 *                    + log Pr(MAT)           (prior on model architecture)
 *                    
 * In other words, we calculate the cost S(j) for each
 * possible previous MAT-assigned column. To do this, we
 * must calculate two emission tables and one transition
 * table. We calculate a MATCH emission vector for column
 * j; an INSERT emission vector for all columns i+1..j-1;
 * and a transition matrix between i and j. We calculate
 * the necessary probabilities by summing individual log prob's
 * for each subsequence i..j in the alignment.
 * 
 * We initialize S(0) = 0 (log 1).
 * 
 **********************************************************
 * Implementation-specific detail:
 * 
 * We can't precalculate much. We count emission vectors for
 * each column. We can't precalculate much transition info,
 * because I-> and ->I transitions can occur across multiple
 * columns in the alignment, so we don't bother.
 * 
 * For this version, all we calculate is an array
 * of M+1 1's or 0's, where a 1 at i means column i is
 * assigned to a match state. We could save
 * more info during the recursion, but it's not clear
 * how much time this might save (or even lose!) at the
 * expense of clarity and memory space.
 * 
 * The algorithm in pseudocode:
 * 
 * - alloc for emission counts   (emcount[1..M])
 *             match scores      (mscore[0..M+1])
 *             traceback indices (tback[0..M+1]) 
 *             emission vectors  (insij, matj)
 *             
 * - count emissions in each column (emcount[1..M])
 * 
 * - for each column j (1..M+1)
 *    - for each possible previous match assignment  0<=i<j
 *       - build transition table for i->j, regularize, logify (transij[3][3])
 *       - build emission vector for INS(i+1..j-1), regularize, logify (insij)
 *       - build emission vector for MAT(j), regularize, logify (matj)
 *       - newscore = S(i)
 *       - for each indiv. sequence
 *           - sum appropriate transij, insij, matj into newscore
 *       - if newscore > S(j)
 *           - S(j) = newscore
 *           - save traceback index, tback[j] = i
 *               
 * - traceback from tback[M+1] to get matassign[] vector of 1's and 0's       
 * 
 ************************************************************************* 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "float.h"

#include "squid.h"
#include "states.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


static int build_transij(char **aseqs, float *weights, int nseq, int M,
			 int i, int j, float transij[3][3]);
static float score_ij_transitions(char **aseqs, int nseq, int M,
				  int  i, int j, float transij[3][3]);

/* Function: Maxmodelmaker()
 * 
 * Purpose:  Construct a maximally likely HMM, given a multiple
 *           sequence alignment.
 *
 *           Obscure feature: Commas ',' in the alignment mark
 *           gaps which don't count. This annotation is produced
 *           by BlockRaggedEdgedAlignment().
 *           
 * Args:     aseqs   - the multiple sequence alignment (flushed)
 *           alen    - length of aseqs alignment
 *           nseq    - number of sequences in aseqs
 *           prior   - priors on parameters to use for model construction
 *           randomseq - random sequence model
 *           mpri    - prior on architecture: probability of new match node
 *           ret_hmm - RETURN: new hmm (counts form)
 *                      (call PriorifyHMM, Renormalize to convert to probabilities)
 *           ret_tr  - RETURN: array of fake tracebacks constructed
 *                      for the sequences. Pass NULL if you don't
 *                      want this array.
 * 
 * Return:   1 on success, 0 on failure; ret_hmm is also returned
 *           and must be free'd by the caller.
 */          
int
Maxmodelmaker(char **aseqs, AINFO *ainfo, int nseq, struct prior_s *prior, float *randomseq,
	      float mpri, struct hmm_struc **ret_hmm, struct trace_s  ***ret_tr)
{
  struct hmm_struc *hmm;	/* RETURN: new hmm       */
  int      i,j;			/* counters for columns  */
  int      idx;			/* counter for sequences */
  int      col;			/* counter over columns  */
  int      x,y;			/* generic counters */
  float  **emcount;             /* precounted emission vectors for all columns */
  float   *mscore;              /* saved S(j) values from recursion */
  float   *inp;			/* structural inputs, for structure priors */
  int     *tback;               /* traceback indices saved during recursion */
  float   *matj;                /* MAT emission vector, log P */
  float    transij[3][3];       /* transition matrix, log P   */
  int     *matassign;           /* MAT state assignments if 1, 0..M */
  float    newscore;
  struct trace_s **tr;          /* fake tracebacks constructed for each seq */
  int      k;
  int      spos;
  int      alen;
  float   *weights;

  /* Make sure we have all the info about the alignment
   * that we need.
   */
  alen = (ainfo->flags & AINFO_ALEN) ? ainfo->alen : strlen(aseqs[0]);
  if ((weights = (float *) malloc (sizeof(float) * nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    weights[idx] = (ainfo->sqinfo[idx].flags & SQINFO_WGT) ? ainfo->sqinfo[idx].weight : 1.0;
  if (prior->strategy == PRI_STRUCT)
    AnnotateAlignment(aseqs, nseq, ainfo, &inp);

  /* Allocations 
   */
  if ((emcount   = (float **) malloc (sizeof(float *) * (alen+1)))  == NULL ||
      (mscore    = (float *)  malloc (sizeof(float)   * (alen+2))) == NULL ||
      (tback     = (int *)    malloc (sizeof(int)     * (alen+2))) == NULL ||
      (matj      = (float *)  malloc (sizeof(float) * Alphabet_size)) == NULL ||
      (matassign = (int *)    malloc (sizeof(int)   * (alen+1))) == NULL)
    Die("malloc failed");
  for (i = 0; i <= alen; i++)
    if ((emcount[i] = (float *) calloc (Alphabet_size, sizeof(float))) == NULL)
      Die("calloc failed");

  /* Precalculations: count emission vectors for each column
   */
  for (i = 1; i <= alen; i++)
    for (idx = 0; idx < nseq; idx++)
      if (aseqs[idx][i-1] != ',' && !isgap(aseqs[idx][i-1]))
	CountSymbol(aseqs[idx][i-1], weights[idx], emcount[i]);
  mpri = LOG2(mpri);
  
  /* Initialization
   */
  mscore[0] = 0.0;
  tback[0]  = 0;

  /* Main recursion
   */
  for (j = 1; j <= alen+1; j++)
    {
				/* build emission vector for this
				 * MAT-assigned column j */
      if (j != alen+1)
	{
	  for (x = 0; x < Alphabet_size; x++)
	    matj[x] = emcount[j][x];
	  if (prior->strategy == PRI_STRUCT) 
	    StructurePerceptron(prior, inp + j*NINPUTS);
	  PriorifyMatchVector(matj, prior);
	  FNorm(matj, Alphabet_size);
				/* log odds */
	  for (x = 0; x < Alphabet_size; x++)
	    matj[x] = LOG2(matj[x] / randomseq[x]);
	}

      mscore[j] = -FLT_MAX;
      for (i = j-1; i >= 0; i--)
	{
	  if (mscore[i] < mscore[j]) continue;

				/* build transition matrix for this 
				 * column pair i,j */
	  if (! build_transij(aseqs, weights, nseq, alen, i, j, transij))
	    Die("build_transij failed");
	  if (prior->strategy == PRI_STRUCT) 
	    StructurePerceptron(prior, inp + i*NINPUTS);
	  PriorifyTransitionVectors(transij[MATCH], transij[INSERT], transij[DELETE], prior);
	  FNorm(transij[MATCH], 3);
	  FNorm(transij[INSERT], 3);
	  FNorm(transij[DELETE], 3);
	  for (x = 0; x < 3; x++)
	    for (y = 0; y < 3; y++)
	      transij[x][y] = LOG2(transij[x][y]);

				/* Record the score if it's the max so far */
	  newscore = mscore[i] + score_ij_transitions(aseqs, nseq, alen, i, j, transij);
	  if (newscore > mscore[j])
	    {
	      mscore[j] = newscore;
	      tback[j]  = i;
	    }
	} /* end loop over possible start points i */

				/* add in the score for emission of column j */
      if (j < alen+1)
	for (x = 0; x < Alphabet_size; x++)
	  mscore[j] += emcount[j][x] * matj[x];

				/* add in the architecture prior */
				/* note ad hoc weighting by nseq */
      if (j < alen+1)
	mscore[j] += mpri * nseq;

    } /* end recursion over all possible match assignments j */
  

  /* Traceback to determine MAT/INS assignment of columns
   */
  for (col = 0; col < alen+1; col++)
    matassign[col] = 0;

  col = alen+1;
  while (col != 0)
    {
      matassign[tback[col]] = 1;
      col = tback[col];
    }

  /* Now construct fake tracebacks for each sequence
   */
  if ((tr = (struct trace_s **) malloc (sizeof(struct trace_s *) * nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    {
      AllocTrace(alen+2, &tr[idx]); /* allow room for BEGIN, END dummies */

      k    = 0;
      spos = 0;
				/* BEGIN */
				/* if there's a blocked left exterior gap, no BEGIN */
      if (aseqs[idx][0] != ',')
	{
	  tr[idx]->nodeidx[spos]   = 0;
	  tr[idx]->statetype[spos] = MATCH;
	  tr[idx]->rpos[spos]      = -1;
	  spos++;
	}

      for (col = 1; col <= alen; col++)
	{
	  if (aseqs[idx][col-1] == ',')
	    {			/* ignore "blocked" exterior gaps */
	      if (matassign[col]) k++;
	      continue;		
	    }

	  else if (matassign[col] && ! isgap(aseqs[idx][col-1]))
	    {			/* MATCH */
	      k++;
	      tr[idx]->nodeidx[spos]   = k;
	      tr[idx]->statetype[spos] = MATCH;
	      tr[idx]->rpos[spos]      = col-1;
	      spos++;
	    }
	  else if (matassign[col])
	    {			/* DELETE */
	      k++;
	      tr[idx]->nodeidx[spos]   = k;
	      tr[idx]->statetype[spos] = DELETE;
	      tr[idx]->rpos[spos]      = -1;
	      spos++;
	    }
	  else if (! isgap(aseqs[idx][col-1]))
	    {			/* INSERT */
	      tr[idx]->nodeidx[spos]   = k;
	      tr[idx]->statetype[spos] = INSERT;
	      tr[idx]->rpos[spos]      = col-1;
	      spos++;
	    }
	}
				/* END */
				/* no END if there's a blocked right exterior gap */
      if (aseqs[idx][alen-1] != ',')
	{
	  tr[idx]->nodeidx[spos]   = k+1;
	  tr[idx]->statetype[spos] = MATCH;
	  tr[idx]->rpos[spos]      = alen;
	  spos++;
	}
      tr[idx]->tlen = spos;

    }

  /* And build the model from those tracebacks.
   * k happens to contain the number of match-assigned columns
   */
  if ((hmm = AllocHMM(k)) == NULL)
    Die("Failed to allocate HMM in maxmodelmaker()");
  for (idx = 0; idx < nseq; idx++)
    TraceCount(hmm, aseqs[idx], weights[idx], tr[idx]);

  /* Annotate the model with ref, cs lines
   */
  if (ainfo->flags & AINFO_RF)
    {
      hmm->ref[0] = ' ';
      for (spos = 1, col = 1; col <= alen; col++)
	if (matassign[col])
	  hmm->ref[spos++] = ainfo->rf[col-1]; 
      hmm->ref[spos] = '\0';
      hmm->flags |= HMM_REF;
    }
  if (ainfo->flags & AINFO_CS)
    {
      hmm->cs[0] = ' ';
      for (spos = 1, col = 1; col <= alen; col++)
	if (matassign[col])
	  hmm->cs[spos++] = ainfo->cs[col-1];
      hmm->cs[spos] = '\0';
      hmm->flags |= HMM_CS;
    }

  /* Set #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from ainfo->rf.
   */
  if (ainfo->flags & AINFO_RF) free(ainfo->rf);
  ainfo->flags |= AINFO_RF;
  ainfo->rf = (char *) MallocOrDie (sizeof(char) * (alen + 1));
  for (i = 0; i < alen; i++)
    ainfo->rf[i] = matassign[i+1] ? 'x' : '.';
  ainfo->rf[alen] = '\0';

  /* Structure annotation
   */
  if (prior->strategy == PRI_STRUCT)
    {
      memcpy(hmm->xray, inp, sizeof(float) * NINPUTS);  /* BEGIN annotation special cased */
      for (spos = 1, col = 1; col <= alen; col++)
	if (matassign[col])
	  {
	    memcpy(hmm->xray + spos*NINPUTS, inp + col*NINPUTS, sizeof(float)*NINPUTS);
	    spos++;
	  }
      hmm->flags |= HMM_XRAY;
    }


  /* Garbage collection and return.
   */
  Free2DArray(emcount, alen+1);
  free(mscore);
  free(tback);
  free(matj);
  free(matassign);
  free(weights);
  if (prior->strategy == PRI_STRUCT) free(inp);

  if (ret_tr) 
    *ret_tr = tr;
  else 
    { for (idx = 0; idx < nseq; idx++) FreeTrace(tr[idx]); free(tr); }
  *ret_hmm = hmm;
  return 1;
}
	      



/* Function: build_transij()
 * 
 * Purpose:  Given a choice for i and j MAT-assigned columns
 *           in the alignment, count transitions into transij[][].
 *           
 * Args:     aseqs    - flushed sequence alignment (seqs are 0..M-1)     
 *           weights  - weights assigned to sequences
 *           nseq     - number of sequences in aseqs           
 *           M        - number of columns in aseqs
 *           i        - starting MAT-assigned column
 *           j        - ending MAT-assigned column
 *           transij  - RETURN; new transition matrix
 *           
 * Return:   1 on success, 0 on failure.          
 */
static int
build_transij(char   **aseqs,
	      float   *weights,
	      int      nseq,
	      int      M,
	      int      i,
	      int      j,
	      float   transij[3][3])
{
  int    x,y;			/* counters for transij */
  int    idx;			/* counter for sequences */
  int    inscount;		/* count of inserts between i and j for a seq */
  int    pos;			/* position in alignment */


  /* Zero transij.
   */
  for (x = 0; x < 3; x++)
    for (y = 0; y < 3; y++)
      transij[x][y] = 0.0;


  /* For each sequence, count transitions.
   */
  for (idx = 0; idx < nseq; idx++)
    {
				/* if either i or j are blocked gaps, ignore
				   transitions involving them */
      if (i > 0   && aseqs[idx][i-1] == ',') continue;
      if (j < M+1 && aseqs[idx][j-1] == ',') continue;

      inscount = 0; 
      for (pos = i+1; pos < j; pos++)
	if (aseqs[idx][pos-1] != ',' && !isgap(aseqs[idx][pos-1])) inscount++;
      
      if (inscount > 0)
	{
	  if (i > 0 && isgap(aseqs[idx][i-1])) 
	    transij[DELETE][INSERT] += weights[idx];
	  else
	    transij[MATCH][INSERT]  += weights[idx];
	    
	  if  (j < M+1 && isgap(aseqs[idx][j-1])) 
	    transij[INSERT][DELETE] += weights[idx];
	  else
	    transij[INSERT][MATCH]  += weights[idx];

	  transij[INSERT][INSERT]   += weights[idx] * (inscount-1);
	}

      else
	{
	  if (i > 0 && isgap(aseqs[idx][i-1]))
	    {
	      if (j < M+1 && isgap(aseqs[idx][j-1])) 
		transij[DELETE][DELETE] += weights[idx];
	      else
		transij[DELETE][MATCH]  += weights[idx];
	    }
	  else
	    {
	      if (j < M+1 && isgap(aseqs[idx][j-1])) 
		transij[MATCH][DELETE] += weights[idx];
	      else
		transij[MATCH][MATCH]  += weights[idx];
	    }
	}
    }

  return 1;
}



/* Function: score_ij_transitions()
 * 
 * Purpose:  Given a choice for i and j MAT-assigned columns in
 *           the alignment, calculate the score that will be added
 *           to S(i) as a result of transitions from i to j.
 *           
 * Args:     aseqs     - flushed sequence alignment (seqs are 0..M-1)     
 *           nseq      - number of seqs in aseqs
 *           M         - number of columns in alignment
 *           i         - starting MAT-assigned column
 *           j         - ending MAT-assigned column
 *           transij   - transition matrix to use
 *           
 * Return:   the cost of making all transitions from i to j.
 */
static float
score_ij_transitions(char   **aseqs,
		     int      nseq,
		     int      M,
		     int      i,
		     int      j,
		     float   transij[3][3])
{
  float score;			/* RETURN: the score to add for transitions */
  int    idx;			/* counter for sequences */
  int    inscount;		/* count of inserts between i and j for a seq */
  int    pos;			/* position in alignment */

  score = 0.0;
  for (idx = 0; idx < nseq; idx++)
    {
      			/* if either i or j are blocked gaps, ignore
			   transitions involving them */
      if (i > 0   && aseqs[idx][i-1] == ',') continue;
      if (j < M+1 && aseqs[idx][j-1] == ',') continue;

      inscount = 0; 
      for (pos = i+1; pos < j; pos++)
	if (!isgap(aseqs[idx][pos-1])) inscount++;
      
      if (inscount > 0)
	{
	  score += (i > 0 && isgap(aseqs[idx][i-1])) ?
	    transij[DELETE][INSERT] : transij[MATCH][INSERT];
	    
	  score += (j < M+1 && isgap(aseqs[idx][j-1])) ?
	    transij[INSERT][DELETE] : transij[INSERT][MATCH];

	  score += transij[INSERT][INSERT] * (inscount-1);
	}

      else
	{
	  if (i > 0 && isgap(aseqs[idx][i-1]))
	    score += (j < M+1 && isgap(aseqs[idx][j-1])) ?
	      transij[DELETE][DELETE] : transij[DELETE][MATCH];
	  else
	    score += (j < M+1 && isgap(aseqs[idx][j-1])) ?
	      transij[MATCH][DELETE] : transij[MATCH][MATCH];
	}
    }

  return score;
}




