/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/
/* forback.c
 * added to 1.7.7 Tue Jan 24 13:35:42 1995
 * Previous forback.c, which has been obsolete for some time, moved to Archive/.
 * 
 * Forward-Backward algorithm for calculating likelihood of an HMM
 * given a sequence. Used for full Baum-Welch expectation maximization
 * and for evaulating confidence of positions in Viterbi alignments.
 * 
 * Derived from saviterbi.c. Related to the whole family of alignment
 * algorithms. Uses sahmm's, like simulated annealing does.
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

/* Function: Forward()
 *
 * Purpose:  Forward calculation.
 * 
 *           The fwd matrix is 0..L+1 rows by 0..M+1 columns for a
 *           sequence of length 1..L and a model of length 1..M.
 *           Each entry represents the scaled summed probability of getting
 *           into this state, starting from (0,0), inclusive of
 *           symbol emission.
 *           
 * Args:     hmm     - model, probability-form, single precision
 *           s       - sequence 0..L-1, upper case             
 *           ret_fwd - RETURN: forward calculation matrix
 *           ret_scl - RETURN: 0..L array of scaling factors
 *           
 * Return:   (void)
 *           fwd and scl are alloc'ed here, and must be free'd by caller:          
 *               Free2DArray(fwd, L+2) 
 *               free(scl) 
 */
void
Forward(struct hmm_struc *hmm, char *s, struct forback_s ***ret_fwd, float **ret_scl)
{
  char  *seq;                   /* sequence, 1..L */
  int    L;			/* length of seq  */
  struct forback_s **fwd;       /* the calculation grid */
  int    i;			/* counter for sequence position: 0,1..L,L+1 */
  int    k;			/* counter for model position: 0,1..M,M+1    */
  int    i_symidx;		/* hint for symbol index in alphabet */
  char  *alphptr;		/* ptr into alphabet (optimization) */
  float *scl;     		/* scaling factors */
  struct forback_s *thisrow;
  struct forback_s *nextrow;

  /********************************************
   * Initial setup and allocations
   ********************************************/
			/* convert sequence to 1..L for ease of indexing */
  L = strlen(s);
  seq = (char *) MallocOrDie((L+2) * sizeof(char));
  strcpy(seq+1, s);
  seq[0] = ' ';

			/* allocate the calculation matrix */
  fwd = (struct forback_s **) MallocOrDie (sizeof(struct forback_s *) * (L+2));
  for (i = 0; i <= L+1; i++)
    fwd[i] = (struct forback_s *) malloc (sizeof(struct forback_s) * (hmm->M + 2));

			/* allocate the scalefactors */
  scl = (float *) MallocOrDie ((L+1) * sizeof(float));

  /********************************************
   * Initialization
   ********************************************/

				/* set up the 0,0 cell */
  fwd[0][0].score_m = 1.0;
  fwd[0][0].score_d = 0.0;
  fwd[0][0].score_i = 0.0;
				/* initialize the top row */
  for (k = 1; k <= hmm->M; k++)
    {
      fwd[0][k].score_m = 0.0;
      fwd[0][k].score_i = 0.0;
    }

  /********************************************
   * Recursion: fill in the fwd matrix
   ********************************************/

  for (i = 0; i <= L; i++)
    {
      				/* get ptrs into current and next row. */
      thisrow = fwd[i];
      nextrow = fwd[i+1];
				/* initialize in the next row */
      nextrow[0].score_m = 0.0;
      nextrow[0].score_d = 0.0;
				/* optimization: get an index for this symbol
				   while we're in the outer loop. */
      if ((alphptr = strchr(Alphabet, seq[i])) != NULL)
	i_symidx = alphptr - Alphabet;
      else
	i_symidx = -1;		/* too bad; it's a degenerate symbol */

				/* First pass across all states.
				   Emission scores, delete states */
      for (k = 0; k <= hmm->M; k++)
	{
	  			/* add emission scores to the current cell. */
	  if (i_symidx >= 0)
	    {
	      thisrow[k].score_m *= hmm->mat[k].p[i_symidx];
	      thisrow[k].score_i *= hmm->ins[k].p[i_symidx];
	    }
	  else if (i > 0) /* watch out for the "gap" at position 0 */
	    {
	      thisrow[k].score_m *= ForbackSymscore(seq[i], hmm->mat[k].p, TRUE);
	      thisrow[k].score_i *= ForbackSymscore(seq[i], hmm->ins[k].p, TRUE);
	    }

				/* deal with transitions to delete state */
	  thisrow[k+1].score_d = thisrow[k].score_d * hmm->del[k].t[DELETE] +
	                         thisrow[k].score_i * hmm->ins[k].t[DELETE] +
                                 thisrow[k].score_m * hmm->mat[k].t[DELETE];
	}


				/* In preparation for the operations that affect
				   the next row, scale the current row. */
      scl[i] = 0.0;
      for (k = 0; k <= hmm->M; k++)
	{
	  if (thisrow[k].score_m > scl[i])  scl[i] = thisrow[k].score_m;
	  if (thisrow[k].score_d > scl[i])  scl[i] = thisrow[k].score_d;
	  if (thisrow[k].score_i > scl[i])  scl[i] = thisrow[k].score_i;
	}
      scl[i] = 1.0 / scl[i];
      
      for (k = 0; k <= hmm->M; k++)
	{
	  thisrow[k].score_m *= scl[i];
	  thisrow[k].score_d *= scl[i];
	  thisrow[k].score_i *= scl[i];
	}

				/* Now deal with insert and match transitions,
				   affecting the next row. */
      for (k = 0; k <= hmm->M; k++)
	{ 
				/* deal with transitions to insert state */
	  nextrow[k].score_i =  thisrow[k].score_d * hmm->del[k].t[INSERT] +
	                        thisrow[k].score_i * hmm->ins[k].t[INSERT] +
				thisrow[k].score_m * hmm->mat[k].t[INSERT];

				/* deal with transitions to match state */
	  nextrow[k+1].score_m = thisrow[k].score_d * hmm->del[k].t[MATCH] +
	                         thisrow[k].score_i * hmm->ins[k].t[MATCH] +
			         thisrow[k].score_m * hmm->mat[k].t[MATCH];
	}

    }

  /******************************************
   * Debugging information
   ******************************************/

#ifdef EXTREME_DEBUG
  DumpForbackMatrix(fwd, L, hmm->M, scl);
#endif /* EXTREME_DEBUG */

  /********************************************
   * Garbage collection and return: caller is responsible for free'ing fwd and scl!
   ********************************************/

  if (ret_scl != NULL) *ret_scl = scl;
  if (ret_fwd != NULL) *ret_fwd = fwd;
  free(seq);
}



/* Function: Backward()
 *
 * Purpose:  Backward calculation.
 * 
 *           The bck matrix is 0..L+1 rows by 0..M+1 columns for a
 *           sequence of length 1..L and a model of length 1..M.
 *           Each entry represents the scaled summed probability of getting
 *           into this state, starting from (0,0), inclusive of
 *           symbol emission. 
 *           
 *           The scalefactors used must be the same as the ones used for
 *           the forward calculation. This means Backward() can only
 *           be called subsequent to a Forward() call.
 *           
 * Args:     hmm     - model, probability form, single precision
 *           s       - sequence 0..L-1, upper case             
 *           scl     - 0..L array of scalefactors from the forward calculation.
 *           ret_bck - RETURN: forward calculation matrix
 *           
 * Return:   (void)
 *           bck and scl are alloc'ed here, and must be free'd by caller:
 *              Free2DArray(bck, L+2)
 *              free(scl)          
 */
void
Backward(struct hmm_struc *hmm, char *s, float *scl, struct forback_s ***ret_bck)
{
  char  *seq;                   /* sequence, 1..L */
  int    L;			/* length of seq  */
  struct forback_s **bck;       /* the calculation grid */
  int    i;			/* counter for sequence position: 0,1..L,L+1 */
  int    k;			/* counter for model position: 0,1..M,M+1    */
  int    i_symidx;		/* hint for symbol index in alphabet */
  char  *alphptr;		/* ptr into alphabet (optimization) */
  struct forback_s *thisrow;
  struct forback_s *nextrow;

  /********************************************
   * Initial setup and allocations
   ********************************************/
			/* convert sequence to 1..L for ease of indexing */
  L = strlen(s);
  seq = (char *) MallocOrDie((L+3) * sizeof(char));
  strcpy(seq+1, s);
  seq[0] = seq[L+1] = ' ';
  seq[L+2] = '\0';

			/* allocate the calculation matrix */
  bck = (struct forback_s **) MallocOrDie (sizeof(struct forback_s *) * (L+2));
  for (i = 0; i <= L+1; i++)
    bck[i] = (struct forback_s *) MallocOrDie (sizeof(struct forback_s) * (hmm->M + 2));

  /********************************************
   * Initialization
   ********************************************/

				/* set up the L+1,M+1 cell */
  bck[L+1][hmm->M+1].score_m = 1.0;
  bck[L+1][hmm->M+1].score_d = 0.0;
  bck[L+1][hmm->M+1].score_i = 0.0;
				/* initialize the last row */
  for (k = 0; k <= hmm->M; k++)
    {
      bck[L+1][k].score_m = 0.0;
      bck[L+1][k].score_i = 0.0;
    }

  /********************************************
   * Recursion: fill in the bck matrix
   ********************************************/

  for (i = L; i >= 0; i--)
    {
      				/* get ptrs into current and next row. */
      thisrow = bck[i];
      nextrow = bck[i+1];
				/* initialize in this row */
      thisrow[hmm->M+1].score_m = 0.0;
      thisrow[hmm->M+1].score_d = 0.0;
      thisrow[hmm->M+1].score_i = 0.0;

				/* optimization: get an index for symbol in next row 
				   while we're in the outer loop. */
      if ((alphptr = strchr(Alphabet, seq[i])) != NULL)
	i_symidx = alphptr - Alphabet;
      else
	i_symidx = -1;		/* too bad; it's a degenerate symbol */

				/* recursion, calculation of backward variable */
      for (k = hmm->M; k >= 0; k--)
	{
	  /* Transitions.
	   */
	  thisrow[k].score_m = nextrow[k+1].score_m * hmm->mat[k].t[MATCH] +
	                       nextrow[k].score_i   * hmm->mat[k].t[INSERT] +
                               thisrow[k+1].score_d * hmm->mat[k].t[DELETE];

	  thisrow[k].score_d = nextrow[k+1].score_m * hmm->del[k].t[MATCH] +
	                       nextrow[k].score_i   * hmm->del[k].t[INSERT] +
                               thisrow[k+1].score_d * hmm->del[k].t[DELETE];

	  thisrow[k].score_i = nextrow[k+1].score_m * hmm->ins[k].t[MATCH] +
	                       nextrow[k].score_i   * hmm->ins[k].t[INSERT] +
                               thisrow[k+1].score_d * hmm->ins[k].t[DELETE];

	  /* Emissions.
	   */
	  if (i_symidx >= 0) {
	    thisrow[k].score_m *= hmm->mat[k].p[i_symidx];
	    thisrow[k].score_i *= hmm->ins[k].p[i_symidx];
	  }
	  else if (i > 0) {
	    thisrow[k].score_m *= ForbackSymscore(seq[i], hmm->mat[k].p, TRUE);
	    thisrow[k].score_i *= ForbackSymscore(seq[i], hmm->ins[k].p, TRUE);
	  }
	}
				/* scale the current row */

      for (k = 0; k <= hmm->M; k++)
	{
	  thisrow[k].score_m *= scl[i];
	  thisrow[k].score_d *= scl[i];
	  thisrow[k].score_i *= scl[i];
	}
    }

  /******************************************
   * Debugging information
   ******************************************/

#ifdef EXTREME_DEBUG
  DumpForbackMatrix(bck, L, hmm->M, scl);
#endif /* EXTREME_DEBUG */

  /********************************************
   * Garbage collection and return: caller is responsible for free'ing bck and scl!
   ********************************************/

  if (ret_bck != NULL) *ret_bck = bck;
  free(seq);
}



/* Function: AlignmentConfidence()
 * 
 * Purpose:  Determines the probability distribution for a symbol i
 *           being aligned to 2M+1 possible insert plus match states.
 *           
 *           These probabilities are stored in a forback_s matrix.
 *           The score_d entries in this matrix are unused.
 *           
 * Args:     hmm      - model, probability form, single precision
 *           L        - length of sequence
 *           fwd      - calculated forward matrix
 *           bck      - calculated backward matrix
 *           scl      - scale factors used in fwd and bck
 *           ret_conf - RETURN: confidence probabilities for each symbol
 *           
 * Return:   void.
 *           conf is allocated here, must be free'd by caller.          
 */
void
AlignmentConfidence(struct hmm_struc *hmm, int L, 
		    struct forback_s **fwd, struct forback_s **bck, float *scl,
		    struct forback_s ***ret_conf)
{
  struct forback_s **conf;
  int i, k;
  float norm;

				/* allocate the calculation matrix */
  conf = (struct forback_s **) MallocOrDie (sizeof(struct forback_s *) * (L+2));
  for (i = 0; i <= L+1; i++)
    conf[i] = (struct forback_s *) MallocOrDie (sizeof(struct forback_s) * (hmm->M + 2));

  /* Sum.
   * Note the removal of a scalefactor from the delete path, to avoid double
   * counting it.
   */
  for (i = 0; i <= L; i++)
    for (k = 0; k <= hmm->M; k++)
      {
	conf[i][k].score_m = bck[i+1][k+1].score_m * hmm->mat[k].t[MATCH] +
                             bck[i][k+1].score_d   * hmm->mat[k].t[DELETE] / scl[i] +
                             bck[i+1][k].score_i   * hmm->mat[k].t[INSERT];
	conf[i][k].score_m *= fwd[i][k].score_m;
	
	conf[i][k].score_i = bck[i+1][k+1].score_m * hmm->ins[k].t[MATCH] +
	                     bck[i][k+1].score_d   * hmm->ins[k].t[DELETE] / scl[i] +
			     bck[i+1][k].score_i   * hmm->ins[k].t[INSERT];
	conf[i][k].score_i *= fwd[i][k].score_i;
      }

  /* Normalize across all match + insert states in row.
   */
  for (i = 0; i <= L; i++)
    {
      norm = 0.0;
      for (k = 0; k <= hmm->M; k++)
	{
	  norm += conf[i][k].score_m;
	  norm += conf[i][k].score_i;
	}
      for (k = 0; k <= hmm->M; k++)
	{
	  conf[i][k].score_m /= norm;
	  conf[i][k].score_i /= norm;
	}
    }

  /* Return.
   */
  *ret_conf = conf;
}


/* Function: TraceConfidence()
 * 
 * Purpose:  Find the confidence estimates on each aligned sequence symbol
 *           for a particular traceback, using the full matrix constructed
 *           by AlignmentConfidence().
 *           
 * Args:     cmx      - full confidence matrix produced by AlignmentConfidence()
 *           L        - length of seq (L+2 rows in cmx)
 *           tr       - traceback to evaulate symbol alignment confidences on
 *           ret_conf - RETURN: 0..L-1 array of confidence estimates 
 *                       for each position in original sequence.              
 */
void 
TraceConfidence(struct forback_s **cmx, int L, struct trace_s *tr, float **ret_conf)
{
  int    tpos;			/* position in trace */
  float *conf;                  /* 0..L-1 array of symbol alignment confidence values */

  conf = (float *) MallocOrDie (L * sizeof(float));

  for (tpos = 1; tpos < tr->tlen-1; tpos++)
    {
      switch (tr->statetype[tpos]) {
      case MATCH:
	conf[tr->rpos[tpos]] = cmx[tr->rpos[tpos]+1][tr->nodeidx[tpos]].score_m;
	break;
      case DELETE:
	break;
      case INSERT:
	conf[tr->rpos[tpos]] = cmx[tr->rpos[tpos]+1][tr->nodeidx[tpos]].score_i;
	break;
      default:
	Die("Bogus statetype %d in TraceConfidence()", tr->statetype[tpos]);
      }
    }
  *ret_conf = conf;
}


/* Function: ForbackCount()
 * 
 * Purpose:  Count probabilities of state transitions and symbol emissions from a
 *           forward-backward calculation on a single sequence. Results are
 *           stored in a counts-based HMM which will be normalized later.
 *           The full Baum-Welch equivalent to TraceCount() for Viterbi. 
 *           
 * Args:     hmm   - model. Probability form. 
 *           prior - Used only for alphabet information.
 *           seq   - sequence, 0..L-1
 *           L     - length of sequence
 *           fwd   - forward matrix produced by Forward()
 *           bck   - backward matrix produced by Backward()
 *           scl   - scale factors used by forward and backward
 *           count - counts-based HMM to accumulate counts in.
 *                   Caller is responsible for allocating this structure
 *                   and initializing it properly.
 *                   
 * Return:   (void)
 *           Some counts are added to the count-based HMM.
 */
void
ForbackCount(struct hmm_struc *hmm, char *seq, int L, float weight,
	     struct forback_s **fwd, struct forback_s **bck,
	     float *scl, struct hmm_struc *count)
{
  int i, k;
  float tmp_m1, tmp_d1, tmp_i0;
  float prob;			/* overall scaled probability of the seq | model */
  float wt;
  int   i_symidx;
  char *alphptr;

  /* "prob" is the total P(seq | model). We have to normalize by this.
   * We want to know, how much of the total alignment probability 
   * flows through a particular state or state transition?
   * 
   * Note the removal of a scale factor from the delete paths, to prevent
   * double-counting it in both the fwd and bck variable.
   */
  prob = fwd[L+1][hmm->M+1].score_m / weight;
  if (prob == 0.0) return;	/* silent! we have precision problems. */
  for (i = 0; i <= L; i++)
    {
				/* optimization: get an index for this symbol
				   while we're in the outer loop. */
      if (i > 0)
	{
	  if ((alphptr = strchr(Alphabet, seq[i-1])) != NULL)
	    i_symidx = alphptr - Alphabet;
	  else
	    i_symidx = -1;		/* too bad; it's a degenerate symbol */
	}

      for (k = 0; k <= hmm->M; k++)
      {
	tmp_m1 = fwd[i][k].score_m * hmm->mat[k].t[MATCH] * bck[i+1][k+1].score_m;
	tmp_d1 = fwd[i][k].score_m * hmm->mat[k].t[DELETE] * bck[i][k+1].score_d / scl[i];
	tmp_i0 = fwd[i][k].score_m * hmm->mat[k].t[INSERT] * bck[i+1][k].score_i;
	count->mat[k].t[MATCH] += tmp_m1 / prob;
	count->mat[k].t[DELETE] += tmp_d1 / prob;
	count->mat[k].t[INSERT] += tmp_i0 / prob;
	if (i > 0)
	  {
	    wt = ((tmp_m1 + tmp_d1 + tmp_i0) / prob);
	    if (i_symidx >= 0)
	      count->mat[k].p[i_symidx] += wt;
	    else
	      CountSymbol(seq[i-1], wt, count->mat[k].p);
	  }

	count->del[k].t[MATCH] += fwd[i][k].score_d * hmm->del[k].t[MATCH] * bck[i+1][k+1].score_m / prob;
	count->del[k].t[DELETE] += fwd[i][k].score_d * hmm->del[k].t[DELETE] * bck[i][k+1].score_d / (scl[i] * prob);
	count->del[k].t[INSERT] += fwd[i][k].score_d * hmm->del[k].t[INSERT] * bck[i+1][k].score_i / prob;
	
	tmp_m1 = fwd[i][k].score_i * hmm->ins[k].t[MATCH] * bck[i+1][k+1].score_m;
	tmp_d1 = fwd[i][k].score_i * hmm->ins[k].t[DELETE] * bck[i][k+1].score_d / scl[i];
	tmp_i0 = fwd[i][k].score_i * hmm->ins[k].t[INSERT] * bck[i+1][k].score_i;
	count->ins[k].t[MATCH] += tmp_m1 / prob;
	count->ins[k].t[DELETE] += tmp_d1 / prob;
	count->ins[k].t[INSERT] += tmp_i0 / prob;
	if (i > 0)
	  {
	    wt = ((tmp_m1 + tmp_d1 + tmp_i0) / prob);
	    if (i_symidx >= 0)
	      count->ins[k].p[i_symidx] += wt;
	    else
	      CountSymbol(seq[i-1], wt, count->ins[k].p);
	  }
      }
    }
}

/* Function: ForwardScore()
 * 
 * Purpose:  Given a forward matrix and the scalefactors,
 *           return the log likeihood log P(S | M) in base 2 (bits)
 */
float
ForwardScore(struct forback_s **fwd, int L, int M, float *scl)
{
  int    i;
  float  score;

  score = LOG2(fwd[L+1][M+1].score_m);
  for (i = 0; i <= L; i++)
    score -= LOG2(scl[i]);
  return score;
}


/* Function: BackwardScore()
 * 
 * Purpose:  Given a backward matrix and the scalefactors,
 *           return the log likeihood log P(S | M) in base 2 (bits)
 */
float
BackwardScore(struct forback_s **bck, int L, float *scl)
{
  int    i;
  float score;

  score = LOG2(bck[0][0].score_m);
  for (i = 0; i <= L; i++)
    score -= LOG2(scl[i]);
  return score;
}


/* Function: RandomScore()
 * 
 * Purpose:  Return the log P(S | R), the log likelihood of the
 *           sequence given the random model.
 */
float 
RandomScore(float *randomseq, char *seq)
{
  float  ll;			/* log likelihood */

  ll = 0.0;
  for (; *seq != '\0'; seq++)
    ll += LOG2(ForbackSymscore(*seq, randomseq, TRUE));

  return ll;
}
    

/* Function: ForbackSymscore()
 *
 * Look up the score of sequence character c, given a table of values
 * from the P(x | y) table for an HMM state. The table must be in
 * precisely the same order as the alphabet. For nucleic acids, both
 * the alphabet and the table must be in the order "ACGT".
 *
 * If hyperbayes is TRUE, it returns the summed probability of degenerate
 * symbols. This is formally correct but produces aberrantly high scores
 * on garbage degenerate sequences in database searches.
 * 
 * Returns the looked-up score from the table.
 */
float
ForbackSymscore(char x, float *scores, int hyperbayes)
{
  float result;
  int    count; 
				/* simple case: x is in the alphabet */
  if (strchr(Alphabet, x) != NULL) return scores[SYMIDX(x)];

  result = 0.0;
  if (Alphabet_type == kAmino)
    {
      switch (x) {
      case 'B': 
	result += scores[SYMIDX('N')];
	result += scores[SYMIDX('D')];
	count = 2;
	break;
      case 'Z':
	result += scores[SYMIDX('Q')];
	result += scores[SYMIDX('E')];
	count = 2;
	break;
      default:
	Warn("Unrecognized character %c (%d) in sequence", x, (int) x);
				/* break thru to case 'X' */
      case 'X':
	result = 1.0;
	count  = 20;
	break;
      }
    }
  else if (Alphabet_type == kDNA || Alphabet_type == kRNA)
    {
      switch (x) {
      case 'B':	result = scores[1] + scores[2] + scores[3]; count = 3; break;
      case 'D':	result = scores[0] + scores[2] + scores[3]; count = 3; break;
      case 'H': result = scores[0] + scores[1] + scores[3]; count = 3; break;
      case 'K': result = scores[2] + scores[3];             count = 2; break;
      case 'M': result = scores[0] + scores[1];             count = 2; break;
      case 'R': result = scores[0] + scores[2];             count = 2; break;
      case 'S': result = scores[1] + scores[2];             count = 2; break;
      case 'T': result = scores[3];                         count = 1; break;
      case 'U': result = scores[3];                         count = 1; break;
      case 'V': result = scores[0] + scores[1] + scores[2]; count = 3; break;
      case 'W': result = scores[0] + scores[3];             count = 2; break;      
      case 'Y': result = scores[1] + scores[3];             count = 2; break;
      default:
	Warn("unrecognized character %c (%d) in sequence", x, (int) x);
				/* break through to case 'N' */
      case 'N':
	result = 1.0; count = 4; break;
      }
    }
  else
    {
      Warn("unrecognized character %c (%d) in sequence\n", x, (int) x);
      result = 1.0;
      count  = Alphabet_size;
    }

  /* If you want correct probabilities, set hyperbayes TRUE and
   * use summed probabilities. But if you want it to work, and not
   * give high probabilities on garbage poly-X sequences, set
   * hyperbayes FALSE and use a simple expected probability.
   */

  if (! hyperbayes) result /= (float) count;
  return result;
}


/* Function: DumpForbackMatrix()
 * 
 * Purpose:  Debugging output; dump the contents of a matrix
 *           to stdout. Because it prints everything, it's
 *           only useful on small test matrices.
 *           
 * Returns:  (void)
 */
void
DumpForbackMatrix(struct forback_s **mx, int L, int M, float *scalefactors)
{
  int i;			/* counter for rows (sequence symbols) */
  int k;			/* counter for columns (states)        */
  
  for (i = 0; i <= L+1; i++)
    {
				/* print match states */
      for (k = 0; k <= M+1; k++) printf("%6.4g ", mx[i][k].score_m);
      printf(":: %.4g\n", scalefactors[i]);

				/* print delete states */
      for (k = 0; k <= M+1; k++) printf("%6.4g ", mx[i][k].score_d);
      putchar('\n');
				/* print insert states */
      for (k = 0; k <= M+1; k++) printf("%6.4g ", mx[i][k].score_i);
      putchar('\n');
      
      putchar('\n');
    }
}

