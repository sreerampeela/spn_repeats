/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* saviterbi.c: created Sun Dec 20 11:39:34 1992 SRE
 * deleted and completely revised, Tue Mar  2 08:43:55 1993 SRE
 * (modified from the core routines in viterbi.c)
 * theory due to Graeme Mitchison
 *
 * Implementation of a simulated annealing variant of the
 * state sequence calculation.
 *
 * Choose a state sequence somewhat randomly, rather than
 * according to the Viterbi maximum likelihood criterion.
 * The degree of randomness is controlled by a "temperature"
 * parameter, kT. At kT == 0 ("absolute zero"), the algorithm
 * reverts to Viterbi behaviour and freezes into choosing the
 * single most likely state sequence. At kT == oo, all paths
 * are equally likely. At kT == 1, paths are chosen with
 * a probability equal to their actual calculated probability.
 *
 * There are a number of pitfalls in trying to devise such
 * a procedure, and my first implementation fell into several.
 * One cannot simply randomize the choices of
 * path during a traceback through a matrix calculated
 * by the Viterbi criterion or the forward criterion of a full-blown
 * Baum-Welch.
 *
 * Graeme nonchalantly suggested the following piece of
 * brilliance. Whereas the forward calculation is:
 * 		_             _
 *     	       | _N_ 	       |
 *    S   (j) =| \   S (i) a   | b (O	)
 *     t+1     | /__  t	    ij |  j  t+1
 * 	       |_i=1 	      _|
 *
 * and the Viterbi approximation is:
 *  	                     _	      _
 *     	           	    | 	       |
 *    S   (j) =   argmax    | S (i) a  |  b (O   )
 *     t+1          1<=i<=N |  t     ij|   j  t+1
 * 	           	    |_        _|
 *
 * Graeme's version is:
 *
 *
 *  		_      _       _ 1/kT   _ kT
 *     	       | _N_  |	        |        |
 *    S   (j) =| \    |S (i) a  |        |    b (O   )
 *     t+1     | /__  | t     ij|        |     j  t+1
 * 	       |_i=1  |_       _|       _|
 *
 * This becomes the forward calculation for kT = 1, and becomes the 
 * Viterbi calculation in the limit kT = 0.
 *
 * The chief difficulty is that this introduces appreciable
 * scaling problems. Non-careful choice of extreme kT's can
 * drive the intermediates in this calculation out of the
 * dynamic range of any computer with ease.
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


/* Function: SaFill()
 *
 * Fill stage of the simulated annealing dynamic programming
 * calculation. Similar to the forward calculation of the Baum-Welch
 * algorithm.
 *
 * ret_mx is filled for the caller, which will typically call SaTrace()
 * (the same traceback function ViterbiFill-generated matrices would call)
 * to extract a state path.
 *
 * Returns 1 on success, 0 on failure.
 */
int
SaFill(struct sa_hmm_s *sahmm,    /* model, in simulated annealing form  */
       char   *s,                 /* sequence 0..L-1                     */
       struct sa_s ***ret_mx)     /* RETURN: the calc'ed matrix          */
{
  char  *seq;                   /* sequence, 1..L */
  int    L;			/* length of seq  */
  struct sa_s **mx;             /* the calculation grid */
  int    i;			/* counter for sequence position: 0,1..L */
  int    k;			/* counter for model position: 0,1..M    */
  int    i_symidx;		/* hint for symbol index in alphabet */
  int    next_k;		/* next k = k+1 (optimization) */
  char  *alphptr;		/* ptr into alphabet (optimization) */
  double *scalefactor;		/* scaling factors */
  struct sa_s *thisrow;
  struct sa_s *nextrow;

  /********************************************
   * Initial setup and allocations
   ********************************************/
				/* convert sequence to 1..L for ease of 
				   indexing, and make sure of upper case */
  L = strlen(s);
  if ((seq = (char *) malloc (L + 2)) == NULL) return 0;
  strcpy(seq+1, s);
  s2upper(seq+1);
  seq[0] = ' ';
				/* allocate the calculation matrix,
				   which is 0..L+1 rows by 0..M+1 cols */
  if (( mx = (struct sa_s **) malloc
       (sizeof(struct sa_s *) * (L+2) )) == NULL)
    Die("memory failure allocating sa matrix");

  for (i = 0; i <= L+1; i++)
    if ((mx[i] = (struct sa_s *) malloc	 (sizeof(struct sa_s) * (sahmm->M + 2) )) == NULL)
      Die("memory failure allocating viterbi matrix, row %d", i);


				/* allocate the scalefactors */
  if ((scalefactor = (double *) calloc (L+1, sizeof(double))) == NULL)
    { fprintf(stderr, "error: calloc() failed\n");  return 0; }

  /********************************************
   * Initialization
   ********************************************/

				/* set up the 0,0 cell */
  mx[0][0].score_m = 1.0;
  mx[0][0].score_d = 0.0;
  mx[0][0].score_i = 0.0;
				/* initialize the top row */
  for (k = 1; k <= sahmm->M; k++)
    {
      mx[0][k].score_m = 0.0;
      mx[0][k].score_i = 0.0;
    }

  /********************************************
   * Recursion: fill in the mx matrix
   ********************************************/

  for (i = 0; i <= L; i++)
    {
      				/* get ptrs into current and next row. */
      thisrow = mx[i];
      nextrow = mx[i+1];
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
      for (k = 0; k <= sahmm->M; k++)
	{
	  next_k = k+1;

	  			/* add emission scores to the current cell. */
	  if (i_symidx >= 0)
	    {
	      thisrow[k].score_m *= sahmm->mat[k].p[i_symidx];
	      thisrow[k].score_i *= sahmm->ins[k].p[i_symidx];
	    }
	  else if (i > 0) /* watch out for the "gap" at position 0 */
	    {
	      thisrow[k].score_m *= SaSymscore(seq[i], sahmm->mat[k].p, TRUE);
	      thisrow[k].score_i *= SaSymscore(seq[i], sahmm->ins[k].p, TRUE);
	    }

				/* deal with transitions to delete state */
	  thisrow[next_k].score_d = thisrow[k].score_d * sahmm->del[k].t[DELETE] +
	                            thisrow[k].score_i * sahmm->ins[k].t[DELETE] +
				    thisrow[k].score_m * sahmm->mat[k].t[DELETE];
	}


				/* In preparation for the operations that affect
				   the next row, scale the current row. */
      scalefactor[i] = 0.0;
      for (k = 0; k <= sahmm->M; k++)
	{
	  if (thisrow[k].score_m > scalefactor[i])  scalefactor[i] = thisrow[k].score_m;
	  if (thisrow[k].score_d > scalefactor[i])  scalefactor[i] = thisrow[k].score_d;
	  if (thisrow[k].score_i > scalefactor[i])  scalefactor[i] = thisrow[k].score_i;
	}
      scalefactor[i] = 1.0 / scalefactor[i];
      
      for (k = 0; k <= sahmm->M; k++)
	{
	  thisrow[k].score_m *= scalefactor[i];
	  thisrow[k].score_d *= scalefactor[i];
	  thisrow[k].score_i *= scalefactor[i];
	}



				/* Now deal with insert and match transitions,
				   affecting the next row. */
      for (k = 0; k <= sahmm->M; k++)
	{ 
	  next_k = k+1;
				/* deal with transitions to insert state */
	  nextrow[k].score_i =  thisrow[k].score_d * sahmm->del[k].t[INSERT] +
	                        thisrow[k].score_i * sahmm->ins[k].t[INSERT] +
				thisrow[k].score_m * sahmm->mat[k].t[INSERT];

				/* deal with transitions to match state */
	  nextrow[next_k].score_m = thisrow[k].score_d * sahmm->del[k].t[MATCH] +
	                            thisrow[k].score_i * sahmm->ins[k].t[MATCH] +
				    thisrow[k].score_m * sahmm->mat[k].t[MATCH];
	}

    }

  /******************************************
   * Debugging information
   ******************************************/

#ifdef EXTREME_DEBUG
  DumpSaMatrix(mx, L, sahmm->M, scalefactor);
#endif /* EXTREME_DEBUG */

  /********************************************
   * Garbage collection and return: caller is responsible for free'ing mx!
   ********************************************/

  free(scalefactor);
  free(seq);
  *ret_mx = mx;
  return 1;
}




/* Function: SaTrace()
 * 
 * Purpose:  Probabilistic traceback through a matrix constructed by SaFill().
 * 
 * Arguments:
 *       mx     - calculation lattice filled by SaFill()
 *       L      - length of sequence. mx has L+2 rows, 0..L+1
 *       sahmm  - the model in simulated annealing form        
 *       ret_tr - RETURN: traceback
 *       
 * Return:  1 on success, 0 on failure. 
 *          ret_tr must be freed by the caller, using FreeTrace().
 */
int
SaTrace(struct sa_s **mx, int L, struct sa_hmm_s *sahmm, struct trace_s **ret_tr)
{
  int   i;			/* counter for rows, 0..L+1 */
  int   k;			/* counter for cols, 0..M+1 */
  int   j;			/* counter for tmp state sequence, 0..L+M+1 */
  double score[3];		/* temp variables for scores */
  struct trace_s *tr;           /* traceback */
  int   N;			/* length of optimal state seq, for return */

  /**************************************************
   * Allocations of temporary space
   **************************************************/

				/* we know N <= L+M+1, so alloc accordingly */
				/* leaving space for dummy BEGIN and END  */
  AllocTrace(L+sahmm->M+3, &tr);

  /**************************************************
   * Initialization
   **************************************************/

				/* start in dummy state at END, aligned to L+1 */
  tr->nodeidx[0]   = sahmm->M+1;
  tr->statetype[0] = MATCH;
  tr->rpos[0]      = -1;
  i = L+1;		
  k = sahmm->M+1;
  j = 1;

 /**************************************************
   * Traceback
   **************************************************/

  while (i > 0 || k > 0)
    {
				/* if we look back one position in tmp_statetype,
				   we find out what substate we're supposed to use
				   here in i,k */
      switch (tr->statetype[j-1]) {

				/* if we're using the MATCH substate here in i,k,
				   we came from some substate in i-1, k-1. */
      case MATCH:		
	if (i == 0 || k == 0) return 0;	/* trace failed! */

	score[FROM_MATCH]  = mx[i-1][k-1].score_m * sahmm->mat[k-1].t[MATCH];
	score[FROM_DELETE] = mx[i-1][k-1].score_d * sahmm->del[k-1].t[MATCH];
	score[FROM_INSERT] = mx[i-1][k-1].score_i * sahmm->ins[k-1].t[MATCH];

	DNorm(score,3);
	if ((tr->statetype[j] = DChoose(score, 3)) == 3)
	  {
	    Warn("Oops: at M %d,%d, failed to choose from %f %f %f\n", 
		 i, k, score[0], score[1], score[2]);
	    return 0;
	  }

	tr->nodeidx[j]    = k-1;
	i--;			/* move to i-1, k-1 cell */
	k--;
	break;

      case DELETE:	/* we use DELETE at i,k , and therefore came from i,k-1 */
	if (k == 0) return 0;	/* trace failed! */

	score[FROM_MATCH]  = mx[i][k-1].score_m * sahmm->mat[k-1].t[DELETE];
	score[FROM_DELETE] = mx[i][k-1].score_d * sahmm->del[k-1].t[DELETE];
	score[FROM_INSERT] = mx[i][k-1].score_i * sahmm->ins[k-1].t[DELETE];

	DNorm(score,3);
	if ((tr->statetype[j]   = DChoose(score, 3)) == 3)
	  {
	    Warn("Oops: at D %d,%d, failed to choose from %f %f %f\n", 
		 i, k, score[0], score[1], score[2]);
	    return 0;
	  }

	tr->nodeidx[j]  = k-1;
	k--;			/* move to i, k-1 cell */
	break;

      case INSERT:	/* we use INSERT at i,k , and therefore came from i-1,k */
	if (i == 0) return 0;	/* trace failed! */

	score[FROM_MATCH]  = mx[i-1][k].score_m * sahmm->mat[k].t[INSERT];
	score[FROM_DELETE] = mx[i-1][k].score_d * sahmm->del[k].t[INSERT];
	score[FROM_INSERT] = mx[i-1][k].score_i * sahmm->ins[k].t[INSERT];

	DNorm(score,3);
	if ((tr->statetype[j]   = DChoose(score, 3)) == 3)
	  {
	    Warn("Oops: at I %d,%d, failed to choose from %f %f %f\n", 
		 i, k, score[0], score[1], score[2]);
	    return 0;
	  }  

	tr->nodeidx[j]  = k;
	i--;
	break;

      default:
	Die("Error: no such state type %d in traceback, i=%d k=%d j=%d\n",
	    tr->statetype[j-1], i, k, j);
      }

      if (tr->statetype[j] == INSERT || tr->statetype[j] == MATCH)
	tr->rpos[j] = i-1;
      else
	tr->rpos[j] = -1;
      j++;
    }

  
  /* Now we know the length of the optimal state sequence,
   * and we can reverse the traceback to make it 0..N-1 the way it
   * should be.
   */
  N = j;			/* length of optimal state seq */
  ReverseTrace(tr, N);

  /**************************************************
   * Ready for return. Dump debugging information
   **************************************************/

#ifdef EXTREME_DEBUG
  PrintTrace(tr);
#endif /* EXTREME_DEBUG */

  *ret_tr = tr;
  return 1;
}




/* Function: DumpSaMatrix()
 * 
 * Purpose:  Debugging output; dump the contents of a matrix
 *           to stdout. Because it prints everything, it's
 *           only useful on small test matrices.
 *           
 * Returns:  (void)
 */
void
DumpSaMatrix(struct sa_s **mx,          /* the matrix       */
	     int           L,  
	     int           M,	        /* dimensions of mx */
	     double *scalefactors)
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



/* Function: CreateSahmm()
 * 
 * Purpose:  Make a copy of a normal hmm for use by the simulated annealing
 *           algorithms. 
 *       
 * Returns:  pointer to allocated sa_hmm structure, or NULL on failure.
 *           Caller is responsible for free'ing the space.
 */          
struct sa_hmm_s *
CreateSahmm(struct hmm_struc *hmm, float kT)
{
  struct sa_hmm_s *sahmm;       /* RETURN: simulated annealing version of hmm */
  int idx, k, ts;		/* counters */

  if ((sahmm = (struct sa_hmm_s *) malloc (sizeof(struct sa_hmm_s))) == NULL)
    return NULL;
  sahmm->M = hmm->M;

  sahmm->ins = (struct sa_state_s *) malloc (sizeof(struct sa_state_s) * (hmm->M + 2));
  sahmm->del = (struct sa_state_s *) malloc (sizeof(struct sa_state_s) * (hmm->M + 2));
  sahmm->mat = (struct sa_state_s *) malloc (sizeof(struct sa_state_s) * (hmm->M + 2));
  if (sahmm->ins == NULL || sahmm->del == NULL || sahmm->mat == NULL) return NULL;

  
  for (k = 0; k <= hmm->M; k++)
    {
				/* state transitions */
      for (ts = 0; ts < 3; ts++)
	{
	  sahmm->mat[k].t[ts] = pow( (double) hmm->mat[k].t[ts], (double) (1.0 / kT) );
	  sahmm->ins[k].t[ts] = pow( (double) hmm->ins[k].t[ts], (double) (1.0 / kT) );
	  if (k > 0)
	    sahmm->del[k].t[ts] = pow( (double) hmm->del[k].t[ts], (double) (1.0 / kT) );
	  else 
	    sahmm->del[k].t[ts] = 0.0;
	}

				/* symbol emission probabilities */
      for (idx = 0; idx < Alphabet_size; idx++)
	{
	  if (k > 0)
	    sahmm->mat[k].p[idx]   = pow((double) hmm->mat[k].p[idx],(double) (1.0 / kT));
	  else
	    sahmm->mat[k].p[idx] = 0.0;
	  sahmm->ins[k].p[idx]   = pow((double) hmm->ins[k].p[idx],(double) (1.0 / kT));
	}

    }

  return sahmm;
}


/* Function: DestroySahmm()
 * 
 * Free the space from a simulated annealing copy of an hmm.
 * 
 * No return value (void).
 */
void
DestroySahmm(struct sa_hmm_s *sahmm)
{
  if (sahmm == NULL) return;
  if (sahmm->mat != NULL) free(sahmm->mat);
  if (sahmm->ins != NULL) free(sahmm->ins);
  if (sahmm->del != NULL) free(sahmm->del);
  free(sahmm);
}





/* Function: SaSymscore()
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
 *  
 *  
 */
double
SaSymscore(char x, double *scores, int hyperbayes)
{
  double result;
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

  if (! hyperbayes) result /= (double) count;
  return result;
}


