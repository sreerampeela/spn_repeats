/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* pam-prior.c
 * Tue Jun 21 14:43:38 1994
 * 
 * Evolving HMM's with substitution matrices; experimental.
 * 
 * Substitution matrices are assumed to be in log probability
 * form, log P(child | root), 20x20, indexed smx[root][child].
 *
 * Additional documentation: See NOTES/pam-prior.tex
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"
#include "externs.h"

static void   normalize_loglikely_vector(double *vec, int num);
static double sum_loglikely_vector(double *vec, int num);


/* Function: EvolveHMM()
 * 
 * Purpose:  Evolve an HMM, either winding it backwards in 
 *           evolutionary time or forwards. Substitution matrix
 *           (or matrices) are used as the evolutionary model.
 *
 *           Only match state emission vectors are affected, for now.
 *           
 * Algorithm: 1. assign each match state emission vector to a maximum
 *               likelihood PAM distance (and matrix, if a mixture
 *               model)
 *               
 *            2. Predict root distribution for each match emission vector
 *            
 *            3. Evolve forwards again by scale * PAM
 *            
 * Args:     hmm   - hmm to evolve
 *           bmx   - array of different base substitution matrices,
 *                   probability form
 *           num   - number of different smx's
 *           scale - evolution scale factor; 1.0 comes close to
 *                   reporoducing original HMM distribution
 *           ret_hmm - RETURN: new evolved HMM        
 *           
 * Return:   (void)
 *           ret_hmm is alloc'ed here and must be free'd by caller.
 */
void
EvolveHMM(struct hmm_struc *hmm, double (*bmx)[20][20], int num, 
	  double scale, struct hmm_struc **ret_hmm)
{
  struct hmm_struc *newhmm;
  int          pamdist;    	/* best PAM distance */
  int          whatmx;          /* idx of best matrix here */
  int          i,k;
  int          dist;            /* tmp for PAM distance */
  double       ll;		/* log likelihood score of a matrix */
  double       best_ll;         /* best log likelihood so far  */     
  double       rd[20];          /* p distribution at root           */
  double       smx[20][20];     /* extrapolated substitution matrix */
  double       p[20];		/* tmp for probability dist's */

  PrintSubstitutionMatrix(bmx[0]);


  newhmm = HMMDup(hmm);
  for (k = 1; k <= hmm->M; k++)
    {
      printf("Working on position %d, boss\n", k);

				/* construct count vector (copy) */
      for (i = 0; i < 20; i++)
	p[i] = (double) hmm->mat[k].p[i] * 100.0;

      /* Step 1. Assign a maximum likelihood choice of matrix and
       *         PAM distance.
       */
      LogifySubstitutionMatrix(bmx[0]);      
      BestMatrixDistance(bmx[0], aafq, p, &pamdist, &best_ll);
      ProbifySubstitutionMatrix(bmx[0]);
      whatmx = 0;

      for (i = 1; i < num; i++)
	{
	  LogifySubstitutionMatrix(bmx[i]);
	  BestMatrixDistance(bmx[i], aafq, p, &dist, &ll);
	  ProbifySubstitutionMatrix(bmx[i]);
	  if (ll > best_ll)
	    {
	      pamdist = dist;
	      whatmx  = i;
	    }
	}
      printf("Assigned matrix %d, PAM distance %d\n", whatmx, pamdist);

      /* Step 2. Predict the root distribution.
       */
      ExpSubstitutionMatrix(bmx[whatmx], pamdist, smx);
      LogifySubstitutionMatrix(smx);
      RootSequenceDistribution(smx, aafq, p, rd);

      /* Step 3. Evolve the root distribution forwards. Store
       *         in new HMM.
       */
      ExpSubstitutionMatrix(bmx[whatmx], (int) (pamdist * scale), smx);
      EvolveCounts(rd, smx, p);

				/* copy p back into new hmm */
      for (i = 0; i < 20; i++)
	newhmm->mat[k].p[i] = (float) p[i];

      
      printf("AA:     ");
      for (i = 0; i < 20; i++)
	printf("  %c   ", Alphabet[i]);

      printf("\nInitial: ");
      for (i = 0; i < 20; i++)
	printf("%5.3f ", hmm->mat[k].p[i]);
      
      printf("\nRoot:    ");
      for (i = 0; i < 20; i++)
	printf("%5.3f ", rd[i]);

      printf("\nFinal:   ");
      for (i = 0; i < 20; i++)
	printf("%5.3f ", p[i]);
      
      printf("\nDone.\n\n");
    }

  *ret_hmm = newhmm;
}  


/* Function: ParseSubstitutionMatrix()
 * 
 * Purpose:  Given a pointer to an open file containing a substitution matrix,
 *           parse the file and allocate and fill a 2D array of
 *           doubles containing the matrix. The file is
 *           assumed to be in root=rows, child=columns 20x20 format.
 *           
 * Returns:  1 on success, 0 if parsing fails.
 *           ret_smx is alloc'ed here and must be free'd by caller
 */
int
ParseSubstitutionMatrix(FILE *fp, double smx[20][20])
{
  char     buffer[512];		/* input buffer from fp                  */
  char    *sptr;
  int      row, col;
				/* parse */
  for (row = 0; row < 20; row++)
    {
				/* skip comments and blank lines */
      do {
	if (fgets(buffer, 512, fp) == NULL) return 0; 
      } while ((sptr = strtok(buffer, " \t\n")) == NULL || *sptr == '#');

				/* parse a line */
      for (col = 0; col < 20; col++)
	{
	  if (sptr == NULL) return 0; 
	  smx [row][col] = atof(sptr);
	  sptr = strtok(NULL, " \t\n");
	}
    }
  return 1;
}


/* Function: LogifySubstitutionMatrix()
 * 
 * Purpose:  Sometimes we read conditional probabilities directly
 *           (for example, the Overington matrices); convert
 *           to log P's.
 *           
 * Args:     smx - substitution matrix, giving conditional 
 *                 probabilities, P(child | root). 20x20;
 *                 indexed as smx[root][child]        
 *                 
 * Return:   (void)
 *           smx is converted to log probabilities                
 */
void
LogifySubstitutionMatrix(double smx[20][20])
{
  int child, root;

  for (root = 0; root < 20; root++)
    for (child = 0; child < 20; child++)
      smx[root][child] = log(smx[root][child]);
}


/* Function: ProbifySubstitutionMatrix()
 * 
 * Purpose:  Convert log P's back to P's.
 *           
 * Args:     smx - substitution matrix, giving log conditional 
 *                 probabilities, log P(child | root). 20x20;
 *                 indexed as smx[root][child]        
 *                 
 * Return:   (void)
 *           smx is converted to probabilities                
 */
void
ProbifySubstitutionMatrix(double smx[20][20])
{
  int child, root;

  for (root = 0; root < 20; root++)
    for (child = 0; child < 20; child++)
      smx[root][child] = exp(smx[root][child]);
}


/* Function: PAM2SubstitutionMatrix()
 * 
 * Purpose:  Convert a PAM scoring matrix (scaled integer log odds)
 *           to a substitution matrix.
 *           
 * Args:     pam   -  scaled integer log odds PAM matrix 27x27, symmetric
 *           scale -  multiplier to convert pam to ln() form
 *                    most PAM matrices are in 1/2 bits now, so
 *                    scale would be 1/2 * ln(2) = 0.34657359
 *           pri   -  prior probabilities of amino acids (initial estimate)
 *           smx   -  substitution matrix, already allocated
 *                    20x20, indexed smx[root][child]
 *                    
 * Return:   (void)
 *           smx is calculated and filled in.
 */
void
PAM2SubstitutionMatrix(int **pam, double scale, double *pri, double smx[20][20])
{
  int c, r;			/* index of child, root */
  double norm;
  
  for (r = 0; r < 20; r++)
    for (c = 0; c < 20; c++)
      smx[r][c] = pri[c] * exp((double) pam[aa_index[r]][aa_index[c]] * scale);

				/* normalize to remove rough edges */
				/* sum over all c for a given r == 1.0 */
  for (r = 0; r < 20; r++)
    {
      for (norm = 0.0, c = 0; c < 20; c++)
	norm += smx[r][c];
      for (c = 0; c < 20; c++)
	smx[r][c] /= norm;
    }
}


/* Function: RootSequenceDistribution()
 * 
 * Purpose:  Given a substitution matrix smx and a count vector
 *           ct, find a maximum likelihood distribution for
 *           the sequence character at the root node.
 *           
 *           Assumes a star topology for the tree. If your
 *           data are related by a tree (say, if you're a
 *           biologist, for instance <- sarcasm intended)
 *           this isn't the function for you and you should
 *           be doing something a la Joe Felsenstein. 
 *           The cheap solution is to "correct" the count
 *           data using a weighting rule, which is what
 *           I'll be doing eventually to keep the computations
 *           simple.
 *           
 * Args:     smx - substitution matrix, giving log conditional 
 *                 probabilities, log P(child | root). 20x20;
 *                 indexed as smx[root][child]
 *                 
 *           pri - prior probabilities of 20 amino acids, 0..19     
 *
 *           ct  - count vector; observed counts in one aligned
 *                 column. 0..19, indexed in alphabetical
 *                 amino acid order.
 *
 *           rd  - RETURN: root distribution for 20 possible
 *                 symbols at this column. 0..19, alphabetical 
 *                 order. Pre-allocated (rd[20]).
 *                 
 * Returns   (void)
 */                
void
RootSequenceDistribution(double smx[20][20], double *pri, double *ct, double *rd)
{
  int    child;                 /* index of child amino acid  */
  int    root;			/* index of root amino acid   */

  for (root = 0; root < 20; root++)
    {
      rd[root] = log(pri[root]);
      for (child = 0; child < 20; child++)
	rd[root] += ct[child] * smx[root][child];
    }
				/* normalize */
  normalize_loglikely_vector(rd, 20);
}



/* Function: MatrixLikelihood()
 * 
 * Purpose:  Given a substitution matrix and an observed count
 *           vector and a calculated root probability distribution,
 *           calculate the log likelihood of the matrix.
 *           
 * Args:     smx  - substitution matrix, giving log conditional 
 *                  probabilities, log P(child | root). 20x20;
 *                  indexed as smx[root][child]
 *                 
 *           pri  - prior probabilities of 20 amino acids
 *                  at root, 0..19     
 *
 *           ct   - count vector; observed counts in one aligned
 *                  column. 0..19, indexed in alphabetical
 *                  amino acid order.
 *
 *           ll   - RETURN: log likelihood of matrix
 *
 * Return:   (void)
 *           log likelihood returned thru passed ptr ll
 */
void
MatrixLikelihood(double smx[20][20], double *pri, double *ct, double *ll)
{
  double term[20];
  int    root, child;
  
  for (root = 0; root < 20; root++)
    {
      term[root] = log(pri[root]);
      for (child = 0; child < 20; child++)
	term[root] += ct[child] * smx[root][child];
    }      
  *ll = sum_loglikely_vector(term, 20);
}



/* Function: BestMatrixDistance()
 * 
 * Purpose:  Given a base substitution matrix (say, PAM2),
 *           find the maximum likelihood extrapolation factor
 *           ("distance") to account for some given count vector.
 *           Return the factor and the log likelihood.
 *           
 *           A much faster method would be to do a bi- or golden-section
 *           style search; bracket the maximum likelihood point
 *           and narrow in on it.
 *           
 * Args:     bmx      - base subsitution matrix log P(c|r), 20x20, smx[r][c]
 *           pri      - prior aa frequencies, 0..19
 *           ct       - observed count vector 0..19
 *           ret_dist - RETURN: ML distance factor
 *           ret_ll   - RETURN: best log likelihood       
 *                      
 * Return:   (void)
 *           ret_dist and ret_ll are passed ptrs.                     
 */
void
BestMatrixDistance(double bmx[20][20], double *pri, double *ct, 
		   int *ret_dist, double *ret_ll)
{
  double smx[20][20];	/* current substitution matrix     */
  double nmx[20][20];	/* new (exponeniated) subst matrix */
  int c,r,i;
  int    dist;
  double ll;
  double best_ll;
  
  /* Start with the base matrix
   */
  for (r = 0; r < 20; r++)
    for (c = 0; c < 20; c++)
      smx[r][c] = bmx[r][c];

  best_ll = -9999999.0;
  dist = 0;
  while (1)
    {
      dist++;
				/* likelihood of current matrix */
      MatrixLikelihood(smx, pri, ct, &ll);
      if (ll > best_ll) best_ll = ll;
      else
	{
	  *ret_ll   = best_ll;
	  *ret_dist = dist-1;
	  return;
	}

				/* extrapolate to next matrix */
      for (r = 0; r < 20; r++)
	for (c = 0; c < 20; c++)
	  {
	    nmx[r][c] = 0.0;
	    for (i = 0; i < 20; i++)
	      nmx[r][c] += exp(smx[r][i] + bmx[i][c]);
	    nmx[r][c] = log(nmx[r][c]);
	  }
      for (r = 0; r < 20; r++)
	for (c = 0; c < 20; c++)
	  smx[r][c] = nmx[r][c];
    }
}


/* Function: EvolveCounts()
 * 
 * Purpose:  Given a count vector, simulate its evolution according
 *           to some substitution matrix (probability form)
 *           
 * Args:     ct     - original count vector
 *           smx    - substitution matrix (probabilities)               
 *           nct    - RETURN: evolved count vector (alloced for 20)
 *           
 * Return:   (void)
 *           nct is calculated and returned          
 */
void
EvolveCounts(double *ct, double smx[20][20], double *nct)
{
  int r, c;

  for (c = 0; c < 20; c++)
    nct[c] = 0.0;
  for (r = 0; r < 20; r++)
    for (c = 0; c < 20; c++)
      nct[c] += smx[r][c] * ct[r];
}

/* Function: ExpSubstitutionMatrix()
 * 
 * Purpose:  Extrapolate a substitution matrix to some higher 
 *           power (distance)
 *           
 * Args:     bmx  - base matrix (PAM1, say) (prob form)
 *           dist - power to raise it to (integer)
 *           smx  - RETURN: exponentiated matrix (prob form)
 *           
 * Return:   (void)
 *           smx is calculated.
 */
void
ExpSubstitutionMatrix(double bmx[20][20], int dist, double smx[20][20])  
{
  double xmx[11][20][20];	/* bmx^0, bmx^2, bmx^4, ..., bmx^1024*/
  double tmx[20][20];		/* tmp matrix */
  int r,c,i;
  int exponent;			/* highest power of 2 < dist */
  int chunk;
  int first;

  /* Calculate powers of two of bmx, up to PAM1024
   */
  for (r = 0; r < 20; r++)
    for (c = 0; c < 20; c++)
      xmx[0][r][c] = bmx[r][c];
  for (exponent = 1; exponent <11; exponent++)
    {
      for (r = 0; r < 20; r++)
	for (c = 0; c < 20; c++)
	  {
	    xmx[exponent][r][c] = 0.0;
	    for (i = 0; i < 20; i++)
	      xmx[exponent][r][c] += xmx[exponent-1][r][i] * xmx[exponent-1][i][c];
	  }
    }

  /* Calculate smx using that family of exponentiated matrices
   */
  first = 1;
  for (exponent = 10, chunk = 1024; chunk >= 1; chunk /= 2, exponent--)
    {
      while (dist >= chunk)
	{
	  dist -= chunk;
	  if (first)
	    {			/* if first, copy into smx */
	      for (r = 0; r < 20; r++)
		for (c = 0; c < 20; c++)
		  smx[r][c] = xmx[exponent][r][c];
	      first = 0;
	    }	      
	  else
	    {			/* else multiply previous smx */
	      for (r = 0; r < 20; r++)
		for (c = 0; c < 20; c++)
		  {
		    tmx[r][c] = 0.0;
		    for (i = 0; i < 20; i++)
		      tmx[r][c] += smx[r][i] * xmx[exponent][i][c];
		  }
	      for (r = 0; r < 20; r++)
		for (c = 0; c < 20; c++)
		  smx[r][c] = tmx[r][c];
	    }
	}
    }
}


/* Function: PrintSubstitutionMatrix()
 * 
 * Purpose:  For debugging: output a matrix to stdout.
 */
void
PrintSubstitutionMatrix(double smx[20][20])
{
  int r,c;

  for (c = 0; c < 20; c++)
    printf("  %c   ", Alphabet[c]);
  puts("");
  for (r = 0; r < 20; r++)
    {
      for (c = 0; c < 20; c++)
	printf("%3.3f ", smx[r][c]);
      puts("");
    }
  puts("");
}



/* Function: sum_loglikely_vector()
 * 
 * Purpose:  Take a vector of log likelihoods, sum them,
 *           and return the summed log likelihood.
 */
static double
sum_loglikely_vector(double *vec, int num)
{
  double scalefactor;           /* for worrying about numerical precision */
  int    idx;			/* index into vec    */
  double sum;

  /* Find the biggest log P term in the vector for scaling
   */
  scalefactor = vec[0];
  for (idx = 1; idx < num; idx++)
    if (vec[idx] > scalefactor) scalefactor = vec[idx];
  
  for (idx = 0; idx < num; idx++) 
    {
      vec[idx] -= scalefactor;	/* rescale all terms */
      vec[idx] = exp(vec[idx]);	/* convert to probs */
    }
  sum = 0.0;
  for (idx = 0; idx < num; idx++) 
    sum += vec[idx];		/* do the sum */
  
  sum = log(sum);		/* back to log */
  sum += scalefactor;		/* rescale */
  return sum;
}


/* Function: normalize_loglikely_vector()
 * 
 * Purpose:  Take a vector of log likelihoods and
 *           convert it to a normalized probability vector.
 *           
 * Args:     vec - log likelihoods (usually negative)
 *           num - dimensionality of vec
 *
 * Return:   (void) 
 *           vec is converted.
 */
static void
normalize_loglikely_vector(double *vec, int num)
{
  double scalefactor;           /* for worrying about numerical precision */
  int    idx;			/* index into vec    */
  double norm;                  /* normalization sum */

  /* Find the biggest log P term in the vector;
   * we'll use it to scale.
   */
  scalefactor = vec[0];
  for (idx = 1; idx < num; idx++)
    if (vec[idx] > scalefactor) scalefactor = vec[idx];
  
  /* Rescale by dividing everything by a constant so 
   * the biggest number becomes 1.0. 
   */
  for (idx = 0; idx < num; idx++)
    {
      vec[idx] -= scalefactor;
      vec[idx] = exp(vec[idx]);
    }

  /* Normalize
   */
  norm = 0.0;
  for (idx = 0; idx < num; idx++)
    norm += vec[idx];
  for (idx = 0; idx < num; idx++)
    vec[idx] /= norm;
}



