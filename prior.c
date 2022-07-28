/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* prior.c
 * SRE, Thu Jul 15 11:50:58 1993
 * 
 * Configuration of default priors; reading customized priors
 * from disk; allocation and free'ing of prior structure;
 * application of priors as pseudocounts to counts-based HMMs.
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

#define TYPE_INT    0
#define TYPE_STRING 1
#define TYPE_REAL   2

static float ln_pnbar(float *cvec, int n, float *alpha);
static void  ln_norm(float *vec, int n);
static char *getword(FILE *fp, int type);


/* Function: DefaultSimplePrior()
 * 
 * Purpose:  Create a simple default prior. 
 *           Global alphabet info is assumed to have been set already.
 *
 * Args:     alphabet_size - how many symbols are in alphabet (4 or 20)
 *           ret_prior     - RETURN: constructed prior
 *
 * Return:   1 on success, 0 on failure.
 */ 
void
DefaultSimplePrior(struct prior_s  **ret_prior)
{
  struct prior_s *pri;
  int             x;
  
  if ((pri = AllocPrior()) == NULL)
    Die("Memory allocation failed for prior");

  pri->strategy = PRI_SIMPLE;
  pri->tnum     = 1;
  pri->tq[0]    = 1.0;
  pri->tm[0][MATCH]  = 0.970 * 3.;  /* match  --> match  */
  pri->tm[0][DELETE] = 0.010 * 3.;  /* match  --> delete */
  pri->tm[0][INSERT] = 0.020 * 3.;  /* match  --> insert */
  pri->td[0][MATCH]  = 0.930 * 3.;  /* delete --> match  */
  pri->td[0][DELETE] = 0.050 * 3.;  /* delete --> delete */
  pri->td[0][INSERT] = 0.020 * 3.;  /* delete --> insert */
  pri->ti[0][MATCH]  = 0.590 * 3.;  /* insert --> match  */
  pri->ti[0][DELETE] = 0.020 * 3.;  /* insert --> delete */
  pri->ti[0][INSERT] = 0.390 * 3.;  /* insert --> insert */
  
  pri->mnum = pri->inum = 1;
  pri->mq[0] = pri->iq[0] = 1.0;
  if (Alphabet_type == kAmino)	/* use real database fi's for protein */
    for (x = 0; x < Alphabet_size; x++)
      {
	pri->mat[0][x] = aafq[x] * 20.;
	pri->ins[0][x] = aafq[x] * 100000.; /* hardwired ins emissions */
      }
  else
    for (x = 0; x < Alphabet_size; x++)
      {
	pri->mat[0][x] = 1.0;
	pri->ins[0][x] = 100000.0; /* hardwired ins emissions */
      }
  *ret_prior = pri;
}

  
/* Function: ReadPrior()
 * 
 * Purpose:  Read a prior from a file.
 */
void
ReadPrior(char *pfile, struct prior_s **ret_prior)
{
  struct prior_s *prior;
  FILE *fp;
  int   x, q;
  char *sptr;

  if ((fp = fopen(pfile, "r")) == NULL)
    Die("Failed to open mixture prior file %s", pfile);
  if ((prior = AllocPrior()) == NULL)
    Die("allocation for prior failed");

  /* First entry is the strategy: 
   * SIMPLE, PAM, MIX, or STRUCT
   */
  sptr = getword(fp, TYPE_STRING);
  s2upper(sptr);
  if      (strcmp(sptr, "SIMPLE") == 0) prior->strategy = PRI_SIMPLE;
  else if (strcmp(sptr, "PAM")    == 0) prior->strategy = PRI_PAM;
  else if (strcmp(sptr, "MIX")    == 0) prior->strategy = PRI_MIX;
  else if (strcmp(sptr, "STRUCT") == 0) prior->strategy = PRI_STRUCT;
  else Die("No such prior strategy %s; failed to parse file %s", sptr, pfile);

  /* Second entry is the alphabet type:
   * Amino or Nucleic
   */
  sptr = getword(fp, TYPE_STRING);
  s2upper(sptr);
  if (strcmp(sptr, "AMINO") == 0 && Alphabet_type != kAmino)
    Die("Prior is for protein, but HMM and/or sequences are nucleic acid");
  if (strcmp(sptr, "NUCLEIC") == 0 && Alphabet_type != kDNA && Alphabet_type != kRNA)
    Die("Prior is for DNA/RNA, but HMM and/or sequences are protein");

  /* State transition priors:
   * # of mixtures.
   * then for each mixture:
   *    prior (if not structural prior)
   *    weight matrix, NINPUTS terms (if structural prior)
   *    Dirichlet terms for Tmm, Tmi, Tmd, Tim, Tii, Tid, Tdm, Tdi, Tdd
   */
  prior->tnum = atoi(getword(fp, TYPE_INT));
  if (prior->tnum < 0)
    Die("%d is bad; need at least one state transition mixture component", prior->tnum);
  if (prior->tnum > MAXDCHLET)
    Die("%d is bad, too many transition components (MAXDCHLET = %d)\n", MAXDCHLET);
  if ((prior->strategy == PRI_SIMPLE || prior->strategy == PRI_PAM) && prior->tnum != 1)
    Die("%d is bad, simple or PAM priors must have only one transition component", prior->tnum);
  for (q = 0; q < prior->tnum; q++)
    {
      if (prior->strategy == PRI_STRUCT)
	for (x = 0; x < NINPUTS; x++) prior->tw[q][x] = (float) atof(getword(fp, TYPE_REAL));
      else                            prior->tq[q]    = (float) atof(getword(fp, TYPE_REAL));
      for (x = 0; x < 3; x++)         prior->tm[q][x] = (float) atof(getword(fp, TYPE_REAL));
      for (x = 0; x < 3; x++)         prior->ti[q][x] = (float) atof(getword(fp, TYPE_REAL));
      for (x = 0; x < 3; x++)         prior->td[q][x] = (float) atof(getword(fp, TYPE_REAL));
    }

  /* Match emission priors:
   * # of mixtures.
   * then for each mixture:
   *    prior (if not structural prior)
   *    weight matrix, NINPUTS terms (if structural prior)
   *    Dirichlet terms for Alphabet_size symbols in Alphabet
   */
  prior->mnum = atoi(getword(fp, TYPE_INT));
  if (prior->mnum < 0)
    Die("%d is bad; need at least one match emission mixture component", prior->mnum);
  if (prior->mnum > MAXDCHLET)
    Die("%d is bad; too many match components (MAXDCHLET = %d)\n", prior->mnum, MAXDCHLET);
  if (prior->strategy == PRI_SIMPLE && prior->mnum != 1)
    Die("%d is bad; simple prior must have one mat emission component.", prior->mnum);
  if (prior->strategy == PRI_PAM && prior->mnum != Alphabet_size)
    Die("%d is bad; PAM prior needs %d match emission components", prior->mnum, Alphabet_size);

  for (q = 0; q < prior->mnum; q++)
    {
      if (prior->strategy == PRI_STRUCT)
	for (x = 0; x < NINPUTS; x++)     prior->mw[q][x]  = (float) atof(getword(fp, TYPE_REAL));
      else                                prior->mq[q]     = (float) atof(getword(fp, TYPE_REAL));
      for (x = 0; x < Alphabet_size; x++) prior->mat[q][x] = (float) atof(getword(fp, TYPE_REAL));
    }

  /* Insert emission priors:
   * # of mixtures.
   * then for each mixture component:
   *    prior (if not structural prior)
   *    weight matrix, NINPUTS terms (if structural prior)
   *    Dirichlet terms for Alphabet_size symbols in Alphabet
   */
  prior->inum = atoi(getword(fp, TYPE_INT));
  if (prior->inum < 0)
    Die("%d is bad; need at least one insert emission mixture component", prior->inum);
  if (prior->inum > MAXDCHLET)
    Die("%d is bad; too many insert components (MAXDCHLET = %d)\n", prior->inum, MAXDCHLET);
  if ((prior->strategy == PRI_SIMPLE || prior->strategy == PRI_PAM) && prior->inum != 1)
    Die("%d is bad; simple or PAM prior must have one ins emission component.", prior->inum);
  for (q = 0; q < prior->inum; q++)
    {
      if (prior->strategy == PRI_STRUCT)
	for (x = 0; x < NINPUTS; x++)     prior->iw[q][x]  = (float) atof(getword(fp, TYPE_REAL));
      else                                prior->iq[q]     = (float) atof(getword(fp, TYPE_REAL));
      for (x = 0; x < Alphabet_size; x++) prior->ins[q][x] = (float) atof(getword(fp, TYPE_REAL));
    }
  *ret_prior = prior;
}
  

/* Function: AllocPrior()
 * 
 * Purpose:  Allocate a struct prior_s, structure for holding information
 *           about prior expectations on probabilities.
 * 
 * Args:     none
 *                  
 * Return:   Pointer to the new prior structure, or NULL on failure.
 */
struct prior_s *
AllocPrior(void)
{
  struct prior_s *pri;

  if ((pri    = (struct prior_s *) malloc (sizeof(struct prior_s))) == NULL)
    return NULL;
  pri->strategy = -1;
  pri->tnum     = 0;
  pri->mnum     = 0;
  pri->inum     = 0;
  return pri;
}


/* Function: FreePrior()
 * 
 * Purpose:  Free space allocated to a prior structure.
 * 
 * Args:     prior - the prior structure to be free'd
 *                   
 * Return:   (void)
 */                  
void
FreePrior(struct prior_s *prior)
{
  free(prior);
}


/* Function: ToPAMPrior()
 * 
 * Purpose:  Convert a simple prior to an ad hoc PAM-based mixture
 *           prior. Insert emission and state transition components
 *           remain unchanged. Match emission prior becomes a 20
 *           component mixture, one per amino acid.
 *           
 * Args:     pam     - 27x27 PAM matrix
 *           scale   - scale on PAM file (units relative to ln) 
 *           wt      - weight on prior (akin to alpha; 20 seems good ?)
 *           prior   - alloc'ed prior to initialize
 *           
 * Return:   (void)
 *           Existing simple prior is revised.
 *           If anything goes wrong, we die here.
 */          
void
ToPAMPrior(int **pam, float scale, float wt, struct prior_s *pri)
{
  float sum;
  int   xi, xj;
  int   idx1, idx2;

  pri->strategy = PRI_PAM;
  pri->mnum     = 20;
  
  /* Convert PAM entries back to conditional prob's P(xj | xi)
   */
  for (xi = 0; xi < Alphabet_size; xi++)
    for (xj = 0; xj < Alphabet_size; xj++)
      {
	idx1 = Alphabet[xi] - 'A';
	idx2 = Alphabet[xj] - 'A';
	pri->mat[xi][xj] = aafq[xj] * exp((float) pam[idx1][idx2] * scale);
      }

  /* Normalize so that rows add up to wt.
   * i.e. Sum(xj) mat[xi][xj] = 1.0 for every row xi
   */
  for (xi = 0; xi < Alphabet_size; xi++)
    {
      pri->mq[xi] = 1. / Alphabet_size;
      sum = 0.0;
      for (xj = 0; xj < Alphabet_size; xj++)
	sum += pri->mat[xi][xj];
      for (xj = 0; xj < Alphabet_size; xj++)
	pri->mat[xi][xj] *= wt / sum;
    }
}



void
PrintPAMPrior(struct prior_s *pri)
{
  int xi, xj;
  float sum;

  printf("    ");
  for (xj = 0; xj < Alphabet_size; xj++)
    printf("  %c    ", Alphabet[xj]);
  puts("");

  for (xi = 0; xi < Alphabet_size; xi++)
    {
      sum = 0.0;
      printf(" %c  ", Alphabet[xi]);
      for (xj = 0; xj < Alphabet_size; xj++)
	{
	  printf("%6.4f ", pri->mat[xi][xj]);
	  sum += pri->mat[xi][xj];
	}
      printf("sum %6.4f\n", sum);
    }

  printf("sum ");
  for (xj = 0; xj < Alphabet_size; xj++)
    {
      sum = 0.0;
      for (xi = 0; xi < Alphabet_size; xi++)
	sum += pri->mat[xi][xj];
      printf("%6.4f ", sum);
    }
  puts("");
}

/* Function: DefaultRandomModel()
 * 
 * Purpose:  Set up a default random sequence model, using
 *           global aafq[]'s for protein or 0.25 for nucleic 
 *           acid. randomseq is alloc'ed in caller. Alphabet information
 *           must already be known.
 */
void
DefaultRandomModel(float *randomseq)
{
  int x;
  if (Alphabet_type == kAmino)
    for (x = 0; x < Alphabet_size; x++)
      randomseq[x] = aafq[x];
  else if (Alphabet_type == kDNA || Alphabet_type == kRNA)
    for (x = 0; x < Alphabet_size; x++)
      randomseq[x] = 0.25;
  else
    Die("No support for non-protein, non-nucleic acid alphabets.");
}


/* Function: ReadRandomModel()
 * 
 * Purpose:  Get a random model from a file.
 *           First non-comment field is "Nucleic" or "Amino".
 *           Remainder of file consists of 4 or 20 floats, 
 *           representing random model emissions.
 *           Alphabet_size and Alphabet_type must already be set.
 */
void
ReadRandomModel(char *fname, float *randomseq)
{
  FILE *fp;
  char *sptr;
  int   x;

  if ((fp = fopen(fname, "r")) == NULL)
    Die("Failed to open random model parameter file %s\n", fname);
  
  /* First entry is Amino or Nucleic
   */
  if ((sptr = getword(fp, TYPE_STRING)) == NULL)
    Die("Random model file %s appears empty", fname);
  s2upper(sptr);
  if (strcmp(sptr, "AMINO") == 0 && Alphabet_type != kAmino)
    Die("Random model is for protein, but HMM and/or sequences are nucleic acid");
  if (strcmp(sptr, "NUCLEIC") == 0 && Alphabet_type != kDNA && Alphabet_type != kRNA)
    Die("Random model is for DNA/RNA, but HMM and/or sequences are protein");
  if (strcmp(sptr, "AMINO") != 0 && strcmp(sptr, "NUCLEIC") != 0)
    Die("Random model does not specify Nucleic or Amino as first field");
  
  /* Remainder are parameters
   */
  for (x = 0; x < Alphabet_size; x++)
    randomseq[x] = (float) atof(getword(fp, TYPE_REAL));
}


/* Function: PriorifyMatchVector()
 * 
 * Purpose:  Add prior pseudocounts to an observed match
 *           emission count vector.
 * 
 * Args:     vec     - the 4 or 20-long vector of counts to modify
 *           prior   - prior information
 *                   
 * Return:   (void)
 *           The counts in vec are changed.
 */                  
void
PriorifyMatchVector(float *vec, struct prior_s *prior)
{
  int   x;			/* counter over vec                     */
  int   q;			/* counter over mixtures                */
  float mix[MAXDCHLET];		/* posterior distribution over mixtures */

  if (prior->strategy == PRI_SIMPLE)
    mix[0] = 1.0;
  else if (prior->strategy == PRI_PAM)
    {
      for (q = 0; q < Alphabet_size; q++) 
	mix[q] = vec[q];
      FNorm(mix, Alphabet_size);
    }
  else if (prior->strategy == PRI_MIX || prior->strategy == PRI_STRUCT)
    {
      for (q = 0; q < prior->mnum; q++)
	{
	  mix[q] =  prior->mq[q] > 0.0 ? log(prior->mq[q]) : -999.;
	  mix[q] += ln_pnbar(vec, Alphabet_size, prior->mat[q]);
	}
      ln_norm(mix, prior->mnum);      
    }

  for (q = 0; q < prior->mnum; q++)
    for (x = 0; x < Alphabet_size; x++)
      vec[x] +=  mix[q] * prior->mat[q][x];
}


/* Function: PriorifyInsertVector()
 * 
 * Purpose:  Add prior pseudocounts to an observed insert
 *           emission count vector.
 * 
 * Args:     vec     - the 4 or 20-long vector of counts to modify
 *           prior   - prior information
 *                   
 * Return:   (void)
 *           The counts in vec are changed.
 */                  
void
PriorifyInsertVector(float *vec, struct prior_s *prior)
{
  int   x;			/* counter over vec                     */
  int   q;			/* counter over mixtures                */
  float mix[MAXDCHLET];		/* posterior distribution over mixtures */

  if (prior->strategy == PRI_SIMPLE || prior->strategy == PRI_PAM)
    mix[0] = 1.0;
  else if (prior->strategy == PRI_MIX || prior->strategy == PRI_STRUCT)
    {
      for (q = 0; q < prior->inum; q++)
	{
	  mix[q] =  prior->iq[q] > 0.0 ? log(prior->iq[q]) : -999.;
	  mix[q] += ln_pnbar(vec, Alphabet_size, prior->ins[q]);
	}
      ln_norm(mix, prior->inum);      
    }

  for (q = 0; q < prior->inum; q++)
    for (x = 0; x < Alphabet_size; x++)
      vec[x] +=  mix[q] * prior->ins[q][x];
}


 
/* Function: PriorifyTransitionVectors()
 * 
 * Purpose:  Add prior pseudocounts to three vectors of 
 *           three state transition observed counts.
 *          
 *           At k=0 in model, there are no valid transitions
 *           from delete. Pass NULL for td.
 *           
 * Args:     tm     - state transition vector from match (counts)
 *           ti     - state transition vector from insert (counts)
 *           td     - state transition vector from delete (counts)         
 *           prior  - prior information for regularization
 *           
 * Return:   void
 *           The tm components are changed.
 */          
void
PriorifyTransitionVectors(float *tm, float *ti, float *td, struct prior_s *prior)
{
  int   ts;
  int   q;
  float mix[MAXDCHLET];

  if (prior->strategy == PRI_SIMPLE || prior->strategy == PRI_PAM)
    mix[0] = 1.0;
  else if (prior->strategy == PRI_MIX || prior->strategy == PRI_STRUCT)
    {
      for (q = 0; q < prior->tnum; q++)
	{
	  mix[q] =  prior->tq[q] > 0.0 ? log(prior->tq[q]) : -999.;
	  mix[q] += ln_pnbar(tm, 3, prior->tm[q]);
	  mix[q] += ln_pnbar(ti, 3, prior->ti[q]);
	  if (td != NULL) mix[q] += ln_pnbar(td, 3, prior->td[q]);
	}
      ln_norm(mix, prior->tnum);
    }

  for (q = 0; q < prior->tnum; q++)
    for (ts = 0; ts < 3; ts++)  
      {
	tm[ts] += mix[q] * prior->tm[q][ts]; 
	ti[ts] += mix[q] * prior->ti[q][ts]; 
	if (td != NULL) td[ts] += mix[q] * prior->td[q][ts]; 
      }
}


/* Function: PriorifyHMM()
 * 
 * Purpose:  Add prior "pseudocounts" to a count-based HMM.
 * 
 * Args:     hmm     - the HMM to convert. Contains structural annotation.
 *           prior   - a prior structure, containing regularization info
 *                    
 * Return:   1 on success, 0 on failure.                   
 */
void
PriorifyHMM(struct hmm_struc *hmm, struct prior_s  *prior)
{
  int   k;			/* counter for states                        */
  
  /* k = 0 BEGIN is special case, no delete state
   */
  if (prior->strategy == PRI_STRUCT) 
    StructurePerceptron(prior, hmm->xray);
  PriorifyTransitionVectors(hmm->mat[0].t, hmm->ins[0].t, NULL, prior);
  PriorifyInsertVector(hmm->ins[0].p, prior);

  /* remainder of model is straightforward
   */
  for (k = 1; k <= hmm->M ; k++)
    {
      if (prior->strategy == PRI_STRUCT) 
	StructurePerceptron(prior, hmm->xray + k*NINPUTS);
      PriorifyTransitionVectors(hmm->mat[k].t, hmm->ins[k].t, hmm->del[k].t, prior);
      PriorifyMatchVector(hmm->mat[k].p, prior);
      PriorifyInsertVector(hmm->ins[k].p, prior);
    }

  /* We've written some crap into state hmm->M, which can only transit
   * to match+1 (the dummy END state) or ins+0. 
   */
  hmm->mat[hmm->M].t[DELETE] = 0.0;
  hmm->ins[hmm->M].t[DELETE] = 0.0;
  hmm->del[hmm->M].t[DELETE] = 0.0;
}

/* Function: StructurePerceptron()
 * 
 * Purpose:  Three-dimensional structure-based prior. Convert an
 *           input vector of structural annotation, xray, to
 *           prior probability distributions over possible
 *           mixtures. Write this distribution directly into the prior
 *           structure.
 *           
 * Args:     prior  - mixture prior data
 *           xray   - structural annotation, 0..NINPUTS-1.         
 */
void
StructurePerceptron(struct prior_s *prior, float *xray)
{
  int q;			/* counter over mixtures */
  int i;			/* counter over inputs   */

				/* Transition mixtures */
  for (q = 0; q < prior->tnum; q++)
    for (prior->tq[q] = 0.0, i = 0; i < NINPUTS; i++)
      prior->tq[q] += prior->tw[q][i] * xray[i];
  ln_norm(prior->tq, prior->tnum); /* softmax */

				/* Match emission mixtures */
  for (q = 0; q < prior->mnum; q++)
    for (prior->mq[q] = 0.0, i = 0; i < NINPUTS; i++)
      prior->mq[q] += prior->mw[q][i] * xray[i];
  ln_norm(prior->mq, prior->mnum); /* softmax */

				/* Insert emission mixtures */
  for (q = 0; q < prior->inum; q++)
    for (prior->iq[q] = 0.0, i = 0; i < NINPUTS; i++)
      prior->iq[q] += prior->iw[q][i] * xray[i];
  ln_norm(prior->iq, prior->inum); /* softmax */
}

/* Function: AnnotateAlignment()
 * 
 * Purpose:  Calculate structural inputs, given an alignment.
 * 
 * Return:   void
 *           ret_inputs is malloced here. free(inputs). 
 */
void
AnnotateAlignment(char **aseq, int nseq, AINFO *ainfo, float **ret_inputs)
{
  int   *rpos;                  /* array of raw sequence positions in aseq */
  int    idx;			/* counter over seqs                       */
  int    apos;			/* position in alignment                   */
  float *inputs;		/* [0..alen][0..NINPUTS-1]                 */
  float  efrac, hfrac, denom;	/* intermediates in calculations           */

  /* inputs is 1..alen, plus zero column for BEGIN;
   * however, alignment is 0..alen-1, so watch for
   * offset stuff
   */
  inputs = (float *) MallocOrDie ((ainfo->alen+1) * NINPUTS * sizeof(float));
  rpos   = (int *) MallocOrDie (nseq * sizeof(int));
  for (idx = 0; idx < nseq; idx++)
    rpos[idx] = 0;

  /* column zero in inputs is special, BEGIN
   */
  inputs[XRAY_bias] = 1.0;	/* always one      */
  inputs[XRAY_E]    = 0.0;	/* not in sheet    */
  inputs[XRAY_H]    = 0.0;	/* not in helix    */
  inputs[XRAY_SA]   = 15.0;	/* very accessible */

  for (apos = 0; apos < ainfo->alen; apos++)
    {
				/* bias is always 1 */
      inputs[(apos+1)*NINPUTS + XRAY_bias] = 1.0; 

				/* fraction in secondary structure */
      efrac = hfrac = denom = 0.0;
      for (idx = 0; idx < nseq; idx++)
	if (ainfo->sqinfo[idx].flags & SQINFO_SS)
	  {
	    if (!isgap(aseq[idx][apos]))
	      {
		if (ainfo->sqinfo[idx].ss[rpos[idx]] == 'E') efrac += 1.0;
		if (ainfo->sqinfo[idx].ss[rpos[idx]] == 'H') hfrac += 1.0;
	      }
	    denom += 1.0;
	  }
      inputs[(apos+1)*NINPUTS + XRAY_E] = (denom > 0.0) ? efrac / denom : 0.33;
      inputs[(apos+1)*NINPUTS + XRAY_H] = (denom > 0.0) ? hfrac / denom : 0.33;

				/* surface accessibility */
      efrac = denom = 0.0;
      for (idx = 0; idx < nseq; idx++)
	if (ainfo->sqinfo[idx].flags & SQINFO_SA)
	  {
	    if (!isgap(aseq[idx][apos]))
	      {
		denom += 1.0;
		switch(ainfo->sqinfo[idx].sa[rpos[idx]]) {
		case '0': efrac += 0.0; break;
		case '1': efrac += 1.0; break;
		case '2': efrac += 2.0; break;
		case '3': efrac += 3.0; break;
		case '4': efrac += 4.0; break;
		case '5': efrac += 5.0; break;
		case '6': efrac += 6.0; break;
		case '7': efrac += 7.0; break;
		case '8': efrac += 8.0; break;
		case '9': efrac += 9.0; break;
		case 'a': efrac +=10.0; break;
		case 'b': efrac +=11.0; break;
		case 'c': efrac +=12.0; break;
		case 'd': efrac +=13.0; break;
		case 'e': efrac +=14.0; break;
		case 'f': efrac +=15.0; break;
		default: Die("bogus surface accessibility %c\n", 
			     ainfo->sqinfo[idx].sa[rpos[idx]]);
		}
	      }
	  }
      inputs[(apos+1)*NINPUTS + XRAY_SA] = (denom > 0.0) ? efrac / denom : 7.5;

				/* increment rpos if necessary */
      for (idx = 0; idx < nseq; idx++)
	if (! isgap(aseq[idx][apos]))
	  rpos[idx]++;
    } /* end loop over all aligned columns. */

  free(rpos);
  *ret_inputs = inputs;
}



/* Function: ln_pnbar()
 * 
 * Calculates ln P(n|dirichlet), the log probability of a count vector
 * given a Dirichlet distribution. Adapted from an implementation
 * by Graeme Mitchison.
 */
static float
ln_pnbar(float *cvec, int n, float *alpha)
{
  float lnp;			/* log likelihood of P(cvec | Dirichlet) */
  float sum1, sum2, sum3;
  int   x;

  sum1 = sum2 = sum3 = lnp = 0.0;
  for (x = 0; x < n; x++)
    {
      sum1 += cvec[x] + alpha[x];
      sum2 += alpha[x]; 
      sum3 += cvec[x];
      lnp  += (float) Gammln((double)(alpha[x] + cvec[x]));
      lnp  -= (float) Gammln((double)(cvec[x] + 1.)); 
      lnp  -= (float) Gammln((double)alpha[x]);
    }
  lnp -= (float) Gammln((double) sum1);
  lnp += (float) Gammln((double) sum2);
  lnp += (float) Gammln((double) (sum3 + 1.));
  return lnp;
}
    
/* Function: ln_norm()
 * 
 * Normalize a vector of log likelihoods, changing it
 * to a probability vector. Be careful of overflowing exp().
 * Implementation adapted from Graeme Mitchison.
 */
static void
ln_norm(float *vec, int n)
{
  int   x;
  float max   = -1.0e30;
  float denom = 0.;

  for (x = 0; x < n; x++)
    if (vec[x] > max) max = vec[x];
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      denom += exp(vec[x] - max);
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      vec[x] = exp(vec[x] - max) / denom;
    else
      vec[x] = 0.0;
}
  
/* Function: getword()
 * 
 * Purpose:  little function used by ReadPrior() to parse
 *           next valid field out of an open file, ignoring
 *           comments. '#' marks the beginning of a comment.
 */
static char *
getword(FILE *fp, int type)
{
  static char buffer[512];
  static char *sptr = NULL;
  
  if (sptr != NULL) sptr = strtok(NULL, " \t\n");

  while (sptr == NULL)
    {
      if ((sptr = fgets(buffer, 512, fp)) == NULL) return NULL;
      if ((sptr = strchr(buffer, '#')) != NULL) *sptr = '\0';
      sptr = strtok(buffer, " \t\n");
    }

  switch (type) {
  case TYPE_STRING: 
    if (strlen(sptr) == 0) Die("Parse of prior failed: expected string, got nothing");
    break;
  case TYPE_INT:    
    if (!IsInt(sptr)) Die("Parse of prior failed: expected integer, got %s", sptr);
    break;
  case TYPE_REAL:
    if (!IsReal(sptr)) Die("Parse of prior failed: expected real value, got %s", sptr); 
    break;
  }

  return sptr;
}



