/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* misc.c
 * SRE, Thu Jul 15 18:49:19 1993
 * (from cove)
 * 
 * Functions that I don't know quite where to put yet.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "states.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: DetermineAlphabet()
 * 
 * Purpose:  From a set of sequences (raw or aligned), make a good
 *           guess whether they're RNA, DNA, protein, or something
 *           else, and set alphabet accordingly.
 */
int
DetermineAlphabet(char  **rseqs, int  nseq)
{
  int idx;
  int othercount, dnacount, rnacount, aminocount;
  
  othercount = dnacount = rnacount = aminocount = 0;
  for (idx = 0; idx < nseq; idx++)
    {
      switch (Seqtype(rseqs[idx])) {
      case kRNA:      rnacount++;   break;
      case kDNA:      dnacount++;   break;
      case kAmino:    aminocount++; break;
      case kOtherSeq: othercount++; break;
      default: Die("No such alphabet type");
      }
    }
  if      (dnacount == nseq)   Alphabet_type = kDNA;
  else if (rnacount == nseq)   Alphabet_type = kRNA;
  else if (aminocount == nseq) Alphabet_type = kAmino;
  else if (dnacount > rnacount && dnacount > aminocount && dnacount > othercount) {
    Warn("Looks like DNA sequence, hope that's right");
    Alphabet_type = kDNA;
  }
  else if (rnacount > dnacount && rnacount > aminocount && rnacount > othercount) {
    Warn("Looks like RNA sequence, hope that's right");
    Alphabet_type = kRNA;
  }
  else if (aminocount > dnacount && aminocount > rnacount && aminocount > othercount) {
    Warn("Looks like protein sequence, hope that's right");
    Alphabet_type = kAmino;
  }
  else Die("Sorry, I can't tell if that's protein or DNA"); 

  switch(Alphabet_type) {
  case kRNA:   strncpy(Alphabet, RNA_ALPHABET, 4);   Alphabet_size = 4;  break;
  case kDNA:   strncpy(Alphabet, DNA_ALPHABET, 4);   Alphabet_size = 4;  break;
  case kAmino: strncpy(Alphabet, AMINO_ALPHABET, 20); Alphabet_size = 20; break;
  default: Die("No support for non-DNA/RNA or protein alphabets");  
  }

  return 1;
}


/* Function: BlockRaggedEdgedAlignment()
 * 
 * Purpose:  A brutal hack for ignoring exterior gaps on an
 *           alignment in Maxmodelmaker(). Convert all
 *           exterior gaps to the symbol ',' and hope to
 *           God nobody ever uses commas to mean anything
 *           in an alignment. 
 *           
 * Args:     aseqs  - [0..nseq-1][0..alen-1] alignment to block
 *           nseq   - number of seqs in the alignment
 *           alen   - width of alignment, columns
 *           
 * Return:   (void). Data in aseqs is changed.
 */
void
BlockRaggedEdgedAlignment(char **aseqs, int nseq, int alen)
{
  int  idx, pos;

  for (idx = 0; idx < nseq; idx++)
    {
      for (pos = 0; pos < alen; pos++)
	{
	  if (isgap(aseqs[idx][pos])) aseqs[idx][pos] = ',';
	  else break;
	}
      for (pos = alen-1; pos >= 0; pos--)
	{
	  if (isgap(aseqs[idx][pos])) aseqs[idx][pos] = ',';
	  else break;
	}
    }
}


/* Functions: DChoose(), FChoose()
 * 
 * Purpose:   Make a random choice from a (possibly
 *            unnormalized) probability distribution.
 *            DChoose() is for double-precision vectors;
 *            FChoose() is for single-precision float vectors.
 *            Returns the number of the choice, or N if
 *            it has failed.
 */
int
DChoose(double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = sre_random();
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }
  return (int) (sre_random() * N);   /* bulletproof */
}
int
FChoose(float *p, int N)
{
  float roll;                   /* random fraction */
  float sum;                    /* integrated prob */
  int   i;                      /* counter over the probs */

  roll    = sre_random();
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      sum += p[i];
      if (roll < sum) return i;
    }
  return (int) (sre_random() * N);    /* bulletproof */
}


/* Function: AlignmentTooBig()
 * 
 * Purpose:  Return TRUE if a full dynamic programming alignment
 *           of a sequence of length L against a model of length M
 *           is likely to exceed the available memory of the machine.
 *           Used to determine switching into the slower, linear-memory
 *           cost alignment algorithms.
 */
int
AlignmentTooBig(int L, int M)
{
  float ram;

  /* DP cells hold three ints; hence 4 bytes * 3 ints = 12 bytes/cell
   */
  ram = 12.0 * (float) (L+2) * (float) (M+2) / 1.0e6;
  if (ram > (float) RAMLIMIT)
    return TRUE;
  else
    return FALSE;
}
