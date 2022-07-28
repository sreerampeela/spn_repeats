/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* profiles.c
 * SRE, Mon Nov 28 16:08:22 1994
 * 
 * Code for mimicking the GCG PROFILE package, both in HMM output
 * as .prf files and HMM input from .prf files.
 * How well this will work, given that HMMs and profiles are
 * substantially different in detail, remains to be seen.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifdef DEBUG
#include <assert.h>
#endif 

#include "states.h"
#include "externs.h"
#include "squid.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: WriteProfile()
 * 
 * Purpose:  Given an HMM, write a profile .prf file as
 *           output. Based on reverse engineering GCG 8.0.
 *           
 *           The profile for symbol i at model position k is:
 *
 *           PROF(i,k) = log P(i,k) + log T[m(k-1)->m(k)]
 *           m_open(k) = -log T[m(k)->i(k)] - log T[i(k)->m(k+1)] 
 *                       + log T[m(k)->m(k+1)]
 *           m_extend(k) =  -log T[i(k)->i(k)]
 *           
 *           This is the simplest possible mapping of HMM to
 *           profile. It ignores the delete transitions entirely,
 *           fitting only the match and inserts. Note the somewhat
 *           tricky method by which the M->M transition score is
 *           smuggled into the profile. It lets the profile do 
 *           its worst on deletes.
 *
 *           A better approach would be to train the HMM optimally
 *           under profile-ish constraints.
 *           
 * Args:     fp        - open file to write to (or stdout, possibly)
 *           hmm       - hmm to write, probability form
 *           shmm      - hmm to write, integer log-odds form
 *           randomseq - random sequence model (symbol frequencies)
 *           
 * Return:   1 on success, 0 on failure.
 */
int
WriteProfile(FILE *fp, struct hmm_struc *hmm, struct shmm_s *shmm, float *randomseq)
{
  int x, k;
  int score, max_score;
  int sym;

  /* Write the header. Apparently, GCG only looks at the sequence type 
   * and the length.
   */
  if (Alphabet_type == kAmino)
    fprintf(fp, "(Peptide) ");
  else if (Alphabet_type == kDNA || Alphabet_type == kRNA)
    fprintf(fp, "(Nucleotide) ");
  else
    { Warn("You can't write profiles of non-biological sequences."); return 0; }
  fprintf(fp, "PROFILEMAKE v0.00 of: foo  Length: %d\n", hmm->M);

  /* Write some commentary. As long as we avoid looking like the
   * first line of the profile (which contains a "..") we should be OK.
   */
  fprintf(fp, "   Ha, fooled you. Actually, this is from HMM %s.\n", RELEASE);
  fprintf(fp, "   Reverse-engineered from GCG 8.0 to fake a profile.\n");
  fprintf(fp, "\n");
  
  /* Write the first line of the profile.
   */
  fprintf(fp, "Cons");
  for (x = 0; x < Alphabet_size; x++)
    fprintf(fp, " %c   ", Alphabet[x]);
  fprintf(fp, "Gap  Len  ..\n");
  
  /* For each consensus match position in the HMM,
   * write a line of profile.
   */
  for (k = 1; k <= shmm->M; k++)
    {
				/* find consensus. */
      max_score = -999;
      for (x = 0; x < 26; x++)
	if (shmm->m_emit[x][k] > max_score)
	  { 
	    max_score = shmm->m_emit[x][k];
	    sym       = 'A' + x;
	  }
      fprintf(fp, " %c ", (char) sym);
				/* emission profile */
      for (x = 0; x < Alphabet_size; x++)
	{
	  sym   = Alphabet[x] - 'A';
	  score = (shmm->m_emit[sym][k] + shmm->t[(k-1)*9 + Tmm]) / 10;
	  fprintf(fp, "%4d ", score);
	}
				/* gap open */
      score = -1 * (shmm->t[k*9 + Tmi] + shmm->t[k*9 + Tim] - shmm->t[k*9 + Tmm]) / 10;
      fprintf(fp, "%4d ", score);
				/* gap extend */
      score = -1 * shmm->t[k*9 + Tii] / 10;
      fprintf(fp, "%4d\n", score);
    }

  /* The final line of the profile contains the base composition.
   */
  fprintf(fp, " * ");
  for (x = 0; x < Alphabet_size; x++)
    fprintf(fp, "%4d ", (int) (1000. * randomseq[x]));
  fprintf(fp, "\n");

  return 1;
}
