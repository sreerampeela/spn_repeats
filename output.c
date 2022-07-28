/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* output.c
 * Tue Jan 18 09:23:13 1994
 * 
 * Fancy output of matches to a model.
 * 
 * 
 */

#define CPL    50
#define CUTOFF 2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <ctype.h>

#include "states.h"
#include "externs.h"

#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
#ifdef DEBUG
#include <assert.h>
#endif

/* Function: PrintFancyTrace()
 * 
 * Purpose:  Print an alignment of an HMM to a sequence, given a traceback.
 *           Somewhat inspired by the output style of BLAST, except that
 *           we're aligning to a complicated model that's difficult to
 *           represent compactly. 
 *
 * Arguments: ofp       - where to print it (open FILE for writing, or stdout)
 *            shmm      - log-odds form HMM
 *            tr        - traceback from ViterbiTrace()
 *            seq       - sequence that is aligned
 *            seqname   - name of seq to print in left margin
 *            from_pos  - first position in seq that aligns (0..seqlen-1)
 *
 * Returns:  (void)
 */
void
PrintFancyTrace(FILE             *ofp,
		struct shmm_s    *shmm,
		struct trace_s   *tr,
		char             *seq,
		char             *seqname,
		int               from_pos)
{
  char *model;                  /* display of model                */
  char *mline;			/* display of match/mismatch       */
  char *rfline;			/* display of reference seq        */
  char *csline;			/* display of consensus struct     */
  char *aseq;                   /* display of aligned sequence     */
  int   rpos;                   /* position in raw seq             */
  int   apos;                   /* position in traceback/alignment */
  float score;			/* score for position              */
  float max_score;		/* best score for position         */
  char  bestsym;		/* best match sym at position      */
  int   idx;			/* counter for alphabet            */
  int   len;			/* current length of display printed */
  char  buffer[CPL+1];          /* buffer for lines of display     */ 
  int   startpos, endpos;

  /* Memory allocation.
   */
  if ((rfline = (char *) malloc (sizeof(char) * (tr->tlen + 1))) == NULL ||
      (csline = (char *) malloc (sizeof(char) * (tr->tlen + 1))) == NULL ||
      (model  = (char *) malloc (sizeof(char) * (tr->tlen + 1))) == NULL ||
      (mline  = (char *) malloc (sizeof(char) * (tr->tlen + 1))) == NULL ||
      (aseq   = (char *) malloc (sizeof(char) * (tr->tlen + 1))) == NULL)
    Die("memory allocation failed at %s:%d", __FILE__, __LINE__);
  memset(rfline, ' ', tr->tlen);
  memset(csline, ' ', tr->tlen);
  memset(model,  ' ', tr->tlen);
  memset(mline,  ' ', tr->tlen);
  memset(aseq,   ' ', tr->tlen);

  /* Create the displays of model and aligned sequence.
   * Ignore BEGIN (apos == 0) and END (apos == N-1) in the traceback.
   */
  rpos = from_pos;
  for (apos = 1; apos < tr->tlen-1; apos++)
    {
				/* find best sym at this model position */
      if (tr->statetype[apos] != INSERT)
	{
	  max_score = -999;
	  for (idx = 0; idx < 26; idx++)
	    if (shmm->m_emit[idx][tr->nodeidx[apos]] > max_score)
	      { 
		max_score = shmm->m_emit[idx][tr->nodeidx[apos]];
		bestsym   = (char) ('A' + idx);
	      }
	  if (max_score > (int)(CUTOFF * INTSCALE))
	    model[apos] = toupper((int) bestsym);
	  else
	    model[apos] = tolower((int) bestsym);
	}
      else
	model[apos] = '.';

				/* construct mline (match/mismatch display), rfline, and csline */
      switch (tr->statetype[apos]) {
      case MATCH:
	score = shmm->m_emit[seq[rpos]-'A'][tr->nodeidx[apos]];
	if (seq[rpos] == bestsym) 
	  mline[apos] = bestsym;
	else if (score > 0)       
	  mline[apos] = '+';
	aseq[apos] = seq[rpos];
	if (shmm->flags & HMM_REF) rfline[apos] = shmm->ref[tr->nodeidx[apos]];
	if (shmm->flags & HMM_CS)  csline[apos] = shmm->cs[tr->nodeidx[apos]];
	rpos++;
	break;

      case INSERT:
	aseq[apos] = seq[rpos];
	rpos++;
	break;

      case DELETE:
	aseq[apos] = '-';
	if (shmm->flags & HMM_REF) rfline[apos] = shmm->ref[tr->nodeidx[apos]];
	if (shmm->flags & HMM_CS)  csline[apos] = shmm->cs[tr->nodeidx[apos]];
	break;

      default: Die("Unrecognized statetype %d at %d in traceback", 
		   tr->statetype[apos], apos);
      }
    }
  /* Null terminate, and tack on asterisks to represent BEGIN and
   * END dummy states in model.
   */
  model[0]          = '*';	/* begin */
  model[tr->tlen-1] = '*';      /* end   */
  model[tr->tlen]   = '\0';
  aseq[tr->tlen]    = '\0';
  mline[tr->tlen]   = '\0';
  csline[tr->tlen]   = '\0';
  rfline[tr->tlen]   = '\0';

  /* Print out the display.
   */
  fprintf(ofp, "  Alignment to HMM consensus:\n");
  buffer[CPL] = '\0';
  len = 0; 
  rpos     = from_pos + 1;
  while (len < tr->tlen)	
    {
      startpos = rpos;
				/* rf line reference coord line */
      if (shmm->flags & HMM_REF)
	{
	  strncpy(buffer, rfline+len, CPL);
	  fprintf(ofp, "               REF %s\n", buffer);
	}
				/* cs consensus structure line */
      if (shmm->flags & HMM_CS)
	{
	  strncpy(buffer, csline+len, CPL);
	  fprintf(ofp, "                CS %s\n", buffer);
	}
				/* model */
      strncpy(buffer, model+len, CPL);
      fprintf(ofp, "                   %s\n", buffer);
				/* mline */
      strncpy(buffer, mline+len, CPL);
      fprintf(ofp, "                   %s\n", buffer);
				/* get coords of this aseq block */
      for (apos = len; aseq[apos] != '\0' && apos < len + CPL; apos++)
	if (! isgap(aseq[apos]))
	  rpos++;
      endpos = rpos-1;

				/* aligned sequence */
      strncpy(buffer, aseq+len, CPL);
      fprintf(ofp, "  %10.10s %5d %s %5d\n", seqname, startpos, buffer, endpos);          

      len += CPL;
      fprintf(ofp, "\n");
    }

  /* Done. Free memory and return.
   */
  fflush(ofp);
  free(model);
  free(aseq);
  free(mline);
  free(rfline);
  free(csline);
  return;
}
