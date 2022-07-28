/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1995 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* msf.c
 * SRE, Sun Jul 11 16:17:32 1993
 * 
 * Import/export of GCG MSF multiple sequence alignment
 * formatted files.
 * 
 *************************************************************
 * Specification of an MSF-formatted file:
 *   (non-official. Empirically derived by examination!)
 *   
 * - The header apparently can consist of arbitrary information,
 *   specific to the program that generated the .msf file.
 *   
 * - After the header, a line appears like this:
 *    
 *        picorna.msf  MSF: 100  Type: P  January 17, 1991  17:53  Check: 541
 *        ..
 *   or        MSF:  171  Type: P    Check:  4694   ..
 *   or   pileup.msf  MSF: 617  Type: P  October 7, 1992  12:14  Check: 1  ..
 *
 *   To check for MSF formatted files, we check for "MSF:" && "Type:" 
 *   && "Check:".
 *   
 *   To determine end of header, we check for "..".
 *   
 * - Then follows a list of N lines, containing name, length, 
 *   checksum, and weight information; i.e.:
 *   
 *      Name: GLB2_MORMR      oo  Len:  171  Check:  6522  Weight: 2.7687
 *
 * - Then a separator:
 *       //
 *       
 * - Then the sequences, in blocks like SELEX format, 50 symbols per line
 *   in groups of 10; there may or may not be coordinates written above each line.
 *   "."'s are gaps. We assume the sequences come in the same order their
 *   names did.
 *       
 *          1                                                   50
 *     Cb3  ...gpvedai .......t.. aaigr..vad tvgtgptnse aipaltaaet
 *       E  gvenae.kgv tentna.tad fvaqpvylpe .nqt...... kv.affynrs
 *******************************************************************
 *
 * Functions provided:
 * 
 * int
 * ReadMSF(char *filename, char ***ret_aseqs, char ***ret_names, 
 *         int *num, struct aliinfo_s *ainfo)
 *         
 * int 
 * WriteMSF(FILE *fp, char **aseqs, char **names, int *weights, int num)
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: ReadMSF()
 * 
 * Purpose:  Read multiple aligned sequences from the file seqfile.
 *           Returns aseqs (aligned sequences), names, weights,
 *           and num (number of sequences).
 *           
 *           Memory is allocated for aseqs, names, and weights;
 *           must be free'd by caller.
 *           
 * Return:   1 on success, 0 on failure.  
 */
int
ReadMSF(char    *seqfile,       /* file to read seqs from      */
	char  ***ret_aseqs,     /* RETURN: aligned seqs        */
	int     *ret_num,       /* RETURN: number of seqs      */
	struct aliinfo_s *ainfo)
{
  FILE    *fp;                  /* ptr to opened seqfile        */
  char   **aseqs;               /* aligned seqs                 */
  int      num;			/* number of seqs read          */
  char     buffer[LINEBUFLEN];	/* input buffer for lines       */
  int      blocknum;		/* number of blocks in file     */
  char    *sptr;                /* ptr into sequence on line    */
  int      currblock;		/* index for blocks             */
  int      i;			/* loop counter                 */
  int      idx;			/* counter for seqs             */
  int      grp;			/* counter for 5 seq groups per line */
  int      pos;			/* counter for position in a seq */
  int      count;
  int      inblock;
  
  /***************************************************
   * First pass across file.
   * Verify that it looks like an MSF file.
   * Count names (Name: Len: Check: Weight:)
   * Count blocks (50 * blocks is the maximum seqlen we'll have to deal with.)
   ***************************************************/
			/* open the file for reading */
  fp = fopen(seqfile, "r");
  if (fp == NULL) { squid_errno = SQERR_NOFILE; return 0; }

				/* look for MSF: header line */
  do
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
    }  while ((strstr(buffer, "MSF:") == NULL) &&
	      (strstr(buffer, "Type:") == NULL) &&
	      (strstr(buffer, "Check:") == NULL));

				/* count the names */
  num = 0;
  do
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
      if (strstr(buffer, "Name:")   != NULL &&
	  strstr(buffer, "Len:")    != NULL &&
	  strstr(buffer, "Check:")  != NULL &&
	  strstr(buffer, "Weight:") != NULL)
	num++;
    } while ((strstr(buffer, "//") == NULL));

				/* count blocks of sequence */
  blocknum = 0;
  inblock  = FALSE;
  while (fgets(buffer, LINEBUFLEN, fp) != NULL)
    if (strtok(buffer, WHITESPACE) == NULL)
      inblock = FALSE;
    else if (! inblock) 
      { inblock = TRUE; blocknum++; }

#ifdef SRE_REMOVED
				/* if seqlines isn't evenly divisible
				   by sequences, we have a problem */
  if (seqlines % num && seqlines % (num+1))
    { squid_errno = SQERR_FORMAT; return 0; }
  blocknum = seqlines / num;
#endif

  /***************************************************
   * Rewind file for second pass; skip header; position on first Name line
   ***************************************************/
  rewind(fp);
  do
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
    }  while ((strstr(buffer, "MSF:") == NULL) ||
	      (strstr(buffer, "Type:") == NULL) ||
	      (strstr(buffer, "Check:") == NULL));

  do 
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
    } while (strstr(buffer, "Name:")   == NULL ||
	     strstr(buffer, "Len:")    == NULL ||
	     strstr(buffer, "Check:")  == NULL ||
	     strstr(buffer, "Weight:") == NULL);

  /***************************************************
   * Parse the name lines, to get name and weight
   ***************************************************/

				/* allocations for num sequences */
  if ((aseqs         = (char **)            malloc 
       (num * sizeof(char *)))           == NULL ||
      (ainfo->sqinfo = (struct seqinfo_s *) malloc 
       (num * sizeof(struct seqinfo_s))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  ainfo->flags = 0;

  idx = 0;
  while (strstr(buffer, "Name:") != NULL)
    {
      ainfo->sqinfo[idx].flags = 0;

				/* Name: foo */
      sptr = strtok(buffer, WHITESPACE);
      sptr = strtok(NULL,   WHITESPACE);
      SetSeqinfoString(&(ainfo->sqinfo[idx]), sptr, SQINFO_NAME);

				/* Len: 000 */
				/* Check: 0000 */
				/* Weight: 1.00 */
      while (strcmp(sptr, "Weight:") != 0)
	sptr = strtok(NULL, WHITESPACE);
      sptr = strtok(NULL, WHITESPACE);
      SetSeqinfoString(&(ainfo->sqinfo[idx]), sptr, SQINFO_WGT);
      
      idx++;
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
    }

				/* skip to separator */
  do
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
    } while ((strstr(buffer, "//") == NULL));
	
  /***************************************************
   * Read the sequences.
   ***************************************************/
      
				/* allocations for sequences */
  for (idx = 0; idx < num; idx++)
    if ((aseqs[idx] = (char *) malloc ((blocknum * 50 + 1) * sizeof(char))) == NULL)
      { squid_errno = SQERR_MEM; return 0; }

				/* for each block of sequences: */
  for (currblock = 0; currblock < blocknum; currblock++)
    {
      
				/* skip blank lines, until name field
				   matches names[0] */
      do
	{
	  if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	    { squid_errno = SQERR_NODATA; return 0; }
	} while ((sptr = strtok(buffer, WHITESPACE)) == NULL ||
		 strcmp(ainfo->sqinfo[0].name, sptr) != 0);

				/* for each sequence */
      for (idx = 0; idx < num; idx++)
	{
				/* for each of 5 groups of 10 */
	  for (grp = 0; grp < 5; grp++)
	    {
	      if ((sptr = strtok(NULL, WHITESPACE)) == NULL)
		break;
	      
	      pos = (currblock * 50) + (grp * 10);
	      for (i = 0; i < 10; i++)
		{
		  if (strchr(WHITESPACE, sptr[i]) != NULL) 
		    { aseqs[idx][pos + i] = '\0'; break; }
		  if (isgap((int) sptr[i])) aseqs[idx][pos+i] = '.';
		  else                      aseqs[idx][pos+i] = sptr[i];
		}
	    }
				/* get next line, if we expect one */
	  if (idx < num-1)
	    {
	      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
		{ squid_errno = SQERR_NODATA; return 0; }
	      if ((sptr = strtok(buffer, WHITESPACE)) == NULL)
		{squid_errno = SQERR_NODATA; return 0; }
	    }
	}
    }
				/* guarantee NULL termination all sequences */
  for (idx = 0; idx < num; idx++)
    aseqs[idx][blocknum * 50] = '\0';

				/* make sure alignment is flushed right */
  FlushAlignment(aseqs, num, &(ainfo->alen));
  ainfo->flags |= AINFO_ALEN;

				/* find raw sequence lengths for sqinfo */
  for (idx = 0; idx < num; idx++)
    {
      count = 0;
      for (sptr = aseqs[idx]; *sptr != '\0'; sptr++)
	if (!isgap(*sptr)) count++;
      ainfo->sqinfo[idx].len    = count;
      ainfo->sqinfo[idx].flags |= SQINFO_LEN;
    }

  fclose(fp);
  *ret_num     = num;
  *ret_aseqs   = aseqs;
  return 1;
}


/* Function: WriteMSF()
 * 
 * Purpose:  Write aseqs, names, weights to an open fp,
 *           in GCG MSF format. The alignment must
 *           be flushed (all aseqs the same length, padded
 *           with gaps)
 * 
 * Returns 1 on success. Returns 0 on failure, and sets
 * squid_errno to indicate the cause.
 */
int
WriteMSF(FILE   *fp,            /* open fp for writing           */
	 char  **aseqs,         /* aligned sequences             */
	 int     num,
	 struct aliinfo_s *ainfo)
{
  int    still_going;		/* True if writing another block */
  int    idx;			/* counter for sequences         */
  int    pos;			/* position counter              */
  int    namelen;		/* maximum name length used      */
  int    len;			/* tmp variable for name lengths */
  char   buffer[51];		/* buffer for writing seq        */
  char **sqptr;                 /* ptrs into each sequence       */
  int    charcount;		/* num. symbols we're writing    */
  double weight;

				/* allocate seq pointers that we'll
				   move across each sequence */
  if ((sqptr = (char **) malloc (num * sizeof(char *))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }

				/* set sqptrs to start of each seq */
  for (idx = 0; idx < num; idx++)
    sqptr[idx] = aseqs[idx];
				/* calculate max namelen used */
  namelen = 0;
  for (idx = 0; idx < num; idx++)
    if ((len = strlen(ainfo->sqinfo[idx].name)) > namelen) 
      namelen = len;

  /*****************************************************
   * Write the title line
   *****************************************************/
  fprintf(fp, "\n");
				/* ack! we're writing bullshit here */
  fprintf(fp, "    MSF:  000  Type: X  Check: 0000  ..\n");
  fprintf(fp, "\n");

  /*****************************************************
   * Write the names
   *****************************************************/

  for (idx = 0; idx < num; idx++)
    {
      weight = 1.0;
      if (ainfo->sqinfo[idx].flags & SQINFO_WGT)
	weight = ainfo->sqinfo[idx].weight;
      fprintf(fp, "  Name: %-*.*s  Len:  %5d  Check:  %5d  Weight: %.4f\n",
	      namelen, namelen,
	      ainfo->sqinfo[idx].name,
	      ainfo->alen,
	      GCGchecksum(aseqs[idx], ainfo->alen),
	      weight);
    }
  fprintf(fp, "\n");
  fprintf(fp, "//\n");
  fprintf(fp, "\n");

  /*****************************************************
   * Write the sequences
   *****************************************************/

  still_going = 1;
  while (still_going)
    {
      still_going = 0;
      for (idx = 0; idx < num; idx++)
	{
	  fprintf(fp, "%-*.*s  ", namelen, namelen, 
		  ainfo->sqinfo[idx].name);

				/* get next line's worth of 50 from seq */
	  strncpy(buffer, sqptr[idx], 50);
	  buffer[50] = '\0';
	  charcount = strlen(buffer);

				/* is there still more to go? */
	  if (charcount == 50 && sqptr[idx][50] != '\0')
	    still_going = 1;

				/* shift the seq ptr by a line */
	  sqptr[idx] += charcount;

				/* draw the sequence line */
	  pos = 0; 
	  while (pos < charcount)
	    {
	      if (isgap(buffer[pos])) fputc('.', fp);
	      else fputc(buffer[pos], fp);
	      pos++;
	      if (!(pos % 10)) fputc(' ', fp);
	    }
	  fputc('\n', fp);
	}
				/* put blank line between blocks */
      fputc('\n', fp);
    }

  free(sqptr);
  return 1;
}



void
FlushAlignment(char **aseqs, int num, int *ret_alen)
{
  int len, alen;
  int idx;
  int apos;

  alen = strlen(aseqs[0]);
  for (idx = 1; idx < num; idx++)
    if ((len = strlen(aseqs[idx])) > alen)
      alen = len;
  
  for (idx = 0; idx < num; idx++)
    {
      if ((aseqs[idx] = (char *) realloc (aseqs[idx], sizeof(char) * (alen+1))) == NULL)
	Die("realloc failed");
      for (apos = strlen(aseqs[idx]); apos < alen; apos++)
	aseqs[idx][apos] = '.';
      aseqs[idx][apos] = '\0';
    }
  *ret_alen = alen;
}
