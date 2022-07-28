/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmio.c
 * modified SRE, Thu Jul 15 12:21:49 1993: v1.1
 * SRE, Wed Sep 21 14:25:34 1994: v1.7
 * 
 * Input/output of hmm models.
 *
 * HMM's can either be saved in a semi-readable, semi-documented
 * ASCII format, or a binary format. The ASCII format is guaranteed
 * to be portable across platforms. The binary format takes up
 * about 6-fold less disk space and is faster to read and write,
 * but is not guaranteed to be portable across machines.
 * 
 * As multiple versions of this program will probably arise
 * eventually, magic numbers (both for the ASCII and binary
 * save formats) are used to label save files with a major
 * version number. ReadHMM() is written in such a way that
 * it determines whether a save file is ASCII or binary, and
 * what version it comes from, and then hands the file to
 * the appropriate parser. This simplifies the task of
 * backwards compatibility as new versions of the program
 * are created.
 *
 * V1.0: original implementation
 * V1.1: regularizers removed from model structure
 * V1.7: ref and cs annotation lines added from alignment, one 
 *       char per match state 1..M
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h> /* to get SEEK_CUR definition on silly Suns */
#include "squid.h"
#include "states.h"
#include "externs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define HMM_VERS10  1
#define HMM_VERS11  2
#define HMM_VERS17  3

static unsigned int v10magic = 0xe8ededb1; /* v1.0 binary files: "hmm1" + 0x80808080 */
static unsigned int v10swap  = 0xb1edede8; /* byteswapped v1.0 */
static unsigned int v11magic = 0xe8ededb2; /* v1.1 binary: "hmm2" + 0x80808080 */
static unsigned int v11swap  = 0xb2edede8; /* byteswapped v1.1 */
static unsigned int v17magic = 0xe8ededb3; /* v1.7 binary: "hmm3" + 0x80808080 */
static unsigned int v17swap  = 0xb3edede8; /* byteswapped v1.7 */

static struct hmm_struc *read_hmm(FILE *fp, int version); 
static struct hmm_struc *read_binhmm(FILE *fp, int version, int swapped);
static void              byteswap(char *swap, int nbytes);

/* Function: WriteHMM()
 * 
 * Write an HMM to disk in ASCII format.
 * Identify it as a major version 1.7 save file.
 * Return 1 on success, 0 on failure.
 */
int 
WriteHMM(FILE *fp, struct hmm_struc *hmm)
{
  int k;			/* counter for nodes */
  int i;			/* counter for symbols */

				/* we write a "magic" header w/ version info */
  fprintf(fp, "# HMM v1.7\n");
				/* write M, the length of the model, and alphabet_size.
				   Because we  can't allocate for an HMM until we know 
				   these numbers, these must be the first pieces of data in 
				   the file. */
  fprintf(fp, "%d\t# M -- length of model\n", hmm->M);
  fprintf(fp, "%d\t# alphabet length\n", Alphabet_size);

				/* write info on the alphabet this HMM is used for */
  fprintf(fp, "%d\t# alphabet type\n", Alphabet_type);
  fprintf(fp, "%s\t# alphabet\n", Alphabet);

				/* version 1.7: optional info: ref and cs lines */
  fprintf(fp, "%s\t# Optional reference line annotation?\n",
	  (hmm->flags & HMM_REF) ? "yes" : "no");
  fprintf(fp, "%s\t# Optional consensus structure annotation?\n",
	  (hmm->flags & HMM_CS) ? "yes" : "no");
  
				/* write every state */
  for (k = 0; k <= hmm->M; k++)
    {
				/* match substate */
      fprintf(fp, "###MATCH_STATE %d (%c) (%c)\n", k, 
	      (k > 0 && (hmm->flags & HMM_REF)) ? hmm->ref[k] : (char) ' ',
	      (k > 0 && (hmm->flags & HMM_CS)) ? hmm->cs[k] : (char) ' ');
      fprintf(fp, "%f\t# t_m1\n", hmm->mat[k].t[MATCH]);
      fprintf(fp, "%f\t# t_d1\n", hmm->mat[k].t[DELETE]);
      fprintf(fp, "%f\t# t_i0\n", hmm->mat[k].t[INSERT]);
      for (i = 0; i < Alphabet_size; i++)
	fprintf(fp, "%f\t# Symbol %c probability\n",
		hmm->mat[k].p[i], Alphabet[i]);
      
      /* delete substate */
      fprintf(fp, "###DELETE_STATE %d\n", k);
      fprintf(fp, "%f\t# t_m1\n", hmm->del[k].t[MATCH]);
      fprintf(fp, "%f\t# t_d1\n", hmm->del[k].t[DELETE]);
      fprintf(fp, "%f\t# t_i0\n", hmm->del[k].t[INSERT]);
      
      /* insert substate */
      fprintf(fp, "###INSERT_STATE %d\n", k);
      fprintf(fp, "%f\t# t_m1\n", hmm->ins[k].t[MATCH]);
      fprintf(fp, "%f\t# t_d1\n", hmm->ins[k].t[DELETE]);
      fprintf(fp, "%f\t# t_i0\n", hmm->ins[k].t[INSERT]);
      for (i = 0; i < Alphabet_size; i++)
	fprintf(fp, "%f\t# Symbol %c probability\n",
		hmm->ins[k].p[i], Alphabet[i]);
    }
  return 1;
}



/* Function: ReadHMM()
 * 
 * Read an HMM from a save file. Automatically
 * determines whether the save file is binary
 * or ASCII, and what version it came from.
 * 
 * Returns a pointer to the (allocated) new HMM,
 * or NULL on failure.
 */
struct hmm_struc *
ReadHMM(FILE *fp)
{
  char buffer[512];
  char *s;
  unsigned int  magic_number;

  /* get header and find version */
  
  /* look to see if it's a binary file */
  fread((char *) &magic_number, 4, 1, fp);
  if      (magic_number == v10magic) return(read_binhmm(fp, HMM_VERS10, FALSE));
  else if (magic_number == v10swap)  return(read_binhmm(fp, HMM_VERS10, TRUE));
  else if (magic_number == v11magic) return(read_binhmm(fp, HMM_VERS11, FALSE));
  else if (magic_number == v11swap)  return(read_binhmm(fp, HMM_VERS11, TRUE));
  else if (magic_number == v17magic) return(read_binhmm(fp, HMM_VERS17, FALSE));
  else if (magic_number == v17swap)  return(read_binhmm(fp, HMM_VERS17, TRUE));

  /* If it looks like a binary but we don't recognize its magic, demur */
  if (magic_number & 0x80000000)
    {
      Warn("Binary HMM file newer than current installed executables!?\n");
      Warn("It would be unsafe to try to read that file.\n");
      Warn("Update your HMM software.\n");
      return 0;
    }

  /* nope, must be an ASCII file */
  rewind(fp);
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  s = strtok(buffer, "\t\n ");
  s = strtok((char *) NULL, "\t\n ");
  if (strcmp(s, "HMM") != 0) return NULL;
  s = strtok((char *) NULL, "\t\n ");
  
  /* hand file to appropriate reader */
  if      (strcmp(s, "v1.0") == 0) return (read_hmm(fp, HMM_VERS10));
  else if (strcmp(s, "v1.1") == 0) return (read_hmm(fp, HMM_VERS11));
  else if (strcmp(s, "v1.7") == 0) return (read_hmm(fp, HMM_VERS17));

  /* else: we failed to parse the file. */
  return NULL;
}





/* Function: read_hmm()
 * 
 * Purpose:  Read ASCII-format save files.
 *           V1.0 contained sympvec and regularizers; these are ignored
 *                in V1.1 and later
 *           V1.7 and later contain ref and cs annotation.
 *
 * Args:     fp      - open save file, header has been read already
 *           version - HMM_VERS17, for instance
 *
 * Returns ptr to the (allocated) new HMM on success,
 * or NULL on failure.
 */
static struct hmm_struc *
read_hmm(FILE *fp, int version)
{
  struct hmm_struc *hmm;
  int   M;			/* length of model  */
  char buffer[512];
  char *statetype;
  char *s;
  int   k;			/* state number  */
  int   i;			/* symbol number */
  
				/* read M from first line */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit(*s)) return NULL;
  M = atoi(s);
				/* read alphabet_length */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit(*s)) return NULL;
  Alphabet_size = atoi(s);
				/* now, create space for hmm */
  if ((hmm = AllocHMM(M)) == NULL)
    Die("malloc failed for reading hmm in\n");
  
				/* read alphabet_type */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (!isdigit(*s)) return NULL;
  Alphabet_type = atoi(s);
				/* read alphabet */
  if (fgets(buffer, 512, fp) == NULL) return NULL;
  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
  if (strlen(s) != Alphabet_size) return NULL;
  strncpy(Alphabet, s, Alphabet_size);
				
  /* skip the random symbol frequencies in V1.0 files. now unused */
  if (version == HMM_VERS10)
    for (i = 0; i < Alphabet_size; i++)
      if (fgets(buffer, 512, fp) == NULL) return NULL;

  /* V1.7 has lines for whether we have valid ref, cs info
   */
  if (version == HMM_VERS17)
    {
      if (fgets(buffer, 512, fp) == NULL) return NULL;
      if (strncmp(buffer, "yes", 3) == 0) hmm->flags |= HMM_REF;
      if (fgets(buffer, 512, fp) == NULL) return NULL;
      if (strncmp(buffer, "yes", 3) == 0) hmm->flags |= HMM_CS;
    }

				/* everything else is states */
  while (fgets(buffer, 512, fp) != NULL)
    {
				/* get state type and index info */
      if ((statetype = strtok(buffer, " \t\n")) == NULL) return NULL;
      if ((s = strtok((char *) NULL, " \t\n")) == NULL) return NULL;
      if (!isdigit(*s)) return NULL;
      k = atoi(s);
      if (k < 0 || k > hmm->M+1) return NULL;
      
      if (strcmp(statetype, "###MATCH_STATE") == 0)
	{
				/* V1.7: get ref, cs info:   */
	                        /* ###MATCH_STATE 16 (x) (H) */
	  if (version == HMM_VERS17)
	    {
	      s = strtok(NULL, "\n");
	      while (*s != '(' && *s != '\0') s++;
	      if (*s != '(') return NULL;
	      hmm->ref[k] = *(s+1);
	      while (*s != '(' && *s != '\0') s++;
	      if (*s != '(') return NULL;
	      hmm->cs[k] = *(s+1);
	    }

	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->mat[k].t[INSERT] = (float) atof(s);
	  
	  for (i = 0; i < Alphabet_size; i++)
	    {
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	      if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	      hmm->mat[k].p[i] = (float) atof(s);
	    }

				/* Skip all regularizer info for V1.0 */
	  if (version == HMM_VERS10)
	    for (i = 0; i < Alphabet_size + 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;

	}
      else if (strcmp(statetype, "###INSERT_STATE") == 0)
	{
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->ins[k].t[INSERT] = (float) atof(s);
	  
	  for (i = 0; i < Alphabet_size; i++)
	    {
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	      if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	      hmm->ins[k].p[i] = (float) atof(s);
	    }
	  
	  /* Skip all regularizer info in V1.0 files */
	  if (version == HMM_VERS10)
	    for (i = 0; i < Alphabet_size + 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;

	}
      else if (strcmp(statetype, "###DELETE_STATE") == 0)
	{
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[MATCH] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[DELETE] = (float) atof(s);
	  
	  if (fgets(buffer, 512, fp) == NULL) return NULL;
	  if ((s = strtok(buffer, " \t\n")) == NULL) return NULL;
	  hmm->del[k].t[INSERT] = (float) atof(s);
	  
	  /* Skip all regularizer info in V1.0 files*/
	  if (version == HMM_VERS10)
	    for (i = 0; i < 3; i++)
	      if (fgets(buffer, 512, fp) == NULL) return NULL;
	}
      else
	return NULL;
    }
  Renormalize(hmm);
  return hmm;
}




/* Function: WriteBinaryHMM()
 * 
 * Save an HMM in more compact, but less portable, binary format. 
 * 
 * Returns 1 on success, 0 on failure.
 */
int 
WriteBinaryHMM(FILE *fp, struct hmm_struc *hmm)
{
  int k;			/* counter for nodes */
  
  /* Write the four-byte magic number. It identifies the file
   * as binary (because the high bits are set), and identifies
   * the major version of the program. It can be used to
   * identify and possibly humor byte-swapped architectures,
   * although we don't bother for now.
   */
  fwrite((char *) &v17magic, 4, 1, fp);

  /* write M, the length of the model, and alphabet_size.
   * Because we  can't allocate for an HMM until we know 
   * these numbers, these must be the first pieces of data 
   * in the file. 
   */
  fwrite((char *) &(hmm->M),        sizeof(int), 1, fp);
  fwrite((char *) &(Alphabet_size), sizeof(int), 1, fp);
  
  /* write info on the alphabet this HMM is used for */
  fwrite((char *) &(Alphabet_type), sizeof(int), 1, fp);
  fwrite((char *) Alphabet, sizeof(char), Alphabet_size, fp);
  
  /* write the optional info, flags first */
  fwrite((char *) &(hmm->flags), sizeof(int), 1, fp);
  if (hmm->flags & HMM_REF) fwrite((char *) hmm->ref, sizeof(char), hmm->M+1, fp);
  if (hmm->flags & HMM_CS)  fwrite((char *) hmm->cs,  sizeof(char), hmm->M+1, fp);

  /* write every state */
  for (k = 0; k <= hmm->M; k++)
    {
      /* match substate */
      fwrite((char *) &(hmm->mat[k].t[MATCH]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->mat[k].t[DELETE]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->mat[k].t[INSERT]), sizeof(float), 1, fp);
      fwrite((char *) hmm->mat[k].p, sizeof(float), Alphabet_size, fp);
      
      /* delete substate */
      fwrite((char *) &(hmm->del[k].t[MATCH]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->del[k].t[DELETE]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->del[k].t[INSERT]), sizeof(float), 1, fp);
      
      /* insert substate */
      fwrite((char *) &(hmm->ins[k].t[MATCH]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->ins[k].t[DELETE]), sizeof(float), 1, fp);
      fwrite((char *) &(hmm->ins[k].t[INSERT]), sizeof(float), 1, fp);
      fwrite((char *) hmm->ins[k].p, sizeof(float), Alphabet_size, fp);
    }
  return 1;
}




/* Function: read_binhmm()
 * 
 * Read binary HMM save files.
 * V1.0 saved regularizer and sympvec info, which V1.1 ignores.
 * V1.7 and later may include optional ref, cs annotation lines.
 * 
 * Returns pointer to the HMM on success; NULL
 * on failure.
 */
static struct hmm_struc *
read_binhmm(FILE *fp, int version, int swapped)
{
  struct hmm_struc *hmm;
  int   M;			/* length of model  */
  int   k;			/* state number  */
  int   x;			/* symbol or transition number */
  
  /* read M and alphabet_size */
  if (! fread((char *) &(M), sizeof(int), 1, fp))  return NULL;
  if (! fread((char *) &Alphabet_size, sizeof(int), 1, fp)) return NULL;
  if (swapped) { 
    byteswap((char *) &M, sizeof(int));
    byteswap((char *) &Alphabet_size, sizeof(int));
  }
  
  /* now, create space for hmm */
  if ((hmm = AllocHMM(M)) == NULL)
    Die("malloc failed for reading hmm in\n");
  
  /* read alphabet_type and alphabet*/
  if (! fread((char *) &Alphabet_type, sizeof(int), 1, fp)) return NULL;
  if (swapped) byteswap((char *) &Alphabet_type, sizeof(int));
  if (! fread((char *) Alphabet, sizeof(char), Alphabet_size, fp)) return NULL;
  /* skip the random symbol frequencies in V1.0 */
  if (version == HMM_VERS10)
    fseek(fp, (long) (sizeof(float) * Alphabet_size), SEEK_CUR);
  
  /* Get optional info in V1.7 and later
   */
  if (version == HMM_VERS17)
    {
      if (! fread((char *) &(hmm->flags), sizeof(int), 1, fp)) return NULL;
      if (swapped) byteswap((char *) &hmm->flags, sizeof(int));
      if ((hmm->flags & HMM_REF) &&
	  ! fread((char *) hmm->ref, sizeof(char), hmm->M+1, fp)) return NULL;
      hmm->ref[hmm->M+1] = '\0';
      if ((hmm->flags & HMM_CS) &&
	  ! fread((char *) hmm->cs,  sizeof(char), hmm->M+1, fp)) return NULL;
      hmm->cs[hmm->M+1]  = '\0';
    }

  /* everything else is states */
  for (k = 0; k <= hmm->M; k++)
    {
      /* get match state info */
      if (! fread((char *) &(hmm->mat[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->mat[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->mat[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) hmm->mat[k].p, sizeof(float), Alphabet_size, fp)) return NULL;
      if (swapped) {
	byteswap((char *) &(hmm->mat[k].t[MATCH]),  sizeof(float));
	byteswap((char *) &(hmm->mat[k].t[DELETE]), sizeof(float));
	byteswap((char *) &(hmm->mat[k].t[INSERT]), sizeof(float));
	for (x = 0; x < Alphabet_size; x++)
	  byteswap((char *) &(hmm->mat[k].p[x]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMM_VERS10)
	fseek(fp, (long)(sizeof(float) * (3 + Alphabet_size)), SEEK_CUR);
      
      /* get delete state info */
      if (! fread((char *) &(hmm->del[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->del[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->del[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (swapped) {
	byteswap((char *) &(hmm->del[k].t[MATCH]),  sizeof(float));
	byteswap((char *) &(hmm->del[k].t[DELETE]), sizeof(float));
	byteswap((char *) &(hmm->del[k].t[INSERT]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMM_VERS10)
	fseek(fp, (long)(sizeof(float) * 3), SEEK_CUR);
      
      /* get insert state info */
      if (! fread((char *) &(hmm->ins[k].t[MATCH]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->ins[k].t[DELETE]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) &(hmm->ins[k].t[INSERT]), sizeof(float), 1, fp)) return NULL;
      if (! fread((char *) hmm->ins[k].p, sizeof(float), Alphabet_size, fp)) return NULL;
      if (swapped) {
	byteswap((char *) &(hmm->ins[k].t[MATCH]),  sizeof(float));
	byteswap((char *) &(hmm->ins[k].t[DELETE]), sizeof(float));
	byteswap((char *) &(hmm->ins[k].t[INSERT]), sizeof(float));
	for (x = 0; x < Alphabet_size; x++)
	  byteswap((char *) &(hmm->ins[k].p[x]), sizeof(float));
      }
      
      /* skip the regularizer info in V1.0 */
      if (version == HMM_VERS10)
	fseek(fp, (long)(sizeof(float) * (3 + Alphabet_size)), SEEK_CUR);
    }
  Renormalize(hmm);
  return hmm;
}


/* Function: byteswap()
 * 
 * Purpose:  Swap between big-endian and little-endian.
 *           For example:
 *               int foo = 0x12345678;
 *               byteswap((char *) &foo, sizeof(int));
 *               printf("%x\n", foo)
 *           gives 78563412.
 *           
 *           I don't fully understand byte-swapping issues.
 *           However, I have tested this on chars through doubles,
 *           on various machines:
 *               SGI IRIX 4.0.5, SunOS 4.1.3, DEC Alpha OSF/1, Alliant
 *
 * Date: Sun Feb 12 10:26:22 1995              
 */
static void
byteswap(char *swap, int nbytes)
{
  int  x;
  char byte;
  
  for (x = 0; x < nbytes / 2; x++)
    {
      byte = swap[nbytes - x - 1];
      swap[nbytes - x - 1] = swap[x];
      swap[x] = byte;
    }
}
