/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* convert_main.c
 * Thu Jun 30 10:23:28 1994
 * 
 * Convert an HMM file to flat file format, or binary.
 * Useful for moving HMM's between architectures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif


#include "states.h"
#include "externs.h" 
#include "squid.h" 
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char Alphabet[MAXABET];		/* ACGT, for instance     */
int  Alphabet_size;		/* 4 or 20                */
int  Alphabet_type;		/* kDNA, kRNA, or kAmino  */

#define OPTIONS "bhP"

static char usage[]  = "\
Usage: hmm_convert [-options] <hmmfile input> <hmmfile output>\n\
where options are:\n\
   -b : write output in efficient binary format (default is portable ASCII)\n\
   -h : show short help, version, usage info\n\
  EXPERIMENTAL OPTIONS:\n\
   -P : write a GCG 8.0 compatible profile in .prf format\n";

static char banner[] = "hmm_convert: hidden Markov model conversion, for model portability";

int
main(int argc, char **argv)
{ 
  struct hmm_struc *hmm;        /* hidden Markov model     */ 
  char  *infile;                /* name of input HMM file  */
  char  *outfile;               /* name of output HMM file */
  FILE  *hmmfp;                 /* open HMM file pointer   */
  float  randomseq[MAXABET];	/* random background model */

  int    do_binary;		/* TRUE to write binary    */
  int    do_prf;		/* TRUE to write a profile */

  int          optc;
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /*********************************************** 
   * Parse command line
   ***********************************************/

				/* set defaults on training parameter configuration */
  do_binary = FALSE;
  do_prf    = FALSE;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {

    case 'b': do_binary = TRUE; break;
    case 'P': do_prf    = TRUE; break;

    case 'h': 
      printf("%s\n version %s, %s\n\n%s\n", banner, RELEASE, RELEASEDATE, usage);
      exit(0);
    default:
      Die("Error: unrecognized option %c\n", optc);
    }

  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);
  infile  = argv[argc-2];
  outfile = argv[argc-1]; 

  /*********************************************** 
   * Ready to start. Print banner information
   ***********************************************/
  printf("%s\n     version %s, %s\n", banner, RELEASE, RELEASEDATE);
  printf(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Do the conversion
   ***********************************************/

  /* Read the model in. Also sets alphabet.
   */
  if ((hmmfp = fopen(infile, "rb")) == NULL)
    Die("Error: failed to open HMM file %s for reading.\n", infile);
  if ((hmm = ReadHMM(hmmfp)) == NULL)
    Die("Error: failed to parse HMM file %s\n", infile);
  fclose(hmmfp);
  DefaultRandomModel(randomseq);

  /* Write the model out.
   */
  if ((hmmfp = fopen(outfile, "wb")) == NULL) 
    Die("Error: failed to open HMM file %s for writing.\n", outfile);
  if (do_prf)
    {
      struct shmm_s  *shmm;

      shmm = AllocSearchHMM(hmm->M);
      MakeSearchHMM(hmm, randomseq, shmm);

      if (! WriteProfile(hmmfp, hmm, shmm, randomseq))
	Die("Failed to convert HMM %s to profile %s\n", infile, outfile);
      printf("Converted %s to GCG PROFILE compatible file %s\n", infile, outfile);
      
      FreeSearchHMM(shmm);
    }
  else if (do_binary)
    {
      if (! WriteBinaryHMM(hmmfp, hmm))
	Die("Failed to save hmm to %s\n", outfile);
      printf("Converted %s to a binary HMM file %s\n", infile, outfile);
    }
  else
    {
      if (!WriteHMM(hmmfp, hmm))
	Die("Failed to save hmm to %s\n", outfile); 
      printf("Converted %s into a portable ASCII HMM file %s\n", infile, outfile);
    }
  fclose(hmmfp);

  free(hmm);
  return (0); 
}

