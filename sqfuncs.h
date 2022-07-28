/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1995 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

#ifndef SQFUNCSH_INCLUDED
#define SQFUNCSH_INCLUDED
/* sqfuncs.h
 * SRE, Mon Jul 12 12:20:00 1993
 * 
 * Prototypes for squid library functions;
 * also makes a good reference list for what the package contains.
 * Slowly being added to as I ANSI-fy squid.
 */

/* 
 * from aligneval.c
 */
extern float ComparePairAlignments(char *known1, char *known2, char *calc1, char *calc2);
extern float CompareRefPairAlignments(int *ref, char *known1, char *known2, char *calc1, char *calc2);
extern float CompareMultAlignments(char **kseqs, char **tseqs, int    N);
extern float CompareRefMultAlignments(int *ref, char **kseqs, char **tseqs, int    N);
extern void  PairwiseIdentity(char *s1, char *s2, int allow_ragged,
			      float *ret_id1, float *ret_id2, float *ret_idtot);

/* 
 * from alignio.c
 */
extern void FreeAlignment(char **aseqs, int nseq, struct aliinfo_s *ainfo);
extern int  MakeAlignedString(char *aseq, int alen, char *ss, char **ret_s);
extern int  MakeDealignedString(char *aseq, int alen, char *ss, char **ret_s);
extern int  ReadAlignment(char *seqfile, int format, char ***ret_aseqs, int *ret_num, 
			  struct aliinfo_s *ainfo);
extern int  WritePairwiseAlignment(FILE *ofp, char *aseq1, char *name1, int spos1,
				   char *aseq2, char *name2, int spos2,
				   int **pam, int indent);
extern int  MingapAlignment(char **aseqs, int num, struct aliinfo_s *ainfo);
extern int  RandomAlignment(char **rseqs, SQINFO *sqinfo, int nseq, double pop, double pex,
			    char ***ret_aseqs, AINFO *ainfo);

/* from cluster.c
 */
extern int Cluster(float **mx, int N, enum clust_strategy mode, struct phylo_s **ret_tree);
extern struct phylo_s *AllocPhylo(int N);
extern void FreePhylo(struct phylo_s *tree, int N);
extern int  MakeDiffMx(char **aseqs, int num, int alen, int allow_ragged, 
		       float ***ret_dmx);
extern int  CompareSeqs(char *s1, char *s2, int alen, int allow_ragged_ends, 
			float *ret_diff);

/* 
 * from dayhoff.c
 */
extern int  ParsePAMFile(FILE *fp, int ***ret_pam, float *ret_scale);
extern void ScalePAM(int **pam, int scale);


/* 
 * from sqio.c
 */
extern void  FreeSequence(char *seq, SQINFO *sqinfo);
extern int   SetSeqinfoString(SQINFO *sqinfo, char *sptr, int flag);
extern void  SeqinfoCopy(SQINFO *sq1, SQINFO *sq2);
extern void  ToDNA(char *seq);
extern void  ToRNA(char *seq);
extern int   ReadMultipleRseqs(char *seqfile, int fformat, char ***ret_rseqs, 
			       SQINFO **ret_sqinfo, int *ret_num);
extern SQFILE *SeqfileOpen(char *filename, int format, char *env);
extern void    SeqfilePosition(SQFILE *sqfp, long offset);
extern void    SeqfileClose(SQFILE *sqfp);
extern int   ReadSeq(SQFILE *fp, int format, char **ret_seq, SQINFO *sqinfo);
extern int   GCGchecksum(char *seq, int seqlen);
extern int   GCGMultchecksum(char **seqs, int nseq);
extern int   SeqfileFormat(char *filename, int  *ret_format, char *env);
extern int   WriteSeq(FILE *outf, int outfmt, char *seq, SQINFO *sqinfo);
extern int   Seqtype(char *seq);
extern char *SeqFormatString(int code);
extern GSIFILE *GSIOpen(char *gsifile);
extern int   GSIGetOffset(GSIFILE *gsi, char *key, char *sqfile, long *ret_offset);
extern void  GSIClose(GSIFILE *gsi);


/* 
 * from msf.c
 */
extern int  ReadMSF(char *seqfile, char ***ret_aseqs, int *ret_num, 
		   struct aliinfo_s *ret_ainfo);
extern int  WriteMSF(FILE *fp, char **aseqs, int num, struct aliinfo_s *ainfo);
extern void FlushAlignment(char **aseqs, int num, int *ret_alen);

/* from revcomp.c
 */
extern char *revcomp(char *comp, char *seq);

/* 
 * from selex.c
 */
extern int ReadSELEX(char *seqfile, char ***ret_aseqs, int *ret_num, 
		     struct aliinfo_s *ret_aliinfo);
extern int WriteSELEX(FILE *fp, char **aseqs, int num, struct aliinfo_s *ainfo, int cpl);
extern int DealignAseqs(char **aseqs, int num, char ***ret_rseqs);
extern int IsSELEXFormat(char *filename);
extern int TruncateNames(char **names, int N); /* OBSOLETE? */
extern FILE *EnvFileOpen(char *fname, char *env);

/* 
 * from sqerror.c
 */
extern int Die(char *format, ...);
extern int Warn(char *format, ...);

/* from sre_ctype.c
 */
extern int sre_tolower(int c);
extern int sre_toupper(int c);

/* from sre_math.c
 */
extern float  Gaussrandom(float mean, float stddev);
extern int    Linefit(double *x, double *y, int N, 
		      double *ret_a, double *ret_b, double *ret_r);
extern double Gammln(double xx);
extern int    DNorm(double *vec, int n);
extern int    FNorm(float *vec, int n);
extern void   DScale(double *vec, int n, double scale);
extern void   FScale(float *vec, int n, float scale);
extern void   DSet(double *vec, int n, double value);
extern void   FSet(float *vec, int n, float value);
extern double DSum(double *vec, int n);
extern float  FSum(float *vec, int n);
extern float  sre_random(void);
extern void   sre_srandom(int seed);

/* from sre_string.c
 */
#ifdef NOSTR
extern char *strstr(char *s, char *subs);
#endif
#ifdef NO_STRDUP
extern char *strdup(char *s);
#endif
extern int strinsert(char *s1, char c, int pos);
extern int strdelete(char *s1, int pos);
extern void s2lower(char *s);
extern void s2upper(char *s);
extern void *MallocOrDie(size_t size);
extern void *ReallocOrDie(void *p, size_t size);

/* from stack.c
 */
extern struct intstack_s *InitIntStack(void);
extern void PushIntStack(struct intstack_s *stack, int data);
extern int  PopIntStack(struct intstack_s  *stack, int *ret_data);
extern void ReverseIntStack(struct intstack_s *stack);
extern int  FreeIntStack( struct intstack_s *stack );

/* 
 * from translate.c
 */
extern char *Translate(char *seq, char **code);

/* 
 * from types.c
 */
extern int IsInt(char *s);
extern int IsReal(char *s);

/* 
 * from weights.c
 */
extern int SonnhammerWeights(char **aseq, int nseq, int alen, float **ret_weights);
extern int VoronoiWeights(char **aseq, int nseq, int alen, float **ret_weights);

/***********************************************************
 * Function declarations: (stuff I haven't ansified yet)
 ***********************************************************/

extern int			/**** Function declarations:  */
  seqcmp(),			/* compare two encoded sequences */
  seqncmp(),			/* compare n bases of two encoded seqs */
  seqencode(),			/* create an encoded seq from char string */
  coded_revcomp(),		/* create encoded revcomp of an encoded seq */
  seqdecode(),			/* decode encoded seq to a char string */
  seqndecode();			/* decode n bases of encoded seq to char string */


#endif /* SQFUNCSH_INCLUDED */
