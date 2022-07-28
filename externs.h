/* 
 * externs.h - declarations of all external functions in the package
 * SRE, Tue Mar 16 11:53:00 1993
 * V1.1 SRE, Thu Jul 15 18:55:46 1993
 */            
/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/


#include "states.h"
#include "squid.h"

/* 
 * from align.c
 * Make multiple alignments from individual tracebacks
 */
extern int Traces2Alignment(char **rseqs, SQINFO *sqinfo, int nseq, int M, struct trace_s **tr,
			    int matchonly, char ***ret_aseqs, AINFO *ret_ainfo);

/* 
 * from dbviterbi.c
 * Alignment algorithm: scanning version (hmmls)
 */
extern int DBViterbi(struct shmm_s *shmm, char *seq, int L, int window, 
		     float thresh, int (*gotone_f)(struct vit_s **,int,int,int,int));


/* 
 * from emit.c
 * Emit sequences from a model.
 */
extern int EmitSequence(struct hmm_struc *hmm, char **ret_seq, 
			struct trace_s **ret_tr);
extern int EmitBestSequence(struct hmm_struc *hmm, char **ret_seq, 
			    struct trace_s **ret_tr);
extern int RD_EmitBestSequence(struct hmm_struc *hmm, float *randomseq, char **ret_seq, 
			       struct trace_s **ret_tr);
extern double HMMInfoContent(struct hmm_struc *hmm, struct shmm_s *shmm);


/* from forback.c
 * Forward-backward algorithm
 */
extern void  Forward(struct hmm_struc *hmm, char *s, struct forback_s ***ret_fwd, float **ret_scl);
extern void  Backward(struct hmm_struc *hmm, char *s, float *scl, struct forback_s ***ret_bck);
extern void  AlignmentConfidence(struct hmm_struc *hmm, int L, 
				 struct forback_s **fwd, struct forback_s **bck, float *scl,
				 struct forback_s ***ret_conf);
extern void  TraceConfidence(struct forback_s **cmx, int L, struct trace_s *tr, float **ret_conf);
extern void  ForbackCount(struct hmm_struc *hmm, char *seq, int L, float weight,
			  struct forback_s **fwd, struct forback_s **bck,
			  float *scl, struct hmm_struc *count);
extern float ForwardScore(struct forback_s **fwd, int L, int M, float *scl);
extern float BackwardScore(struct forback_s **bck, int L, float *scl);
extern float RandomScore(float *randomseq, char *seq);
extern float ForbackSymscore(char x, float *scores, int hyperbayes);
extern void  DumpForbackMatrix(struct forback_s **mx, int L, int M, float *scalefactors);


/* from fragviterbi.c
 * alignment algorithm: multiple-hit Smith/Waterman scanning version (hmmfs)
 */
extern int FragViterbi(struct shmm_s *shmm, char *seq, int L, int singlehit,
		       float P1, float P2, float P3, float thresh,
		       int (*gotone_f)(struct shmm_s *, char *, int,int,int,float));

/* from hmmio.c
 * Input/output (saving/reading) of models
 */
extern int WriteHMM(FILE *fp, struct hmm_struc *hmm);
extern int WriteBinaryHMM(FILE *fp, struct hmm_struc *hmm);
extern struct hmm_struc *ReadHMM(FILE *fp);

/* 
 * from maxmodelmaker.c
 * Construction of models from multiple alignments
 */
extern int Maxmodelmaker(char **aseqs, AINFO *ainfo, int nseq,
			 struct prior_s *prior, float *randomseq, float mpri, 
			 struct hmm_struc **ret_hmm, struct trace_s ***tr);


/* 
 * from misc.c
 * Miscellaneous functions with no home
 */
extern int   DetermineAlphabet(char **rseqs, int nseq);
extern void  BlockRaggedEdgedAlignment(char **aseqs, int nseq, int alen);
extern int   DChoose(double *probs, int N);
extern int   FChoose(float *probs, int N);
extern int   AlignmentTooBig(int L, int M);

/* from output.c
 * "User-friendly" output
 */
extern void PrintFancyTrace(FILE *ofp, struct shmm_s *shmm, 
			    struct trace_s *tr,
			    char *seq, char *seqname, int from_pos);

/* from pam_prior.c
 * Experimental: HMM evolution using PAM-like substitution matrices
 */
extern void EvolveHMM(struct hmm_struc *hmm, double (*bmx)[20][20], int num, 
		      double scale, struct hmm_struc **ret_hmm);
extern int  ParseSubstitutionMatrix(FILE *fp, double smx[20][20]);
extern void LogifySubstitutionMatrix(double smx[20][20]);
extern void ProbifySubstitutionMatrix(double smx[20][20]);
extern void PAM2SubstitutionMatrix(int **pam, double scale, 
				   double *pri, double smx[20][20]);
extern void RootSequenceDistribution(double smx[20][20], double *pri, 
				     double *ct, double *rd);
extern void MatrixLikelihood(double smx[20][20], double *pri, double *ct, double *ll);
extern void BestMatrixDistance(double bmx[20][20], double *pri, double *ct, 
		   int *ret_dist, double *ret_ll);
extern void EvolveCounts(double *ct, double smx[20][20], double *nct);
extern void ExpSubstitutionMatrix(double bmx[20][20], int dist, double smx[20][20]);
extern void PrintSubstitutionMatrix(double smx[20][20]);

/* 
 * from prior.c
 * Dirichlet priors
 */
extern void DefaultSimplePrior(struct prior_s **ret_prior);
extern void ReadPrior(char *pfile, struct prior_s **ret_prior);
extern void ToPAMPrior(int **pam, float scale, float wt, struct prior_s *prior);
extern void PrintPAMPrior(struct prior_s *pri);
extern struct prior_s *AllocPrior(void);
extern void FreePrior(struct prior_s *prior);
extern void DefaultRandomModel(float *randomseq);
extern void ReadRandomModel(char *fname, float *randomseq);
extern void PriorifyMatchVector(float *vec, struct prior_s *prior);
extern void PriorifyInsertVector(float *vec, struct prior_s *prior);
extern void PriorifyTransitionVectors(float *tm, float *ti, float *td, struct prior_s *prior);
extern void PriorifyHMM(struct hmm_struc *hmm, struct prior_s *prior);
extern void StructurePerceptron(struct prior_s *prior, float *xray);
extern void AnnotateAlignment(char **aseq, int nseq, AINFO *ainfo, float **ret_inputs);

/* from profiles.c
 * Emulation of GCG PROFILE package by Michael Gribskov
 */
extern int WriteProfile(FILE *fp, struct hmm_struc *hmm, struct shmm_s *shmm, 
			float *randomseq);


/* 
 * from saviterbi.c
 * Alignment algorithm: simulated annealing version
 */
extern int    SaFill(struct sa_hmm_s *sahmm, char *s, struct sa_s ***ret_mx);
extern int    SaTrace(struct sa_s **mx, int L, struct sa_hmm_s *sahmm, 
		      struct trace_s **ret_tr);
extern void   DumpSaMatrix(struct sa_s **mx, int L, int M, double *scalefactors);
extern struct sa_hmm_s *CreateSahmm(struct hmm_struc *hmm, float kT);
extern void   DestroySahmm(struct sa_hmm_s *sahmm);
extern double SaSymscore(char x, double *scores, int hyperbayes);


/* 
 * from scorestack.c and new_scorestack.c
 * Overlap filtering routines used by hmmsw and hmmls
 */
extern int ReportScanHit(int left, int right, double score, int (*print_hit)(int,int,double));
extern int Simple_ReportScanHit(int left, int right, double score, int (*print_hit)(int,int,double));

/* 
 * from states.c
 * Support for the basic data structures
 */
extern struct hmm_struc *AllocHMM(int M);
extern int  WriteFlatPriorHMM(struct hmm_struc *hmm, struct prior_s *prior);
extern struct hmm_struc *HMMDup(struct hmm_struc *hmm);
extern int  FreeHMM(struct hmm_struc *hmm);
extern int  CountSymbol(char sym, double wt, float *counters);
extern float HMMDistance(struct hmm_struc *newhmm, struct hmm_struc *oldhmm);
extern void VerifyHMM(struct hmm_struc *hmm);
extern void Renormalize(struct hmm_struc *hmm);
extern void HybridizeHMMs(struct hmm_struc *newhmm, struct hmm_struc *oldhmm, 
			  double damp_factor);


/* 
 * from swviterbi.c
 * Alignment algorithm: single-hit Smith/Waterman version
 */
extern int SWViterbi(struct shmm_s *shmm, char *seq, int L, 
		     float P1, float P2, float P3, int fullseq,
		     int *ret_i, int *ret_j, int *ret_kstart, int *ret_kend, 
		     float *ret_score, struct trace_s **ret_tr);


/* from trace.c
 * Support for the alignment traceback data structure
 */
extern void AllocTrace(int tlen, struct trace_s **ret_tr);
extern void ReallocTrace(struct trace_s *tr, int tlen);
extern void FreeTrace(struct trace_s *tr);
extern void ReverseTrace(struct trace_s *tr, int tlen);
extern void PrintTrace(struct trace_s *tr);
extern void TraceCount(struct hmm_struc *hmm, char *seq, float wt, struct trace_s *tr);
extern int  TraceScore(struct shmm_s *shmm, char *seq, struct trace_s *tr,
		       float *ret_score);
extern void DealignTrace(struct trace_s *tr, char *aseq, int alen);


/* from viterbi.c
 * Alignment algorithm: global alignment (hmms, hmma, hmmt)
 */
extern int   Symscore(char x, float *scores, float *priors);
extern void  MakeSearchHMM(struct hmm_struc *hmm, float *randomseq, struct shmm_s *shmm);
extern void  PrintSearchHMM(FILE *fp, struct shmm_s *shmm);
extern struct shmm_s *AllocSearchHMM(int M);
extern void  FreeSearchHMM(struct shmm_s *shmm);
extern void  PrepareSequence(char *s, char **ret_seq, int *ret_L);
extern int   ViterbiFill(struct shmm_s *shmm, char *seq, int L, 
			 struct vit_s ***ret_mx, float *ret_score);
extern int   ViterbiTrace(struct vit_s **mx, struct shmm_s *shmm, char *seq, int window,
			  int end_i, int end_k, 
			  struct trace_s **ret_tr, int *ret_i, int *ret_k);
extern void  FreeViterbiMatrix(struct vit_s **mx, int L);	
extern void  PrintViterbiMatrix(struct vit_s **mx, char *seq1, int L, int M);
extern void  ViterbiAlignAlignment(struct shmm_s *shmm, char **aseqs, int alen, int nseq,
				   struct trace_s ***ret_tr, float *ret_sc);
extern void PrintFragViterbiMatrix(struct fvit_s **mx, int L, int M);



/* from weeviterbi.c
 * Linear memory alignment for Needleman/Wunsch and full seq Smith Waterman
 */
extern void WeeViterbi(struct shmm_s *shmm, char *seq, int L, 
		       int smithwaterman, float P2, float P3,
		       struct trace_s **ret_tr, float *ret_sc);



