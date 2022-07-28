/************************************************************
 * HMMER - Biological sequence analysis with HMMs
 * Copyright 1992-1995 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

#ifndef STATESH_INCLUDED
#define STATESH_INCLUDED
/* states.h
 * 
 * in which we declare various kinds of HMM states, and a
 * structure to organize a whole model in.
 */


/* RAMLIMIT determines the point at which we switch from fast,
 * full dynamic programming to slow, linear-memory divide and conquer
 * dynamic programming algorithms. It is the minimum amount of available
 * RAM on the systems the package will run on. It can be overridden
 * from the Makefile. By default, we assume we have 32 Mb RAM available.
 */
#ifndef RAMLIMIT
#define RAMLIMIT 32
#endif

/* an idiom for determining a symbol's position in the array
 * by pointer arithmetic.
 * does no error checking, so caller must already be damned sure x is
 * valid in the alphabet!
 */
extern char  Alphabet[];        /* ACGT, for instance            */ 
extern int   Alphabet_type;	/* kDNA, for instance, or kAmino */
extern int   Alphabet_size;	/* 4 or 20                       */

#define SYMIDX(x)  (strchr(Alphabet, (x)) - Alphabet)
#define LOG2(x)   ((x) > 0 ? log(x) * 1.44269504 : -999.)
#define EXP2(x)   (exp((x) * 0.69314718 )) 
#define SQR(x)    ((x) * (x))
#define INTSCALE  1000.0
#define MAXABET   20		/* maximum size of alphabet, used for some arrays  */
#define MAXDCHLET 20		/* maximim # Dirichlet components in mixture prior */
#define NINPUTS   4		/* number of inputs into structural prior */


/* We define a "basic" state, which covers the basic match, insert, and
 * delete states from the Haussler paper. Numbers are stored as
 * pre-calculated negative logs.
 */
struct basic_state {
  float t[3];			/* state transitions to +1 M, +0 I, +1 D */
  float p[MAXABET];            	/* symbol emission probabilities         */
};

/* A complete hidden Markov model
 */
struct hmm_struc {
  int    M;			/* length of the model                */
  struct basic_state *ins;      /* insert states 0..M+1               */
  struct basic_state *mat;      /* match 0..M+1; 0 = BEGIN, M+1 = END */
  struct basic_state *del;      /* delete 0..M+1                      */

  /* Optional annotation on the HMM, taken from alignment
   */
  char  *ref;			/* reference coords and annotation        */
  char  *cs;                    /* consensus structure annotation         */         
  int    flags;			/* flags for what optional info is in HMM */

  /* Structural annotation: xray[0..M+1][NINPUTS], indexed manually 
   */
  float *xray;	
};

/* Flags for optional info in an HMM structure
 */
#define HMM_REF   (1<<0) 
#define HMM_CS    (1<<1)
#define HMM_XRAY  (1<<2)

/* Array indices for structural inputs in HMM 
 */
#define XRAY_bias  0		/* bias: always 1 */
#define XRAY_E     1		/* P(sheet), 0..1 */
#define XRAY_H     2		/* P(helix), 0..1 */
#define XRAY_SA    3		/* relative side chain solv. access. (0..15) */


/* We have to increase the dynamic range of the calculations to doubles
 * when we do simulated annealing.
 * Rather than sacrificing speed in the rest of the code, we define
 * new HMM structures specifically for simulated annealing use. 
 */
struct sa_state_s {
  double t[3];			/* transition probability to M, I, D */
  double p[MAXABET];		/* symbol emission probabilities     */
};

struct sa_hmm_s {
  int    M;			/* length of the model            */
  struct sa_state_s *ins;       /* insert states 0..M+1           */
  struct sa_state_s *mat;       /* match 0..M+1; 0 = impossible   */
  struct sa_state_s *del;       /* delete 0..M+1; 0 = BEGIN state */
}; 


/* Search HMM structure:
 * We use a special structure to store HMM's in integer log odds
 * form, for use in the various Viterbi alignment procedures.
 * The format is very compact and a bit fragile because this
 * is optimized for speed & memory access patterns.
 */
struct shmm_s {
  int  M;                       /* length of the model           */
  int *m_emit[26];		/* 26 x M+1 emission scores      */
  int *i_emit[26];		/* 26 x M+1 insertion scores     */
  int *t;                       /* 9 x M+1 state transitions: 
				   order is dd,di,dm,id,ii,im,md,mi,mm */
  /* plus extra annotation:
   */
  int   flags;
  char *ref;
  char *cs;
};
/* Order of transition probabilities in shmm_s
 */
#define Tdd 0
#define Tdi 1
#define Tdm 2
#define Tid 3
#define Tii 4
#define Tim 5
#define Tmd 6
#define Tmi 7
#define Tmm 8


/* Prior information: expectations for state transition probabilities,
 * emission probabilities, and alphas (regularizers) describing the 
 * strengths of belief in the prior. alphas of 1 result in Laplace
 * small sample size corrections; very high alphas result in "hard-wiring"
 * of probabilities to the prior expectation.
 */
struct prior_s {
  int   strategy;		/* PRI_SIMPLE, etc.                          */

  int   tnum;			/* number of transition Dirichlet mixtures   */
  float tw[MAXDCHLET][NINPUTS];	/* weight matrix for struct prior perceptron */
  float tq[MAXDCHLET];		/* probabilities of tnum components          */
  float tm[MAXDCHLET][3];	/* mat transition terms per mix component    */
  float td[MAXDCHLET][3];	/* del transition terms per mix component    */
  float ti[MAXDCHLET][3];	/* ins transition terms per mix component    */

  int   mnum;			/* number of mat emission Dirichlet mixtures */
  float mw[MAXDCHLET][NINPUTS];	/* weight matrix for struct prior perceptron */
  float mq[MAXDCHLET];		/* probabilities of mnum components          */
  float mat[MAXDCHLET][MAXABET];/* match emission terms per mix component    */

  int   inum;			/* number of ins emission Dirichlet mixtures */
  float iw[MAXDCHLET][NINPUTS];	/* weight matrix for struct prior perceptron */
  float iq[MAXDCHLET];		/* probabilities of inum components          */
  float ins[MAXDCHLET][MAXABET];/* insert emission m_x per mix component     */
};

#define PRI_SIMPLE 0		/* simple single-component Dirichlets (orig) */
#define PRI_PAM    1		/* ad hoc PAM-based mixture prior            */
#define PRI_MIX    2		/* mixture prior, a la Brown/Haussler        */
#define PRI_STRUCT 3		/* structure annotation based priors         */



/* Traceback structure for alignments of model to sequence.
 * Each array in a trace_s is 0..tlen-1.
 * Element 0 and tlen-1 are dummy "alignments" to BEGIN and END.
 */
struct trace_s {
  int   tlen;			/* length of traceback                       */
  int  *nodeidx;                /* index of aligned node, 0..hmm->M+1        */
  char *statetype;              /* state type used for alignment             */
  int  *rpos;                   /* position in raw sequence, 0..slen-1 or -1 */ 
};


/* Bookkeeping for individual cells of the dynamic programming calculations
 */
				/* used by most Viterbi procedures */
struct vit_s {
  int  score_m;			/* score to match here  */
  int  score_d;			/* score to delete here */
  int  score_i;			/* score to insert here */
};

				/* used by simulated annealing */
struct sa_s {
  double score_m;
  double score_d;
  double score_i;
};

				/* used by forward-backwards */
struct forback_s {
  float score_m;
  float score_d;
  float score_i;
};

				/* used in fragviterbi.c (hmmfs) */
struct fvit_s {
  int  score_m;			/* score to match here   */
  int  score_d;			/* score to delete here  */
  int  score_i;			/* score to insert here  */
  int  tback_m;			/* traceback from match  */
  int  tback_d;			/* traceback from delete */
  int  tback_i;			/* traceback from insert */
};


/* some define's for storing backtracking info
 */
#define FROM_NOWHERE 127
#define FROM_MATCH   0
#define FROM_INSERT  1
#define FROM_DELETE  2

/* same numbers, different names -- just for readability, makes
 * more semantic sense after the tracebacks
 */
#define MATCH  FROM_MATCH
#define INSERT FROM_INSERT
#define DELETE FROM_DELETE

#define BEGIN  MATCH
#define END    MATCH

#endif /* STATESH_INCLUDED */
