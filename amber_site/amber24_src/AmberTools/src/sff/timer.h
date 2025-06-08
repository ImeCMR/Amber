#ifndef TIMER_H
#define TIMER_H

extern REAL_T *tinit;
/* MD timers */
extern REAL_T *tmd, *tmdOther, *tmdIO;
/* NMODE timers */
extern REAL_T *tnmode, *tnmodeOther, *tnmodeFact, *tnmodeEigenvec, 
  *tnmodeInvdiag, *tnmodeAmat, *tnmodeSort, *tnmodeHessian, 
  *tnmodeDSYEV, *tnmodeDSYEVD, *tnmodeDiagA, *tnmodeNorm, 
  *tnmodeLan;
/* conjgrad timers */
extern REAL_T *tconjgrad, *tconjgradMME, *tconjgradOther;
/* newton timers */
extern REAL_T *tnewton, *tnewtonLevel, *tnewtonCholesky, *tnewtonMME, 
  *tnewtonMME2, *tnewtonOther, *tnewtonDSYEV;
/* lmod timers */ /* lmod_opt.lmod_time lmod_opt.aux_time see sff.h */

/* xmin timers */ /* xmin_opt.xmin_time see sff.h */

/* mme timers */
extern REAL_T *tmme, *tmmeCons, *tmmeNonb, *tmmePair, *tmmeBond, *tmmeAngl,
  *tmmePhi, *tmmeBorn, *tmmePB, *tmmeRism, *tmmeOther;
/* mme2 timers */
extern REAL_T *tmme2Cons, *tmme2Nonb, *tmme2Pair, *tmme2Bond, *tmme2Angl,
  *tmme2Phi, *tmme2Born, *tmme2, *tmme2Other;
/* GB timers */
extern REAL_T *tgb2, *tgb2dgemm1, *tgb2dgemm2, *tgb2dgemm3, *tgb2Other;
/* MPI timers */
extern REAL_T *treduceegb, *treducemme;

void init_timers(void);
int mme_timer(void);

#endif
