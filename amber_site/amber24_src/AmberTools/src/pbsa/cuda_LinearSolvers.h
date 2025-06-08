//
// CUDA function header used by linear PB solvers
// Ruxi Qi @ UC Irvine, 2017-2018
//
#ifndef CUDA_LINEARSOLVERS_H_
#define CUDA_LINEARSOLVERS_H_

#if defined(AMBER_PLATFORM_AMD)
#include <hip/hip_runtime.h>
#endif

extern "C" void init_param_c_(int *nx, int *ny, int *nz, int *p_maxitn, int *p_bcopt, float *p_accept, float *p_pbkappa, float *p_epsout, float *p_h, float *p_wsor);
extern "C" void allocate_array_cuda_(int *solvopt);
extern "C" void deallocate_array_cuda_();
extern "C" void init_array_cuda_(int *solvopt, float *epsx, float *epsy, float *epsz, float *p_bv, float *p_iv, float *p_xs);
extern "C" void pb_mg_cuda_(float *phi_f, float *xs_f);
__host__ void init_vector(float *vec, int m);
__host__ void restrict_eps_map(float *epsx, float *epsy, float *epsz, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr);
__host__ void set_am_ad(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn, int yn, int zn, float lfactor, float epsout);
__host__ void restrict_v(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__host__ void restrict_cuda(int level, float *vf, float *vt, int *index, int coef, float divider);
__host__ void interpolate(int level);
__host__ void ipl_chain_h(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__host__ void ipl_chain2_h(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__host__ void interpolate_cuda(int level);
__host__ void relax(int level, int ncyc);
__host__ void relax_cuda(int level, int ncyc);
__host__ void cg_cuda(int level, int ncyc);
__host__ void jacobi_dia_solver(float *x, float *values, int *dia_offsets, float *rhs, int bwidth, int xm, int ym, int zm, int maxitn, float acpt, float *lnorm);
__host__ void FCycle(int level);
__host__ void VCycle(int level);
#endif // CUDA_LINEARSOLVERS_H_
