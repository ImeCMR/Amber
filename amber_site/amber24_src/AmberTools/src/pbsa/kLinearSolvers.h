//
// CUDA kernel header used by linear PB solvers
// Ruxi Qi @ UC Irvine, 2017-2018
// 
#ifndef KLINEARSOLVERS_H_
#define KLINEARSOLVERS_H_

#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime.h>
#endif

__global__ void init_vector_kernel(float *vec, int m);
__global__ void copy_vector_kernel(float *vec, float *vec_f, int m);
__global__ void inv_vector_kernel(float *vec, float *inv, int m);
__host__ __device__ float r_map_exp_x(float *epsxmp, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float r_map_exp_y(float *epsymp, int i2, int j2, int k2, int xn, int yn); 
__host__ __device__ float r_map_exp_z(float *epszmp, int i2, int j2, int k2, int xn, int yn);
__host__ __device__ float hmav(float a, float b);
__global__ void restrict_eps_map_kernel(float *epsxf, float *epsyf, float *epszf, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr); 
__global__ void feedepsintoam_kernel(int lxm, int lym, int lzm, float *am1, float *am2, float *am3, float *eps1, float *eps2, float *eps3);
__global__ void set_am_ad_kernel_head(float *epsx, float *epsy, float *epsz, float *lam1, float *lam2, float *lam3, int xnynzn);
__global__ void set_am_ad_kernel_body(float *lam1, float *lam2, float *lam3, float *lad, float *lbz, float *iv, int xn, int yn, int zn, float lfactor, float epsout);
__global__ void set_am_ad_kernel_tail(float *lam1, float *lam2, float *lam3, int xn, int yn, int zn);
__global__ void restrict_v_kernel(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr);
__host__ __device__ int f_id(int i, int j, int k, int nx, int ny);
__global__ void interpolate_kernel_head(int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float epsout);
__global__ void interpolate_kernel_body(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz);
__device__ void ipl_chain_d(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn); 
__device__ void ipl_chain2_d(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn);
__host__ __device__ float ipl_comp1(float v, int l, float *lbz, float *am_1, int xnyn, int xnynzn, int shift_1);
__host__ __device__ float ipl_comp2(float v, int l, float *lbz, float *am_1, float *am_2, int xnyn, int xnynzn, int shift_1, int shift_2);
__host__ __device__ float ipl_comp3(float v, int l, float *lbz, float *am_1, float *am_2, float *am_3, int xnyn, int xnynzn, int shift_1, int shift_2, int shift_3);
__global__ void solver_black_kernel(float *phi, float *epsi, float *epsj, float *epsk,float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void solver_red_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm);
__global__ void residue_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *epsc, float *rho, int xm, int ym, int zm, float* res);

#endif //KLINEARSOLVERS_H_
