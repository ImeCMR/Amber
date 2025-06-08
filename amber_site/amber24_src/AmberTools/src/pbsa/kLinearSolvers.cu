//
// CUDA kernel routines used by linear PB solvers
// Ruxi Qi @ UC Irvine, 2017-2018
// 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <sys/time.h>
#include <signal.h> // For gdb
#include "helper_cuda.h" // For error handling and device pickup

#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime_api.h>
#  include "hip_definitions.h"
#endif

#include "kLinearSolvers.h"

__global__
void init_vector_kernel(float *vec, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) vec[i] = 0.0;
}

__global__
void copy_vector_kernel(float *vec, float *vec_f, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) vec[i] = vec_f[i];
}

__global__
void inv_vector_kernel(float *vec, float *inv, int m) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) inv[i] = 1.0/vec[i];
}

__host__ __device__
float r_map_exp_x(float *epsxmp, int i2, int j2, int k2, int xn, int yn) {
    float exp_x = hmav(epsxmp[f_id(i2  , j2  , k2  , xn, yn)], epsxmp[f_id(i2+1, j2  , k2  , xn, yn)]) / 4.0 +
                 (hmav(epsxmp[f_id(i2  , j2-1, k2  , xn, yn)], epsxmp[f_id(i2+1, j2-1, k2  , xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2  , xn, yn)], epsxmp[f_id(i2+1, j2+1, k2  , xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2  , k2-1, xn, yn)], epsxmp[f_id(i2+1, j2  , k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2  , k2+1, xn, yn)], epsxmp[f_id(i2+1, j2  , k2+1, xn, yn)])) / 8.0 +
                 (hmav(epsxmp[f_id(i2  , j2-1, k2-1, xn, yn)], epsxmp[f_id(i2+1, j2-1, k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2-1, xn, yn)], epsxmp[f_id(i2+1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2-1, k2+1, xn, yn)], epsxmp[f_id(i2+1, j2-1, k2+1, xn, yn)]) +
                  hmav(epsxmp[f_id(i2  , j2+1, k2+1, xn, yn)], epsxmp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;
    return exp_x;
}

__host__ __device__
float r_map_exp_y(float *epsymp, int i2, int j2, int k2, int xn, int yn) {
    float exp_y = hmav(epsymp[f_id(i2  , j2  , k2  , xn, yn)], epsymp[f_id(i2  , j2+1, k2  , xn, yn)]) / 4.0 +
                 (hmav(epsymp[f_id(i2-1, j2  , k2  , xn, yn)], epsymp[f_id(i2-1, j2+1, k2  , xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2  , xn, yn)], epsymp[f_id(i2+1, j2+1, k2  , xn, yn)]) +
                  hmav(epsymp[f_id(i2  , j2  , k2-1, xn, yn)], epsymp[f_id(i2  , j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2  , j2  , k2+1, xn, yn)], epsymp[f_id(i2  , j2+1, k2+1, xn, yn)])) / 8.0 +
                 (hmav(epsymp[f_id(i2-1, j2  , k2-1, xn, yn)], epsymp[f_id(i2-1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2-1, xn, yn)], epsymp[f_id(i2+1, j2+1, k2-1, xn, yn)]) +
                  hmav(epsymp[f_id(i2-1, j2  , k2+1, xn, yn)], epsymp[f_id(i2-1, j2+1, k2+1, xn, yn)]) +
                  hmav(epsymp[f_id(i2+1, j2  , k2+1, xn, yn)], epsymp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;
    return exp_y;
}

__host__ __device__
float r_map_exp_z(float *epszmp, int i2, int j2, int k2, int xn, int yn) {
    float exp_z = hmav(epszmp[f_id(i2  , j2  , k2  , xn, yn)], epszmp[f_id(i2  , j2  , k2+1, xn, yn)]) / 4.0 +
                 (hmav(epszmp[f_id(i2  , j2-1, k2  , xn, yn)], epszmp[f_id(i2  , j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2  , j2+1, k2  , xn, yn)], epszmp[f_id(i2  , j2+1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2-1, j2  , k2  , xn, yn)], epszmp[f_id(i2-1, j2  , k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2  , k2  , xn, yn)], epszmp[f_id(i2+1, j2  , k2+1, xn, yn)])) / 8.0 +
                 (hmav(epszmp[f_id(i2-1, j2-1, k2  , xn, yn)], epszmp[f_id(i2-1, j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2-1, j2+1, k2  , xn, yn)], epszmp[f_id(i2-1, j2+1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2-1, k2  , xn, yn)], epszmp[f_id(i2+1, j2-1, k2+1, xn, yn)]) +
                  hmav(epszmp[f_id(i2+1, j2+1, k2  , xn, yn)], epszmp[f_id(i2+1, j2+1, k2+1, xn, yn)])) / 16.0;
    return exp_z;
}

__host__ __device__
float hmav(float a, float b) {
    return 2.0 * a * b / (a + b);
}

__global__
void restrict_eps_map_kernel(float *epsxf, float *epsyf, float *epszf, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int i2 = 2 * i + 1;
    int j2 = 2 * j + 1;
    int k2 = 2 * k + 1;
    if (i < xnr  && j < ynr  && k < znr){
        int flatid = i + xnr * j + xnr * ynr * k;
        epsxr[flatid] = r_map_exp_x(epsxf, i2, j2, k2, xn, yn);
        epsyr[flatid] = r_map_exp_y(epsyf, i2, j2, k2, xn, yn);
        epszr[flatid] = r_map_exp_z(epszf, i2, j2, k2, xn, yn);
    }
}

__global__
void feedepsintoam_kernel(int lxm, int lym, int lzm, float *am1, float *am2, float *am3, float *eps1, float *eps2, float *eps3) {
    // eps1 has (0:xm, ym, zm) dimension (0:ym/0:zm for eps2/eps3), which contains am1 (xm, ym, zm),
    // so need to recheck the passed array index.
    // Edit: Solved, using 3-D mapping
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < lxm && j < lym && k < lzm) {
        int flatid = i + lxm * j + lxm * lym * k;
        // Warning of bug, x/y dimensions are incorrect, that's why only lepsz is good
        //int flatid_x_plus1 = (i + 1) + lxm * j + lxm * lym * k;
        //int flatid_y_plus1 = i + lxm * (j + 1) + lxm * lym * k;
        int flatid_x_plus1 = (i + 1) + (lxm + 1) * j + (lxm + 1) * lym * k;
        int flatid_y_plus1 = i + lxm * (j + 1) + lxm * (lym + 1) * k;
        int flatid_z_plus1 = i + lxm * j + lxm * lym * (k + 1);

        // eps* causes CUDA Exception 14 error. Solved
        am1[flatid] = eps1[flatid_x_plus1];
        am2[flatid] = eps2[flatid_y_plus1];
        am3[flatid] = eps3[flatid_z_plus1];
    }
}

__global__
void set_am_ad_kernel_head(float *epsx, float *epsy, float *epsz, float *lam1, float *lam2, float *lam3, int xnynzn) {
    // 1-D, xnynzn
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < xnynzn) {
        lam1[i] = epsx[i];
        lam2[i] = epsy[i];
        lam3[i] = epsz[i];
    }
}

__global__
void set_am_ad_kernel_body(float *lam1, float *lam2, float *lam3, float *lad, float *lbz, float *iv, int xn, int yn, int zn, float lfactor, float epsout) {
    // 3-D, (xn, yn, zn)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < xn && j < yn && k < zn) {
        int flatid = i + xn * j + xn * yn * k;
        int flatid_x_minus1 = (i - 1) + xn * j + xn * yn * k;
        int flatid_y_minus1 = i + xn * (j - 1) + xn * yn * k;
        int flatid_z_minus1 = i + xn * j + xn * yn * (k - 1);

        lad[flatid] = lam1[flatid] + lam2[flatid] + lam3[flatid];
        if (i == 0) lad[flatid] += epsout; else lad[flatid] += lam1[flatid_x_minus1];
        if (j == 0) lad[flatid] += epsout; else lad[flatid] += lam2[flatid_y_minus1];
        if (k == 0) lad[flatid] += epsout; else lad[flatid] += lam3[flatid_z_minus1];

        lbz[flatid] = lfactor * iv[flatid];
        lad[flatid] += lbz[flatid];
    }
}

__global__
void set_am_ad_kernel_tail(float *lam1, float *lam2, float *lam3, int xn, int yn, int zn) {
    // 2-D, (i, j) -> (max_xy, max_yz)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < yn && j < zn) {
        lam1[xn - 1 + xn * i + xn * yn * j] = 0.0;
    }

    if (i < xn && j < zn) {
        lam2[i + xn * (yn - 1) + xn * yn * j] = 0.0;
    }

    if (i < xn && j < yn) {
        // Bug here!!! Detected from DIA padding in Jacobi
        //lam2[i + xn * j + xn * yn * (zn - 1)] = 0.0;
        lam3[i + xn * j + xn * yn * (zn - 1)] = 0.0;
    }
}

__global__
void restrict_v_kernel(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int i2 = 2 * i + 1;
    int j2 = 2 * j + 1;
    int k2 = 2 * k + 1;
    if (i < nxr  && j < nyr  && k < nzr){
        int b_id = i + nxr * j + nxr * nyr * k;

        bvr[b_id] = ( bvf[f_id(i2-1, j2-1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2-1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2  , k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2-1, nx, ny)] + bvf[f_id(i2+1, j2  , k2-1, nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2+1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2-1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2-1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2  , nx, ny)] + bvf[f_id(i2+1, j2-1, k2  , nx, ny)] ) +
                4 * ( bvf[f_id(i2-1, j2  , k2  , nx, ny)] + 2 * bvf[f_id(i2, j2  , k2  , nx, ny)] + bvf[f_id(i2+1, j2  , k2  , nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2+1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2  , nx, ny)] + bvf[f_id(i2+1, j2+1, k2  , nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2-1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2+1, nx, ny)] ) +
                2 * ( bvf[f_id(i2-1, j2  , k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2+1, nx, ny)] + bvf[f_id(i2+1, j2  , k2+1, nx, ny)] ) +
                    ( bvf[f_id(i2-1, j2+1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2+1, nx, ny)] );
        bvr[b_id] /= divider;
    }
}

__host__ __device__
int f_id(int i, int j, int k, int nx, int ny) {
    return i + nx * j + nx * ny * k;
}

__global__
void interpolate_kernel_head(int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float epsout) {
    int xniyni = xni * yni;
    // 2-D, [i, j] -> [max(xni, yni), max(yni, zni)]
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Padding offset <xniyni> included
    if (i < yni && j < zni) {
        lam1[-1 + (i + 1) * xni + (j + 1) * xniyni] = epsout;
    }

    if (i < xni && j < zni) {
        lam2[i - xni + (j + 2) * xniyni] = epsout;
    }

    if (i < xni && j < yni) {
        lam3[i + j * xni + xni * yni * zni] = epsout;
        // 0:xniyni-1
        lam1[i + j * xni] = epsout;
        lam2[i + j * xni] = epsout;
        lam3[i + j * xni] = epsout;
    }
}

__global__
void interpolate_kernel_body(float *v, int xn, int yn, int zn, float *vi, int xni, int yni, int zni, float *lam1, float *lam2, float *lam3, float *lbz) {
    int xniyni = xni * yni;
    int xniynizni = xni * yni * zni;

    // 3-D, (xn, yn, zn)
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < xn  && j < yn  && k < zn){
        int flatid_c = i + xn * j + xn * yn * k; // Coarse grid
        int flatid_f = (2 * i + 1) + (2 * j + 1) * xni + (2 * k + 1) * xniyni; // Fine grid

        vi[flatid_f] += v[flatid_c];

        // offset for offsetting the index
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      -1, lam2, xni, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      +1, lam2, xni, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    -xni, lam1,   1, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    +xni, lam1,   1, lam3, xniyni, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, -xniyni, lam2, xni, lam1,      1, xn, yn, zn);
        ipl_chain_d(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, +xniyni, lam2, xni, lam1,      1, xn, yn, zn);
    }

}

__device__
void ipl_chain_d(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    // Here v is value of coarse grid; vi is array of fine grid
    float v1 = ipl_comp1(v, l, lbz, am_1, xnyn, xnynzn, shift_1);
    int l1 = l + shift_1;
    atomicAdd(&vi[l1], v1);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2, -shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2,  shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3, -shift_3, am_2, shift_2, xn, yn, zn);
    ipl_chain2_d(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3,  shift_3, am_2, shift_2, xn, yn, zn);
}

__device__
void ipl_chain2_d(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    float v2 = ipl_comp2(v1, l1, lbz, am_1, am_2, xnyn, xnynzn, shift_1, shift_2);
    int l2 = l1 + shift_2;
    atomicAdd(&vi[l2], v2);
    float v3 = ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, -shift_3);
    atomicAdd(&vi[l2 - shift_3], v3);
    float v4 = ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, +shift_3);
    atomicAdd(&vi[l2 + shift_3], v4);
}

__host__ __device__
float ipl_comp1(float v, int l, float *lbz, float *am_1, int xnyn, int xnynzn, int shift_1) {
    float ipl_comp1_v;
    // For offsetting
    int bz_l = l;
    l += xnyn;

    if (shift_1 < 0)
        ipl_comp1_v = v * am_1[l + shift_1] / ( lbz[bz_l + shift_1] + am_1[l + 2 * shift_1] + am_1[l + shift_1] );
    else
        ipl_comp1_v = v * am_1[l] / ( lbz[bz_l + shift_1] + am_1[l] + am_1[l + shift_1] );
    return ipl_comp1_v;
}

__host__ __device__
float ipl_comp2(float v, int l, float *lbz, float *am_1, float *am_2, int xnyn, int xnynzn, int shift_1, int shift_2) {
    // For offsetting
    int bz_l = l;
    l += xnyn;

    float lad = am_1[l + shift_2] + am_1[l + shift_2 - abs(shift_1)] + lbz[bz_l + shift_2];
    float ipl_comp2_v;

    if (shift_2 < 0)
        ipl_comp2_v = v * am_2[l + shift_2] / ( am_2[l + 2 * shift_2] + am_2[l + shift_2] + lad);
    else
        ipl_comp2_v = v * am_2[l] / ( am_2[l] + am_2[l + shift_2] + lad );
    return ipl_comp2_v;
}

__host__ __device__
float ipl_comp3(float v, int l, float *lbz, float *am_1, float *am_2, float *am_3, int xnyn, int xnynzn, int shift_1, int shift_2, int shift_3) {
    // For offsetting
    int bz_l = l;
    l += xnyn;

    float lad = am_1[l + shift_3] + am_1[l + shift_3 - abs(shift_1)] + am_2[l + shift_3] + am_2[l + shift_3 - abs(shift_2)] + lbz[bz_l + shift_3];
    float ipl_comp3_v;

    if (shift_3 < 0)
        ipl_comp3_v = v * am_3[l + shift_3] / ( am_3[l + 2 * shift_3] + am_3[l + shift_3] + lad);
    else
        ipl_comp3_v = v * am_3[l] / ( am_3[l] + am_3[l + shift_3] + lad );
    return ipl_comp3_v;
}

__global__
void solver_black_kernel(float *phi, float *epsi, float *epsj, float *epsk,float *repsc, float *rho, float wsor, int xm, int ym, int zm) {

    int xmym = xm * ym;
    int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x) + xmym;

    if (i < xm * ym * zm + xmym) {
        phi[i] -= wsor * (phi[i] - (epsi[i - 1   ] * phi[i - 1   ] + epsi[i      ] * phi[i + 1   ] +
                                    epsj[i - xm  ] * phi[i - xm  ] + epsj[i      ] * phi[i + xm  ] +
                                    epsk[i - xmym] * phi[i - xmym] + epsk[i      ] * phi[i + xmym] + rho[i - xmym]) * repsc[i - xmym]);
    }

}

__global__
void solver_red_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *repsc, float *rho, float wsor, int xm, int ym, int zm) {

    int xmym = xm * ym;
    int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x) + 1 + xmym;

    if (i < xm * ym * zm + xmym) {
        phi[i] -= wsor * (phi[i] - (epsi[i - 1   ] * phi[i - 1   ] + epsi[i      ] * phi[i + 1   ] +
                                    epsj[i - xm  ] * phi[i - xm  ] + epsj[i      ] * phi[i + xm  ] +
                                    epsk[i - xmym] * phi[i - xmym] + epsk[i      ] * phi[i + xmym] + rho[i - xmym]) * repsc[i - xmym]);
    }

}

__global__
void residue_kernel(float *phi, float *epsi,float *epsj, float *epsk, float *epsc, float *rho, int xm, int ym, int zm, float* res) {

    int xmym = xm * ym;
    int i = (blockIdx.x * blockDim.x + threadIdx.x) + xmym;

    if (i < xm * ym * zm + xmym) {
        res[i - xmym] = epsi[i - 1   ] * phi[i - 1   ] + epsi[i       ] * phi[i + 1   ] +
                        epsj[i - xm  ] * phi[i - xm  ] + epsj[i       ] * phi[i + xm  ] +
                        epsk[i - xmym] * phi[i - xmym] + epsk[i       ] * phi[i + xmym] + rho[i - xmym] - epsc[i - xmym] * phi[i];
    }

}

