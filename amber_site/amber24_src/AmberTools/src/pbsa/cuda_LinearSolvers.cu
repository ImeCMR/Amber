//
// CUDA linear solvers with Unified Memory
// Ruxi Qi @ UC Irvine, 2017-2018
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <signal.h> // For gdb
#include <sys/time.h> // For timing
#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime.h>
#  include <hip/hip_runtime_api.h>
#  include <hipblas.h>
#  include "hip_definitions.h"
#else
#  include <cuda_runtime.h> // CUDA Runtime
#  include "cublas_v2.h" // CUDA BLAS Library
#endif
#include "helper_cuda.h" // For error handling and device pickup

// Jacobi-PCG with DIA format
#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/krylov/cg.h>
#include <cusp/precond/diagonal.h>
#include <thrust/device_ptr.h>

// local files
#include "cuda_LinearSolvers.h"
#include "kLinearSolvers.h"

// Global vairables
int l_xm;
int l_ym;
int l_zm;
int l_xmym;
int l_xmymzm;
int l_maxitn;
// For initialization
int l_l;
int l_m;
int l_n;
#if defined(AMBER_PLATFORM_AMD)
__device__ int l_bcopt;
#else
__device__ __managed__ int l_bcopt;
#endif
float l_accept;
//int mg_nlevel;
int ncyc_before;
int ncyc_after;
float l_pbkappa;
float l_epsout;
float l_h;
float l_wsor;
int l_itn;
float l_inorm;
float l_norm;

int threshold;

#define MG_NLEVEL 4
int mg_index[MG_NLEVEL + 1];
int mg_index_ext[MG_NLEVEL + 1];
int mg_x_idx[MG_NLEVEL + 1];
int mg_size[MG_NLEVEL][3];
float mg_onorm[MG_NLEVEL];

float *l_zv;
float *l_ad;
float *l_bv;
float *l_rv;
float *l_iv;
float *l_bz;
float *l_am1;
float *l_am2;
float *l_am3;
float *l_xv;

int devThreadsPerBlock = 8;

extern "C"
void init_param_c_(int *nx, int *ny, int *nz, int *p_maxitn, int *p_bcopt, float *p_accept, float *p_pbkappa, float *p_epsout, float *p_h, float *p_wsor) {
    l_xm = *nx;
    l_ym = *ny;
    l_zm = *nz;
    l_xmym = *nx * *ny;
    l_xmymzm = *nx * *ny * *nz;
    l_maxitn = *p_maxitn;
#ifdef AMBER_PLATFORM_AMD
    cudaMemcpyToSymbol(l_bcopt, p_bcopt, sizeof(int));
#else
    l_bcopt = *p_bcopt;
#endif
    l_accept = *p_accept;
    ncyc_before = 10;
    ncyc_after = 10;
    l_pbkappa = *p_pbkappa;
    l_epsout = *p_epsout;
    l_h = *p_h;
    l_wsor = *p_wsor;

    threshold = 4;
}

extern "C"
void allocate_array_cuda_(int *solvopt) {
    if (!(*solvopt == 2 || *solvopt == 4)) {
        printf("Error: Only MG/SOR is supported now.\n");
        exit(2);
    }

    int m, l, n;

    // set indices for the finest level for all solvers
    mg_index_ext[0] = 0;
    mg_index[0] = 0;
    mg_x_idx[0] = 0;
    mg_size[0][0] = l_xm;
    mg_size[0][1] = l_ym;
    mg_size[0][2] = l_zm;
    m = l_xmymzm;
    l = m + l_xmym;
    n = l + l_xmym;

    // set indices for all other levels for MG only
    if (*solvopt == 2) {
        for (int i = 1; i < MG_NLEVEL; i++) {
            mg_index_ext[i] = l;
            mg_index[i] = m;
            mg_x_idx[i] = n;

            //l_bcopt != 10 for now
            for (int j = 0; j < 3; j++) {
                mg_size[i][j] = mg_size[i - 1][j] / 2;
            }
            m += mg_size[i][0] * mg_size[i][1] * mg_size[i][2];
            l += mg_size[i][0] * mg_size[i][1] * mg_size[i][2] + mg_size[i][0] * mg_size[i][1];
            n += mg_size[i][0] * mg_size[i][1] * mg_size[i][2] + 2 * mg_size[i][0] * mg_size[i][1];
        }

        mg_index_ext[MG_NLEVEL] = l;
        mg_index[MG_NLEVEL] = m;
        mg_x_idx[MG_NLEVEL] = n;
    }

    // Now for all arrays
    // Try __managed__ declaration later for performance tuning
    //__device__ __managed__ l_xv[n];
    // Note in Fortran these arrays index from 1, not 1-xmym etc.
    cudaErrorCheck(cudaMallocManaged(&l_zv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_ad, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_bv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_rv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_iv, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_bz, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&l_am1, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_am2, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_am3, sizeof(float) * l));
    cudaErrorCheck(cudaMallocManaged(&l_xv, sizeof(float) * n));

    l_l = l;
    l_m = m;
    l_n = n;
}

extern "C"
void deallocate_array_cuda_() {
    cudaFree(l_zv);
    cudaFree(l_ad);
    cudaFree(l_bv);
    cudaFree(l_rv);
    cudaFree(l_iv);
    cudaFree(l_bz);
    cudaFree(l_am1);
    cudaFree(l_am2);
    cudaFree(l_am3);
    cudaFree(l_xv);
    cudaDeviceReset();
}

extern "C"
void init_array_cuda_(int *solvopt, float *epsx, float *epsy, float *epsz, float *p_bv, float *p_iv, float *p_xs) {
    if (!(*solvopt == 2 || *solvopt == 4) ) {
        printf("Error: Only MG/SOR is supported now.\n");
        exit(2);
    }

    // Initialize arrays l_ad, l_am*, l_*v to 0 on device. Use 1D thread block.
    int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;

    //!!!!!!!!!!Mark to replace!!!!!!!!!! page faults, do not initialize on levels >= threshold, where cpu version shoule be used
    // DO: split init into GPU, CPU two parts, using threshold as condition
    // -----------------------------

    cudaDeviceProp p;
    int deviceId;
    cudaGetDevice(&deviceId);
    cudaGetDeviceProperties(&p, deviceId);

    // HIP-TODO: Remove #if check when https://github.com/ROCm-Developer-Tools/HIP/issues/1696
    // is fixed
#if !defined(AMBER_PLATFORM_AMD)
    if (p.concurrentManagedAccess) {
        size_t const mbytes = l_m * sizeof(float);
        size_t const lbytes = l_l * sizeof(float);
        size_t const nbytes = l_n * sizeof(float);

        cudaMemPrefetchAsync(l_zv, mbytes, deviceId);
        cudaMemPrefetchAsync(l_ad, mbytes, deviceId);
        cudaMemPrefetchAsync(l_bv, mbytes, deviceId);
        cudaMemPrefetchAsync(l_rv, mbytes, deviceId);
        cudaMemPrefetchAsync(l_iv, mbytes, deviceId);

        cudaMemPrefetchAsync(l_am1, lbytes, deviceId);
        cudaMemPrefetchAsync(l_am2, lbytes, deviceId);
        cudaMemPrefetchAsync(l_am3, lbytes, deviceId);

        cudaMemPrefetchAsync(l_xv, nbytes, deviceId);
    }
#endif


    // -----------------------------
    /** All Levels **/
    // m
    int nblocks = (l_m - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_zv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_ad, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_bv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_rv, l_m);
    init_vector_kernel<<<nblocks, blocksize>>>(l_iv, l_m);
    cudaLaunchErrorCheck();

    // l
    nblocks = (l_l - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_am1, l_l);
    init_vector_kernel<<<nblocks, blocksize>>>(l_am2, l_l);
    init_vector_kernel<<<nblocks, blocksize>>>(l_am3, l_l);
    cudaLaunchErrorCheck();

    // n
    nblocks = (l_n - 1) / blocksize + 1;
    init_vector_kernel<<<nblocks, blocksize>>>(l_xv, l_n);
    cudaLaunchErrorCheck();

    /** Level == 0 **/

    // Ready to array assignment: do the easiest arrays first, directly copying from the caller
    // starting index passed from F90: l_xv(1:n) which covers 2*xmym buffer; p_xs(1:l_xmymzm) just core part
    // so l_xv(1+xmym:l_xmymzm+xmym), p_xs(1:l_xmymzm)
    // Eidt: copy to device memory to limit CPU page faults
    float *lp_bv, *lp_iv, *lp_xs;
    cudaErrorCheck(cudaMalloc(&lp_bv, sizeof(float) * l_xmymzm));
    cudaErrorCheck(cudaMalloc(&lp_iv, sizeof(float) * l_xmymzm));
    cudaErrorCheck(cudaMalloc(&lp_xs, sizeof(float) * l_xmymzm));
    cudaMemcpy(lp_bv, p_bv, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lp_iv, p_iv, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lp_xs, p_xs, l_xmymzm * sizeof(float), cudaMemcpyHostToDevice);


    // ---------------------------
    // HIP-TODO: Remove #if check when https://github.com/ROCm-Developer-Tools/HIP/issues/1696
    // is fixed
#if !defined(AMBER_PLATFORM_AMD)
    if (p.concurrentManagedAccess) {
        // Seems no need of this two...TBC
        cudaMemPrefetchAsync(l_xv + l_xmym, l_xmymzm * sizeof(float), deviceId);
        cudaMemPrefetchAsync(l_bv, l_xmymzm * sizeof(float), deviceId);
    }
#endif
    // ---------------------------
    nblocks = (l_xmymzm  - 1) / blocksize + 1;
    copy_vector_kernel<<<nblocks, blocksize>>>(l_xv + l_xmym, lp_xs, l_xmymzm);
    cudaLaunchErrorCheck();

    copy_vector_kernel<<<nblocks, blocksize>>>(l_bv, lp_bv, l_xmymzm);
    cudaLaunchErrorCheck();

    // Salt term
    copy_vector_kernel<<<nblocks, blocksize>>>(l_iv, lp_iv, l_xmymzm);
    cudaLaunchErrorCheck();

    // Immediate memory release
    cudaFree(lp_bv);
    cudaFree(lp_iv);
    cudaFree(lp_xs);

    // Set up local eps arrays for data assignment on the kernel
    int m = 0;
    for (int i = 0; i < MG_NLEVEL; i++) {
        m += (mg_size[i][0] + 1) * (mg_size[i][1] + 1) * (mg_size[i][2] + 1);
    }

    /** lepsx/yz/ - All Levels **/
    float *lepsx, *lepsy, *lepsz;
    cudaErrorCheck(cudaMallocManaged(&lepsx, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&lepsy, sizeof(float) * m));
    cudaErrorCheck(cudaMallocManaged(&lepsz, sizeof(float) * m));
    float *epsx_f, *epsy_f, *epsz_f;
    cudaErrorCheck(cudaMalloc(&epsx_f, (l_xmymzm + l_ym * l_zm) * sizeof(float)));
    cudaErrorCheck(cudaMalloc(&epsy_f, (l_xmymzm + l_xm * l_zm) * sizeof(float)));
    cudaErrorCheck(cudaMalloc(&epsz_f, (l_xmymzm + l_xm * l_ym) * sizeof(float)));

    // Copy passed array to UM
    // Wanring: This can make managed epsx_f[i] accessable from CPU, but not from kernel, which will cause
    // CUDA_EXCEPTION_1/14 errors in the kernel.
    //epsx_f = epsx;
    //epsy_f = epsy;
    //epsz_f = epsz;
    cudaMemcpy(epsx_f, epsx, (l_xmymzm + l_ym * l_zm) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(epsy_f, epsy, (l_xmymzm + l_xm * l_zm) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(epsz_f, epsz, (l_xmymzm + l_xm * l_ym) * sizeof(float), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock);
    dim3 blocks((l_xm - 1)/devThreadsPerBlock + 1, (l_ym - 1) / devThreadsPerBlock + 1, (l_zm - 1) / devThreadsPerBlock + 1);
    //!!!!!!!!!!Mark to replace!!!!!!!!!! This is too time-consuming... improve it
    // --------------------------- temporary try. E: it works

    // HIP-TODO: Remove #if check when https://github.com/ROCm-Developer-Tools/HIP/issues/1696
    // is fixed
#if !defined(AMBER_PLATFORM_AMD)
    if (p.concurrentManagedAccess) {
        cudaMemPrefetchAsync(lepsx, m * sizeof(float), deviceId);
        cudaMemPrefetchAsync(lepsy, m * sizeof(float), deviceId);
        cudaMemPrefetchAsync(lepsz, m * sizeof(float), deviceId);
    }
#endif
    // ---------------------------
    feedepsintoam_kernel<<<blocks, threadsPerBlock>>>(l_xm, l_ym, l_zm, lepsx, lepsy, lepsz, epsx_f, epsy_f, epsz_f);
    cudaLaunchErrorCheck();

    //!!!!!!!!!!Mark to replace!!!!!!!!!! Rethink all cudaDeviceSynchronize() - blocks everything, too time-consuming... reduce them where possible
    cudaDeviceSynchronize();

    float lfactor = l_epsout * (l_h * l_pbkappa) * (l_h * l_pbkappa);

    //!!!!!!!!!!Mark to replace!!!!!!!!!! same story, page faults, all inter-threshold-level initializations are inappropriate

    // Finally we are ready to set up the A matrix
    // set up am/ad arrays at the finest level for all solvers
    // so only 1_xmymzm elements of leps* are initialized
    int j = 0;
    m = mg_index[j]; // m == 0 here
    int n = mg_index_ext[j];
    int lxmym = mg_size[j][0] * mg_size[j][1];
    // 1-D grid
    int lxmymzm = mg_size[j][0] * mg_size[j][1] * mg_size[j][2];
    dim3 h_threadsPerBlock (512);
    dim3 h_blocks((lxmymzm - 1) / 512 + 1);
    set_am_ad_kernel_head<<<h_blocks, h_threadsPerBlock>>>(lepsx + m, lepsy + m, lepsz + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, lxmymzm);
    cudaLaunchErrorCheck();
    // 3-D grid
    set_am_ad_kernel_body<<<blocks, threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
    cudaLaunchErrorCheck();
    // 2-D grid
    dim3 t_threadsPerBlock(16, 16);
    dim3 t_blocks((max(mg_size[j][0], mg_size[j][1]) - 1) / 16 + 1, (max(mg_size[j][1], mg_size[j][2]) - 1) / 16 + 1);
    set_am_ad_kernel_tail<<<t_blocks, t_threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
    cudaLaunchErrorCheck();

    cudaDeviceSynchronize();

    if (*solvopt == 2) {
        for (j = 1; j < MG_NLEVEL; j++) {
            int l = mg_index[j-1];
            m = mg_index[j];
            n = mg_index_ext[j];

            lfactor *= 4;
            lxmym = mg_size[j][0] * mg_size[j][1];
            lxmymzm = mg_size[j][0] * mg_size[j][1] * mg_size[j][2];

            if (j < threshold) {
                // Resize
                dim3 blocks((mg_size[j][0] - 1)/devThreadsPerBlock + 1, (mg_size[j][1] - 1) / devThreadsPerBlock + 1, (mg_size[j][2] - 1) / devThreadsPerBlock + 1);
                restrict_eps_map_kernel<<<blocks, threadsPerBlock>>>(lepsx + l, lepsy + l, lepsz + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], lepsx + m, lepsy + m, lepsz+ m, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                cudaLaunchErrorCheck();
                restrict_v_kernel<<<blocks, threadsPerBlock>>>(64.0, l_iv + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]); //iv
                cudaLaunchErrorCheck();

                // 1-D grid
                dim3 h_blocks ((lxmymzm - 1) / 512 + 1);
                set_am_ad_kernel_head<<<h_blocks, h_threadsPerBlock>>>(lepsx + m, lepsy + m, lepsz + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, lxmymzm);
                cudaLaunchErrorCheck();
                // 3-D grid
                set_am_ad_kernel_body<<<blocks, threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
                cudaLaunchErrorCheck();
                // 2-D grid
                dim3 t_blocks((max(mg_size[j][0], mg_size[j][1]) - 1) / 16 + 1, (max(mg_size[j][1], mg_size[j][2]) - 1) / 16 + 1);
                set_am_ad_kernel_tail<<<t_blocks, t_threadsPerBlock>>>(l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                cudaLaunchErrorCheck();

                cudaDeviceSynchronize();

            } else {
                restrict_eps_map(lepsx + l, lepsy + l, lepsz + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], lepsx + m, lepsy + m, lepsz + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]);
                restrict_v(64.0, l_iv + l, mg_size[j - 1][0], mg_size[j - 1][1], mg_size[j - 1][2], l_iv + m, mg_size[j][0], mg_size[j][1], mg_size[j][2]); //iv
                set_am_ad(lepsx + m, lepsy + m, lepsz + m, l_iv + m, l_am1 + n + lxmym, l_am2 + n + lxmym, l_am3 + n + lxmym, l_ad + m, l_bz + m, mg_size[j][0], mg_size[j][1], mg_size[j][2], lfactor, l_epsout);
            }
        }
    }
    cudaFree(lepsx);
    cudaFree(lepsy);
    cudaFree(lepsz);
    cudaFree(epsx_f);
    cudaFree(epsy_f);
    cudaFree(epsz_f);
}

__host__
void restrict_eps_map(float *epsxf, float *epsyf, float *epszf, int xn, int yn, int zn, float *epsxr, float *epsyr, float *epszr, int xnr, int ynr, int znr) {
#ifdef AMBER_PLATFORM_AMD
    int l_bcopt_h;
    cudaMemcpyFromSymbol(&l_bcopt_h, l_bcopt, sizeof(int));
    if (l_bcopt_h == 10) {
#else
    if (l_bcopt == 10) {
#endif
        printf("Not yet for PBC.\n");
        //Do nothing;
    } else {
        for(int k = 0; k < znr; k++) {
            int k2 = 2 * k + 1;
            for(int j = 0; j < ynr; j++) {
                int j2 = 2 * j + 1;
                for(int i = 0; i < xnr; i++) {
                    int i2 = 2 * i + 1;
                    int flatid = i + xnr * j + xnr * ynr * k;
                    // eps*r causes CUDA Exception 15 error. Solved
                    epsxr[flatid] = r_map_exp_x(epsxf, i2, j2, k2, xn, yn);
                    epsyr[flatid] = r_map_exp_y(epsyf, i2, j2, k2, xn, yn);
                    epszr[flatid] = r_map_exp_z(epszf, i2, j2, k2, xn, yn);
                }
            }
        }
    }
}

__host__
void set_am_ad(float *epsx, float *epsy, float *epsz, float *iv, float *lam1, float *lam2, float *lam3, float *lad, float *lbz, int xn, int yn, int zn, float lfactor, float epsout) {
    for (int i = 0; i < xn * yn * zn; i++) {
        lam1[i] = epsx[i];
        lam2[i] = epsy[i];
        lam3[i] = epsz[i];
    }

    for (int k = 0; k < zn; k++) {
        for (int j = 0; j < yn; j++) {
            for (int i = 0; i < xn; i++) {
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
    }

#ifdef AMBER_PLATFORM_AMD
    int l_bcopt_h;
    cudaMemcpyFromSymbol(&l_bcopt_h, l_bcopt, sizeof(int));
    if (l_bcopt_h != 10) {
#else
    if (l_bcopt != 10) {
#endif
        for (int k = 0; k < zn; k++) {
            for (int j = 0; j < yn; j++) {
                for (int i = 0; i < xn; i++) {
                    int flatid = i + xn * j + xn * yn * k;
                    if (i == xn - 1) lam1[flatid] = 0;
                    if (j == yn - 1) lam2[flatid] = 0;
                    if (k == zn - 1) lam3[flatid] = 0;
                }
            }
        }
    }
}

__host__
void restrict_v(float divider, float *bvf, int nx, int ny, int nz, float *bvr, int nxr, int nyr, int nzr) {
    for(int k = 0; k < nzr; k++) {
        int k2 = 2 * k + 1;
        for(int j = 0; j < nyr; j++) {
            int j2 = 2 * j + 1;
            for(int i = 0; i < nxr; i++) {
                int i2 = 2 * i + 1;
                int flatid = i + nxr * j + nxr * nyr * k;
                bvr[flatid] = ( bvf[f_id(i2-1, j2-1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2-1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2  , k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2-1, nx, ny)] + bvf[f_id(i2+1, j2  , k2-1, nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2+1, k2-1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2-1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2-1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2-1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2  , nx, ny)] + bvf[f_id(i2+1, j2-1, k2  , nx, ny)] ) +
                          4 * ( bvf[f_id(i2-1, j2  , k2  , nx, ny)] + 2 * bvf[f_id(i2, j2  , k2  , nx, ny)] + bvf[f_id(i2+1, j2  , k2  , nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2+1, k2  , nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2  , nx, ny)] + bvf[f_id(i2+1, j2+1, k2  , nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2-1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2-1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2-1, k2+1, nx, ny)] ) +
                          2 * ( bvf[f_id(i2-1, j2  , k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2  , k2+1, nx, ny)] + bvf[f_id(i2+1, j2  , k2+1, nx, ny)] ) +
                              ( bvf[f_id(i2-1, j2+1, k2+1, nx, ny)] + 2 * bvf[f_id(i2, j2+1, k2+1, nx, ny)] + bvf[f_id(i2+1, j2+1, k2+1, nx, ny)] );
                bvr[flatid] /= divider;
            }
        }
    }
}

__host__
void restrict_cuda(int level, float *vf, float *vt, int *index, int coef, float divider) {
    //float div = 16.0; //restrict_bv, either from rv or from bv

    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxr = mg_size[level + 1][0];
    int nyr = mg_size[level + 1][1];
    int nzr = mg_size[level + 1][2];

    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock);
    dim3 blocks((nx - 1)/devThreadsPerBlock + 1, (ny - 1) / devThreadsPerBlock + 1, (nz - 1) / devThreadsPerBlock +1);
    restrict_v_kernel<<<blocks, threadsPerBlock>>>(divider, vf + index[level] + coef * nx * ny, nx, ny, nz, vt + index[level + 1] + coef * nxr * nyr, nxr, nyr, nzr);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();
}

__host__
void interpolate(int level) {
    // Derive all vectors
    int xn = mg_size[level + 1][0];
    int yn = mg_size[level + 1][1];
    int zn = mg_size[level + 1][2];
    int xni = mg_size[level][0];
    int yni = mg_size[level][1];
    int zni = mg_size[level][2];
    int xniyni = xni * yni;
    int xniynizni = xni * yni * zni;

    int p1 = mg_x_idx[level + 1] + mg_size[level + 1][0] * mg_size[level + 1][1];
    int p2 = mg_x_idx[level] + mg_size[level][0] * mg_size[level][1];

    float *v = &l_xv[p1];
    float *vi = &l_xv[p2];
    float *lam1 = &l_am1[mg_index_ext[level]]; // 1-xniyni:~
    float *lam2 = &l_am2[mg_index_ext[level]];
    float *lam3 = &l_am3[mg_index_ext[level]];
    float *lbz = &l_bz[mg_index[level]];
    float epsout = l_epsout;

    if (xn * 2 + 1 != xni || yn * 2 + 1 != yni || zn * 2 + 1 != zni) {
        printf("Interpolation failed because of incorrect dimension (interpolate_host)\n");
        printf("xn %d, yn %d, zn %d\n", xn, yn, zn);
        printf("xni %d, yni %d, zni %d\n", xni, yni, zni);
        exit(2);
    }

    // Initialize to epsout, should be before ipl_chain function calls
    // 1-xniyni:0, outside 3-D loop
    for ( int i = 0; i < xniyni; i++) {
        lam1[i] = epsout;
        lam2[i] = epsout;
        lam3[i] = epsout;
    }

    // Three surfaces
    // Need adding index offset. TBF. Fixed
    for (int k = 0; k < zni; k++) {
        for (int j = 0; j < yni; j++) {
            // [(xni - 1) + j * xni + k * xniyni] + xniyni
            lam1[-1 + (j + 1) * xni + (k + 1) * xniyni] = epsout;
        }
    }
    for (int k = 0; k < zni; k++) {
        for (int i = 0; i < xni; i++) {
            // [i + (yni - 1) * xni + k * xniyni] + xniyni
            lam2[i - xni + (k + 2) * xniyni] = epsout;
        }
    }
    for (int j = 0; j < yni; j++) {
        for (int i = 0; i < xni; i++) {
            lam3[i + j * xni + xniynizni] = epsout;
        }
    }

    for(int k = 0; k < zn; k++) {
        for(int j = 0; j < yn; j++) {
            for(int i = 0; i < xn; i++) {
                // Caution with the indexes! Starting from 0
                int flatid_c = i + xn * j + xn * yn * k; // Coarse grid
                int flatid_f = (2 * i + 1) + (2 * j + 1) * xni + (2 * k + 1) * xniyni; // Fine grid

                vi[flatid_f] += v[flatid_c];

                // offset for offsetting the index
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      -1, lam2, xni, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam1,      +1, lam2, xni, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    -xni, lam1,   1, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam2,    +xni, lam1,   1, lam3, xniyni, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, -xniyni, lam2, xni, lam1,      1, xn, yn, zn);
                ipl_chain_h(vi, xniyni, xniynizni, flatid_f, *(v + flatid_c), lbz, lam3, +xniyni, lam2, xni, lam1,      1, xn, yn, zn);
            }
        }
    }

    // 1-xniyni:0, outside 3-D loop
    for ( int i = 0; i < xniyni; i++) {
        lam1[i] = 0.0;
        lam2[i] = 0.0;
        lam3[i] = 0.0;
    }

    // Three surfaces
    for (int k = 0; k < zni; k++) {
        for (int j = 0; j < yni; j++) {
            lam1[-1 + (j + 1) * xni + (k + 1) * xniyni] = 0.0;
        }
    }
    for (int k = 0; k < zni; k++) {
        for (int i = 0; i < xni; i++) {
            lam2[i - xni + (k + 2) * xniyni] = 0.0;
        }
    }
    for (int j = 0; j < yni; j++) {
        for (int i = 0; i < xni; i++) {
            lam3[i + j * xni + xniynizni] = 0.0;
        }
    }

    // A better way to handel and fix all index issues, using index 1 ~ xnynzn, not
    // 1-nxny:nxnynz+nxny.
    // i.e. grid[level]->xs[0:lxlylz-1] v.s. (arrays + index dereference),
    // but at the expense of struct indrection overhead. Figure out to what degree.
    // May define a general Grid structure with all vector pointers inside, then create an array
    // Grid *grid[MG_NLEVEL],
    // then initialize all matrix vectors via grid[level]->xs = (float *) malloc();
}

// Refine the four types point value assignment using directly thread block, layer by layer,
// pass single array value distributed on single thread. TBD
__host__
void ipl_chain_h(float *vi, int xnyn, int xnynzn, int l, float v, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    // Here v is value of coarse grid; vi is array of fine grid
    float v1 = ipl_comp1(v, l, lbz, am_1, xnyn, xnynzn, shift_1);
    int l1 = l + shift_1;
    vi[l1] += v1;
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2, -shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_2,  shift_2, am_3, shift_3, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3, -shift_3, am_2, shift_2, xn, yn, zn);
    ipl_chain2_h(vi, xnyn, xnynzn, l1, v1, lbz, am_1, shift_1, am_3,  shift_3, am_2, shift_2, xn, yn, zn);
}

__host__
void ipl_chain2_h(float *vi, int xnyn, int xnynzn, int l1, float v1, float *lbz, float *am_1, int shift_1, float *am_2, int shift_2, float *am_3, int shift_3, int xn, int yn, int zn) {
    float v2 = ipl_comp2(v1, l1, lbz, am_1, am_2, xnyn, xnynzn, shift_1, shift_2);
    int l2 = l1 + shift_2;
    vi[l2] += v2;
    vi[l2 - shift_3] += ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, -shift_3);
    vi[l2 + shift_3] += ipl_comp3(v2, l2, lbz, am_1, am_2, am_3, xnyn, xnynzn, shift_1, shift_2, +shift_3);
}

__host__
void interpolate_cuda(int level) {
    int xn = mg_size[level + 1][0];
    int yn = mg_size[level + 1][1];
    int zn = mg_size[level + 1][2];
    int xni = mg_size[level][0];
    int yni = mg_size[level][1];
    int zni = mg_size[level][2];

    int p1 = mg_x_idx[level + 1] + mg_size[level + 1][0] * mg_size[level + 1][1];
    int p2 = mg_x_idx[level] + mg_size[level][0] * mg_size[level][1];
    int p3 = mg_index_ext[level];
    int p4 = mg_index[level];

    if (xn * 2 + 1 != xni || yn * 2 + 1 != yni || zn * 2 + 1 != zni) {
        printf("Interpolation failed because of incorrect dimension (interpolate_cuda)\n");
        printf("xn %d, yn %d, zn %d\n", xn, yn, zn);
        printf("xni %d, yni %d, zni %d\n", xni, yni, zni);
        exit(2);
    }

    // 2-D grid
    dim3 h_threadsPerBlock(16, 16);
    dim3 h_blocks((max(xni, yni) - 1) / 16 + 1, (max(yni, zni) - 1) / 16 + 1);
    interpolate_kernel_head<<<h_blocks, h_threadsPerBlock>>>(xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, l_epsout);
    cudaLaunchErrorCheck();

    // Using (8, 8, 4) grid solved the 'CUDA launch failed: too many resources requested' issue.
    // Optimize the register use, block, grid size later. TBD
    dim3 threadsPerBlock(devThreadsPerBlock, devThreadsPerBlock, devThreadsPerBlock/2);
    dim3 blocks((xn - 1) / devThreadsPerBlock + 1, (yn - 1) / devThreadsPerBlock + 1, 2 * (zn - 1) / devThreadsPerBlock + 1);
    interpolate_kernel_body<<<blocks, threadsPerBlock>>>(l_xv + p1, xn, yn, zn, l_xv + p2, xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, l_bz + p4);
    cudaLaunchErrorCheck();

    interpolate_kernel_head<<<h_blocks, h_threadsPerBlock>>>(xni, yni, zni, l_am1 + p3, l_am2 + p3, l_am3 + p3, 0.0);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();
}

__host__
void relax(int level, int ncyc) {
        // Newage test
        //printf("relax on level %d\n", level);
    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;
    float *xs = l_xv + mg_x_idx[level];
    float *lam1 = l_am1 + mg_index_ext[level];
    float *lam2 = l_am2 + mg_index_ext[level];
    float *lam3 = l_am3 + mg_index_ext[level];
    float *lzv = l_zv + mg_index[level];
    float *lad = l_ad + mg_index[level];
    float *lbv = l_bv + mg_index[level];
    float *lrv = l_rv + mg_index[level];
    float accept = l_accept;

    float onorm = mg_onorm[level];

    int itn_checknorm;
    float wsor;//, wsor1;
    float linorm = 0.0;
    float lnorm;
    int itmax = 20; // Should move to upper layer

    if (ncyc > 0) {
        itn_checknorm = ncyc;
        wsor = 1.0; //1.0; Debug: 1.0 is Gauss-Seidel
    } else {
        itn_checknorm = 10;
        wsor = 1.9; // 1.0; Debug of omega - use GS to solve on coarsest level
    }

    for (int i = 0; i < nxnynz; i++) {
        linorm += abs(lbv[i]);
        lzv[i] = 1.0 / lad[i];
    }

    bool converged = false;
    int litn = 0;
    while (!converged) {
        for (int i = nxny;  i < nxnynz+nxny; i++) {
            xs[i] -= wsor * (xs[i] - (lam1[i - 1   ] * xs[i - 1   ] + lam1[i         ] * xs[i + 1   ] +
                                      lam2[i - nx  ] * xs[i - nx  ] + lam2[i         ] * xs[i + nx  ] +
                                      lam3[i - nxny] * xs[i - nxny] + lam3[i         ] * xs[i + nxny] + lbv[i - nxny]) * lzv[i - nxny]);
        }

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            // residual
            for (int i = nxny; i < nxnynz+nxny; i++) {
                lrv[i - nxny] = lam1[i - 1   ] * xs[i - 1   ] + lam1[i       ] * xs[i + 1   ] +
                                lam2[i - nx  ] * xs[i - nx  ] + lam2[i       ] * xs[i + nx  ] +
                                lam3[i - nxny] * xs[i - nxny] + lam3[i       ] * xs[i + nxny] + lbv[i - nxny] - lad[i - nxny] * xs[i];
                }

            // norm
            lnorm = 0.0;
            for (int i = 0; i < nxnynz; i++) {
                lnorm += abs(lrv[i]);
            }

            // Newage
            //printf("    ncyc %d\t, litn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);

            if (litn >= itmax || (ncyc > 0 && (litn >= ncyc && lnorm < onorm)) || lnorm <= accept * linorm) {
                converged = true;
                if (ncyc > 0 && litn >= ncyc && lnorm > onorm) {
                    printf("PB_MG FAILED: ncyc %d\t, itn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);
                    // Continue to shift to cg_cuda
                    //exit(2);
                    break;
                }

                if (ncyc > 0) mg_onorm[level] = lnorm; // Update global array
                if (litn >= itmax) printf("PB_MG WARNING: SOR maxitn exceeded (relax_host)!\n");
            }

        }
    } // while
}

__host__
void relax_cuda(int level, int ncyc) {

    int threadsPerBlock = 512; // This shouldn't be final, RL

    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;

    float *xs = l_xv + mg_x_idx[level];
    float *lam1 = l_am1 + mg_index_ext[level];
    float *lam2 = l_am2 + mg_index_ext[level];
    float *lam3 = l_am3 + mg_index_ext[level];
    float *lzv = l_zv + mg_index[level];
    float *lad = l_ad + mg_index[level];
    float *lbv = l_bv + mg_index[level];
    float *lrv = l_rv + mg_index[level];
    float accept = l_accept;
    float onorm = mg_onorm[level];

    int itn_checknorm = 10;
    float wsor;//, wsor1;
    int itmax = 100;

    if (ncyc > 0) {
        itn_checknorm = ncyc;
        wsor = 1.0;// 1.0 Test
    } else {
        itn_checknorm = 10;
        wsor = 1.9; // 1.9 Test
    }

    ncyc = 10; // This shouldn't be final, RL

    // Inverse AD for fast processing later
    int blocks = (nxnynz - 1) / threadsPerBlock + 1;
    inv_vector_kernel<<<blocks, threadsPerBlock>>>(lad, lzv, nxnynz);
    cudaLaunchErrorCheck();

    // Initial norm. Create cuBlAS context
    float linorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, nxnynz, lbv, 1, &linorm);
    cublasErrorCheck(status);

    bool converged = false;
    int litn = 0;
    int sblocks = ((int)(nxnynz / 2) + (nxnynz & 1) - 1) / threadsPerBlock + 1;
    while (!converged) {
        // non-periodic
        solver_red_kernel<<<sblocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lzv, lbv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        solver_black_kernel<<<sblocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lzv, lbv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        //cudaDeviceSynchronize(); // Warning of bug: CUDA_EXCEPTION_15

        //printf("    litn %d:", litn);

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            residue_kernel<<<blocks, threadsPerBlock>>>(xs, lam1, lam2, lam3, lad, lbv, nx, ny, nz, lrv);
            cudaLaunchErrorCheck();

            float lnorm = 0.0;
            status = cublasSasum(handle, nxnynz, lrv, 1, &lnorm);
            cublasErrorCheck(status);

            //printf("    ncyc %d\t, litn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);

            if (litn >= itmax || (ncyc > 0 && (litn >= ncyc && lnorm < onorm)) || lnorm <= accept * linorm) {
                converged = true;
                if (ncyc > 0 && litn >= ncyc && lnorm > onorm) {
                    printf("RELAX_CUDA FAILED: ncyc %d\t, litn %d\t, norm %e\t, onorm %e\n", ncyc, litn, lnorm, onorm);
                    // Continue to solve with jacobi
                    //exit(2);
                    break;
                }

                if (ncyc > 0) mg_onorm[level] = lnorm; // Update global array
                if (litn >= itmax) printf("PB_MG WARNING: SOR maxitn exceeded (relax_kernel)!\n");
            }
        } // if
    } // while

    // Destroy context
    cublasDestroy(handle);
}

__host__
void cg_cuda(int level, int ncyc) {
        //printf("    cg_cuda on level %d\n", level);
    int nx = mg_size[level][0];
    int ny = mg_size[level][1];
    int nz = mg_size[level][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;

    float *xs = l_xv + mg_x_idx[level] + nxny;
    float *lbv = l_bv + mg_index[level];

    // Padding DIA matrix
    // Note: the initial padding at the tail (e.g. bottom z face for l_am3) of l_am* should be ZERO, or will incorrectly regonize the padding elements outside of DIA matix as part of matrix. As in the following example, "7" will be counted in the DIA matrix in CUSP.
    // 0   [1 4 0]
    //   0 [0 2 5]
    //     [6 0 3] 7
    // A_values = [lam6 lam5 lam4 lad lam1 lam2 lam3]
    float *lam1 = l_am1 + mg_index_ext[level] + nxny;
    float *lam2 = l_am2 + mg_index_ext[level] + nxny;
    float *lam3 = l_am3 + mg_index_ext[level] + nxny;
    float *lam4 = lam1 - 1;
    float *lam5 = lam2 - nx;
    float *lam6 = lam3 - nxny;
    float *lad = l_ad + mg_index[level];

    int h_dia_offsets[7] = {-nxny, -nx, -1, 0, 1, nx, nxny};
    int *d_dia_offsets;
    cudaErrorCheck(cudaMalloc(&d_dia_offsets, 7 * sizeof(int)));
    cudaMemcpy(d_dia_offsets, h_dia_offsets, 7 * sizeof(int), cudaMemcpyHostToDevice);

    // Copy all padded diagonals to A.values
    float *values;
    cudaErrorCheck(cudaMalloc(&values, 7 * nxnynz * sizeof(float)));
    cudaMemcpy(values             , lam6, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values +     nxnynz, lam5, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values + 2 * nxnynz, lam4, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values + 3 * nxnynz, lad , nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values + 4 * nxnynz, lam1, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values + 5 * nxnynz, lam2, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(values + 6 * nxnynz, lam3, nxnynz * sizeof(float), cudaMemcpyDeviceToDevice);

    // Change the signs of off-diagonal matrix arrays
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    float alpha = -1.0;
    status = cublasSscal(handle, 3 * nxnynz, &alpha, values, 1);
    cublasErrorCheck(status);
    status = cublasSscal(handle, 3 * nxnynz, &alpha, values + 4 * nxnynz, 1);
    cublasErrorCheck(status);
    // Destroy context
    cublasDestroy(handle);

    //int itn_checknorm = 10;
    //float wsor;//, wsor1;
    int bwidth = 7;
    int itmax = 3000;
    float accept = l_accept;
    //int litn = 0;

    // Need to use pointer to return l_norm due to device implementation
    //jacobi_dia_solver(xs, values, d_dia_offsets, lbv, bwidth, nx, ny, nz, itmax, accept, &litn, &l_norm);
    jacobi_dia_solver(xs, values, d_dia_offsets, lbv, bwidth, nx, ny, nz, itmax, accept, &l_norm);

    //printf("    ncyc %d\t, litn %d\t, itmax %d\t, norm %e\n", ncyc, litn, itmax, l_norm);
    //printf("    ncyc %d\t, litn %d\t, itmax %d\t, norm %e\t, onorm %e\n", ncyc, litn, itmax, l_norm, oldnorm);

    cudaFree(values);
    cudaFree(d_dia_offsets);
}

// cusp DIA CG solver with the jacobi preconditioner
// see cusp_LinearSolvers.cu for details
__host__
void jacobi_dia_solver(float *x, float *values, int *dia_offsets, float *rhs, int bwidth, int xm, int ym, int zm, int maxitn, float acpt, float *lnorm) {
	int xmym = xm * ym;
	int N = xm * ym * zm;
	int nnz = N + 2 * (N - 1 + N - xm + N - xmym);

        cusp::dia_matrix<int, float, cusp::device_memory> A(N, N, nnz, bwidth);

        thrust::device_ptr<int> wrap_offsets(dia_offsets);
        thrust::device_ptr<float> wrap_values(values);

        typedef typename cusp::array1d_view< thrust::device_ptr<int> > DeviceIndexArray1dView;
        typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArray1dView;
        DeviceIndexArray1dView view_offsets(wrap_offsets, wrap_offsets + 7);
        DeviceValueArray1dView view_values(wrap_values, wrap_values + 7 * N);

        typedef cusp::array2d_view<DeviceValueArray1dView, cusp::column_major> DeviceValueArray2dView;
        DeviceValueArray2dView view2d_values(N, 7, N, view_values);

        A.diagonal_offsets = view_offsets;
        A.values = view2d_values;

        thrust::device_ptr<float> wrap_x(x), wrap_b(rhs);
        DeviceValueArray1dView cg_x(wrap_x, wrap_x + N), cg_b(wrap_b, wrap_b + N);

	cusp::monitor<float> monitor(cg_b, maxitn, acpt, 0, false);

	cusp::precond::diagonal<float, cusp::device_memory> M(A);

	cusp::krylov::cg(A, cg_x, cg_b, monitor, M);

	*lnorm = monitor.residual_norm();

}

// Recursive F-Cycle
__host__
void FCycle(int level) {

    if (level == MG_NLEVEL - 1) {
        // Solve on coarsest grid
    //!!!!!!!!!!Mark to replace!!!!!!!!!! Use threshold as if-condition, relax_cuda(level, -1)
        //relax_cuda(level, -1);
        if (level < threshold)
            relax_cuda(level, -1);
        else
            relax(level, -1);
    } else {
        // Restrict bv & initialize xv
        //int vnx = mg_size[level][0];
        //int vny = mg_size[level][1];
        //int vnz = mg_size[level][2];
        //int vnxny = vnx * vny;
        //int vnxnynz = vnxny * vnz;
        //int pa = mg_x_idx[level] + vnxny;

        if (level < threshold) {
            // On CUDA
            // Restrict bv
            restrict_cuda(level, l_bv, l_bv, mg_index, 0, 16.0); // May save to another array if too much overhead
            //restrict_cuda(level, l_xv, l_xv, mg_x_idx, 1, 128.0); // Failed. No benefit for re-guessing xv uing more FCycle

            // Reinitialize l_xv on level (not level+1 as in VCycle!) to zero
            // For only 1 round FCycle, no need to reinitialize, as already done in init_array_cuda
            //
            //int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;
            //int nblocks = (vnxnynz - 1) / blocksize + 1;
            //init_vector_kernel<<<nblocks, blocksize>>>(l_xv + pa, vnxnynz);
            //cudaLaunchErrorCheck();
            //cudaDeviceSynchronize();
        } else {
            // On CPU
            // Restrict bv
            restrict_v(16.0, l_bv + mg_index[level], mg_size[level][0], mg_size[level][1], mg_size[level][2], l_bv + mg_index[level + 1], mg_size[level + 1][0], mg_size[level + 1][1], mg_size[level + 1][2]); //bv
            // Restrict xv (for more FCycle)
            //restrict_v(16.0, l_xv + mg_x_idx[level] + vnxny, mg_size[level][0], mg_size[level][1], mg_size[level][2], l_xv + mg_x_idx[level + 1], mg_size[level + 1][0], mg_size[level + 1][1], mg_size[level + 1][2]); //xv

            // Reinitialize l_xv on level (not level+1 as in VCycle!) to zero
            // Same as above
            //int pb = pa + vnxnynz;
            //for (int i = pa; i < pb; i++) {
            //    l_xv[i] = 0.0;
            //}
        }

        // Recursive call
        FCycle(level + 1);

        // non-addition Interpolate & solve with one VCycle
        // On CUDA
        if (level < threshold)
    //!!!!!!!!!!Mark to replace!!!!!!!!!! check if this just works with pinned mem on level == threshold already set
            interpolate_cuda(level);
        else
            // On CPU
            interpolate(level);

        VCycle(level);
    }
}

//*****************************************************
// Recursive V-Cycle
__host__
void VCycle(int level) {

    if (level == MG_NLEVEL - 1) {
        // Solve on coarsest grid
    //!!!!!!!!!!Mark to replace!!!!!!!!!! Use threshold as if-condition, relax_cuda(level, -1)
        //relax_cuda(level, -1);

        if (level < threshold)
            relax_cuda(level, -1);
        else
            relax(level, -1);

    } else {
        // Relax & restrict
        if (level < threshold) {
            // On CUDA
            relax_cuda(level, ncyc_before);

//if (level == 1 || level == 2) {
//    relax_cuda(level, 1); // SUR
//} else {
//    relax_cuda(level, 10);
//}
///*@@ Debug of omega
//if (level == 0) {
//    relax_cuda(level, 0); // GS
//} else if (level == 1) { // level 1
//    relax_cuda(level, 1); // SUR, downgoing
//} else { // level 2
//    relax_cuda(level, 10);
//}*/


    //!!!!!!!!!!Mark to replace!!!!!!!!!! page faults on level+1 == threshold, use page-lock mem on level+1
            restrict_cuda(level, l_rv, l_bv, mg_index, 0, 16.0); // from rv to bv

        } else {
            // On CPU
            relax(level, ncyc_before);

            restrict_v(16.0, l_rv + mg_index[level], mg_size[level][0], mg_size[level][1], mg_size[level][2], l_bv + mg_index[level + 1], mg_size[level + 1][0], mg_size[level + 1][1], mg_size[level + 1][2]); //bv

        }

        // Reinitialize l_xv on level+1 to zero
        int vnx = mg_size[level + 1][0];
        int vny = mg_size[level + 1][1];
        int vnz = mg_size[level + 1][2];
        int vnxny = vnx * vny;
        int vnxnynz = vnxny * vnz;
        int pa = mg_x_idx[level + 1] + vnxny;
        if (level < threshold - 1) {
            int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;
            int nblocks = (vnxnynz - 1) / blocksize + 1;
            init_vector_kernel<<<nblocks, blocksize>>>(l_xv + pa, vnxnynz);
            cudaLaunchErrorCheck();
            cudaDeviceSynchronize();
        } else {
            int pb = pa + vnxnynz;
            for (int i = pa; i < pb; i++) {
                l_xv[i] = 0.0;
            }
        }

        // Recursive call
        VCycle(level + 1);

        // Interpolate & relax
        if (level < threshold) {
            // On CUDA
    //!!!!!!!!!!Mark to replace!!!!!!!!!! check if this just works with pinned mem on level == threshold already set
    /*
    FILE *ptt = fopen("ip_xv.dat", "a");
        int inx = mg_size[level ][0];
        int iny = mg_size[level ][1];
        int inz = mg_size[level ][2];
        int st = mg_x_idx[level ] + inx * iny;
    if (level == 2) {
        fprintf(ptt, "before: interplate_cuda on level %d\n", level);
        for (int i = st; i < st + inx * iny * inz; i++) {
            fprintf(ptt, "l_xv[%d] = %e\n", i, l_xv[i]);
        }
    }
    */

            interpolate_cuda(level);

    /*
    if (level == 2) {
        fprintf(ptt, "after: interplate_cuda on level %d\n", level);
        for (int i = st; i < st + inx * iny * inz; i++) {
            fprintf(ptt, "l_xs[%d] = %e\n", i, l_xv[i]);
        }
    }
    fclose(ptt);
    */
            relax_cuda(level, ncyc_after);

//if (level == 1 || level == 2) {
//    relax_cuda(level, 2); // SOR
//} else {
//    relax_cuda(level, 10);
//}
///*@@ Debug of omega
//if (level == 0) {
//    relax_cuda(level, 0); // GS
//} else if (level == 1) { // level 1
//    relax_cuda(level, 2); // SOR, upgoing
//} else { // level 2
//    relax_cuda(level, 10);
//}*/

        } else {
            // On CPU
            interpolate(level);
            relax(level, ncyc_after);
        }
    }
}

//*****************************************************
// PB/MG CUDA driver
extern "C"
void pb_mg_cuda_(float *phi_f, float *xs_f) {
    l_itn = 0;

    // Initial norm
    l_inorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, l_xmymzm, l_bv, 1, &l_inorm);
    cublasErrorCheck(status);

    // Newage test
    for (int i = 0; i < MG_NLEVEL; i++) {
        mg_onorm[i] = 9.9E99;
    }
    //printf("Preheat with one FCycle # %d, linorm %e\n", l_itn, l_inorm);
    FCycle(0);

    l_norm = l_inorm; // Initialize to l_inorm
    float pre_norm;
    float rate = 1.0;
    bool mgconverged = false;

    while (!mgconverged) {
        l_itn++;

        // Newage test
        //printf("  Before VCycle # %d, linorm %e\n", l_itn, l_inorm);

        pre_norm = l_norm;

        // Norm
        status = cublasSasum(handle, l_xmymzm, l_rv, 1, &l_norm);
        cublasErrorCheck(status);

        rate = l_norm / pre_norm;

        if (l_norm <= l_inorm * l_accept) {
            mgconverged = true;
            break;
        } else if (rate > 0.90) {
            // Continue to solve with cg_cuda
            // Move this info to stdout (log file) of F90. TBD
            printf("PB_MG Notice: Slow to converge, continuing to solve with Jacobi-PCG.\n");
            cg_cuda(0, -1);

            mgconverged = true;
            break;
        } else if (l_itn >= l_maxitn) {
            mgconverged = true;
            printf("PB_MG WARNING: maxitn exceeded (pb_mg_cuda)!\n");
            break;
        }

        // VCycle
        for (int i = 0; i < MG_NLEVEL; i++) {
            mg_onorm[i] = 9.9E99;
        }
        VCycle(0);
    } // while

    // Destroy context
    cublasDestroy(handle);

    cudaMemcpy(xs_f, l_xv + l_xmym, l_xmymzm * sizeof(float), cudaMemcpyDeviceToHost);

#ifdef AMBER_PLATFORM_AMD
    int l_bcopt_h;
    cudaMemcpyFromSymbol(&l_bcopt_h, l_bcopt, sizeof(int));
    if (l_bcopt_h != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
#else
    if (l_bcopt != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
#endif
}

//*****************************************************
// PB/SOR CUDA driver
extern "C"
void pb_sor_cuda_(float *phi_f, float *xs_f) {
    int threadsPerBlock = 512; // This shouldn't be final, RL

    int nx = l_xm; //mg_size[0][0];
    int ny = l_ym; //mg_size[0][1];
    int nz = l_zm; //mg_size[0][2];
    int nxny = nx * ny;
    int nxnynz = nxny * nz;

    int itn_checknorm = 10;
    float wsor = 1.9;//, wsor1;

    // inverse AD for fast processing later
    int blocks = (nxnynz - 1) / threadsPerBlock + 1;
    inv_vector_kernel<<<blocks, threadsPerBlock>>>(l_ad, l_zv, nxnynz);
    cudaLaunchErrorCheck();

    // initial norm
    float linorm = 0.0;
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    cublasErrorCheck(status);
    status = cublasSasum(handle, nxnynz, l_bv, 1, &linorm);
    cublasErrorCheck(status);

    bool converged = false;
    int litn = 0;
    float lnorm = 0.0;
    int sblocks = ((int)(nxnynz / 2) + (nxnynz & 1) - 1) / threadsPerBlock + 1;
    while (!converged) {
        solver_red_kernel<<<sblocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_zv, l_bv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        solver_black_kernel<<<sblocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_zv, l_bv, wsor, nx, ny, nz);
        cudaLaunchErrorCheck();
        //cudaDeviceSynchronize(); // Warning of bug: CUDA_EXCEPTION_15

        litn++;

        // Check convergence
        if (litn % itn_checknorm == 0) {
            // residue
            residue_kernel<<<blocks, threadsPerBlock>>>(l_xv, l_am1, l_am2, l_am3, l_ad, l_bv, nx, ny, nz, l_rv);
            cudaLaunchErrorCheck();

            // norm
            status = cublasSasum(handle, nxnynz, l_rv, 1, &lnorm);
            cublasErrorCheck(status);

            if (litn >= l_maxitn || lnorm <= l_accept * linorm) {
                converged = true;
                if (litn >= l_maxitn) printf("PB_SOR WARNING: maxitn exceeded (kernel)!\n");
            } // if
        }
    } // while

    // Destroy context
    cublasDestroy(handle);

    l_itn = litn;
    l_inorm = linorm;
    l_norm = lnorm;

    cudaMemcpy(xs_f, l_xv + l_xmym, l_xmymzm * sizeof(float), cudaMemcpyDeviceToHost);
#ifdef AMBER_PLATFORM_AMD
    int l_bcopt_h;
    cudaMemcpyFromSymbol(&l_bcopt_h, l_bcopt, sizeof(int));
    if (l_bcopt_h != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
#else
    if (l_bcopt != 10 || l_pbkappa != 0) memcpy(phi_f, xs_f, l_xmymzm * sizeof(float));
#endif
}
//*****************************************************

// Return values
extern "C"
int get_itn_() {
    return l_itn;
}

extern "C"
float get_inorm_() {
    return l_inorm;
}

extern "C"
float get_norm_() {
    return l_norm;
}
