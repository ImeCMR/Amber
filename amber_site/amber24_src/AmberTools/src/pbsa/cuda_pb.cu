/*
 * CUDA kernels for MMPBSA routines
 * Ruxi Qi @ UC Irvine, 2018-2019
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
// For gdb
//#include <signal.h>
// For timing
//#include <sys/time.h>
// For error handling and device pickup
#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime.h>
#  include <hip/hip_runtime_api.h>
#  include "hip_definitions.h"
#else
// CUDA Runtime
#  include <cuda_runtime.h>
#endif
//#include "cublas_v2.h"
//#include <cub/cub.cuh>
#include "helper_cuda.h"
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

/*
 * TODO: Common variables shared by more than one routines, remove redundancy,
 * Redesign PBSA code flows, transport common resources to/from device only once;
 * Change atomicAdd into two splited kernels: calculating one and reducing one.
 *******
   float h;
   int devThreadsPerBlock = 8;
   atmlast == natom
   acrd == gcrd
 *******
 */

/***********************************************************
 * Reaction field energy
 ***********************************************************/

__global__
void epb_kernel(int nbnd, float *iepsav, float *polchrg, int atmlast, float *acrd, float *acrg, float gox, float goy, float goz, float h, float *eel, float dprob, float *radi, float *arccrd) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nbnd) {
        float g1 = gox + h * iepsav[4 * i    ];
        float g2 = goy + h * iepsav[4 * i + 1];
        float g3 = goz + h * iepsav[4 * i + 2];
        int iatm = iepsav[4 * i + 3];

        float x, y, z;
        float dist, rdist;
        if (iatm > 0) {
            x = acrd[3 * (iatm - 1)    ];
            y = acrd[3 * (iatm - 1) + 1];
            z = acrd[3 * (iatm - 1) + 2];
            dist = radi[iatm - 1];
        } else if (iatm < 0) {
            x = arccrd[3 * (- iatm - 1)    ];
            y = arccrd[3 * (- iatm - 1) + 1];
            z = arccrd[3 * (- iatm - 1) + 2];
            dist = dprob;
        } else {
            x = g1;
            y = g2;
            z = g3;
            dist = 1.0;
        }

        float dx = g1 - x;
        float dy = g2 - y;
        float dz = g3 - z;
        float d2 = dx * dx + dy * dy + dz * dz;

        if (d2 == 0.0)
            rdist = 1.0;
        else
            rdist = dist / sqrt(d2);

        float crd1 = x + dx * rdist;
        float crd2 = y + dy * rdist;
        float crd3 = z + dz * rdist;

        float pc = polchrg[i];

        // Loop over atom
        float eelrf = 0.0;
        for (int j = 0; j < atmlast; j++) {
            dx = crd1 - acrd[3 * j    ];
            dy = crd2 - acrd[3 * j + 1];
            dz = crd3 - acrd[3 * j + 2];

            eelrf += pc * acrg[j] / sqrt(dx * dx + dy * dy + dz * dz);

            /*
            if (i == 180 && j < 3) {
                printf("Inside kernel, iepsav: i [%d] %f %f %f\n", i, iepsav[4 * i], iepsav[4 * i + 1], iepsav[4 * i + 2]);
                printf("Inside kernel, i [%d] j [%d] xi %f yi %f zi %f xj %f yj %f zj %f\n", i, j, crd1, crd2, crd3, acrd[3 * j], acrd[3 * j + 1], acrd[3 * j + 2]);
                printf("Inside kernel, i [%d] j [%d] iatm %d dx %f dy %f dz %f de %e eelrf %e\n", i, j, iatm, dx, dy, dz, de, eelrf);
            }
            */
        }
        eel[i] = eelrf;
    }
}

extern "C"
void cuda_lpbene_(int *nbnd, float *iepsav, float *pol_charge, int *atmlast, float *acrd, float *acrg, float *gox, float *goy, float *goz, float *h, float *eel, float *dprob, float *radi, float *arccrd, int *alen) {
    float *d_iepsav;
    float *d_polchrg;
    float *d_eel;
    float *d_acrd;
    float *d_acrg;
    float *d_radi;
    float *d_arccrd;
    //float *eel_t;

    //printf("Passed alen arccrd : %d %f %f %f\n", *alen, arccrd[3*38372], arccrd[3*38372+1], arccrd[3*38372+2]);

    cudaErrorCheck(cudaMalloc(&d_iepsav, sizeof(float) * 4 * *nbnd));
    cudaErrorCheck(cudaMalloc(&d_polchrg, sizeof(float) * *nbnd));
    cudaErrorCheck(cudaMalloc(&d_eel, sizeof(float) * *nbnd));
    cudaErrorCheck(cudaMalloc(&d_acrd, sizeof(float) * 3 * *atmlast));
    cudaErrorCheck(cudaMalloc(&d_acrg, sizeof(float) * *atmlast));
    cudaErrorCheck(cudaMalloc(&d_radi, sizeof(float) * *atmlast));
    cudaErrorCheck(cudaMalloc(&d_arccrd, sizeof(float) * 3 * *alen));
    //cudaErrorCheck(cudaMalloc(&eel_t, sizeof(float)));

    cudaMemcpy(d_iepsav, iepsav, sizeof(float) * 4 * *nbnd, cudaMemcpyHostToDevice);
    cudaMemcpy(d_polchrg, pol_charge, sizeof(float) * *nbnd, cudaMemcpyHostToDevice);
    cudaMemcpy(d_acrd, acrd, sizeof(float) * 3 * *atmlast, cudaMemcpyHostToDevice);
    cudaMemcpy(d_acrg, acrg, sizeof(float) * *atmlast, cudaMemcpyHostToDevice);
    cudaMemcpy(d_radi, radi, sizeof(float) * *atmlast, cudaMemcpyHostToDevice);
    cudaMemcpy(d_arccrd, arccrd, sizeof(float) * 3 * *alen, cudaMemcpyHostToDevice);

    int devThreadsPerBlock = 8;
    int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;
    int nblocks = (*nbnd - 1) / blocksize + 1;
    epb_kernel<<<nblocks, blocksize>>>(*nbnd, d_iepsav, d_polchrg, *atmlast, d_acrd, d_acrg, *gox, *goy, *goz, *h, d_eel, *dprob, d_radi, d_arccrd);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();

    // Reduction sum by CUB. For developing purpose only, avoiding library redistribution
    /*void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_eel, eel_t, *nbnd);
    cudaErrorCheck(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_eel, eel_t, *nbnd);
    */

    // Use thrust for reduction sum. Note cublas only supports absolute vector sum
    thrust::device_ptr<float> wrap_eel(d_eel);
    *eel = thrust::reduce(wrap_eel, wrap_eel + *nbnd);

    //cudaMemcpy(eel, eel_t, sizeof(float), cudaMemcpyDeviceToHost);

    //printf(" in C, eel = %f\n", *eel);

    cudaFree(d_iepsav);
    cudaFree(d_polchrg);
    cudaFree(d_eel);
    cudaFree(d_acrd);
    cudaFree(d_acrg);
    cudaFree(d_radi);
    cudaFree(d_arccrd);
    //cudaFree(eel_t);

    //cudaDeviceReset(); // Called in MG
}

/***********************************************************
 * Levelset
 ***********************************************************/

__device__
float density_atom(float dist) {
    float dash[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    float spcoef[5][4] = {
        {1.0 , -4.527143, -3.640532 ,  32.631235},
        {0.21, -2.067608, 15.938209 , -35.500854},
        {0.15,  0.047573, -5.362303 ,  13.122180},
        {0.05, -0.522686,  2.5110050,  -4.487867},
        {0.01, -0.056828, -0.181716 ,   1.079289},
    };
    float density = 0.0;

    if (dist <= 0.0)
        density = 1.0 - 4.527143 * dist;
    else if (dist > 0.0 && dist <= 1.0) {
        for (int i = 0; i < 5; i++) {
            if (dist > dash[i] && dist <= dash[i + 1])
                density +=
                    spcoef[i][0] +
                    spcoef[i][1] * (dist - dash[i]) +
                    spcoef[i][2] * (dist - dash[i]) * (dist - dash[i]) +
                    spcoef[i][3] * (dist - dash[i]) * (dist - dash[i]) * (dist - dash[i]);
        }
    }

    return density;
}

__global__
void lvlset_kernel(int natom, float *lvlset, int xm, int ym, int zm, float *radi, float rh, float *gcrd, float cutoffh, float dprobh) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < natom) {
        float range0 = radi[id] * rh;
        if (range0 == 0.0) return;
        float xi = gcrd[3 * id    ];
        float yi = gcrd[3 * id + 1];
        float zi = gcrd[3 * id + 2];
        // Skip eneopt == 4 now, add later
        float range1 = max(cutoffh, range0 + 2.0 * dprobh + 1.0);
        if (zi + range1 < 0.0 || zi - range1 > zm + 1.0) return;
        float lowk = max(1.0, ceil(zi - range1));
        float highk = min((float) zm, floor(zi + range1));
        for (int k = lowk; k <= highk; k++) {
            float range2 = sqrt(range1 * range1 - (zi - k) * (zi - k));
            if (yi + range2 < 0.0 || yi - range2 > ym + 1.0) return;
            float lowj = max(1.0, ceil(yi - range2));
            float highj = min((float) ym, floor(yi + range2));
            for (int j = lowj; j <= highj; j++) {
                float range3 = sqrt(range2 * range2 - (yi - j) * (yi - j));
                if (range3 == 0.0) return;
                float lowi = max(1.0, ceil(xi - range3));
                float highi = min((float) xm, floor(xi + range3));
                for ( int i = lowi; i <= highi; i++) {
                    float dist = sqrt((xi - i) * (xi -i) + (yi -j) * (yi - j) + (zi - k) * (zi - k)) - range0;
                    dist *= 0.5 / dprobh;
                    // Dimension: xm + 2, ym + 2, zm + 2
                    atomicAdd(&lvlset[i + (xm + 2) * j + (xm + 2) * (ym + 2) * k], density_atom(dist));
                    //ndenatm[i, j, k] += 1; //Skip now
                }
            }
        }
    }
}

extern "C"
void cuda_lvlset_(int *natom, float *lvlset, int *xm, int *ym, int *zm, float *radi, float *rh, float *gcrd, float *cutoffh, float *dprobh) {
    float *d_lvlset;
    float *d_radi;
    float *d_gcrd;

    int N = (*xm + 2) * (*ym + 2) * (*zm + 2);
    cudaErrorCheck(cudaMalloc(&d_lvlset, sizeof(float) * N));
    cudaErrorCheck(cudaMalloc(&d_radi, sizeof(float) * *natom));
    cudaErrorCheck(cudaMalloc(&d_gcrd, sizeof(float) * 3 * *natom));
    cudaMemcpy(d_lvlset, lvlset, sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_radi, radi, sizeof(float) * *natom, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gcrd, gcrd, sizeof(float) * 3 * *natom, cudaMemcpyHostToDevice);

    int devThreadsPerBlock = 8;
    int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;
    int nblocks = (*natom - 1) / blocksize + 1;
    lvlset_kernel<<<nblocks, blocksize>>>(*natom, d_lvlset, *xm, *ym, *zm, d_radi, *rh, d_gcrd, *cutoffh, *dprobh);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();

    cudaMemcpy(lvlset, d_lvlset, sizeof(float) * N, cudaMemcpyDeviceToHost);

    cudaFree(d_lvlset);
    cudaFree(d_radi);
    cudaFree(d_gcrd);

    //cudaDeviceReset(); // Called in MG
}

/***********************************************************
 * Non-bonded list | SA
 ***********************************************************/

// Algorithm
// 1. Mark all (i, j) pairs with the flag 1 for paired and 0 for non-paired in ij_pair; [On device]
//    To get the pair indices from the linear thread index of non-repetitive comparison, the
//    following mappings are used:
//        i = n - 2 - floor(sqrt(-8 * id + 4 * n * (n - 1) - 7) / 2.0 - 0.5),
//        j = id + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
//    and the inversed operation:
//        id = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
// 2. Sum all the paired for each atom and save the counts into nb_offset; [On device] [Could be reduction sum]
// 3. Update nb_offset by incremental offset; [On host]
// 4. Allocate memory for nb_list, save all paired j's into nb_list for each atom. [On device]
//
// TODO: use shared memory for natex - copy exclude list of atom i into shared memory first

#ifdef NBLIST
// Update: since the SA part in pb_atmlist is heavily tangled with other nblist/offset, so need to split them first to make this work

__global__
void ijpair_kernel(int n, int *nshrt, int *natex, int *ij_pair, float *acrd, float cutsa) { // n(n-1)/2 threads
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < n * (n - 1) / 2) { // Non-repeatitive (i, j) pairs, i < j
        int i = n - 2 - floor(sqrt(-8 * id + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
        int j = id + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;

        ij_pair[id] = 0;
        for (int k = nshrt[i]; k < nshrt[i + 1]; k++) {
            if (j == natex[k]) {
                // j is excluded by i
                return;
            }
        }

        float dx = acrd[3 * i    ] - acrd[3 * j    ];
        float dy = acrd[3 * i + 1] - acrd[3 * j + 1];
        float dz = acrd[3 * i + 2] - acrd[3 * j + 2];
        float d2 = dx * dx + dy * dy + dz * dz;

        if (d2 <= cutsa) ij_pair[id] = 1; // Paired
    }
}

__global__
void nboffset_kernel(int n, int *ij_pair, int *nb_offset) { // Natom threads
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n - 1) {
        // To retrieve ij_pair beginning and ending indices for atom i. Note i == n - 1 is skipped as i < j
        int id_begin = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2; // From (i, i + 1), including itself
        int id_end = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + n - i - 2; // From (i, n - 1), including itself

        int count = 0;
        for (int k = id_begin; k <= id_end; k++) {
            if (ij_pair[k] == 1) count++;
        }
        nb_offset[i + 1] = count;
    }
    if (i == n - 1) { // Use the idle thread to assign terminal values
        nb_offset[0] = 0;
        nb_offset[n] = 0; // TODO Confirm the offset convention of the last element in Fortran code
    }
}

__host__
void nboffset_incr(int n, int *nb_offset) {
    for (int i = 1; i < n + 1; i++) {
        nb_offset[i] += nb_offset[i - 1];
    }
}

__global__
void nblist_kernel(int n, int *ij_pair, int *nb_offset, int *nb_list) { // Natom threads
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n - 1) {
        int id_begin = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2; // From (i, i + 1), including itself
        int id_end = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + n - i - 2; // From (i, n - 1), including itself

        int count = 0;
        for (int k = id_begin; k <= id_end; k++) {
            if (ij_pair[k] == 1) {
                int j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;
                nb_list[nb_offset[i] + count] = j; // Fill j into nb_list from starting pointer of each atom in nb_offset
                count++;
            }
        }
    }
}

extern "C"
void cuda_atomlist_(int *natom, int *nshrt, int *natex, float *acrd, float *cutsa) {
    float *d_nshrt;
    float *d_natex;
    float *d_ijpair;
    float *d_acrd;
    float *d_nboffset;
    float *d_nblist;
    float *nb_offset;
    float *nb_list;

    int devThreadsPerBlock = 8;
    int blocksize = devThreadsPerBlock * devThreadsPerBlock * devThreadsPerBlock;

    int N = *natom * (*natom - 1) / 2;
    cudaErrorCheck(cudaMalloc(&d_nshrt, sizeof(int) * (*natom + 1)));
    cudaErrorCheck(cudaMalloc(&d_natex, sizeof(int) * i10)); // i10 to be determined
    cudaErrorCheck(cudaMalloc(&d_ijpair, sizeof(int) * N));
    cudaErrorCheck(cudaMalloc(&d_acrd, sizeof(float) * 3 * *natom));
    cudaMemcpy(d_nshrt, nshrt, sizeof(int) * (*natom + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(d_natex, natex, sizeof(int) * i10, cudaMemcpyHostToDevice); // i10 to be determined
    cudaMemcpy(d_acrd, acrd, sizeof(float) * 3 * *natom, cudaMemcpyHostToDevice);

    /** A.1 **/
    int nblocks = (N - 1) / blocksize + 1;
    ijpair_kernel<<<nblocks, blocksize>>>(*natom, d_nshrt, d_natex, d_ijpair, *acrd, *cutsa);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();

    /** A.2 **/
    // Can also use thrust::count
    int nblocks2 = (*natom - 1) / blocksize + 1;
    cudaErrorCheck(cudaMalloc(&d_nboffset, sizeof(int) * (*natom + 1)));
    nboffset_kernel<<<nblocks2, blocksize>>>(*natom, d_ijpair, d_nboffset);
    cudaLaunchErrorCheck();
    cudaDeviceSynchronize();

    /** A.3 **/
    nb_offset = (int *) malloc(sizeof(int) * (*natom + 1));
    cudaMemcpy(nb_offset, d_nboffset, sizeof(int) * (*natom + 1), cudaMemcpyDeviceToHost);
    nboffset_incr(*natom, nb_offset);

    /** A.4 **/
    cudaMemcpy(d_nboffset, nb_offset, sizeof(int) * (*natom + 1), cudaMemcpyHostToDevice);
    cudaErrorCheck(cudaMalloc(&d_nblist, sizeof(int) * nb_offset[*natom]));
    nblist_kernel(*natom, d_ijpair, d_nboffset, d_nblist);

    // Copy nblist to host
    nb_list = (int *) malloc(sizeof(int) * nb_offset[*natom]);
    cudaMemcpy(nb_list, d_nblist, sizeof(int) * nb_offset[*natom], cudaMemcpyDeviceToHost);

    cudaFree(d_nshrt);
    cudaFree(d_natex);
    cudaFree(d_ijpair);
    cudaFree(d_acrd);
    cudaFree(d_nboffset);
    cudaFree(d_nblist);
    free(nb_offset);
    //free(nb_list); // Should pass back to Fortran
}
#endif
