//
// CUDA Linear Solvers with the CUSPARSE Library
// Only CSR matrices are supported.
//
// Aaron Luo & Ruxi Qi @ UC Irvine, Jul 2018
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <sys/time.h> // For timing
#if defined(AMBER_PLATFORM_AMD)
#  include <hipblas.h>
#  include <hipsparse.h>
#  include "hip_definitions.h"
#else
#  include <cuda_runtime.h> // CUDA Runtime
#  include <cusparse_v2.h> // Using updated (v2) interfaces for CUBLAS and CUSPARSE
#  include <cublas_v2.h>
#endif
#include "helper_cuda.h" // For error handling and device pickup

//
// Standard conjugate gradient method without preconditioner
// Matrix in the CSR format. Coded with the CUSPARSE library
// Both PBC and Non-PBC can be handled as long as the matrix
// is correctly set up by the CSR routine
//
extern "C" void cusparse_cg_wrapper_(float *x, float *b, int *I, int *J, float *val, int *xmymzm, int *nz_num, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol = *accept; // accept is 1-norm in pbsa

    // input CSR matrix parameters
    int N = *xmymzm;
    int nz = *nz_num;

    // Create CUBLAS context
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus = cublasCreate(&cublasHandle);
    cublasErrorCheck(cublasStatus);

    // Create CUSPARSE context 
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus = cusparseCreate(&cusparseHandle);
    cusparseErrorCheck(cusparseStatus);

    // Set up required arrays on device
    int *d_col, *d_row;
    float *d_val;
    float *d_x;
    float *d_r, *d_p, *d_q;
    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_q, N*sizeof(float)));

    // Get data into working device arrays
    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b, N*sizeof(float), cudaMemcpyHostToDevice);

    const float constONE = 1.0;
    const float constZERO = 0.0;    
#if __CUDACC_VER_MAJOR__ >= 11
    cusparseSpMatDescr_t descr_A;
    cusparseDnVecDescr_t descr_p, descr_q;
    void *spmvBuffer = NULL;
    size_t spmvBufferSize = 0;
    cusparseCreateCsr(&descr_A, N, N, nz, d_row, d_col, d_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
    cusparseCreateDnVec(&descr_p, N, d_p, CUDA_R_32F);
    cusparseCreateDnVec(&descr_q, N, d_q, CUDA_R_32F);
    cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &constONE, descr_A, descr_p, &constZERO, descr_q,
                            CUDA_R_32F, 
#if __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ <= 2
                            CUSPARSE_MV_ALG_DEFAULT,
#else
                            CUSPARSE_SPMV_ALG_DEFAULT,
#endif
                            &spmvBufferSize);
    cudaMalloc(&spmvBuffer, spmvBufferSize);
#else    
    // Description of the A matrix as in Ax = b    
    cusparseMatDescr_t descrA = 0;
    cusparseStatus = cusparseCreateMatDescr(&descrA);
    cusparseErrorCheck(cusparseStatus);
    // Define the properties of the A matrix
    cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO);
#endif
    
    // CG iteration
    float r0, r1; // initial 2-norm
    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r0);
    r1 = r0;
    float b0; // initial 1-norm
    cublasSasum(cublasHandle, N, d_r, 1, &b0);
    float alpha, beta; // CG initial parameters
    float dot, nalpha;
    int k;
    for (k = 1; k < maxiter; k++) { // iteration starts
        // compute beta & p
	beta = r1/r0;
	cublasSscal(cublasHandle, N, &beta, d_p, 1);
	cublasSaxpy(cublasHandle, N, &constONE, d_r, 1, d_p, 1) ;
        // compute Ap & alpha
#if __CUDACC_VER_MAJOR__ >= 11
        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &constONE, descr_A, descr_p, &constZERO, descr_q,
                     CUDA_R_32F, 
#if __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ <= 2
                     CUSPARSE_MV_ALG_DEFAULT, 
#else
                     CUSPARSE_SPMV_ALG_DEFAULT,
#endif
                     spmvBuffer);
#else        
	cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &constONE, descrA, d_val, d_row, d_col, d_p, &constZERO, d_q);
#endif        
	cublasSdot(cublasHandle, N, d_p, 1, d_q, 1, &dot);
	alpha = r1/dot;
        // update x
	cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
	nalpha = -alpha;
        // update r
	cublasSaxpy(cublasHandle, N, &nalpha, d_q, 1, d_r, 1);
        // update norm and check convergence
	r0 = r1;
	cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        //printf("itn %d residue %e tol %e init norm %e\n", k, tol, b0);
        if (r1 < tol*b0){
            break;
        }
    }

    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
    *itn = k;
    *residual = r1;

#if __CUDACC_VER_MAJOR__ >= 11
    cusparseDestroySpMat(descr_A);
    cusparseDestroyDnVec(descr_p);
    cusparseDestroyDnVec(descr_q);
    cudaFree(spmvBuffer);
#endif    
    
    // Destroy contexts
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    // Free device memory
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_q);

    // clean up all state, flush all profile data
    cudaDeviceReset();

}

//
// Bi-conjugate gradient (BICGSTAB) method without preconditioner
// Matrix in the CSR format. Coded with the CUSPARSE library
// Both PBC and Non-PBC can be handled as long as the matrix
// is correctly set up by the CSR routine
//
extern "C" void cusparse_bicg_wrapper_(float *x, float *b, int *I, int *J, float *val, int *xmymzm, int *nz_num, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol2 = *accept * *accept; // accept is 1-norm in pbsa

    // input CSR matrix parameters
    int N = *xmymzm;
    int nz = *nz_num;

    // Create CUBLAS context
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus = cublasCreate(&cublasHandle);
    cublasErrorCheck(cublasStatus);

    // Create CUSPARSE context 
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus = cusparseCreate(&cusparseHandle);
    cusparseErrorCheck(cusparseStatus);

    // Set up required arrays on device
    int *d_col, *d_row;
    float *d_val;
    float *d_x, *d_p, *d_q;
    float *d_r, *d_r0, *d_t;
    cudaErrorCheck(cudaMalloc((void **)&d_col, nz*sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **)&d_val, nz*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_p, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_q, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_r, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_r0, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_t, N*sizeof(float)));

    // Get data into working device arrays, assuming x0 = 0.0
    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, b, N*sizeof(float), cudaMemcpyHostToDevice);

    const float zero = 0.0;
    const float one = 1.0;
#if __CUDACC_VER_MAJOR__ >= 11
    cusparseSpMatDescr_t descr_A;
    cusparseDnVecDescr_t descr_p, descr_q, descr_r, descr_t;
    void *spmvBuffer = NULL;
    size_t spmvBufferSize = 0;
    cusparseCreateCsr(&descr_A, N, N, nz, d_row, d_col, d_val,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
    cusparseCreateDnVec(&descr_p, N, d_p, CUDA_R_32F);
    cusparseCreateDnVec(&descr_q, N, d_q, CUDA_R_32F);
    cusparseCreateDnVec(&descr_r, N, d_r, CUDA_R_32F);
    cusparseCreateDnVec(&descr_t, N, d_t, CUDA_R_32F);    
    cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &one, descr_A, descr_p, &zero, descr_q,
                            CUDA_R_32F, 
#if __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ <= 2
                            CUSPARSE_MV_ALG_DEFAULT,
#else
                            CUSPARSE_SPMV_ALG_DEFAULT,
#endif
                            &spmvBufferSize);
    cudaMalloc(&spmvBuffer, spmvBufferSize);
#else        
    // Description of the A matrix as in Ax = b
    cusparseMatDescr_t descrA = 0;
    cusparseStatus = cusparseCreateMatDescr(&descrA);
    cusparseErrorCheck(cusparseStatus);
    // Define the properties of the A matrix
    cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO);
#endif
    
    // CG iteration
    // initialize norm, p, r
    float nrmr0, nrmr;
    cublasScopy(cublasHandle, N, d_r, 1, d_r0, 1);
    cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
    cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &nrmr0);
    nrmr = nrmr0;
    float rhop, rho = nrmr0;
    int k;
    for (k = 0; k < maxiter; k++) {
        float alpha, beta, omega;
        float temp1, temp2;

        // compute q=Ap & alpha
#if __CUDACC_VER_MAJOR__ >= 11
        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &one, descr_A, descr_p, &zero, descr_q,
                     CUDA_R_32F, 
#if __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ <= 2
                     CUSPARSE_MV_ALG_DEFAULT, 
#else
                     CUSPARSE_SPMV_ALG_DEFAULT,                     
#endif
                     spmvBuffer);
#else                
	cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,N,N,nz,&one,descrA,d_val,d_row,d_col,d_p,&zero,d_q);
#endif        
	cublasSdot(cublasHandle, N, d_r0, 1, d_q, 1, &temp1);
	alpha = rho/temp1;
        // compute s = r - \alpha q, note s is just r
        float nalpha = -alpha;
        cublasSaxpy(cublasHandle, N, &nalpha, d_q, 1, d_r, 1);
        // compute omega = (t^{T} s) / (t^{T} t), t = As
#if __CUDACC_VER_MAJOR__ >= 11
        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &one, descr_A, descr_r, &zero, descr_t,
                     CUDA_R_32F, 
#if __CUDACC_VER_MAJOR__ == 11 && __CUDACC_VER_MINOR__ <= 2
                     CUSPARSE_MV_ALG_DEFAULT, 
#else
                     CUSPARSE_SPMV_ALG_DEFAULT,
#endif
                     spmvBuffer);
#else                
        cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,N,N,nz,&one,descrA,d_val,d_row,d_col,d_r,&zero,d_t);
#endif        
        cublasSdot(cublasHandle, N, d_t, 1, d_r, 1, &temp1);
        cublasSdot(cublasHandle, N, d_t, 1, d_t, 1, &temp2);
        omega = temp1/temp2;
        float nomega = -omega;
        // update x = x + alpha p + omega s
        // update r = s - omega t, t = As
	cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1);
	cublasSaxpy(cublasHandle, N, &omega, d_r, 1, d_x, 1);
        cublasSaxpy(cublasHandle, N, &nomega, d_t, 1, d_r, 1);
        // compute new norm
        cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &nrmr);
        //printf("itn %d residue %e rel res %e accept %e\n", k, nrmr, sqrt(nrmr/nrmr0), sqrt(tol2));
        if (nrmr < tol2*nrmr0){
            break;
        }
        // compute beta = (rho_i/rho_i-1) (alpha/omega)
        rhop = rho;
        cublasSdot(cublasHandle, N, d_r0, 1, d_r, 1, &rho);
        beta = (rho/rhop)*(alpha/omega);
        // compute p = r + beta (p - omega q), q = Ap
        cublasSaxpy(cublasHandle, N, &nomega, d_q, 1, d_p, 1);
        cublasSscal(cublasHandle, N, &beta, d_p, 1);
        cublasSaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1);
    }

    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
    *itn = k;
    *residual = sqrt(nrmr); // sqrt it to be consistent with 1 norm in pbsa

#if __CUDACC_VER_MAJOR__ >= 11
    cusparseDestroySpMat(descr_A);
    cusparseDestroyDnVec(descr_p);
    cusparseDestroyDnVec(descr_q);
    cusparseDestroyDnVec(descr_r);
    cusparseDestroyDnVec(descr_t);    
    cudaFree(spmvBuffer);
#endif    
    
    // Destroy contexts
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    // Free device memory
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_p);
    cudaFree(d_q);
    cudaFree(d_r);
    cudaFree(d_r0);
    cudaFree(d_t);

    // clean up all state, flush all profile data
    cudaDeviceReset();

}
