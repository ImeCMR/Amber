// 
// CUDA Linear System Solvers with the CUSP Library
//
// Aaron Luo & Ruxi Qi @ UC Irvine, Jul 2019
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <sys/time.h> // For timing
#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime_api.h>
#  include "hip_definitions.h"
#else 
#include <cuda_runtime.h> // CUDA Runtime
#endif
#include "helper_cuda.h" // For error handling and device pickup

// CUSP libraries
#include <cusp/csr_matrix.h>
#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/precond/diagonal.h>
#include <cusp/print.h>
#include <thrust/device_ptr.h>

//
// Conjugate gradient method with and without preconditioner
// Matrix input in the CSR format, but can be converted
// to DIA for the preconditioned method. Coded with the CUSP library.
// Both PBC and Non-PBC can be handled as long as the matrix
// is correctly set up by the CSR routines.
//
#ifdef CSR
extern "C" void cusp_cg_wrapper_(float *x, float *b, int *I, int *J, float *val, int *xmymzm, int *nz_num, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol = *accept; // accept is 1-norm in pbsa and cusp also uses 1-norm

    // Input CSR matrix parameters
    int N = *xmymzm;
    int nz = *nz_num;

    // Initialize CSR vectors on device (otherwise cannot assign values to A)
    thrust::device_vector<int> d_row(I, I + N + 1), d_col(J, J + nz);
    thrust::device_vector<float> d_val(val, val + nz);

    // Initialize cusp matrix A on device
    cusp::csr_matrix<int, float, cusp::device_memory> A(N, N, nz);
    A.row_offsets = d_row;
    A.column_indices = d_col;
    A.values = d_val;

    typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArray1dView;
#endif
#ifdef DIA
extern "C" void cusp_cg_wrapper_(float *x, float *b, int *I, float *val, int *xm, int *ym, int *zm, int *nz_num, int *nbnd, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol = *accept; // accept is 1-norm in pbsa and cusp also uses 1-norm
    const int nz = *nz_num;
    const int nb = *nbnd;

    // Input DIA matrix parameters
    int xmym = *xm * *ym;
    int N = xmym * *zm;
    //int nz = N + 2 * (N - 1 + N - *xm + N - xmym); ! Different for PBC, passed by Fortran

    // Initialize matrix values on device
    // This implementation is to use off and val in the raw memory
    int *d_off;
    cudaErrorCheck(cudaMalloc((void **)&d_off, nb * sizeof(int)));
    cudaMemcpy(d_off, I, nb * sizeof(int), cudaMemcpyHostToDevice);
    thrust::device_ptr<int> wrap_d_off(d_off);
    float *d_val;
    cudaErrorCheck(cudaMalloc((void **)&d_val, nb * N * sizeof(float)));
    cudaMemcpy(d_val, val, nb * N * sizeof(float), cudaMemcpyHostToDevice);
    thrust::device_ptr<float> wrap_d_val(d_val);

    typedef typename cusp::array1d_view< thrust::device_ptr<int> > DeviceIndexArray1dView;
    typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArray1dView;
    typedef cusp::array2d_view<DeviceValueArray1dView, cusp::column_major> DeviceValueArray2dView;
    DeviceIndexArray1dView d1_off(wrap_d_off, wrap_d_off + nb);
    DeviceValueArray1dView d1_val(wrap_d_val, wrap_d_val + nb * N);
    DeviceValueArray2dView d2_val(N, nb, N, d1_val);

    // Initialize cusp matrix A on device
    cusp::dia_matrix<int, float, cusp::device_memory> cg_A(N, N, nz, nb);
    cg_A.diagonal_offsets = d1_off;
    cg_A.values = d2_val;
#endif

    // Set up vectors x and b on device
    // This implementation is to use x and b in the raw memory
    float *d_x, *d_b;
    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_b, N*sizeof(float)));
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, N*sizeof(float), cudaMemcpyHostToDevice);
    // Raw device pointers are wrapped before using
    thrust::device_ptr<float> wrap_d_x(d_x);
    thrust::device_ptr<float> wrap_d_b(d_b);
    // Use array1d_view to wrap the individual arrays
    DeviceValueArray1dView cg_x (wrap_d_x, wrap_d_x + N);
    DeviceValueArray1dView cg_b (wrap_d_b, wrap_d_b + N);

    // Iteration control
    cusp::monitor<float> monitor(cg_b, maxiter, tol, 0, false);

    // Solve Ax = b
#ifdef CSR
    //cusp::krylov::cg(A, cg_x, cg_b, monitor);
    cusp::precond::diagonal<float, cusp::device_memory> M(A);
    cusp::krylov::cg(A, cg_x, cg_b, monitor, M);
#endif
#ifdef DIA
    cusp::precond::diagonal<float, cusp::device_memory> M(cg_A);
    cusp::krylov::cg(cg_A, cg_x, cg_b, monitor, M);
#endif

    // Returning
    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
    *itn = monitor.iteration_count();
    *residual = monitor.residual_norm();

}

//
// Bi-conjugate gradient (BICGSTAB) method with and without preconditioner
// Matrix input in the CSR format, but can be converted
// to DIA for the preconditioned method. Coded with the CUSP library.
// Both PBC and Non-PBC can be handled as long as the matrix
// is correctly set up by the CSR routine
//
#ifdef CSR
extern "C" void cusp_bicg_wrapper_(float *x, float *b, int *I, int *J, float *val, int *xmymzm, int *nz_num, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol = *accept; // accept is 1-norm in pbsa but we are using cusp's 2-norm monitor as in all bicg calls

    // Input CSR matrix parameters
    int N = *xmymzm;
    int nz = *nz_num;

    // Initialize CSR vectors on device (or cannot assign values to A)
    thrust::device_vector<int> d_row(I, I + N + 1), d_col(J, J + nz);
    thrust::device_vector<float> d_val(val, val + nz);

    // Initialize cusp matrix A on device
    cusp::csr_matrix<int, float, cusp::device_memory> cg_A(N, N, nz);
    cg_A.row_offsets = d_row;
    cg_A.column_indices = d_col;
    cg_A.values = d_val;

    typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArray1dView;
#endif
#ifdef DIA
extern "C" void cusp_bicg_wrapper_(float *x, float *b, int *I, float *val, int *xm, int *ym, int *zm, int *maxitn, float *accept, int *itn, float *residual)
{
    const int maxiter = *maxitn;
    const float tol = *accept; // accept is 1-norm in pbsa but we are using cusp's 2-norm monitor as in all bicg calls

    // Input DIA matrix parameters
    int xmym = *xm * *ym;
    int N = xmym * *zm;
    int nz = N + 2 * (N - 1 + N - *xm + N - xmym); // Be careful when implementing bicg solver for PBC

    // Initialize matrix values on device
    // This implementation is to use off and val in the raw memory
    int *d_off;
    cudaErrorCheck(cudaMalloc((void **)&d_off, 7*sizeof(int)));
    cudaMemcpy(d_off, I, 7*sizeof(int), cudaMemcpyHostToDevice);
    thrust::device_ptr<int> wrap_d_off(d_off);
    float *d_val;
    cudaErrorCheck(cudaMalloc((void **)&d_val, 7*N*sizeof(float)));
    cudaMemcpy(d_val, val, 7*N*sizeof(float), cudaMemcpyHostToDevice);
    thrust::device_ptr<float> wrap_d_val(d_val);

    typedef typename cusp::array1d_view< thrust::device_ptr<int> > DeviceIndexArray1dView;
    typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArray1dView;
    typedef cusp::array2d_view<DeviceValueArray1dView, cusp::column_major> DeviceValueArray2dView;
    DeviceIndexArray1dView d1_off(wrap_d_off, wrap_d_off + 7);
    DeviceValueArray1dView d1_val(wrap_d_val, wrap_d_val + 7*N);
    DeviceValueArray2dView d2_val(N, 7, N, d1_val);

    // Initialize cusp matrix A on device
    cusp::dia_matrix<int, float, cusp::device_memory> cg_A(N, N, nz, 7);
    cg_A.diagonal_offsets = d1_off;
    cg_A.values = d2_val;
#endif

    // Set up vectors x and b on device
    // This implementation is to use x and b in the raw memory
    float *d_x, *d_b;
    cudaErrorCheck(cudaMalloc((void **)&d_x, N*sizeof(float)));
    cudaErrorCheck(cudaMalloc((void **)&d_b, N*sizeof(float)));
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, N*sizeof(float), cudaMemcpyHostToDevice);
    // Raw device pointers need to be wrapped before using
    // Use array1d_view to wrap the individual arrays
    thrust::device_ptr<float> wrap_d_x(d_x);
    thrust::device_ptr<float> wrap_d_b(d_b);
    DeviceValueArray1dView cg_x (wrap_d_x, wrap_d_x + N);
    DeviceValueArray1dView cg_b (wrap_d_b, wrap_d_b + N);

    // Iteration control
    cusp::monitor_l2<float> monitor(cg_b, maxiter, tol, 0, false);

    // Set preconditioner (identity)
    //cusp::identity_operator<float, cusp::device_memory> M(cg_A.num_rows, cg_A.num_rows);
    cusp::precond::diagonal<float, cusp::device_memory> M(cg_A);

    // Solve Ax = b
    cusp::krylov::bicgstab(cg_A, cg_x, cg_b, monitor, M);

    // Returning
    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
    *itn = monitor.iteration_count();
    *residual = monitor.residual_norm();

}
