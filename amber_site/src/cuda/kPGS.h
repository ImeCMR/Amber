#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
{
  const int GRADSUMGRID     = 8;
  const PMEFloat ONEHALF    = 0.5;
  const PMEFloat ONETHIRD   = (1.0 / 3.0);
  const PMEFloat ONE        = 1.0;
  const PMEFloat TWO        = 2.0;
  const PMEFloat THREE      = 3.0;
  struct GradSumAtomData
  {
    PMEFloat charge;
    int ix;
    int iy;
    int iz;
    PMEFloat tx[4];
    PMEFloat ty[4];
    PMEFloat tz0;
    PMEFloat tz1;
    PMEFloat tz2;
    PMEFloat tz3;
    PMEFloat dtx[4];
    PMEFloat dty[4];
    PMEFloat dtz0;
    PMEFloat dtz1;
    PMEFloat dtz2;
    PMEFloat dtz3;
    PMEFloat fx;
    PMEFloat fy;
    PMEFloat fz;
#if defined(SC_REGION_1) || defined(SC_REGION_2)
    unsigned int TIRegion;
#endif
  };

  __shared__ VOLATILE GradSumAtomData sAtom[GRADSUMLOADSIZE];
#ifdef PME_VIRIAL
  __shared__ PMEFloat sRecipf[9];
  if (threadIdx.x < 9) {
    sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
  }
  __syncthreads();
#endif

  // Determine warp constants
  unsigned int tgx = threadIdx.x & (GRADSUMGRID - 1);

  // Determine grid offsets
  const int tOffsetX = tgx & 0x03;
  const int tOffsetY = (tgx & 0x07) >> 2;
  const int iOffsetX = tOffsetX;
  const int iOffsetY = tOffsetY;
  unsigned int pos   = (blockIdx.x + offset) * GRADSUMLOADSIZE;

  // Read batch of atoms and procedurally generate spline weights
  unsigned int maxatom = min(pos + GRADSUMLOADSIZE, cSim.atoms);
  unsigned int pos1    = pos + threadIdx.x;
  if (pos1 < maxatom) {
    PMEFloat charge = cSim.pAtomChargeSP[pos1];
#if defined(SC_REGION_1) || defined(SC_REGION_2)
    unsigned int TIRegion = cSim.pImageTIRegion[pos1];
#endif
    PMEFloat fx = cSim.pFractX[pos1];
    PMEFloat fy = cSim.pFractY[pos1];
    PMEFloat fz = cSim.pFractZ[pos1];
    int ix      = int(fx);
    int iy      = int(fy);
    int iz      = int(fz);
    fx         -= ix;
    fy         -= iy;
    fz         -= iz;
    ix         -= cSim.orderMinusOne;
    iy         -= cSim.orderMinusOne;
    iz         -= cSim.orderMinusOne;
    if (ix < 0) {
      ix += cSim.nfft1;
    }
    if (iy < 0) {
      iy += cSim.nfft2;
    }
    if (iz < 0) {
      iz += cSim.nfft3;
    }

    // Order 2 B-spline accumulation
    PMEFloat4 tx;
    tx.x = ONE - fx;
    tx.y = fx;

    // Order 3 B-spline accumulation
    tx.z = ONEHALF * fx * tx.y;
    tx.y = ONEHALF * ((fx + ONE) * tx.x + (TWO - fx) * tx.y);
    tx.x = ONEHALF * (ONE - fx)  * tx.x;
    sAtom[threadIdx.x].dtx[0] = -tx.x;
    sAtom[threadIdx.x].dtx[1] =  tx.x - tx.y;
    sAtom[threadIdx.x].dtx[2] =  tx.y - tx.z;
    sAtom[threadIdx.x].dtx[3] =  tx.z;

    // Order 4 B-spline accumulation
    tx.w = ONETHIRD * fx * tx.z;
    tx.z = ONETHIRD * ((fx + ONE) * tx.y + (THREE - fx) * tx.z);
    tx.y = ONETHIRD * ((fx + TWO) * tx.x + (TWO - fx) * tx.y);
    tx.x = ONETHIRD *  (ONE - fx) * tx.x;

    // Order 2 B-spline accumulation
    PMEFloat4 ty;
    ty.x = ONE - fy;
    ty.y = fy;

    // Order 3 B-spline accumulation
    ty.z = ONEHALF * fy * ty.y;
    ty.y = ONEHALF * ((fy + ONE) * ty.x + (TWO - fy) * ty.y);
    ty.x = ONEHALF * (ONE - fy)  * ty.x;
    sAtom[threadIdx.x].dty[0] = -ty.x;
    sAtom[threadIdx.x].dty[1] =  ty.x - ty.y;
    sAtom[threadIdx.x].dty[2] =  ty.y - ty.z;
    sAtom[threadIdx.x].dty[3] =  ty.z;

    // Order 4 B-spline accumulation
    ty.w = ONETHIRD * fy * ty.z;
    ty.z = ONETHIRD * ((fy + ONE) * ty.y + (THREE - fy) * ty.z);
    ty.y = ONETHIRD * ((fy + TWO) * ty.x + (TWO - fy) * ty.y);
    ty.x = ONETHIRD *  (ONE - fy) * ty.x;

    // Order 2 B-spline accumulation
    PMEFloat4 tz;
    tz.x = ONE - fz;
    tz.y = fz;

    // Order 3
    tz.z = ONEHALF * fz * tz.y;
    tz.y = ONEHALF * ((fz + ONE) * tz.x + (TWO - fz) * tz.y);
    tz.x = ONEHALF * (ONE - fz)  * tz.x;
    sAtom[threadIdx.x].dtz0 = -tz.x;
    sAtom[threadIdx.x].dtz1 =  tz.x - tz.y;
    sAtom[threadIdx.x].dtz2 =  tz.y - tz.z;
    sAtom[threadIdx.x].dtz3 =  tz.z;

    // Order 4
    tz.w = ONETHIRD * fz * tz.z;
    tz.z = ONETHIRD * ((fz + ONE) * tz.y + (THREE - fz) * tz.z);
    tz.y = ONETHIRD * ((fz + TWO) * tz.x + (TWO - fz) * tz.y);
    tz.x = ONETHIRD *  (ONE - fz) * tz.x;

    sAtom[threadIdx.x].charge = charge;
#if defined(SC_REGION_1) || defined(SC_REGION_2)
    sAtom[threadIdx.x].TIRegion = TIRegion;
#endif
    sAtom[threadIdx.x].ix    = ix;
    sAtom[threadIdx.x].iy    = iy;
    sAtom[threadIdx.x].iz    = iz;
    sAtom[threadIdx.x].tx[0] = tx.x;
    sAtom[threadIdx.x].tx[1] = tx.y;
    sAtom[threadIdx.x].tx[2] = tx.z;
    sAtom[threadIdx.x].tx[3] = tx.w;
    sAtom[threadIdx.x].ty[0] = ty.x;
    sAtom[threadIdx.x].ty[1] = ty.y;
    sAtom[threadIdx.x].ty[2] = ty.z;
    sAtom[threadIdx.x].ty[3] = ty.w;
    sAtom[threadIdx.x].tz0   = tz.x;
    sAtom[threadIdx.x].tz1   = tz.y;
    sAtom[threadIdx.x].tz2   = tz.z;
    sAtom[threadIdx.x].tz3   = tz.w;
  }
  __syncthreads();

  // Process batch
  pos1 = threadIdx.x / GRADSUMGRID;
  unsigned int lastAtom = min(GRADSUMLOADSIZE, cSim.atoms - pos);
  PMEMask mask1 = __BALLOT(WARP_MASK, pos1 < lastAtom);
  while (pos1 < lastAtom) {

    // Calculate values
    int ix  = sAtom[pos1].ix + iOffsetX;
    int iy0 = sAtom[pos1].iy + iOffsetY;
    int iy1 = iy0 + 2;
    int iz0 = sAtom[pos1].iz;
    int iz1 = iz0 + 1;
    int iz2 = iz0 + 2;
    int iz3 = iz0 + 3;

    // Insure coordinates stay in bounds
    ix  -= (ix >= cSim.nfft1)  * cSim.nfft1;
    iy0 -= (iy0 >= cSim.nfft2) * cSim.nfft2;
    iy1 -= (iy1 >= cSim.nfft2) * cSim.nfft2;
    iz0 -= (iz0 >= cSim.nfft3) * cSim.nfft3;
    iz1 -= (iz1 >= cSim.nfft3) * cSim.nfft3;
    iz2 -= (iz2 >= cSim.nfft3) * cSim.nfft3;
    iz3 -= (iz3 >= cSim.nfft3) * cSim.nfft3;

    // Calculate interpolation values and destinations
#ifdef use_DPFP
    int2 i2term0 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz0 * cSim.nfft2 + iy0) * cSim.nfft1 + ix));
    int2 i2term1 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz1 * cSim.nfft2 + iy0) * cSim.nfft1 + ix));
    int2 i2term2 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz2 * cSim.nfft2 + iy0) * cSim.nfft1 + ix));
    int2 i2term3 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz3 * cSim.nfft2 + iy0) * cSim.nfft1 + ix));
    int2 i2term4 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz0 * cSim.nfft2 + iy1) * cSim.nfft1 + ix));
    int2 i2term5 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz1 * cSim.nfft2 + iy1) * cSim.nfft1 + ix));
    int2 i2term6 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz2 * cSim.nfft2 + iy1) * cSim.nfft1 + ix));
    int2 i2term7 = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz3 * cSim.nfft2 + iy1) * cSim.nfft1 + ix));
    PMEFloat qterm0 = __hiloint2double(i2term0.y, i2term0.x);
    PMEFloat qterm1 = __hiloint2double(i2term1.y, i2term1.x);
    PMEFloat qterm2 = __hiloint2double(i2term2.y, i2term2.x);
    PMEFloat qterm3 = __hiloint2double(i2term3.y, i2term3.x);
    PMEFloat qterm4 = __hiloint2double(i2term4.y, i2term4.x);
    PMEFloat qterm5 = __hiloint2double(i2term5.y, i2term5.x);
    PMEFloat qterm6 = __hiloint2double(i2term6.y, i2term6.x);
    PMEFloat qterm7 = __hiloint2double(i2term7.y, i2term7.x);
#else
    PMEFloat qterm0 = tex1Dfetch<float>(cSim.texXYZ_q, (iz0 * cSim.nfft2 + iy0) * cSim.nfft1 + ix);
    PMEFloat qterm1 = tex1Dfetch<float>(cSim.texXYZ_q, (iz1 * cSim.nfft2 + iy0) * cSim.nfft1 + ix);
    PMEFloat qterm2 = tex1Dfetch<float>(cSim.texXYZ_q, (iz2 * cSim.nfft2 + iy0) * cSim.nfft1 + ix);
    PMEFloat qterm3 = tex1Dfetch<float>(cSim.texXYZ_q, (iz3 * cSim.nfft2 + iy0) * cSim.nfft1 + ix);
    PMEFloat qterm4 = tex1Dfetch<float>(cSim.texXYZ_q, (iz0 * cSim.nfft2 + iy1) * cSim.nfft1 + ix);
    PMEFloat qterm5 = tex1Dfetch<float>(cSim.texXYZ_q, (iz1 * cSim.nfft2 + iy1) * cSim.nfft1 + ix);
    PMEFloat qterm6 = tex1Dfetch<float>(cSim.texXYZ_q, (iz2 * cSim.nfft2 + iy1) * cSim.nfft1 + ix);
    PMEFloat qterm7 = tex1Dfetch<float>(cSim.texXYZ_q, (iz3 * cSim.nfft2 + iy1) * cSim.nfft1 + ix);
#endif
    PMEFloat wx       = sAtom[pos1].tx[tOffsetX];
    PMEFloat wy0      = sAtom[pos1].ty[tOffsetY];
    PMEFloat wy1      = sAtom[pos1].ty[tOffsetY + 2];
    PMEFloat dwx      = sAtom[pos1].dtx[tOffsetX];
    PMEFloat dwy0     = sAtom[pos1].dty[tOffsetY];
    PMEFloat dwy1     = sAtom[pos1].dty[tOffsetY + 2];
    PMEFloat tz0      = sAtom[pos1].tz0;
    PMEFloat dtz0     = sAtom[pos1].dtz0;
    PMEFloat tz1      = sAtom[pos1].tz1;
    PMEFloat dtz1     = sAtom[pos1].dtz1;
    PMEFloat tz2      = sAtom[pos1].tz2;
    PMEFloat dtz2     = sAtom[pos1].dtz2;
    PMEFloat tz3      = sAtom[pos1].tz3;
    PMEFloat dtz3     = sAtom[pos1].dtz3;
    PMEFloat dwxwy0   = dwx * wy0;
    PMEFloat wxdwy0   = wx  * dwy0;
    PMEFloat wxwy0    = wx  * wy0;
    PMEFloat qterm0a  = qterm0 * tz0;
    PMEFloat qterm1a  = qterm1 * tz1;
    PMEFloat qterm2a  = qterm2 * tz2;
    PMEFloat qterm3a  = qterm3 * tz3;
    PMEFloat f1       = -qterm0a * dwxwy0;
    PMEFloat f2       = -qterm0a * wxdwy0;
    PMEFloat f3       = -qterm0  * wxwy0  * dtz0;
    f1 -=  qterm1a * dwxwy0;
    f2 -=  qterm1a * wxdwy0;
    f3 -=  qterm1  * wxwy0  * dtz1;
    f1 -=  qterm2a * dwxwy0;
    f2 -=  qterm2a * wxdwy0;
    f3 -=  qterm2  * wxwy0  * dtz2;
    f1 -=  qterm3a * dwxwy0;
    f2 -=  qterm3a * wxdwy0;
    f3 -=  qterm3  * wxwy0  * dtz3;
    PMEFloat dwxwy1   = dwx * wy1;
    PMEFloat wxdwy1   = wx  * dwy1;
    PMEFloat wxwy1    = wx  * wy1;
    PMEFloat qterm4a  = qterm4 * tz0;
    PMEFloat qterm5a  = qterm5 * tz1;
    PMEFloat qterm6a  = qterm6 * tz2;
    PMEFloat qterm7a  = qterm7 * tz3;
    f1 -=  qterm4a * dwxwy1;
    f2 -=  qterm4a * wxdwy1;
    f3 -=  qterm4  * wxwy1  * dtz0;
    f1 -=  qterm5a * dwxwy1;
    f2 -=  qterm5a * wxdwy1;
    f3 -=  qterm5  * wxwy1  * dtz1;
    f1 -=  qterm6a * dwxwy1;
    f2 -=  qterm6a * wxdwy1;
    f3 -=  qterm6  * wxwy1  * dtz2;
    f1 -=  qterm7a * dwxwy1;
    f2 -=  qterm7a * wxdwy1;
    f3 -=  qterm7  * wxwy1  * dtz3;
    unsigned int tx = threadIdx.x & (warpSize - 1);
    f1 += __SHFL(mask1, f1, tx ^ 4);
    f2 += __SHFL(mask1, f2, tx ^ 4);
    f3 += __SHFL(mask1, f3, tx ^ 4);
    f1 += __SHFL(mask1, f1, tx ^ 2);
    f2 += __SHFL(mask1, f2, tx ^ 2);
    f3 += __SHFL(mask1, f3, tx ^ 2);
    f1 += __SHFL(mask1, f1, tx ^ 1);
    f2 += __SHFL(mask1, f2, tx ^ 1);
    f3 += __SHFL(mask1, f3, tx ^ 1);
    if (tgx == 0) {
      sAtom[pos1].fx = f1;
      sAtom[pos1].fy = f2;
      sAtom[pos1].fz = f3;
    }
    pos1 += GRADSUMTHREADS / GRADSUMGRID;
    mask1 = __BALLOT(mask1, pos1 < lastAtom);
  }
  __syncthreads();

  // Set up one store for results
  if (threadIdx.x < lastAtom) {
#ifdef use_DPFP
    PMEFloat charge = sAtom[threadIdx.x].charge * FORCESCALE;
#  ifdef SC_REGION_1
    // For Alchemical Free Energy simulations, we need to mask contributions from the
    // other region.  We can do this by zeroing the charges from the other region.
    if (sAtom[threadIdx.x].TIRegion & 4) {
      charge = (PMEFloat)0.0;
    }
#  endif
#  ifdef SC_REGION_2
    if (sAtom[threadIdx.x].TIRegion & 2) {
      charge = (PMEFloat)0.0;
    }
#  endif
#else  // use_DPFP
    PMEFloat charge = sAtom[threadIdx.x].charge * FORCESCALEF;
#endif // use_DPFP
    PMEFloat fx = cSim.nfft1 * charge * sAtom[threadIdx.x].fx;
    PMEFloat fy = cSim.nfft2 * charge * sAtom[threadIdx.x].fy;
    PMEFloat fz = cSim.nfft3 * charge * sAtom[threadIdx.x].fz;
    unsigned int pos1 = pos + threadIdx.x;
#ifdef PME_VIRIAL
    fz = sRecipf[6] * fx + sRecipf[7] * fy + sRecipf[8] * fz;
    fy = sRecipf[3] * fx + sRecipf[4] * fy;
    fx = sRecipf[0] * fx;
#else
    fz = cSim.recipf[2][0] * fx + cSim.recipf[2][1] * fy + cSim.recipf[2][2] * fz;
    fy = cSim.recipf[1][0] * fx + cSim.recipf[1][1] * fy;
    fx = cSim.recipf[0][0] * fx;
#endif
#ifdef use_DPFP
#  if !defined(SC_REGION_1) && !defined(SC_REGION_2)
    // Only do the regular forces if we're not doing TI
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(llrint(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(llrint(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(llrint(fz)));
#  endif
#  ifdef SC_REGION_1
    fx = cSim.AFElambda[0] * fx;
    fy = cSim.AFElambda[0] * fy;
    fz = cSim.AFElambda[0] * fz;
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(llrint(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(llrint(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(llrint(fz)));
#  endif
#  ifdef SC_REGION_2
    fx = cSim.AFElambda[1] * fx;
    fy = cSim.AFElambda[1] * fy;
    fz = cSim.AFElambda[1] * fz;
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(llrint(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(llrint(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(llrint(fz)));
#  endif
#else  // use_DPFP
#  if !defined(SC_REGION_1) && !defined(SC_REGION_2)
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(fast_llrintf(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(fast_llrintf(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(fast_llrintf(fz)));
#  endif
#  ifdef SC_REGION_1
    fx = cSim.AFElambda[0] * fx;
    fy = cSim.AFElambda[0] * fy;
    fz = cSim.AFElambda[0] * fz;
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(fast_llrintf(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(fast_llrintf(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(fast_llrintf(fz)));
#  endif
#  ifdef SC_REGION_2
    fx = cSim.AFElambda[1] * fx;
    fy = cSim.AFElambda[1] * fy;
    fz = cSim.AFElambda[1] * fz;
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[pos1],
              llitoulli(fast_llrintf(fx)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[pos1],
              llitoulli(fast_llrintf(fy)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[pos1],
              llitoulli(fast_llrintf(fz)));
#  endif
#endif // use_DPFP
  }
}
#undef VOLATILE
