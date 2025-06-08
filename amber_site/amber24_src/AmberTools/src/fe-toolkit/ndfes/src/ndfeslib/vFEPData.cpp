#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "vFEPData.hpp"
#include "MeshUtils.hpp"
#include "PeriodicUtils.hpp"

#include "GetBeta.hpp"

extern "C"
{
  extern void
  dgemv_( char const * TRANS,
	  int const * M,
	  int const * N,
	  double const * alp,
	  double const * A,
	  int const * LDA,
	  double const * x,
	  int const * INCX,
	  double const * bet,
	  double * Y,
	  int const * INCY );

  extern double
  ddot_( int const * N,
	 double const * DX,
	 int const * INCX,
	 double const * DY,
	 int const * INCY );
}



namespace ndfes
{

  void bspline_one_pass( double * c, double const w, int const n );

  void bspline_eval( double const w, int const order, double * array );
  
  void bspline_aperiodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );

  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
}



void ndfes::bspline_one_pass( double * c, double const w, int const n )
{
  int nm1 = n-1;
  double const div = 1. / nm1;
  c[nm1] = div * w * c[nm1 - 1];
  for ( int j=1; j<nm1; ++j )
    c[nm1 - j] = div * ((w + j) * c[nm1 - j - 1] + (n - j - w) * c[nm1 - j]);
  c[0] = div * (1 - w) * c[0];
}


void ndfes::bspline_eval( double const w, int const order, double * array )
{

  array[0] = 1. - w;
  array[1] = w;
  if (order > 2)
    {
      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
      if (order > 3)
	{
	  // One pass to order 4:         
	  double const div = 1./3.;
	  array[3] = div * w * array[2];
	  array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
	  array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
	  array[0] = div * (1. - w) * array[0];
	  // and the rest
	  for ( int k = 5; k < order+1; ++k )
	    ndfes::bspline_one_pass(array, w, k);
	};
    };
}

void ndfes::bspline_aperiodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  //double const ix = x;

  //std::printf("x=%20.14f ",x);
  
  x = (x/lenx)*(nx-1) - (order%2)*0.5;
  int ilo = std::floor(x);

  //std::printf("%20.14f %5i\n",x,ilo);

  
  //std::printf("%23.14e %20.14e ilo %i\n",x-ilo,x,ilo);
  // maybe this should really be
  // if ( ilo + 1 - order/2 == -1 )
  if ( ilo < nx/2-1 )
    {
      if ( x-ilo > 1. - 5.e-14 )
	{
	  //std::printf("boundary case increase ilo %20.10e %i\n",x,ilo);
	  ilo += 1;
	  x = ilo;
	}
    }
  
  else if ( ilo + 1 - (order+1)/2 + (order-1) == nx )
    {
      if ( x-ilo < 5.e-14 )
	{
	  //std::printf("boundary case decrease ilo %20.10e %i\n",x,ilo);
	  ilo -= 1;
	  x = ilo+1;
	}
    }
  
  //else if ( ilo > order/2 and x-ilo
  //std::printf("x,ilo: %26.17e %i\n",x,ilo);
  if ( order == 1 )
    {
      w[0] = 1.;
      gidx[0] = (x-ilo < 0.5) ? ilo : ilo+1;
      //std::printf("x %20.10e %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i (ilo=%i) = %20.10e\n",ix,x-ilo,lenx,nx,order,0,gidx[0],ilo,w[0]);
      //std::cout.flush();

    }
  else
    {
      //std::printf("bspline frac %20.10e %i\n",x-ilo,ilo);
      ndfes::bspline_eval(x-ilo,order,w);
      ilo += 1 - order/2;
      for ( int b=0; b<order; ++b, ++ilo )
	{
	  gidx[b] = ilo;
	  //std::printf("x %20.10e lenx %20.10e nx %4i order %4i gidx[%i] %i = %20.10e\n",x,lenx,nx,order,b,ilo,w[b]);
	  if ( ilo < 0 or ilo > nx-1 )
	    {
	      std::cerr << "ndfes::bspline_aperiodic node index out of bounds " << ilo << " " << nx << "\n";
	      std::exit(1);
	    }
      //gidx[b] = std::min(nx-1,std::max(0,ilo));
      //gidx[b] = ((ilo+nx)%nx+nx)%nx;
      //if ( ilo < 0 or ilo >= nx ) w[b]=0.;
      //std::printf("%4i(%4i) ",gidx[b],nx);
	}
    }
  //std::printf("\n");
}


void ndfes::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  x = (x/lenx)*nx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  if ( order == 1 )
    {
      w[0] = 1.;
      gidx[0] = (x < 0.5) ? ilo : ilo+1;
    }
  else
    {
      ndfes::bspline_eval(x,order,w);
      // the first modulus takes us to [-nx+1,nx-1]
      // then we add nx, so we are in the range [0,2*nx-1]
      // the second modulus then puts us in [0,nx-1]
      // Each ilo in the loop uses a modulus, so we don't need
      // the second modulus on the first calculation of ilo
      ilo   =  (ilo+1-order/2)%nx+nx;
      for ( int b=0; b<order; ++b, ++ilo )
	{
	  gidx[b] = ilo%nx;
	}
    }
  //std::printf("\n");
}


static void GaussLegendreRule
( double const x1,
  double const x2,
  std::size_t const n,
  double * rVecX, 
  double * rVecW )
{
  std::size_t const MAXIT=10;
  double const EPS = 3.0e-14;
  std::size_t const m = (n+1)/2;
  double const xm = 0.5 * ( x2 + x1 );
  double const xl = 0.5 * ( x2 - x1 );  
  double const C1 = M_PI / (n+0.5);
  
  std::vector<double> z(m);
  std::vector<double> z1(m);
  std::vector<double> p1(m);
  std::vector<double> p2(m);
  std::vector<double> p3(m);
  std::vector<double> pp(m);
  std::vector<bool> unfinished(m);

  for ( std::size_t i=0; i<m; ++i )
    z[i] = std::cos( C1 * ( i + 0.75 ) );
  for ( std::size_t i=0; i<m; ++i )
    unfinished[i] = true;
  for ( std::size_t its=0; its < MAXIT; ++its )
    {
      for ( std::size_t i=0; i < m; ++i )
        {
          if ( unfinished[i] )
            {
              p1[i] = 1.0;
              p2[i] = 0.0;
            };
        };
      for ( std::size_t j=0; j < n; ++j )
        {
          for ( std::size_t i=0; i < m; ++i )
            {
              if ( unfinished[i] )
                {
                  p3[i]=p2[i];
                  p2[i]=p1[i];
                  p1[i]=((2.0*j+1.0)*z[i]*p2[i]-j*p3[i])/(j+1.0);
                };
            };
        };

      bool HasUnfinished = false;
      for ( std::size_t i=0; i < m; ++i )
        {
          if ( unfinished[i] )
            {
              pp[i]=n*(z[i]*p1[i]-p2[i])/(z[i]*z[i]-1.0);
              z1[i]=z[i];
              z[i]=z1[i]-p1[i]/pp[i];
              unfinished[i] = ( std::abs(z[i]-z1[i]) > EPS );
              if ( unfinished[i] )
                HasUnfinished = true;
            };
        };
      if ( ! HasUnfinished ){ break; };
      
      if ( its == MAXIT )
        std::cerr << "Too many iterations in GaussLegendreRule" << std::endl;
    };

  for ( std::size_t i=0; i < m; ++i )
    rVecX[i] = xm-xl*z[i];
  for ( std::size_t i=0; i < n-m; ++i )
    rVecX[n-i-1] = xm+xl*z[i];
  for ( std::size_t i=0; i < m; ++i )
    rVecW[i] = 2.0*xl/((1.-z[i]*z[i])*pp[i]*pp[i]);
  for ( std::size_t i=0; i < n-m; ++i )
    rVecW[n-i-1] = rVecW[i];
}



// ndfes::vFEPBin::vFEPBin()
//   : nbspl(0),
//     nsubbins(0),
//     nstate(0),
//     c2sc(NULL),
//     sbin(NULL),
//     ndim(0),
//     ncorner(0),
//     subbinsize(0)
// {
// }




// ndfes::vFEPBin::vFEPBin
// ( int const inbspl,
//   int const insubbins,
//   int const instate,
//   double const * ic2sc,
//   ndfes::SpatialBin const * isbin,
//   std::vector<int> const & glbcidxs )
//   : nbspl(inbspl),
//     nsubbins(insubbins),
//     nstate(instate),
//     c2sc(ic2sc),
//     sbin(isbin),
//     ndim(sbin->bidxs.size()),
//     ncorner(sbin->cidxs.size()),
//     subbinsize(std::pow(insubbins,(int)(sbin->bidxs.size()))),
//     ucidxs(sbin->cidxs.size(),0),
//     zwts(instate*std::pow(insubbins,(int)(sbin->bidxs.size())),0.)
// {  
//   typedef std::vector<int>::const_iterator citer; 
//   for ( int ic=0; ic<ncorner; ++ic )
//     {
//       citer p=std::find(glbcidxs.begin(),glbcidxs.end(),sbin->cidxs[ic]);
//       if ( p != glbcidxs.end() )
// 	{
// 	  ucidxs[ic] = std::distance(glbcidxs.begin(),p);
// 	}
//       else
// 	{
// 	  std::cerr << "Failed to locate corner index "
// 		    << sbin->cirds[ic]
// 		    << " from the global array" << std::endl;
// 	  std::exit(EXIT_FAILURE);
// 	}
//     }
// }


ndfes::vFEPData::vFEPData()
  : ndim(0),
    bsplorder(0),
    nsubbins(0),
    nbspl(0),
    nstates(0),
    nbins(0),
    ncorner(0),
    subbinsize(0)
{
}

ndfes::vFEPData::vFEPData
( ndfes::DimInfo const & iinfo,
  std::size_t const insubbins,
  std::vector< ndfes::State > const & istates,
  std::vector< ndfes::SpatialBin > const & isbins )
{
  reset( iinfo, insubbins, istates, isbins );
}


void ndfes::vFEPData::reset
( ndfes::DimInfo const & info,
  std::size_t const insubbins,
  std::vector< ndfes::State > const & states,
  std::vector< ndfes::SpatialBin > const & sbins )
{
  diminfo    = info;
  ndim       = info.GetNumDims();
  bsplorder  = info.GetBsplOrder();
  nbspl      = bsplorder + (bsplorder%2);
  nsubbins   = insubbins;
  nstates    = states.size();
  nbins      = sbins.size();
  ncorner    = info.GetNumCorners();
  subbinsize = std::pow(nsubbins,ndim);

  ////////////////////////////////
  
  glbcidxs.resize(0);
  {
    for ( std::size_t ibin=0; ibin<nbins; ++ibin )
      {
	for ( std::size_t ic=0; ic<ncorner; ++ic )
	  {
	    glbcidxs.push_back( sbins[ibin].cidxs[ic] );
	  }
      }
    std::sort( glbcidxs.begin(), glbcidxs.end() );
    std::vector<std::size_t>::iterator p =
      std::unique(glbcidxs.begin(), glbcidxs.end());
    glbcidxs.resize( std::distance(glbcidxs.begin(),p) );
  }
  
  ///////////////////////////////
  
  ucidxs.assign( ncorner*nbins, 0 );
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      for ( std::size_t ic=0; ic<ncorner; ++ic )
	{
	  
	  std::vector<std::size_t>::iterator p
	    = std::find( glbcidxs.begin(), glbcidxs.end(),
			 sbins[ibin].cidxs[ic] );
	  
	  if ( p != glbcidxs.end() )
	    {
	      ucidxs[ic+ibin*ncorner] = std::distance(glbcidxs.begin(),p);
	    }
	  else
	    {
	      std::cerr << "Failed to locate corner index "
			<< sbins[ibin].cidxs[ic]
			<< " from the global array" << std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	}
    };

  ///////////////////////////////

  std::vector< std::vector<double> > linqpts(ndim);
  std::vector< std::vector<double> > linqwts(ndim);
  {
    double const * width = diminfo.GetTargetWidths();
    for ( std::size_t dim=0; dim<ndim; ++dim )
      {
	linqpts[dim].resize( nsubbins );
	linqwts[dim].resize( nsubbins );
	GaussLegendreRule( -width[dim]/2.,
			   width[dim]/2.,
			   nsubbins,
			   linqpts[dim].data(),
			   linqwts[dim].data() );
      }
  }

  ///////////////////////////////
  
  {
    double wsize = nstates*subbinsize*nbins * 8. / 1024. / 1024. / 1024.;
    std::cout << "Configuration integral matrix size "
	      << std::fixed << std::setw(8) << std::setprecision(2)
	      << wsize << " GB" << std::endl;
    
    zwts.assign( nstates*subbinsize*nbins, 0. );

    std::vector<double> widths(ndim,0.);
    {
      double const * cmin = diminfo.GetXmin();
      double const * cmax = diminfo.GetXmax();
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  widths[dim] = cmax[dim]-cmin[dim];
	}
    }
      
    std::vector<double> qwts;
    ndfes::LinearWtsToMeshWts( linqwts, qwts );
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for ( std::size_t ibin=0; ibin<nbins; ++ibin )
      {
	for ( std::size_t i=0; i<nstates; ++i )
	  {
	    std::vector< std::vector<double> > linswts(linqwts);
	    double const * Rc = states[i].GetCenter();
	    double const * ks = states[i].GetFConst();

	    for ( std::size_t dim=0; dim<ndim; ++dim )
	      {
		double const hw = widths[dim];
		long double const k = ks[dim];
		long double const x0 = sbins[ibin].center[dim] - Rc[dim];
		for ( std::size_t iq=0, nq=linqpts[dim].size(); iq<nq; ++iq )
		  {
		    double dx = linqpts[dim][iq] + x0;
		    if ( diminfo.IsPeriodic(dim) )
		      {
			dx = wrap(dx,hw);
		      }
		    double ene = states[i].GetBeta() * ( k*dx*dx );
		    linswts[dim][iq] *= std::exp(-ene);
		  }
	      }
	    std::vector<double> twts;
	    ndfes::LinearWtsToMeshWts( linswts, twts );
	    for ( std::size_t isub=0; isub < subbinsize; ++isub )
	      zwts[i+(isub+ibin*subbinsize)*nstates] = twts[isub];
	  };
      }
  }

  ///////////////////////////////
  
  {
    c2sc.assign( ncorner*subbinsize, 0. );

    std::vector<double> subcrds;
    ndfes::LinearSpacingsToMeshgrid( linqpts, subcrds );

    double const * widths = diminfo.GetTargetWidths();
    for ( std::size_t isub=0; isub<subbinsize; ++isub )
      {
	std::vector< std::vector<double> > linsubwts(ndim);
	std::vector<double> gwts(bsplorder,0.);
	std::vector<int> gidx(bsplorder,0);
	for ( std::size_t dim=0; dim<ndim; ++dim )
	  {
	    linsubwts[dim].assign(nbspl,0.);
	    std::fill( gwts.data(), gwts.data() + bsplorder, 0. );
	    double w = (nbspl-1)*widths[dim];
	    double c = w/2.;
	    double dx = subcrds[dim+isub*ndim] + c;
	    //std::printf("bsplorder %i/%i %i %i %12.6f %12.6f\n",isub,subbinsize,bsplorder,nbspl,dx,w);
	    ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
				      gwts.data(), gidx.data() );
	    for ( std::size_t ib=0; ib<bsplorder; ++ib )
	      {
		linsubwts[dim][ gidx[ib] ] = gwts[ib];
	      }
	  }
	std::vector<double> gridsubwts;
	ndfes::LinearWtsToMeshWts(linsubwts,gridsubwts);
	for ( std::size_t icorner=0, nc=ncorner; icorner<nc; ++icorner )
	  {
	    c2sc[icorner+isub*nc] = gridsubwts[icorner];
	  };
      }
  }
  
}





void ndfes::vFEPData::CornerValsToSubCellVals
( std::size_t const ibin,
  double const * cornervals,
  double * subcellvals ) const
{
  //std::fill( subcellvals, subcellvals+subbinsize, 0. );
  std::vector<double> cvals(ncorner,0.);
  for ( std::size_t ic=0; ic<ncorner; ++ic )
    {
      cvals[ic] = cornervals[ ucidxs[ic+ibin*ncorner] ];
    }
  double alp = 1.;
  double bet = 0.;
  int inc = 1;
  int nc = ncorner;
  int sbs = subbinsize;
  dgemv_( "T", &nc, &sbs, &alp, c2sc.data(),
	  &nc, cvals.data(), &inc, &bet, subcellvals, &inc );
}


void ndfes::vFEPData::SubCellGrdsToCornerGrds
( std::size_t const ibin,
  double const * subcellvals,
  double * cornervals ) const
{
  std::vector<double> cvals(ncorner,0.);  
  double alp = 1.;
  double bet = 0.;
  int inc = 1;
  int nc = ncorner;
  int sbs = subbinsize;
  dgemv_( "N", &nc, &sbs, &alp, c2sc.data(),
	  &nc, subcellvals, &inc, &bet, cvals.data(), &inc );
  for ( std::size_t ic=0; ic<ncorner; ++ic )
    {
      cornervals[ ucidxs[ic+ibin*ncorner] ] += cvals[ic];
    }
}



std::vector<double> ndfes::vFEPData::CptLinearTerm
( std::vector< ndfes::Sample > const & samples,
  std::vector< std::vector<std::size_t> > const & sbinsamples,
  std::vector< ndfes::SpatialBin > const & sbins ) const
{

  std::vector<double> cornerhs( glbcidxs.size(), 0. );
  
  std::size_t const nsbins = sbinsamples.size();
  if ( nsbins != nbins )
    {
      std::cerr << "ndfes::vFEPData::CptLinearTerm bin size mismatch "
		<< std::setw(6) << nsbins
		<< " but expected "
		<< std::setw(5) << nbins
		<< std::endl;
      std::exit(EXIT_FAILURE);
    };
  
  double const * cmin = diminfo.GetXmin();
  double const * cmax = diminfo.GetXmax();
  double const * twidths = diminfo.GetTargetWidths();

  std::vector<double> histsums(nstates,0.);
  for ( std::size_t ibin=0; ibin<nbins; ++ibin )
    {
      std::size_t const nobs = sbinsamples[ibin].size()/2;
      for ( std::size_t iobs=0; iobs<nobs; ++iobs )
	{
	  std::size_t const l = sbinsamples[ibin][0+iobs*2];
	  std::size_t const degen = sbinsamples[ibin][1+iobs*2];
	  histsums[ samples[l].GetStateIdx() ] += degen;
	}
    }
  

#ifdef WITH_OPENMP
#pragma omp parallel
  {
    std::vector<double> mych( glbcidxs.size(), 0. );
#pragma omp for
#endif
    for ( std::size_t ibin=0; ibin < nbins; ++ibin )
    {
      std::size_t const * idxs = sbins[ibin].bidxs.data();
      std::vector<double> gwts(nbspl,0.);
      std::vector<int> gidx(nbspl,0);
      std::vector< std::vector<double> > linh(ndim);
      for ( std::size_t dim=0; dim<ndim; ++dim )
	{
	  linh[dim].resize(nbspl);
	}
      std::size_t const nobs = sbinsamples[ibin].size()/2;
      std::vector<double> pt(ndim,0.);
      for ( std::size_t iobs=0; iobs<nobs; ++iobs )
	{
	  std::size_t const sidx = sbinsamples[ibin][0+iobs*2];
	  std::size_t const degen = sbinsamples[ibin][1+iobs*2];
	  ndfes::Sample const * s = samples.data() + sidx;
	  double const * R = s->GetPt();
	  for ( std::size_t dim=ndim; dim-->0; )
	    {
	      double dx = R[dim];

	      //std::printf("iobs %8i dim %2i R %10.6f ",iobs,dim,dx);
	      if ( diminfo.IsPeriodic(dim) )
		{
		  dx = wraptorange( R[dim], cmin[dim], cmax[dim] );
		}
	      std::fill( linh[dim].data(), linh[dim].data() + nbspl, 0. );
	      double w = (nbspl-1)*twidths[dim];
	      //std::printf(" dx %10.6f w %10.6f\n",dx,w);
	      // //int idx = (dx-cmin[dim])/twidths[dim];
	      //std::printf("cmin %20.10e idxs %12i nbspl %6i twidths %12.4f\n",cmin[dim],idxs[dim],nbspl,twidths[dim]);
	      dx = dx - ( cmin[dim] + ((int)idxs[dim]-(nbspl/2-1)) * twidths[dim] );
	      //std::printf("dx2 %10.6f\n",dx);
	      dx = std::max(dx,0.);
	      //std::printf("dx3 %10.6f\n",dx);
	      //std::printf("bspl %i/%i %20.10e %20.10e %i %i\n",iobs,nobs,dx,w,idxs[dim],idx);
	      ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
					gwts.data(), gidx.data() );
	      for ( std::size_t ib=0; ib<bsplorder; ++ib )
		{
		  linh[dim][ gidx[ib] ] = gwts[ib];
		};
	    }
	  
	  std::vector<double> cornerh;
	  ndfes::LinearWtsToMeshWts( linh, cornerh );

	  double wt = degen / histsums[ s->GetStateIdx() ];
	  for ( std::size_t ic=0; ic<ncorner; ++ic )
	    {
#ifdef WITH_OPENMP
	      mych[ ucidxs[ic+ibin*ncorner] ] += cornerh[ic]*wt;
      
#else
	      cornerhs[ ucidxs[ic+ibin*ncorner] ] += cornerh[ic]*wt;
#endif	      
	    }
	  
	}
    }
#ifdef WITH_OPENMP
#pragma omp critical
    {
      for ( std::size_t i=0, n=glbcidxs.size(); i<n; ++i )
	{
	  cornerhs[i] += mych[i];
	}
    }
  }
#endif

  return cornerhs;
}


double ndfes::vFEPData::CptChisq
( std::vector<double> const & x,
  std::vector<double> & g,
  std::vector<double> const & cornerhs ) const
{
  std::size_t const nparam = x.size();
  int const inc = 1;
  double const one = 1.;
  double chisq = 0.;

  
  if ( g.size() > 0 )
    {      
      std::vector<double> Za(nstates,0.);
      std::vector<double> dZdx(nparam*nstates,0.);
#ifdef WITH_OPENMP
#pragma omp parallel
      {
	std::vector<double> myZa(nstates,0.);
	std::vector<double> mydZdx(nparam*nstates,0.);
#endif
	std::vector<double> Fsub(subbinsize,0.);
	std::vector<double> subdZdx(subbinsize*nstates,0.);
#ifdef WITH_OPENMP
#pragma omp for
#endif
	for ( std::size_t ibin=0; ibin<nbins; ++ibin )
	  {
	    CornerValsToSubCellVals( ibin, x.data(), Fsub.data() );
	    std::fill( subdZdx.data(), subdZdx.data() + subdZdx.size(), 0. );
	    for ( std::size_t isub=0; isub<subbinsize; ++isub )
	      {
		double expFsub = std::exp(-Fsub[isub]);
		for ( std::size_t i=0; i<nstates; ++i )
		  {
		    double we = zwts[i+(isub+ibin*subbinsize)*nstates]*expFsub;
		    subdZdx[isub+i*subbinsize] -= we;
#ifdef WITH_OPENMP
		    myZa[i] += we;
#else
		    Za[i] += we;
#endif
		  } // i
	      } // isub
	    for ( std::size_t i=0; i<nstates; ++i )
	      {
#ifdef WITH_OPENMP
		SubCellGrdsToCornerGrds( ibin, subdZdx.data()+i*subbinsize,
					 mydZdx.data()+i*nparam );
#else
		SubCellGrdsToCornerGrds( ibin, subdZdx.data()+i*subbinsize,
					 dZdx.data()+i*nparam );
#endif
	      } // i
	  } // ibin
#ifdef WITH_OPENMP
#pragma omp critical
	{
	  for ( std::size_t i=0; i<nstates; ++i )
	    {
	      Za[i] += myZa[i];
	    }
	  for ( std::size_t i=0, n=nstates*nparam; i<n; ++i )
	    {
	      dZdx[i] += mydZdx[i];
	    }
	} // critical
      } // omp parallel
#endif
      {
	int mynparam = nparam;
	chisq = ddot_(&mynparam,x.data(),&inc,cornerhs.data(),&inc);
      }
      for ( std::size_t i=0; i<nparam; ++i )
	{
	  g[i] = cornerhs[i];
	} // i
      
      for ( std::size_t i=0; i<nstates; ++i )
	{
	  chisq += std::log(Za[i]);
	  Za[i] = 1. / Za[i];
	};

      {
	int mynparam = nparam;
	int mynstates = nstates;
	dgemv_( "N", &mynparam, &mynstates, &one, dZdx.data(), &mynparam,
		Za.data(), &inc, &one, g.data(), &inc );
      }
      
      
    } // g.size() > 0
  else //////////////////////////////////////////////
    { // g.size() == 0
      
      std::vector<double> Za(nstates,0.);
#ifdef WITH_OPENMP
#pragma omp parallel
      {
	std::vector<double> myZa(nstates,0.);
#endif
	std::vector<double> Fsub(subbinsize,0.);
#ifdef WITH_OPENMP
#pragma omp for
#endif
	for ( std::size_t ibin=0; ibin<nbins; ++ibin )
	  {
	    CornerValsToSubCellVals( ibin, x.data(), Fsub.data() );
	    for ( std::size_t isub=0; isub<subbinsize; ++isub )
	      {
		double expFsub = std::exp(-Fsub[isub]);
		for ( std::size_t i=0; i<nstates; ++i )
		  {
		    double we = zwts[i+(isub+ibin*subbinsize)*nstates]*expFsub;
#ifdef WITH_OPENMP
		    myZa[i] += we;
#else
		    Za[i] += we;
#endif
		  } // i
	      } // isub
	  } // ibin
#ifdef WITH_OPENMP
#pragma omp critical
	{
	  for ( std::size_t i=0; i<nstates; ++i )
	    {
	      Za[i] += myZa[i];
	    }
	}
      } // omp parallel
#endif
      {
	int mynparam=nparam;
	chisq = ddot_(&mynparam,x.data(),&inc,cornerhs.data(),&inc);
      }
      for ( std::size_t i=0; i<nstates; ++i )
	{
	  chisq += std::log(Za[i]);
	};
    } // g.size() == 0
  return chisq;
}




void ndfes::vFEPData::InterpvFEP
( double const * R,
  std::vector<double> const & p,
  std::vector<double> const & dp,
  double & val,
  double & err ) const
{
  std::size_t const ndim         = diminfo.GetNumDims();
  int const * dimsizes   = diminfo.GetDimSizes();
  double const * cmin    = diminfo.GetXmin();
  double const * cmax    = diminfo.GetXmax();
  double const * twidths = diminfo.GetTargetWidths();

  // std::printf("eval @ ");
  // for ( int dim=0; dim<ndim; ++dim )
  //   {
  //     std::printf("%8.3f",R[dim]);
  //   }
  // std::printf("\n");
  
  std::vector<double> gwts(bsplorder,0.);
  std::vector<int> gidx(bsplorder,0);
  
  std::vector< std::vector<double> > linh(ndim);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      linh[dim].assign(bsplorder,0.);
    }

  std::vector< std::vector<std::size_t> > lincidxs(ndim);
  for ( std::size_t dim=0; dim<ndim; ++dim )
    {
      lincidxs[dim].assign(bsplorder,0);
    }

  for ( std::size_t dim=ndim; dim-->0; )
    {
      //std::printf("dim %3i  ",dim);
      double dx = R[dim];
      if ( diminfo.IsPeriodic(dim) )
	{
	  dx = wraptorange( R[dim], cmin[dim], cmax[dim] );
	}

      std::size_t binidx = (dx-cmin[dim])/twidths[dim];      
      
      double w = (nbspl-1)*twidths[dim];
      dx = dx - ( cmin[dim] + (binidx-(nbspl/2-1)) * twidths[dim] );
      dx = std::max(dx,0.);

      ndfes::bspline_aperiodic( dx, w, nbspl, bsplorder,
				gwts.data(), gidx.data() );
      
      for ( std::size_t ib=0; ib<bsplorder; ++ib )
	{
	  int lidx = gidx[ib] + binidx - (nbspl/2-1);
	  if ( diminfo.IsPeriodic(dim) )
            {
              lidx = INTWRAP(lidx,dimsizes[dim]);
            }
	  //std::printf("%5i",lidx);
	  linh[dim][ib] = gwts[ib];
	  lincidxs[dim][ib] = (std::size_t)lidx;
	};
      //std::printf("\n");
    }

  
  std::vector<double> cornerh;
  ndfes::LinearWtsToMeshWts( linh, cornerh );
  std::vector<std::size_t> cidxs;
  ndfes::LinearSpacingsToMeshgrid( lincidxs, cidxs );
  std::size_t const bsize = std::pow(bsplorder,ndim);
  
  val = 0.;
  err = 0.;
  for ( std::size_t ic=0; ic<bsize; ++ic )
    {
      std::size_t cidx=0;
      for ( std::size_t dim=ndim; dim-->0; )
        {
          if ( diminfo.IsPeriodic(dim) )
            {
              cidx = cidxs[dim+ic*ndim] + cidx*((std::size_t)dimsizes[dim]);
            }
          else
            {
              cidx = cidxs[dim+ic*ndim] + cidx*((std::size_t)(dimsizes[dim]+1));
            }
        }
      std::vector<std::size_t>::const_iterator pidx = std::find( glbcidxs.begin(), glbcidxs.end(), cidx );
      if ( pidx != glbcidxs.end() )
	{
	  std::size_t const iparam = std::distance( glbcidxs.begin(), pidx );
	  if ( iparam >= p.size() )
	    {
	      std::cout << "ndfes::vFEPData::InterpvFEP Parameter index "
			<< "out of bounds " << iparam << " " << p.size()
			<< "\n";
	      std::exit(EXIT_FAILURE);
	    }
	  double x = p[iparam];
	  double dx = dp[iparam];
	  double w = cornerh[ic];
	  val += w*x;
	  err += std::pow( w*dx, 2 );

	  //double beta = ndfes::GetBeta(298.);
	  //std::printf("gidx %6i %20.10e %20.10e\n",cidx,x/beta,w);

	}
      else
	{
	  std::cout << "ndfes::vFEPData::InterpvFEP Failed to find "
		    << "parameter index of global corner index "
		    << cidx << "\n";
	  std::exit(EXIT_FAILURE);
	}
    };

  if ( err > 0 )
    {
      err = std::sqrt(err);
    }
  else
    {
      err = 0.;
    };
}
