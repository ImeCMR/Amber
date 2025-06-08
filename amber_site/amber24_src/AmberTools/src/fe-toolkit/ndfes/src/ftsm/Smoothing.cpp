#include "Smoothing.hpp"


std::vector<double> ndfes::GetWinAvg
( int const N,
  double const * xs,
  int wlen )
{
  wlen = std::min(N,wlen);
  if ( wlen % 2 == 0 )
    {
      --wlen;
    }

  std::vector<double> avgs(N,0);
  int disp = (wlen-1)/2;
  for ( int i=0; i<N; ++i )
    {
      for ( int j=i-disp; j<i+disp+1; ++j )
	{
	  double x = 0;
	  if ( j < 0 )
	    {
	      x = 2*xs[0] - xs[-j];
	    }
	  else if ( j > N-1 )
	    {
	      x = 2*xs[N-1] - xs[2*N-j-2];
	    }
	  else
	    {
	      x = xs[j];
	    };
	  avgs[i] += x/wlen;
	}
    };
  return avgs;
}
  
std::vector<double> ndfes::GetIterWinAvg
( int const N,
  double const * xs,
  int wlen,
  int wmax )
{
  std::vector<double> xnew(N,0);
  double fact = 1;
  int wlim = std::max(wmax,wlen);
  for ( int it=0, w=wlen; w<wlim+1; w += 2, ++it )
    {
      
      std::vector<double> dx(N,0);
      for ( int i=0; i<N; ++i )
	{
	  dx[i] = fact * (xs[i]-xnew[i]);
	}

      std::vector<double> odx = ndfes::GetWinAvg(N,dx.data(),w);

      for ( int i=0; i<N; ++i )
	{
	  xnew[i] += odx[i];
	}
      
      if ( it > 1 )
	{
	  fact *= 0.5;
	}
    }
  return xnew;
}

std::vector<double> ndfes::GetIterWinAvg
( int const ndim,
  int const n,
  double const * xs,
  int const wlen,
  int const wmax )
{
  std::vector<double> odata(ndim*n,0);
  for ( int k=0; k<ndim; ++k )
    {
      std::vector<double> data(n,0);
      for ( int i=0; i<n; ++i )
	{
	  data[i] = xs[k+i*ndim];
	}

      data = GetIterWinAvg(n,data.data(),wlen,wmax);

      for ( int i=0; i<n; ++i )
	{
	  odata[k+i*ndim] = data[i];
	}
    }
  return odata;
}
