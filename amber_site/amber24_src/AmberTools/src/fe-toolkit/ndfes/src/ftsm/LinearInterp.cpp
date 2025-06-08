#include <utility>
#include <algorithm>
#include <cmath>

#include "LinearInterp.hpp"


ndfes::LinearInterp::LinearInterp
( int const n,  double const * x, double const * y )
  : mX(x,x+n),
    mY(y,y+n),
    mB(n,0)
{

  std::vector< std::pair<double,double> > tpairs;
  for ( int i=0; i<n; ++i )
    {
      tpairs.push_back( std::make_pair( mX[i], mY[i] ) );
    }
  std::sort( tpairs.begin(), tpairs.end() );
  for ( int i=0; i<n; ++i )
    {
      mX[i] = tpairs[i].first;
      mY[i] = tpairs[i].second;
    }
  for ( int i=0; i<n-1; ++i )
    {
      double dx = mX[i+1]-mX[i];
      double dy = mY[i+1]-mY[i];
      mB[i] = dy/dx;
    }
}


void ndfes::LinearInterp::reset
( int const n,  double const * x, double const * y )
{
  mX.assign( x,x+n );
  mY.assign( y,y+n );
  mB.assign( n,0 );
  
  std::vector< std::pair<double,double> > tpairs;
  for ( int i=0; i<n; ++i )
    {
      tpairs.push_back( std::make_pair( mX[i], mY[i] ) );
    }
  std::sort( tpairs.begin(), tpairs.end() );
  for ( int i=0; i<n; ++i )
    {
      mX[i] = tpairs[i].first;
      mY[i] = tpairs[i].second;
    }
  for ( int i=0; i<n-1; ++i )
    {
      double dx = mX[i+1]-mX[i];
      double dy = mY[i+1]-mY[i];
      mB[i] = dy/dx;
    }
}


bool ndfes::LinearInterp::push_back( double const x, double const y )
{
  bool pushed = false;
  bool already_has_x = false;
  for ( std::size_t i=0; i<mX.size(); ++i )
    {
      if ( std::abs(mX[i]-x) < 1.e-8 )
	{
	  already_has_x = true;
	  break;
	}
    }
  if ( ! already_has_x )
    {
      pushed = true;
      std::vector<double> xs(mX);
      std::vector<double> ys(mY);
      xs.push_back(x);
      ys.push_back(y);
      reset( xs.size(), xs.data(), ys.data() );
    }
  
  return pushed;
}


int ndfes::LinearInterp::FindIdx( double const x ) const
{
  int n = mX.size();
  int ind = 0;
  
  if ( x < mX[0] or mX.size() == 1 )
    {
      ind = 0;
    }
  else if ( x > mX[n-1] )
    {
      ind = n-1;
    }
  else
    {
      int low = 0;
      int high = n-1;
      while ( low <= high )
        {
          ind = (low+high)/2;
          if ( x < mX[ind] )
            {
              high = ind-1;
            }
          else if ( x > mX[ind+1] )
            {
              low = ind + 1;
            }
          else
            {
              break;
            }
        }
    }
  return ind;
}


double ndfes::LinearInterp::GetValue( double const x ) const
{
  double v = 0;
  if ( x < mX[0] or mX.size() == 0 )
    {
      v = mY[0];
    }
  else if ( x > mX.back() )
    {
      v = mY.back();
    }
  else
    {
      int ifound = FindIdx(x);
      double dx = x - mX[ifound];
      v = mY[ifound] + dx*mB[ifound];
    }
  return v;
}

double ndfes::LinearInterp::GetDeriv( double const x ) const
{
  double v = 0;
  if ( x < mX[0] or mX.size() == 0 )
    {
      v = mB[0];
    }
  else if ( x > mX.back() )
    {
      v = mB[ mX.size() - 2 ];
    }
  else
    {
      int ifound = FindIdx(x);
      v = mB[ifound];
    }
  return v;
}
