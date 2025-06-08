#ifndef _pertools_hpp_
#define _pertools_hpp_

#include <cmath>
#include <vector>

inline double anint( double const x )
{
  return ((x)>0? std::floor((x)+0.5) : std::ceil((x)-0.5));
}

inline double wrap( double const x, double const l )
{
  return x - l * anint( x/l );
}

double wraptorange( double const x, double const lo, double const hi );

#define INTWRAP(i,n) ( ( (i)%(n) + (n) )%(n) )


// namespace ndfes
// {

//   template<typename T>
//   inline void LinearSpacingsToMeshgrid
//   ( std::vector< std::vector<T> > const & lsp,
//     std::vector<T> & meshgrid );
  
//   void LinearWtsToMeshWts
//   ( std::vector< std::vector<double> > const & lsp,
//     std::vector<double> & meshgrid ); 
// }




// template<typename T>
// inline void ndfes::LinearSpacingsToMeshgrid
// ( std::vector< std::vector<T> > const & lsp,
//   std::vector<T> & meshgrid )
// {
//   int ndim = lsp.size();
  
//   // dim 0
//   meshgrid.resize( ndim*lsp[0].size() );
//   for ( std::size_t i=0, n=lsp[0].size(); i<n; ++i )
//     meshgrid[0+i*ndim] = lsp[0][i];

//   for ( int dim=1; dim<ndim; ++dim )
//     {
//       std::vector<T> blk(meshgrid);
//       int nsize = blk.size();
//       int npts = nsize / ndim;
//       meshgrid.resize( nsize * lsp[dim].size() );
//       int o=0;
//       for ( int i=0, ni = lsp[dim].size(); i<ni; ++i )
//         {
//           for ( int j=0; j<npts; ++j )
//             blk[dim+j*ndim] = lsp[dim][i];
//           for ( int j=0; j<nsize; ++j )
//             meshgrid[j+o] = blk[j];
//           o += nsize;
//         };
//     };
// }


#endif
