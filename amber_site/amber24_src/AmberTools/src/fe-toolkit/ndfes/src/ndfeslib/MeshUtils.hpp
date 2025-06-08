#ifndef _ndfes_MeshUtils_hpp_
#define _ndfes_MeshUtils_hpp_

#include <vector>

extern "C"
{
  extern void
  dger_(const int *m, const int *n, double const *alpha,
        double const *x, const int *incx,
        double const *y, const int *incy,
        double *a, const int *lda);

}


namespace ndfes
{
  template<typename T>
  void LinearSpacingsToMeshgrid
  ( std::vector< std::vector<T> > const & lsp,
    std::vector<T> & meshgrid );
  
  template<typename T>
  void LinearWtsToMeshWts
  ( std::vector< std::vector<T> > const & lsp,
    std::vector<T> & meshgrid );
}



template<typename T>
void ndfes::LinearSpacingsToMeshgrid
( std::vector< std::vector<T> > const & lsp,
  std::vector<T> & meshgrid )
{
  std::size_t ndim = lsp.size();
  
  // dim 0
  meshgrid.resize( ndim*lsp[0].size() );
  for ( std::size_t i=0, n=lsp[0].size(); i<n; ++i )
    meshgrid[0+i*ndim] = lsp[0][i];

  for ( std::size_t dim=1; dim<ndim; ++dim )
    {
      std::vector<T> blk(meshgrid);
      std::size_t nsize = blk.size();
      std::size_t npts = nsize / ndim;
      meshgrid.resize( nsize * lsp[dim].size() );
      std::size_t o=0;
      for ( std::size_t i=0, ni = lsp[dim].size(); i<ni; ++i )
        {
          for ( std::size_t j=0; j<npts; ++j )
            blk[dim+j*ndim] = lsp[dim][i];
          for ( std::size_t j=0; j<nsize; ++j )
            meshgrid[j+o] = blk[j];
          o += nsize;
        };
    };
}

template<typename T>
void ndfes::LinearWtsToMeshWts
( std::vector< std::vector<T> > const & lsp,
  std::vector<T> & Ws )
{
  std::size_t ndim = lsp.size();
  //int nquad_per_dim = lsp[0].size();
  
  Ws = lsp[0];

  int inc=1;
  double alpha=1.;
  for ( std::size_t dim=1; dim < ndim; ++dim )
    {
      std::vector<double> Wfast(Ws);
      std::vector<double> Wslow(lsp[dim]);
      int nf=Wfast.size();
      int ns=Wslow.size();
      Ws.resize( nf*ns );
      std::fill( Ws.data(), Ws.data() + Ws.size(), 0. );
      dger_( &nf, &ns, &alpha, Wfast.data(), &inc, Wslow.data(), &inc, Ws.data(), &nf );
    }
}


#endif
