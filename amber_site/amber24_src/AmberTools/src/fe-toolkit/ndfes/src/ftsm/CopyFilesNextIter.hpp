#ifndef _ndfes_CopyFilesNextIter_hpp_
#define _ndfes_CopyFilesNextIter_hpp_

#include <string>
#include <vector>


namespace ndfes
{
  
  void CopyFilesNextIter
  ( std::size_t const ndim,
    std::size_t const nsim,
    std::string const disang,
    std::string const mdin,
    int const curit,
    std::string const nextdir,
    std::string const prefix,
    std::vector<std::string> const extra_prefixes,
    bool cp_nearest,
    bool not_isapath,
    int firstit,
    int pad,
    double const * rcs,
    double const * fcs );

}

#endif

