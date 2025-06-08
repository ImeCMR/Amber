#ifndef _edgembar_ReadInput_hpp_
#define _edgembar_ReadInput_hpp_

#include <string>
#include <memory>
#include "ReadMode.hpp"
#include "Edge.hpp"

namespace edgembar
{
  std::shared_ptr<edgembar::Edge> ReadInput
  ( std::string fname,
    double const beta,
    edgembar::ReadMode mode,
    double const fstart,
    double const fstop,
    int const stride,
    double const ptol,
    int const autoeqmode );
  
}

#endif
