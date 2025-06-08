#ifndef _edgembar_options_hpp_
#define _edgembar_options_hpp_

#include <string>
#include "../edgembar/ReadMode.hpp"

namespace edgembar
{
  class cli_options
  {
  public:

    cli_options();

    void print_help();

    double temp;
    double tol;
    double btol;
    double ptol;
    int nbootstrap;
    int verbosity;
    bool nullguess;
    edgembar::ReadMode readmode;
    std::string inpfile;
    std::string outfile;
    int ncon;
    double dcon;
    int ntimes;
    bool fwdrev;
    bool halves;
    double fstart;
    double fstop;
    int stride;
    bool autoequil;
    bool uwts;
    int autoeqmode;
  };


  cli_options read_options( int argc, char * argv[] );
}


#endif
