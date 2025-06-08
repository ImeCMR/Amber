#ifndef _ndfes_options_hpp_
#define _ndfes_options_hpp_

#include <string>
#include <vector>


namespace ndfes
{
  class options
  {
  public:

    options();

    void print_help();
    
    double temp;
    int nham;
    std::string metafile;
    std::string chkpt;
    int maxiter;
    int nquad;
    int order;
    int nboot;
    double reltol;
    double btol;
    bool vfep;
    bool mbar;
    bool sdosFullSampling;
    std::string sdosHistPrefix;
    double sdos;
    double startfrac;
    double stopfrac;
    std::vector<double> target_widths;
    std::vector<int> perdims;
    bool uniq;
    int maxuniq;
    bool expert;
  };


  ndfes::options ReadOptions( int argc, char * argv[] );
}


#endif
