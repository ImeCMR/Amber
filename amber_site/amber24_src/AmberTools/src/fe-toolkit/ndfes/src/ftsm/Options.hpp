#ifndef _ndfes_ftsm_options_hpp_
#define _ndfes_ftsm_options_hpp_

#include <string>
#include <vector>

class Options
{
public:
  
  Options();
  void print_help() const;
  
  std::string disang;
  std::string mdin;
  std::string ipath;
  std::string opath;
  std::string chk;
  int curit;
  std::string odir;
  //bool dacc;
  int model;
  bool rbf;
  bool bspl;
  bool hist;
  int arbf;
  double shape;
  int minsize;
  double rbfmaxerr;
  int maxit;
  double ptol;
  double distol;
  double angtol;
  int mlen;
  //bool closest;
  double temp;
  bool fix0;
  bool fix1;
  //double dampfactor;
  //double smooth_frac;
  //int smooth_nmin;
  //int smooth_nmax;
  bool smoothfc;
  //bool scalefc;
  bool predfc;
  double maxfc;
  double minfc;
  double maxfc_ang;
  double minfc_ang;
  double tdispfc;
  //bool pkl;
  //int nquad;
  bool cp_nearest;
  bool boot;
  double booterror;
  int npathpts;
  int nsplpts;
  //double nebk;
  //int conv_layers;
  //int conv_samples;
  //bool tetris;
  //bool tetris_smooth;
  bool acc;
  double acc_oobk;
  int acc_maxit;
  bool buffer;
  //int sim_layers;
  int nlayers;
  bool dry_run;
  int neqit;
  int wavg;
  int wavg_niter;
  std::string prefix;
  std::vector<std::string> extra_prefixes;
  bool smooth;
  bool msm;
  bool gbsm;
  int stype;
  int printfreq;
  int pad;
  //double occpen;
  bool sasm;
};

Options ReadOptions( int argc, char * argv[] );

#endif

