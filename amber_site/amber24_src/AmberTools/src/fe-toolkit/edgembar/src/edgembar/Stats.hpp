#ifndef _edgembar_Stats_hpp_
#define _edgembar_Stats_hpp_

namespace edgembar
{

  // // set tol to larger values to reject test more frequently
  // bool HaveTheSameMean
  // ( int N1, double mu1, double serr1,
  //   int N2, double mu2, double serr2,
  //   double tol );

  // // set tol to larger values to reject test more frequently
  // bool HaveTheSameVar
  // ( int N1, double serr1,
  //   int N2, double serr2,
  //   double tol );
  
  // void CptMeanAndStderr
  // ( int const N,
  //   double const * A,
  //   double & mu,
  //   double & err );

  
  // void CptMeanAndStderr
  // ( int const istart,
  //   int const stride,
  //   int const M,
  //   double const * A,
  //   double & mu,
  //   double & err );
  
  // int CptStatIneff
  // ( int const n,
  //   double const * x );

  
  // int CptSampleStride
  // ( int const n,
  //   double const * x );

  
  // bool CptSampleStartAndStride
  // ( int const n,
  //   double const * x,
  //   double const maxeq,
  //   int & istart,
  //   int & stride );

  // bool CptSampleStartAndStrideByBlock
  // ( int const n,
  //   double const * x,
  //   double const maxeq,
  //   int const nblocks,
  //   int & istart,
  //   int & stride );


  bool CptSampleStartAndStrideByBlock_v3
  ( int const n,
    double const * x,
    double const maxeq,
    int const nblocks,
    int & istart,
    int & stride,
    double const ptol );

  
  class TimeseriesSegment
  {
  public:

    TimeseriesSegment();
    
    TimeseriesSegment( int istart, int istop, int nxs, double const * xs );

    void reset( int istart, int istop, int nxs, double const * xs );

    bool IsSimilarTo( edgembar::TimeseriesSegment const & rhs, double const ptol ) const;
    
  public:
    int istart;
    int istop;
    int s;
    
    int cn;
    double cavg;
    double cvar;
    double cstd;
    double cerr;

    int un;
    double uavg;
    double uvar;
    double ustd;
    double uerr;
    
    std::vector<double> cdata;
    std::vector<double> udata;
  };

  /*
  class TimeseriesHalves
  {
  public:

    TimeseriesHalves();
    
    TimeseriesHalves( int istart, int istop, int nxs, double const * xs );

    void reset( int istart, int istop, int nxs, double const * xs );

    //bool IsSelfSimilar( double const ptol ) const;
    bool HasSameMeans( double const ptol ) const;
    bool HasSameMeansAndDist( double const ptol ) const;
    bool HasSameDist( double const ptol ) const;

  public:
    int istart;
    int istop;
    int s;
    int n;
    int n1;
    int n2;
    double avg1;
    double err1;
    double avg2;
    double err2;
    std::vector<double> data1;
    std::vector<double> data2;
  };
  */
  
}

#endif
