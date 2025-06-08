#ifndef _amber_Disang_hpp_
#define _amber_Disang_hpp_

#include <vector>
#include <string>
#include <ostream>
#include <istream>

namespace amber
{
  
  class Restraint
  {
  public:
    
    Restraint( std::string line );

    void Write( std::ostream & cout ) const;
    
    std::vector<int> iat;
    std::vector<double> rstwt;
    double r1;
    double r2;
    double r3;
    double r4;
    double rk2;
    double rk3;

    int tidx;
    std::vector<std::string> extra;
    
    bool isr12;
    bool isangle;
    bool isdihed;
    
  };


  class Disang
  {
  public:

    Disang() {};
    
    Disang( std::string fname );

    void Read( std::string fname );
    
    void Read( std::istream & fname );
    
    void Write( std::string fname ) const;

    void Write( std::ostream & cout ) const;

    bool ValidateTemplate() const;
    
    std::vector<int> GetTemplateIdxs() const;

    std::vector<bool> GetTemplateIdxIsAngle() const;
    
    std::vector<double> GetTemplateForceConsts() const;
    
    amber::Disang FillTemplateValues( double const * rcs, double const * fcs ) const;

    std::vector<double> GetCenters( std::vector<int> const & idxs ) const;
    
    std::vector<double> GetForceConsts( std::vector<int> const & idxs ) const;

  public:

    std::vector<amber::Restraint> mRestraints;
  };
  
}

std::ostream & operator<<( std::ostream & cout, amber::Restraint const & res );
std::ostream & operator<<( std::ostream & cout, amber::Disang const & dis );


#endif
