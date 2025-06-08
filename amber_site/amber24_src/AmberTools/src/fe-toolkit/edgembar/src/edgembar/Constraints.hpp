#ifndef _edgembar_Constraints_hpp_
#define _edgembar_Constraints_hpp_

#include <vector>

namespace edgembar
{
  class Constraints
  {
  public:
    
    Constraints();
    
    Constraints( int nparam );
    

    void reset( int nparam );
    
    void push_back( std::vector<int> const & pidxs,
		    std::vector<double> const & coefs,
		    double const cval );

    void finalize();

    int GetNumFreeParam() const { return mNumFreeParam; }

    int GetNumParam() const { return mNumParam; }

    int GetNumCon() const { return mNumCon; }

    std::vector<double> GetParamsFromFree( std::vector<double> const & q ) const;

    std::vector<double> GetFreeFromParams( std::vector<double> const & q ) const;

    
    
  private:

    int mNumParam;
    int mNumFreeParam;
    int mNumCon;

    std::vector<double> mConMat;  // fast: mNumParam, slow: mNumCon
    std::vector<double> mConVals; // fast: mNumCon
    std::vector<double> VT; // fast: mNumParam, slow: mNumParam
    std::vector<double> P0; // fast: mNumParam

  };
}

#endif
