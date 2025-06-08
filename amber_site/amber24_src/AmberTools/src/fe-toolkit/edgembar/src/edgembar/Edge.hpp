#ifndef _edgembar_Edge_hpp_
#define _edgembar_Edge_hpp_

#include <memory>
#include <vector>
#include "Env.hpp"
#include "Constraints.hpp"
#include "CalcIdx.hpp"
#include "StateValues.hpp"

namespace edgembar
{
  class Edge
  {
  public:
    
    Edge( std::string const name,
	  edgembar::Env const & complex,
	  edgembar::Env const & solvated );



    std::shared_ptr<edgembar::Edge> ExtractCopy() const;

    std::shared_ptr<edgembar::Edge> ExtractRange
    ( double const start, double const stop ) const;

    std::shared_ptr<edgembar::Edge> ExtractProd() const;
    
    std::shared_ptr<edgembar::Edge> ExtractBootstrap() const;

    
    std::vector<edgembar::Sim *> GetSims();
    std::vector<edgembar::Sim const *> GetConstSims() const;
    
    std::vector<edgembar::Trial *> GetTrials();
    std::vector<edgembar::Trial const *> GetConstTrials() const;
    
    std::string GetName() const { return mName; }

    std::vector<double> GetFreeFromParams
    ( std::vector<double> const & p ) const;
    
    std::vector<double> GetParamsFromFree
    ( std::vector<double> const & p ) const;

    // std::vector<double> GetStateFreeEneFromParams
    // ( std::vector<double> const & p ) const;

    void ApplyConstraint( double const cval );
    
    void ApplyMinimalConstraints();

    void RemoveConstraints();

    int GetNumParam() const { return mNumParam; }

    int GetNumFreeParam() const { return mConstraints.GetNumFreeParam(); }

    void StoreFreeAvgEnes();
    
    double Optimize( std::vector<double> const & pguess,
		     std::vector<double> & popt,
		     std::vector<double> & fopt,
		     double const tol,
		     int const verbosity ) const;

    std::vector<double> MakeExpAvgGuess() const;

    double GetFreeEnergy( std::vector<double> const & fopt ) const;

    std::string GetPython
    ( int argc,
      char * argv[],
      double const beta,
      std::vector<edgembar::CalcIdx> const & calcs,
      std::vector<edgembar::StateValues> const & vals,
      double const ptol ) const;

    void UseUniformWts();
    
  private:
    
    std::string mName;
    int mNumParam;
    std::vector<int> mFreeEneIdxs;
    std::vector<double> mFreeEneCoefs;
    std::vector<double> mFreeShifts;
    std::vector<double> mFreeAvgEnes;
    std::vector<double> mAvgEnes;
    std::vector<double> mShifts;
    edgembar::Env mComplex;
    edgembar::Env mSolvated;
    edgembar::Constraints mConstraints;
  };
}

#endif
