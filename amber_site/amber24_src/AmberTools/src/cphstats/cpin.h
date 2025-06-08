#ifndef CPIN_H
#define CPIN_H

#include <string>
#include <vector>
#include "constants.h"

// External Fortran functions for parsing the cpin

typedef struct {
   int num_states, first_atom,  num_atoms,  first_state, first_charge;
} StateInfo;

extern "C" {
#ifdef REDOX
   void parse_cein_(int*, int*, StateInfo*,
                        char(*)[40], char[FN_LEN], int*, int*,
                        const int*, const int*, const int*);
#else
   void parse_cpin_(int*, int*, StateInfo*,
                        char(*)[40], char[FN_LEN], int*, int*,
                        const int*, const int*, const int*);
#endif
void nmlsrc_from_file_(char[FN_LEN], int*, const char[FN_LEN], int*, int*);
void parse_limits_(const char[FN_LEN], int*, int*, int*, int*, int*);
}

/// A titratable residue
class TitratableResidue {
   public:
      // Constructors (only the implemented/used ones are uncommented)
//    TitratableResidue(const char*);
//    TitratableResidue(std::string const&);
      TitratableResidue(std::string const&, std::vector<int> const&);
//    TitratableResidue(const char*, int*);

      // A way to add a list of states
      void defineStates(std::vector<int>);
      void defineStates(int*);

      // In-line Setters
      void setResname(const char* resname) { resname_ = std::string(resname); }
      void setResname(std::string const& resname) { resname_ = resname; }
      void setResnum(int resnum) { resid_ = resnum; }

      // In-line Getters
      std::string getResname() const { return resname_; }
      int getResnum() const { return resid_; }

      // Determine how many states are present
      int numStates() const { return protonated_.size(); }

      // Determine if a particular state is protonated or not
      bool isProtonated(int state) const { return protonated_[state]; }
      // Determine how many protons are in a specific state
      int numProtons(int state) const { return protcnts_[state]; }
#ifdef REDOX
      // Computing the maximum and the minimum number of electrons at ELECCNT on cein, in order to get the value of v of the Nernst equation
      int getvNernst() const {
       int nelecmin = numProtons(0);
       int nelecmax = numProtons(0);
         for (int j = 1; j < numStates(); j++) {
           nelecmax = (numProtons(j) > nelecmax) ? numProtons(j) : nelecmax;
           nelecmin = (numProtons(j) < nelecmin) ? numProtons(j) : nelecmin;
         }
       return (nelecmax-nelecmin);
      }
#endif

   private:
      /// How many protons are in each state?
      std::vector<int> protcnts_;

      /// Is each state protonated?
      std::vector<bool> protonated_;

      /// Name of the titratable residue
      std::string resname_;

      /// Number of the residue
      int resid_;
};

/// Wrapper for calling the Fortran cpin
class Cpin {
   public:
      // Constructors
      Cpin();

      // Get the data
      std::vector<TitratableResidue> getResidues() const { return residues_; }

      // Provide an iterator over the data
      typedef std::vector<TitratableResidue>::const_iterator ResIterator;
      ResIterator begin() const { return residues_.begin(); }
      ResIterator end()   const { return residues_.end();   }

      // Parse the cpin
      int Parse(const char*);
      int Parse(std::string const&);

      // Get the number of residues
      int getTrescnt() const { return trescnt_; }

      // Get the name of the file
      std::string getFilename() const { return filename_; }

      bool isRead() const { return is_valid_; }

      bool isCpein() const { return is_cpein_ == 1; }

   private:
      /// Limits for array sizes
      int TITR_RES_C;
      int TITR_STATES_C;
      int ATOM_CHRG_C;
      int MAX_H_COUNT;

      /// Is our cpin file valid?
      bool is_valid_;

      /// Is the input file a cpein file? If not, it is a cpin
      int is_cpein_;

      /// Number of titratable residues defined in the cpin
      int trescnt_;

      /// Name of each titratable residue
      std::vector<std::string> resname_;

      /// List of titratable residues
      std::vector<TitratableResidue> residues_;

      /// Name of the file
      std::string filename_;

};
#endif /* CPIN_H */
