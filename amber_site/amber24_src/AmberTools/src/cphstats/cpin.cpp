/** cpin.cpp: Handles the code for parsing cpin files and defining titratable
  * residues.
  */

#include <iostream>
#include <fstream>

#include "constants.h"
#include "cpin.h"
#include "exceptions.h"
#include "string_manip.h"

using namespace std;

/// Cpin constructor from a char array
Cpin::Cpin() :
is_valid_(false),
is_cpein_(1)
{ }

int Cpin::Parse(const char* cpinname) {
   filename_ = string(cpinname);
   if (filename_.size() > FN_LEN)
      throw StringBufferOverflow("File name [[ " + filename_ +
            " ]] too large. Adjust FN_LEN in cpin.h and parse_cpin.F90 and recompile");
   char *my_fname = (char*) filename_.c_str();
   // Declare variables
   int* protcnt_;
   StateInfo* stateinf_;
   char (*resname)[40];
   int ierr;
   // Check if cnstphe namelist is in the input file
   char txt[FN_LEN] = "cnstphe";
   int lentxt = 7;
   nmlsrc_from_file_(txt,&lentxt,cpinname,&is_cpein_,&ierr);
   if (ierr == 0) {
      // Read limits from input file or use default values
      parse_limits_(cpinname,&ierr,&TITR_RES_C,&TITR_STATES_C,&ATOM_CHRG_C,&MAX_H_COUNT);
      if (ierr == 0) {
         // Allocate variables on heap
         resname = new char[TITR_RES_C+1][40];
         stateinf_ = new StateInfo[TITR_RES_C];
         protcnt_ = new int[TITR_STATES_C];
#ifdef REDOX
         parse_cein_(&trescnt_, protcnt_, stateinf_, resname, my_fname, &is_cpein_, &ierr, &TITR_RES_C, &TITR_STATES_C, &ATOM_CHRG_C);
#else
         parse_cpin_(&trescnt_, protcnt_, stateinf_, resname, my_fname, &is_cpein_, &ierr, &TITR_RES_C, &TITR_STATES_C, &ATOM_CHRG_C);
#endif
      }
   }

   // Error catch
   if (ierr != 0) {
      cerr << "Error: Could not open or parse " << filename_ << " for reading!" << endl;
      if (protcnt_ != nullptr) { delete[] protcnt_; }
      if (stateinf_ != nullptr) { delete[] stateinf_; }
      if (resname != nullptr) { delete[] resname; }
      return ierr;
   }

   for (int i = 0; i <= trescnt_; i++) {
      resname_.push_back(string(resname[i]));
   }

   is_valid_ = true;

   // Now define all of the residues
   for (int i = 0; i < trescnt_; i++) {
      int nstates = stateinf_[i].num_states;
      int firststate = stateinf_[i].first_state;
      vector<int> res_protcnt;
      for (int j = 0; j < nstates; j++)
         res_protcnt.push_back( protcnt_[firststate+j] );
      residues_.push_back( TitratableResidue(resname_[i+1], res_protcnt) );
   }

   if (protcnt_ != nullptr) { delete[] protcnt_; }
   if (stateinf_ != nullptr) { delete[] stateinf_; }
   if (resname != nullptr) { delete[] resname; }
   return 0;
}

int Cpin::Parse(string const& cpinname) {
   return Parse(cpinname.c_str());
}

/// Constructor for TitratableResidue
TitratableResidue::TitratableResidue(string const& resname,
                                     vector<int> const& protcnts) :
resid_(0)
{
   protcnts_ = protcnts;

   // Process the resname
   vector<string> words = split(resname);
   resname_ = words[1]; resid_ = StringToInt(words[2]);

   // Now determine which states are "protonated" and which are "deprotonated"
   int max_prots = protcnts_[0];

   for (size_t i = 1; i < protcnts_.size(); i++)
      max_prots = protcnts_[i] > max_prots ? protcnts_[i] : max_prots;

   for (size_t i = 0; i < protcnts_.size(); i++)
      protonated_.push_back(protcnts_[i] == max_prots);
}
