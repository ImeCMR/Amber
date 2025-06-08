#ifndef INC_DSL_H
#define INC_DSL_H
#include <string>
class DataSet {
  public:
    DataSet();
};
/** AmbPDB will not support the DataSetList-related features Cpptraj
  * has for trajectory writes. This class serves as a fill in
  * for the missing functionality.
  */
class DataSetList {
  public:
    DataSetList();
    ~DataSetList();
    DataSet* GetDataSet(std::string const&) const;
};
#endif
