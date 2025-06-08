#include "DSL.h"

DataSet::DataSet() {}

DataSetList::DataSetList() {}

DataSetList::~DataSetList() {}

DataSet* DataSetList::GetDataSet(std::string const&) const { return 0; }
