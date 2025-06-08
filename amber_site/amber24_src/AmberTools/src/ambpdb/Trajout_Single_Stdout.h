#ifndef INC_TRAJOUT_SINGLE_STDOUT_H
#define INC_TRAJOUT_SINGLE_STDOUT_H
#include "Trajout_Single.h"
/// Version of Trajout_Single that writes to stdout
class Trajout_Single_Stdout : private Trajout_Single {
  public:
    Trajout_Single_Stdout();
    int PrepareStdoutTrajWrite(ArgList const&, Topology*, CoordinateInfo const&,
                               int, TrajectoryFile::TrajFormatType);
    using Trajout_Single::WriteSingle;
    using Trajout_Single::EndTraj;
};
#endif
