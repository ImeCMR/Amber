#include "Trajout_Single_Stdout.h"
#include "DSL.h"

Trajout_Single_Stdout::Trajout_Single_Stdout() {}

/** Initialize and set up output trajectory for STDOUT write. */
int Trajout_Single_Stdout::PrepareStdoutTrajWrite(ArgList const& argIn, Topology *tparmIn,
                                           CoordinateInfo const& cInfoIn, int nFrames,
                                           TrajectoryFile::TrajFormatType writeFormatIn)
{
  DataSetList blankDsl;
  if (InitTrajout("", argIn, blankDsl, writeFormatIn)) return 1;
  if (SetupTrajWrite(tparmIn, cInfoIn, nFrames)) return 1;
  return 0;
}
