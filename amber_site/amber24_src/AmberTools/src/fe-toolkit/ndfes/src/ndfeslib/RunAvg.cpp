#define _USE_MATH_DEFINES
#include <cmath>

#include "RunAvg.hpp"

  
void ndfes::RunAvg::push_back(double x)
{
  numpts++;

  // See Knuth TAOCP vol 2, 3rd edition, page 232
  if (numpts == 1)
    {
      oldmu = newmu = x;
      oldS = 0.0;
    }
  else
    {
      newmu = oldmu + (x - oldmu)/numpts;
      newS = oldS + (x - oldmu)*(x - newmu);
      oldmu = newmu;
      oldS = newS;
    }
}

