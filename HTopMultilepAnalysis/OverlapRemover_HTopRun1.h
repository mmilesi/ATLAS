#ifndef HTOPMULTILEPANALYSIS_OVERLAPREMOVER_HTOPRUN1_H
#define HTOPMULTILEPANALYSIS_OVERLAPREMOVER_HTOPRUN1_H

// xAH stuff
#include "xAODAnaHelpers/Algorithm.h"
#include "xAODAnaHelpers/OverlapRemover.h"

class OverlapRemover_HTopRun1 : public OverlapRemover
{
public:

  // this is a standard constructor
  OverlapRemover_HTopRun1 () : OverlapRemover() {};

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode initialize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(OverlapRemover_HTopRun1, 1);
};

#endif
