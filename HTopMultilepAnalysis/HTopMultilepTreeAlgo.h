#ifndef HTopMultilepAnalysis_HTopMultilepTreeAlgo_H
#define HTopMultilepAnalysis_HTopMultilepTreeAlgo_H

// EL include(s)
#include <EventLoop/Algorithm.h>

// package include(s):
#include <xAODAnaHelpers/TreeAlgo.h>
#include <HTopMultilepAnalysis/HTopMultilepTree.h>

class HTopMultilepTreeAlgo : public TreeAlgo
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.

public:

  // additional data members which are not already in xAH::TreeAlgo.h 
  std::string m_lepContainerName;

  // this is a standard constructor
  HTopMultilepTreeAlgo () : TreeAlgo() {};

  // these are the functions inherited from Algorithm

  // overload only the ones that somehow differ from the original methods in TreeAlgo (e.g., the ones manipulating m_helpTree)
  virtual EL::StatusCode treeInitialize ();     
  virtual EL::StatusCode execute ();            

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();                   
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepTreeAlgo, 1);
};

#endif
