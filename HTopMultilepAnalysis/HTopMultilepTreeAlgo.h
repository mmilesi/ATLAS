#ifndef xAODAnaHelpers_HTopMultilepTreeAlgo_H
#define xAODAnaHelpers_HTopMultilepTreeAlgo_H

#include <EventLoop/Algorithm.h>
// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "TTree.h"

// package include(s):
#include "xAODAnaHelpers/TreeAlgo.h" 
#include <HTopMultilepAnalysis/HTopMultilepTree.h>

class HTopMultilepTreeAlgo : public TreeAlgo
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  const std::string m_name;                  
  std::string m_configName;            

  xAOD::TEvent *m_event;               //!
  xAOD::TStore *m_store;               //!

  std::string m_lepContainerName;

private:

  HTopMultilepTree* m_helpTree;        //!

public:
  // this is a standard constructor
  HTopMultilepTreeAlgo ();  //!

  // these are the functions inherited from Algorithm
  
  /* overload only the ones that somehow differ from the original methods in TreeAlgo (e.g., the ones manipulating m_helpTree) */
  
  virtual EL::StatusCode treeInitialize ();                 //!
  virtual EL::StatusCode initialize ();                     //!
  virtual EL::StatusCode execute ();                        //!
  virtual EL::StatusCode finalize ();                       //!

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();                      //!

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepTreeAlgo, 1);                                 //!
};

#endif
