#ifndef HTopMultilepAnalysis_HTopMultilepEventSelector_H
#define HTopMultilepAnalysis_HTopMultilepEventSelector_H

// EL include(s):
#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

// EDM include(s):
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// ROOT include(s):
#include "TH1D.h"

class HTopMultilepEventSelector : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  // configuration variables
  std::string  m_inContainerName_el;	 
  std::string  m_inContainerName_mu;	 
  std::string  m_inContainerName_jets;   
  std::string  m_inContainerName_tau;	 
  std::string  m_outContainerName_lep;

  bool         m_useCutFlow;    

  bool         m_DC14;
  
  bool         m_doMinObjCut;
  bool         m_doMaxObjCut;
  int          m_n_leptons_min;
  int          m_n_leptons_max;
  int          m_n_leptons_with_tau_min;
  int          m_n_jets_min;
  int          m_n_jets_max; 
  int          m_n_bjets_min;
  int          m_n_taus_min;
  float        m_leptons_eta_max;	 
  float        m_leading_lep_pT_min;		       
  float        m_subleading_lep_pT_min; 

  std::string  m_BTag_WP;

  std::string              m_passAuxDecorKeys; 
  std::string              m_failAuxDecorKeys; 
  std::vector<std::string> m_passKeys; 
  std::vector<std::string> m_failKeys; 

private:
  int m_numEvent;           //!
  int m_numObject;          //!
  int m_numEventPass;       //!
  int m_weightNumEventPass; //!
  int m_numObjectPass;      //!

  // cutflow
  TH1D* m_cutflowHist;      //!
  TH1D* m_cutflowHistW;     //!
  int   m_cutflow_bin;      //!


  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist;  //!

  // this is a standard constructor
  HTopMultilepEventSelector ();

  ~HTopMultilepEventSelector();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();

  // added functions not from Algorithm

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepEventSelector, 1);
};

#endif
