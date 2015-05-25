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

// ROOT include(s):
#include "TH1D.h"

namespace TauAnalysisTools{
  class TauSelectionTool;
}

class HTopMultilepEventSelector : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  xAOD::TEvent *m_event;    //!
  xAOD::TStore *m_store;    //!
  int m_numEvent;           //!
  int m_numObject;          //!
  int m_numEventPass;       //!
  int m_weightNumEventPass; //!
  int m_numObjectPass;      //!

  std::string m_name;
  std::string m_configName;

  bool m_debug;                 //!

  // cutflow
  bool m_useCutFlow;            //!
  TH1D* m_cutflowHist;          //!
  TH1D* m_cutflowHistW;         //!
  int   m_cutflow_bin;          //!

private:

  // tools
  TauAnalysisTools::TauSelectionTool         *m_TauSelTool ; //!

  // configuration variables
  std::string    m_inContainerName_el;      // input container name
  std::string    m_inContainerName_mu;      // input container name
  std::string    m_inContainerName_jets;    // input container name
  std::string    m_inContainerName_tau;     // input container name
  
  bool         m_doMinObjCut;
  bool         m_doMaxObjCut;
  unsigned int m_n_leptons_min;
  unsigned int m_n_leptons_max;
  unsigned int m_n_leptons_with_tau_min;
  unsigned int m_n_jets_min;
  unsigned int m_n_jets_max; 
  unsigned int m_n_bjets_min;
  unsigned int m_n_taus_min;
  bool         m_JetBDTLoose;
  bool         m_JetBDTMedium;
  bool         m_JetBDTTight;
  float        m_leptons_eta_max;	 
  float        m_leading_lep_pT_min;		       
  float        m_subleading_lep_pT_min; 
  float	       m_taus_pT_min;	 
  
  std::string              m_passAuxDecorKeys;  //!
  std::string              m_failAuxDecorKeys;  //!
  std::vector<std::string> m_passKeys;  //!
  std::vector<std::string> m_failKeys;  //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  // this is a standard constructor
  HTopMultilepEventSelector ();
  HTopMultilepEventSelector (std::string name, std::string configName);

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
