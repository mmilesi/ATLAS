#ifndef HTopMultilepAnalysis_TagAndProbeRFRateMeasurement_H
#define HTopMultilepAnalysis_TagAndProbeRFRateMeasurement_H

// EL include(s):
#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

// EDM include(s):
#include "xAODTruth/TruthEventContainer.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"


// ROOT include(s):
#include "TH1D.h"

namespace CP{
  class ElectronIsolationSelectionTool;
}
class AsgElectronLikelihoodTool;
class AsgElectronIsEMSelector;

class TagAndProbeRFRateMeasurement : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!
  int m_numEvent;         //!
  int m_numObject;        //!
  int m_numEventPass;     //!
  int m_weightNumEventPass; //!
  int m_numObjectPass;    //!

  std::string m_name;
  std::string m_configName;

  bool m_debug;                 //!

  // cutflow
  bool m_useCutFlow;            //!
  TH1D* m_cutflowHist;          //!
  TH1D* m_cutflowHistW;         //!
  int   m_cutflow_bin;          //!

private:
  
  std::multimap<float, const xAOD::IParticle*> m_lepton_map;

  // tools
  CP::ElectronIsolationSelectionTool *m_electronIsolationSelectionTool; //!
  AsgElectronLikelihoodTool* m_LHToolTight2012;     //!
  AsgElectronLikelihoodTool* m_LHToolVeryTight2012; //!
  AsgElectronIsEMSelector* m_TightPP2012;  //!
  
  // configuration variables
  std::string    m_inContainerName_el;      // input container name
  std::string    m_inContainerName_mu;      // input container name
  std::string    m_inContainerName_jet;     // input container name

  // isolation for Muons
  bool         m_useRelativeIso_Mu;
  std::string  m_CaloBasedIsoType_Mu;
  float        m_CaloBasedIsoCut_Mu;
  std::string  m_TrackBasedIsoType_Mu;
  float        m_TrackBasedIsoCut_Mu;
  // isolation for Electrons
  bool         m_useRelativeIso_El;
  std::string  m_CaloBasedIsoType_El;
  float        m_CaloBasedIsoCut_El;
  std::string  m_TrackBasedIsoType_El;
  float        m_TrackBasedIsoCut_El;
    
  bool m_doLHPIDCut;
  bool m_doCutBasedPIDCut;
  
  std::string              m_passAuxDecorKeys;  //!
  std::string              m_failAuxDecorKeys;  //!
  std::vector<std::string> m_passKeys;  //!
  std::vector<std::string> m_failKeys;  //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  TagAndProbeRFRateMeasurement ();
  TagAndProbeRFRateMeasurement (std::string name, std::string configName);

  ~TagAndProbeRFRateMeasurement();

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

  // this is needed to distribute the algorithm to the workers
  ClassDef(TagAndProbeRFRateMeasurement, 1);
};

#endif
