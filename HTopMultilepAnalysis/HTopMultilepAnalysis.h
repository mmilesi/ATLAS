#ifndef HTopMultilepAnalysis_HTopMultilepAnalysis_H
#define HTopMultilepAnalysis_HTopMultilepAnalysis_H

// EL include(s):
#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
  #include "xAODRootAccess/Init.h"
  #include "xAODRootAccess/TEvent.h"
  #include "xAODRootAccess/TStore.h"

// ROOT include(s):
#include "TH1D.h"

#include "xAODAnaHelpers/JetHists.h"

namespace CP{
  class ElectronIsolationSelectionTool;
}
class AsgElectronLikelihoodTool;
class AsgElectronIsEMSelector;

namespace TauAnalysisTools{
  class TauSelectionTool;
}

class HTopMultilepAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;
  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!
  int m_eventCounter;     //!

  std::string m_name;
  std::string m_configName;
  bool m_debug;           //!
  bool m_useCutFlow;      //!
  TH1D* m_cutflowHist;    //!
  TH1D* m_cutflowHistW;   //!
  TH1D* m_histEventCount; //! 
  
  int m_cutflow_bin;      //!

  /* for Francesco */
  TH1D* m_totalEvents;    //!  
  TH1D* m_totalEventsW;   //!

  // input containers
  /* Muons */
  std::string m_inContainerName_Muons;    
  /* Electrons */
  std::string m_inContainerName_Electrons;
  /* Jets */
  std::string m_inContainerName_Jets;     
    
  // electron ID stuff
  bool m_doLHPIDCut;  
  bool m_doCutBasedPIDCut;  
    
  // isolation for Muons
  bool         m_doIsolation_Mu;
  bool         m_useRelativeIso_Mu;
  std::string  m_CaloBasedIsoType_Mu;
  float        m_CaloBasedIsoCut_Mu;
  std::string  m_TrackBasedIsoType_Mu;
  float        m_TrackBasedIsoCut_Mu;
  // isolation for Electrons
  bool         m_doIsolation_El;
  bool         m_useRelativeIso_El;
  std::string  m_CaloBasedIsoType_El;
  float        m_CaloBasedIsoCut_El;
  std::string  m_TrackBasedIsoType_El;
  float        m_TrackBasedIsoCut_El;
  
private:

  JetHists* m_jetPlots; //!

  // tools
  CP::ElectronIsolationSelectionTool *m_electronIsolationSelectionTool; //!
  AsgElectronLikelihoodTool* m_LHToolTight2012;     //!
  AsgElectronLikelihoodTool* m_LHToolMedium2012;    //!
  AsgElectronLikelihoodTool* m_LHToolVeryTight2012; //!
  AsgElectronIsEMSelector* m_MediumPP2012; //!
  AsgElectronIsEMSelector* m_TightPP2012;  //!
  TauAnalysisTools::TauSelectionTool         *m_TauSelTool; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  HTopMultilepAnalysis ();
  HTopMultilepAnalysis (std::string name, std::string configName);

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

  virtual EL::StatusCode applyTagAndProbeRFRateMeasurement( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  
  virtual EL::StatusCode addChannelDecorations( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  virtual EL::StatusCode fakeWeightCalculator ( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  std::vector<double> calc_el_weights( float pt, float eta, bool isDataDerived=true, bool isFakeLep=true );  /* NB: internally converts MeV into GeV --> pass pT in MeV!!! */
  std::vector<double> calc_mu_weights( float pt, float eta, bool isDataDerived=true, bool isFakeLep=true );  /* NB: internally converts MeV into GeV --> pass pT in MeV!!! */
  double calc_fake_weight( std::string region, double f1, double f2, double r1 = 1.0, double r2 = 1.0 );
  double scaleFactorToRate( double val );

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepAnalysis, 1);
};

#endif
