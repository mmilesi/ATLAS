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

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

namespace TauAnalysisTools{
  class TauSelectionTool;
}

class HTopMultilepAnalysis : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  // configuration variables
  /* Muons */
  std::string m_inContainerName_Muons;    
  /* Electrons */
  std::string m_inContainerName_Electrons;
  /* Leptons */
  std::string m_inContainerName_Leptons;
  /* Jets */
  std::string m_inContainerName_Jets;     
  /* Taus */
  std::string m_inContainerName_Taus;

  bool m_useCutFlow;   

  bool m_useLH_ElPID;  
  bool m_useCutBased_ElPID;  

private:

  bool m_isDerivation;    //!

  int  m_eventCounter;    //!

  TH1D* m_cutflowHist;    //!
  TH1D* m_cutflowHistW;   //!
  TH1D* m_histEventCount; //! 
  
  int m_cutflow_bin;      //!

  /* for Francesco */
  TH1D* m_totalEvents;    //!  
  TH1D* m_totalEventsW;   //!

  JetHists* m_jetPlots;   //!

  // tools
  TauAnalysisTools::TauSelectionTool    *m_TauSelTool; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  HTopMultilepAnalysis ();

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
