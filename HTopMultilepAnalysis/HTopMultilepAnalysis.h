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

  std::string m_inContainerName_PreSelectedElectrons;
  std::string m_inContainerName_PreSelectedMuons;
  std::string m_inContainerName_PreSelectedJets;

  bool m_useCutFlow;   
  
  /* to define "Tight" leptons */
  std::string m_TightElectronPID_WP;
  std::string m_TightElectronIso_WP;
  float       m_TightMuonD0sig_cut;
  std::string m_TightMuonIso_WP;
  
  /* BTag WP to define nbjets*/
  std::string m_BTag_WP;

  bool m_useMCForTagAndProbe; // To define tag and probe leptons for RF rate mesurement using MC truth
                              // NB: use it only for pure MC estimate (e.g. ttbar MM closure test)

private:

  bool m_isMC;            //!
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

  /* MM/FF method stuff */
  
  std::map< std::string, TH1D* > m_el_hist_map; //!
  std::map< std::string, TH1D* > m_mu_hist_map; //!
  
  int m_n_el_bins_eta;        //!
  int m_n_el_bins_pt_fr;      //!
  int m_n_el_bins_pt_rr;      //!
  int m_n_mu_bins_eta;        //!
  int m_n_mu_bins_pt_fr;      //!
  int m_n_mu_bins_pt_rr;      //!  
  float m_el_rr_tot;	      //!
  float m_el_fr_tot;	      //! 
  float m_mu_rr_tot;	      //!
  float m_mu_fr_tot;	      //!    
  
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

  virtual EL::StatusCode defineTagAndProbeRFRateVars( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  virtual EL::StatusCode defineTagAndProbeRFRateVars_MC( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  
  virtual EL::StatusCode addChannelDecorations( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  virtual EL::StatusCode fakeWeightCalculator ( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  
  std::vector<double>  calc_weights( std::map< std::string, TH1D* > &histograms, 
				     float pt,    /* NB: internally converts MeV into GeV --> pass pT in MeV!!! */                        
				     float eta, 
				     bool isFakeLep, 
				     int n_bins_eta,
				     int n_bins_pt_fr,
				     int n_bins_pt_rr,
				     float fr_tot,
				     float rr_tot
			            );  
  double calc_final_event_weight( std::string region, double f1, double f2, double r1 = 1.0, double r2 = 1.0 );
  double scaleRateToFactor( double rate );

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepAnalysis, 1);
};

#endif
