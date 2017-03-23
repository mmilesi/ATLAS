/**
 + * @file   HTopMultilepAnalysis.h
 + * @Author Marco Milesi <marco.milesi@cern.ch>
 + * @brief The actual analysis algorithm. Here the user categorises events, and performs the background estimation.
 + *
 + */

#ifndef HTopMultilepAnalysis_HTopMultilepAnalysis_H
#define HTopMultilepAnalysis_HTopMultilepAnalysis_H

// EL include(s):
#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "AthContainers/ConstDataVector.h"

// EDM include(s):
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"

// ROOT include(s):
#include "TH1D.h"

// external tools include(s):
#include "TauAnalysisTools/TauSelectionTool.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"
#include "xAODAnaHelpers/JetHists.h"


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

  bool m_useOLRDecision;

  std::string m_TightElectronPID_WP;
  std::string m_TightElectronIso_WP;
  float       m_TightElectronD0sig_cut;
  float       m_TightElectronTrkz0sinTheta_cut;

  float       m_TightMuonD0sig_cut;
  float       m_TightMuonTrkz0sinTheta_cut;
  std::string m_TightMuonIso_WP;

  /* to define "Tight" taus */
  std::string m_ConfigPathTightTaus;

  bool m_useLooseAsLoosest;
  bool m_useMediumAsLoosest;
  bool m_vetoMediumNonTight;
  bool m_useMediumAsTightest;
  bool m_useTightAsTightest;

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

  TH1D* m_totalEvents;    //!
  TH1D* m_totalEventsW;   //!

  JetHists* m_jetPlots;   //!

  enum SFType {
    RECO       = 0,
    ISOLATION  = 1,
    ID         = 2,
    TTVA       = 3,
    JVT        = 4,
  };

  // tools
  TauAnalysisTools::TauSelectionTool    *m_TauSelTool; //!

  /* QMisID, MM/FF method stuff */

  std::map< std::string, TH2D* > m_QMisID_hist_map;  //!
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

  EL::StatusCode readQMisIDRates ( const std::string& input_path );
  EL::StatusCode readFakeRates ( const std::string& input_path );

  EL::StatusCode selectTightLeptons( const xAOD::ElectronContainer* electrons, const xAOD::MuonContainer* muons );

  inline  float computeHT( ConstDataVector<xAOD::IParticleContainer>& objects ) {
    float HT(0.0);
    for ( auto obj_itr : objects ) { HT += obj_itr->pt(); }
    return HT;
  };

  virtual EL::StatusCode computeEventLepTrigSF ( const xAOD::EventInfo* eventInfo,
					         const xAOD::IParticleContainer& leptons
					        );
  virtual EL::StatusCode computeEventLepSF( const xAOD::EventInfo* eventInfo,
                                            const xAOD::IParticleContainer& leptons,
					    SFType TYPE
		                          );
  virtual EL::StatusCode computeEventJetSF( const xAOD::EventInfo* eventInfo,
                                            const xAOD::JetContainer* jets,
					    SFType TYPE
		                           );

  virtual EL::StatusCode defineTagAndProbeRFRateVars( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  virtual EL::StatusCode defineTagAndProbeRFRateVars_MC( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );

  virtual EL::StatusCode addChannelDecorations( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );

  virtual EL::StatusCode fakeWeightCalculator ( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons );
  virtual EL::StatusCode QMisIDWeightCalculator ( const xAOD::EventInfo* eventInfo, const xAOD::ElectronContainer* electrons );

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
  double scaleRateToEfficiency( double rate );

  EL::StatusCode calc_QMisID_weights( std::vector<float>& weights, const xAOD::Electron* elA, const xAOD::Electron* elB );
  EL::StatusCode readRatesAndError(TH2D* rate_map, TH1D* proj_X, TH1D* proj_Y,
                                   const float& x, const float& y,
				   float& r, float& r_up, float& r_dn );

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepAnalysis, 1);
};

#endif
