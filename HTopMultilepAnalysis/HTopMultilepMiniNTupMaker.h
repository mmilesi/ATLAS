/**
 * @file   HTopMultilepMiniNTupMaker.h
 * @author Marco Milesi <marco.milesi@cern.ch>
 * @brief  EventLoop algorithm to skim/slim/augment HTop group ntuples
 *
 */

#ifndef HTopMultilepAnalysis_HTopMultilepMiniNTupMaker_H
#define HTopMultilepAnalysis_HTopMultilepMiniNTupMaker_H

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// EL include(s):
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>
#include <EventLoop/Worker.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Job.h>
#include <EventLoop/OutputStream.h>

// ROOT include(s):
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

class eventObj {

public:
  eventObj():
    isMC(0), isSS01(0), isSS12(0),
    dilep(0), trilep(0),
    notightlep(0),
	weight_event(1.0),weight_event_trig(1.0),weight_event_lep(1.0),weight_tag(1.0),weight_probe(1.0)
  { };

  char isMC;
  char isSS01;
  char isSS12;
  char dilep;
  char trilep;
  char notightlep;

  float weight_event;
  float weight_event_trig;
  float weight_event_lep;
  float weight_tag;
  float weight_probe;

};

class leptonObj {

public:
  leptonObj():
    pt(-1.0),eta(-999.0),etaBE2(-999.0),ID(0),flavour(0),charge(-999.0),d0sig(-999.0),z0sintheta(-999.0),
    pid(0),isolated(0),tight(0),trigmatched(0),prompt(0),fake(0),brems(0),tag(0),
    SFIDLoose(1.0),
    SFIDTight(1.0),
    SFTrigLoose(1.0),
    SFTrigTight(1.0),
    EffTrigLoose(0.0),
    EffTrigTight(0.0),
    SFIsoLoose(1.0),
    SFIsoTight(1.0),
    SFReco(1.0),
    SFTTVA(1.0),
    SFObjLoose(1.0),
    SFObjTight(1.0)
  { };

  float pt;
  float eta;
  float etaBE2;
  int ID;
  int flavour;
  float charge;
  float d0sig;
  float z0sintheta;
  char pid;
  char isolated;
  char tight;
  char trigmatched;
  char prompt;
  char fake;
  char brems;
  int  truthType;
  int  truthOrigin;
  char tag;

  float SFIDLoose;
  float SFIDTight;
  float SFTrigLoose;
  float SFTrigTight;
  float EffTrigLoose;
  float EffTrigTight;
  float SFIsoLoose;
  float SFIsoTight;
  float SFReco;
  float SFTTVA;
  float SFObjLoose;
  float SFObjTight;

};


class HTopMultilepMiniNTupMaker : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  /** The name of the output TTree */
  std::string m_outputNTupName;

  std::string m_outputNTupStreamName;

  /** A comma-separated list of input branches to be activated */
  std::string m_inputBranches;

  /** Perform event skimming through an EL::AlgSelect svc, or do it manually */
  bool         m_useAlgSelect;

  /** Add an output stream (aka, directory) for the histogram containing the total number of generated raw/weighted events */
  bool         m_addStreamEventsHist;

private:

  /** Input TTree */

  TTree*          m_inputNTuple;

  /** Input TTree with the sum of weights */

  TTree*          m_sumWeightsTree;

  /** Histogram containing the total number of generated raw/weighted events */

  TH1F*           m_sumGenEventsHist;

  /** Output TTree (svc) */

  EL::NTupleSvc*  m_outputNTuple;

  /** Input TTree branches whcih need to be used by the algorithm */

  ULong64_t       m_EventNumber;
  UInt_t          m_RunNumber;
  Bool_t          m_passEventCleaning;
  UInt_t          m_mc_channel_number; /** for DATA, mc_channel_number=0 */

  Double_t        m_mcWeightOrg;
  Double_t	  m_pileupEventWeight_090;
  Double_t	  m_MV2c20_77_EventWeight;
  Double_t	  m_JVT_EventWeight;

  Int_t 	  m_dilep_type;
  Int_t 	  m_trilep_type;
  Int_t           m_nJets_OR;
  Int_t           m_nJets_OR_MV2c20_77;
  Int_t           m_nJets_OR_T;
  Int_t           m_nJets_OR_T_MV2c20_77;

  UInt_t          m_HLT_mu20_iloose_L1MU15;
  Float_t	  m_HLT_mu20_iloose_L1MU15_PS;
  UInt_t	  m_HLT_mu50;
  Float_t	  m_HLT_mu50_PS;
  UInt_t	  m_HLT_e24_lhmedium_L1EM18VH;
  Float_t	  m_HLT_e24_lhmedium_L1EM18VH_PS;
  UInt_t	  m_HLT_e24_lhmedium_L1EM20VH;
  Float_t	  m_HLT_e24_lhmedium_L1EM20VH_PS;
  UInt_t	  m_HLT_e60_lhmedium;
  Float_t	  m_HLT_e60_lhmedium_PS;
  UInt_t	  m_HLT_e120_lhloose;
  Float_t	  m_HLT_e120_lhloose_PS;

  Float_t	  m_lep_ID_0;
  Float_t	  m_lep_Pt_0;
  Float_t	  m_lep_E_0;
  Float_t	  m_lep_Eta_0;
  Float_t	  m_lep_Phi_0;
  Float_t	  m_lep_EtaBE2_0;
  Float_t	  m_lep_sigd0PV_0;
  Float_t	  m_lep_Z0SinTheta_0;
  Char_t	  m_lep_isTightLH_0;
  Char_t	  m_lep_isMediumLH_0;
  Char_t	  m_lep_isLooseLH_0;
  Char_t	  m_lep_isTight_0;
  Char_t	  m_lep_isMedium_0;
  Char_t	  m_lep_isLoose_0;
  Int_t 	  m_lep_isolationLooseTrackOnly_0;
  Int_t 	  m_lep_isolationLoose_0;
  Int_t 	  m_lep_isolationFixedCutTight_0;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_0;
  Int_t 	  m_lep_isolationFixedCutLoose_0;
  Char_t	  m_lep_isTrigMatch_0;
  Char_t	  m_lep_isPrompt_0;
  Char_t	  m_lep_isBremsElec_0;
  Char_t	  m_lep_isFakeLep_0;
  Int_t	          m_lep_truthType_0;
  Int_t	          m_lep_truthOrigin_0;
  Float_t	  m_lep_SFIDLoose_0;
  Float_t	  m_lep_SFIDTight_0;
  Float_t	  m_lep_SFTrigLoose_0;
  Float_t	  m_lep_SFTrigTight_0;
  Float_t         m_lep_EffTrigLoose_0;
  Float_t         m_lep_EffTrigTight_0;
  Float_t	  m_lep_SFIsoLoose_0;
  Float_t	  m_lep_SFIsoTight_0;
  Float_t	  m_lep_SFReco_0;
  Float_t	  m_lep_SFTTVA_0;
  Float_t	  m_lep_SFObjLoose_0;
  Float_t	  m_lep_SFObjTight_0;

  Float_t	  m_lep_ID_1;
  Float_t	  m_lep_Pt_1;
  Float_t	  m_lep_E_1;
  Float_t	  m_lep_Eta_1;
  Float_t	  m_lep_Phi_1;
  Float_t	  m_lep_EtaBE2_1;
  Float_t	  m_lep_sigd0PV_1;
  Float_t	  m_lep_Z0SinTheta_1;
  Char_t	  m_lep_isTightLH_1;
  Char_t	  m_lep_isMediumLH_1;
  Char_t	  m_lep_isLooseLH_1;
  Char_t	  m_lep_isTight_1;
  Char_t	  m_lep_isMedium_1;
  Char_t	  m_lep_isLoose_1;
  Int_t 	  m_lep_isolationLooseTrackOnly_1;
  Int_t 	  m_lep_isolationLoose_1;
  Int_t 	  m_lep_isolationFixedCutTight_1;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_1;
  Int_t 	  m_lep_isolationFixedCutLoose_1;
  Char_t	  m_lep_isTrigMatch_1;
  Char_t	  m_lep_isPrompt_1;
  Char_t	  m_lep_isBremsElec_1;
  Char_t	  m_lep_isFakeLep_1;
  Int_t	          m_lep_truthType_1;
  Int_t	          m_lep_truthOrigin_1;
  Float_t	  m_lep_SFIDLoose_1;
  Float_t	  m_lep_SFIDTight_1;
  Float_t	  m_lep_SFTrigLoose_1;
  Float_t	  m_lep_SFTrigTight_1;
  Float_t         m_lep_EffTrigLoose_1;
  Float_t         m_lep_EffTrigTight_1;
  Float_t	  m_lep_SFIsoLoose_1;
  Float_t	  m_lep_SFIsoTight_1;
  Float_t	  m_lep_SFReco_1;
  Float_t	  m_lep_SFTTVA_1;
  Float_t	  m_lep_SFObjLoose_1;
  Float_t	  m_lep_SFObjTight_1;

  Float_t	  m_lep_ID_2;
  Float_t	  m_lep_Pt_2;
  Float_t	  m_lep_E_2;
  Float_t	  m_lep_Eta_2;
  Float_t	  m_lep_Phi_2;
  Float_t	  m_lep_EtaBE2_2;
  Float_t	  m_lep_sigd0PV_2;
  Float_t	  m_lep_Z0SinTheta_2;
  Char_t	  m_lep_isTightLH_2;
  Char_t	  m_lep_isMediumLH_2;
  Char_t	  m_lep_isLooseLH_2;
  Char_t	  m_lep_isTight_2;
  Char_t	  m_lep_isMedium_2;
  Char_t	  m_lep_isLoose_2;
  Int_t 	  m_lep_isolationLooseTrackOnly_2;
  Int_t 	  m_lep_isolationLoose_2;
  Int_t 	  m_lep_isolationFixedCutTight_2;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_2;
  Int_t 	  m_lep_isolationFixedCutLoose_2;
  Char_t	  m_lep_isTrigMatch_2;
  Char_t	  m_lep_isPrompt_2;
  Char_t	  m_lep_isBremsElec_2;
  Char_t	  m_lep_isFakeLep_2;
  Int_t	          m_lep_truthType_2;
  Int_t	          m_lep_truthOrigin_2;
  Float_t	  m_lep_SFIDLoose_2;
  Float_t	  m_lep_SFIDTight_2;
  Float_t	  m_lep_SFTrigLoose_2;
  Float_t	  m_lep_SFTrigTight_2;
  Float_t         m_lep_EffTrigLoose_2;
  Float_t         m_lep_EffTrigTight_2;
  Float_t	  m_lep_SFIsoLoose_2;
  Float_t	  m_lep_SFIsoTight_2;
  Float_t	  m_lep_SFReco_2;
  Float_t	  m_lep_SFTTVA_2;
  Float_t	  m_lep_SFObjLoose_2;
  Float_t	  m_lep_SFObjTight_2;

  ULong64_t       m_totalEvents;
  Float_t         m_totalEventsWeighted;

  /** Reco jets BEFORE overlap removal */

  std::vector<float>   *m_jet_pt;
  std::vector<float>   *m_jet_eta;
  std::vector<float>   *m_jet_phi;
  std::vector<float>   *m_jet_E;
  std::vector<int>     *m_jet_flavor_truth_label;
  std::vector<int>     *m_jet_flavor_truth_label_ghost;

  /** Indexes of jets that pass overlap removal */

  std::vector<short>   *selected_jets;
  std::vector<short>   *selected_jets_T;

  /** Truth jets */

  std::vector<float>   *m_truth_jet_pt;
  std::vector<float>   *m_truth_jet_eta;
  std::vector<float>   *m_truth_jet_phi;
  std::vector<float>   *m_truth_jet_e;

  /** Extra branches to be stored in output TTree */

  char      m_isMC;

  float     m_weight_event;
  float     m_weight_event_trig;
  float     m_weight_event_lep;
  float     m_weight_tag;
  float     m_weight_probe;

  char	    m_isSS01;
  char	    m_isSS12;

  char	    m_is_T_T;
  char	    m_is_T_AntiT;
  char	    m_is_AntiT_T;
  char	    m_is_AntiT_AntiT;
  char      m_is_Tel_AntiTmu;
  char      m_is_Tmu_AntiTel;
  char      m_is_AntiTel_Tmu;
  char      m_is_AntiTmu_Tel;

  int       m_nmuons;
  int       m_nelectrons;
  int       m_nleptons;

  char	    m_lep_isTightSelected_0;
  char	    m_lep_isTightSelected_1;
  char	    m_lep_isTightSelected_2;

  float     m_lep_Tag_Pt;
  float     m_lep_Tag_Eta;
  float     m_lep_Tag_EtaBE2;
  float     m_lep_Tag_sigd0PV;
  float     m_lep_Tag_Z0SinTheta;
  float	    m_lep_Tag_ID;
  char	    m_lep_Tag_isTrigMatch;
  char	    m_lep_Tag_isTightSelected;
  char	    m_lep_Tag_isPrompt;
  char	    m_lep_Tag_isBremsElec;
  char	    m_lep_Tag_isFakeLep;
  int	    m_lep_Tag_truthType;
  int	    m_lep_Tag_truthOrigin;

  float     m_lep_Probe_Pt;
  float     m_lep_Probe_Eta;
  float     m_lep_Probe_EtaBE2;
  float     m_lep_Probe_sigd0PV;
  float     m_lep_Probe_Z0SinTheta;
  float	    m_lep_Probe_ID;
  char	    m_lep_Probe_isTrigMatch;
  char	    m_lep_Probe_isTightSelected;
  char	    m_lep_Probe_isPrompt;
  char	    m_lep_Probe_isBremsElec;
  char	    m_lep_Probe_isFakeLep;
  int	    m_lep_Probe_truthType;
  int	    m_lep_Probe_truthOrigin;

  std::vector<float> m_lep_Pt;
  std::vector<float> m_lep_Eta;
  std::vector<float> m_lep_EtaBE2;

  /** Jets AFTER overlap removal */

  std::vector<float> m_jet_OLR_Pt;
  std::vector<float> m_jet_OLR_Eta;
  std::vector<float> m_jet_OLR_Phi;
  std::vector<float> m_jet_OLR_E;
  std::vector<float> m_jet_OLR_truthMatch_Pt;
  std::vector<float> m_jet_OLR_truthMatch_Eta;
  std::vector<float> m_jet_OLR_truthMatch_Phi;
  std::vector<float> m_jet_OLR_truthMatch_E;
  std::vector<char>  m_jet_OLR_truthMatch_isBJet;
  std::vector<char>  m_jet_OLR_truthMatch_isCJet;
  std::vector<char>  m_jet_OLR_truthMatch_isLFJet;
  std::vector<char>  m_jet_OLR_truthMatch_isGluonJet;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

  /** Other private members */

  unsigned int m_numEntry;   //!
  float        m_sumGenEvents;  //!
  float        m_sumGenEventsWeighted; //!

  std::shared_ptr<eventObj>                 m_event;   //!
  std::vector< std::shared_ptr<leptonObj> > m_leptons; //!

public:

  // this is a standard constructor
  HTopMultilepMiniNTupMaker (std::string className = "HTopMultilepMiniNTupMaker");

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

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepMiniNTupMaker, 1);

private:

  EL::StatusCode enableSelectedBranches ();
  EL::StatusCode checkIsTightLep( std::shared_ptr<leptonObj> lep );
  EL::StatusCode decorateEvent ();
  EL::StatusCode decorateWeights ();

  /**
    * @brief  Set which lepton is tag and which is probe for the r/f efficiency measurement
    *
    * The first lepton found that is tight && trigger-matched will be the tag, the other the probe
    * Note that the lepton container is sorted, so in case both are T & TM, the leading will be the tag and the subleading the probe
    * If none satisfies the above, choose the leading as tight and the subleading as probe, but flag the event as "bad"
    *
    * Use the same logic for OS and SS events
    *
  */
  EL::StatusCode defineTagAndProbe ();

  EL::StatusCode setOutputBranches ();

  EL::StatusCode clearBranches ();

  EL::StatusCode jetTruthMatching();

};

#endif
