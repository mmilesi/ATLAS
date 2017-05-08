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
#include "TRandom3.h"

namespace MiniNTupMaker {

  struct Branch_Types {
      Branch_Types() : f(-999.0), c(-1), i(-999) {};
      float f;
      char c;
      int i;
      std::vector<float> vec_f;
      std::vector<char> vec_c;
      std::vector<int> vec_i;
  };

  class eventObj {

  public:
    eventObj():
      isSS01(0), isSS12(0),
      dilep_type(0),
      trilep_type(0),
      quadlep_type(0),
      njets(0),nbjets(0),
      isBadTPEvent_SLT(1), /** Be pesimistic */
      weight_event(1.0),
      weight_event_trig_SLT(1.0),
      weight_event_lep(1.0),
      weight_lep_tag(1.0),weight_lep_probe(1.0),
      weight_trig_tag(1.0),weight_trig_probe(1.0)
    { };

    char isSS01;
    char isSS12;
    int  dilep_type;
    int  trilep_type;
    int  quadlep_type;
    int  njets;
    int  nbjets;

    char isBadTPEvent_SLT; /** No T&TM (SLT) leptons found */

    float weight_event;
    float weight_event_trig_SLT;
    float weight_event_lep;

    float weight_lep_tag;
    float weight_lep_probe;
    float weight_trig_tag;
    float weight_trig_probe;

  };

  class jetObj {
  public:
    jetObj():
      pt(-1.0),eta(-999.0),phi(-999.0),E(-1.0),flavor_weight_MV2c10(-999.0),isbtag(0)
      { };

      float pt;
      float eta;
      float phi;
      float E;
      float flavor_weight_MV2c10;
      char  isbtag;
  };

  class leptonObj {

  public:
    leptonObj()
    {

	props.insert( std::make_pair( "Pt", Branch_Types() ) );	                  props["Pt"].f = -1.0;
	props.insert( std::make_pair( "Eta", Branch_Types() ) );                  props["Eta"].f = -999.0;
	props.insert( std::make_pair( "EtaBE2", Branch_Types() ) );	          props["EtaBE2"].f = -999.0;
	props.insert( std::make_pair( "Phi", Branch_Types() ) );	          props["Phi"].f = -999.0;
	props.insert( std::make_pair( "ID", Branch_Types() ) );	                  props["ID"].f = 0.0;
	props.insert( std::make_pair( "Flavour", Branch_Types() ) );	          props["Flavour"].i = 0;
	props.insert( std::make_pair( "Charge", Branch_Types() ) );	          props["Charge"].f = -999.0;
	props.insert( std::make_pair( "sigd0PV", Branch_Types() ) );	          props["sigd0PV"].f = -999.0;
	props.insert( std::make_pair( "Z0SinTheta", Branch_Types() ) );	          props["Z0SinTheta"].f = -999.0;
	props.insert( std::make_pair( "isTightLH", Branch_Types() ) );	          props["isTightLH"].c = 0;
	props.insert( std::make_pair( "isolationLoose", Branch_Types() ) );	  props["isolationLoose"].i = 0;
	props.insert( std::make_pair( "Isolated", Branch_Types() ) );	          props["Isolated"].c = 0;
	props.insert( std::make_pair( "TrackIsoOverPt", Branch_Types() ) );	  props["TrackIsoOverPt"].f = -1.0;
	props.insert( std::make_pair( "CaloIsoOverPt", Branch_Types() ) );	  props["CaloIsoOverPt"].f = -1.0;
	props.insert( std::make_pair( "ptVarcone20", Branch_Types() ) );	  props["ptVarcone20"].f = -1.0;
	props.insert( std::make_pair( "ptVarcone30", Branch_Types() ) );	  props["ptVarcone30"].f = -1.0;
	props.insert( std::make_pair( "topoEtcone20", Branch_Types() ) );	  props["topoEtcone20"].f = -1.0;
	props.insert( std::make_pair( "promptLeptonIso_TagWeight", Branch_Types() ) ); props["promptLeptonIso_TagWeight"].f = 999.0;
	props.insert( std::make_pair( "chargeIDBDTLoose", Branch_Types() ) );	  props["chargeIDBDTLoose"].f = -999.0;
	props.insert( std::make_pair( "chargeIDBDTMedium", Branch_Types() ) );    props["chargeIDBDTMedium"].f = -999.0;
	props.insert( std::make_pair( "chargeIDBDTTight", Branch_Types() ) );	  props["chargeIDBDTTight"].f = -999.0;
	props.insert( std::make_pair( "isTightSelected", Branch_Types() ) );	  props["isTightSelected"].c = 0;
	props.insert( std::make_pair( "isTightSelectedMVA", Branch_Types() ) );   props["isTightSelectedMVA"].c = 0;
	props.insert( std::make_pair( "isTrigMatch", Branch_Types() ) );	  props["isTrigMatch"].c = 0;
	props.insert( std::make_pair( "isTrigMatchDLT", Branch_Types() ) );	  props["isTrigMatchDLT"].c = 0;
	props.insert( std::make_pair( "Rank3L", Branch_Types() ) );	          props["Rank3L"].i = 999; /** Lepton ranking for 3L:
													       0 : the lepton which is OS wrt. the other two
													       1 : the lepton closets in DR to 0
													       2 : the other lepton
													     */
	props.insert( std::make_pair( "truthType", Branch_Types() ) );	          props["truthType"].i = 0;
	props.insert( std::make_pair( "truthOrigin", Branch_Types() ) );	  props["truthOrigin"].i = 0;
	props.insert( std::make_pair( "truthParentType", Branch_Types() ) );	  props["truthParentType"].i = 0;
	props.insert( std::make_pair( "truthParentOrigin", Branch_Types() ) );	  props["truthParentOrigin"].i = 0;
	props.insert( std::make_pair( "truthParentPdgId", Branch_Types() ) );	  props["truthParentPdgId"].i = 0;
	props.insert( std::make_pair( "isPrompt", Branch_Types() ) );	          props["isPrompt"].c = 0;
	props.insert( std::make_pair( "isBrems", Branch_Types() ) );	          props["isBrems"].c = 0;
	props.insert( std::make_pair( "isQMisID", Branch_Types() ) );	          props["isQMisID"].c = 0;
        /** New definition of truth variables */
	props.insert( std::make_pair( "isPrompt_v2", Branch_Types() ) );	      props["isPrompt_v2"].c = 0;
	props.insert( std::make_pair( "isBremsPrompt_v2", Branch_Types() ) );	      props["isBremsPrompt_v2"].c = 0;
	//
	props.insert( std::make_pair( "isQMisID_v2", Branch_Types() ) );	      props["isQMisID_v2"].c = 0;
	props.insert( std::make_pair( "isPromptQMisID_v2", Branch_Types() ) );	      props["isPromptQMisID_v2"].c = 0;
	props.insert( std::make_pair( "isBremsQMisID_v2", Branch_Types() ) );	      props["isBremsQMisID_v2"].c = 0;
	//
	props.insert( std::make_pair( "isPhotonConv_v2", Branch_Types() ) );	      props["isPhotonConv_v2"].c = 0;
	props.insert( std::make_pair( "isPromptPhotonConv_v2", Branch_Types() ) );    props["isPromptPhotonConv_v2"].c = 0;
	props.insert( std::make_pair( "isRadiationPhotonConv_v2", Branch_Types() ) ); props["isRadiationPhotonConv_v2"].c = 0;
	props.insert( std::make_pair( "isOtherPhotonConv_v2", Branch_Types() ) );     props["isOtherPhotonConv_v2"].c = 0;
	//
	props.insert( std::make_pair( "isNonPrompt_v2", Branch_Types() ) );	      props["isNonPrompt_v2"].c = 0;
	//
	props.insert( std::make_pair( "isOtherFake_v2", Branch_Types() ) );	      props["isOtherFake_v2"].c = 0;
	//
	props.insert( std::make_pair( "isUnknown_v2", Branch_Types() ) );	      props["isUnknown_v2"].c = 0;
	//
	props.insert( std::make_pair( "isTagSLT", Branch_Types() ) );	          props["isTagSLT"].c = 0;
	props.insert( std::make_pair( "isTagDLT", Branch_Types() ) );	          props["isTagDLT"].c = 0;
	props.insert( std::make_pair( "deltaRClosestJet", Branch_Types() ) );	  props["deltaRClosestJet"].f = -1.0;
	props.insert( std::make_pair( "deltaRClosestBJet", Branch_Types() ) );    props["deltaRClosestBJet"].f = -1.0;
	props.insert( std::make_pair( "massClosestJet", Branch_Types() ) );	  props["massClosestJet"].f = -1.0;
	props.insert( std::make_pair( "massClosestBJet", Branch_Types() ) );	  props["massClosestBJet"].f = -1.0;
	props.insert( std::make_pair( "SFIDLoose", Branch_Types() ) );	          props["SFIDLoose"].f = 1.0;
	props.insert( std::make_pair( "SFIDTight", Branch_Types() ) );	          props["SFIDTight"].f = 1.0;
	props.insert( std::make_pair( "SFTrigLoose", Branch_Types() ) );	  props["SFTrigLoose"].f = 1.0;
	props.insert( std::make_pair( "SFTrigTight", Branch_Types() ) );	  props["SFTrigTight"].f = 1.0;
	props.insert( std::make_pair( "EffTrigLoose", Branch_Types() ) );	  props["EffTrigLoose"].f = 1.0;
	props.insert( std::make_pair( "EffTrigTight", Branch_Types() ) );	  props["EffTrigTight"].f = 1.0;
	props.insert( std::make_pair( "SFIsoLoose", Branch_Types() ) );	          props["SFIsoLoose"].f = 1.0;
	props.insert( std::make_pair( "SFIsoTight", Branch_Types() ) );	          props["SFIsoTight"].f = 1.0;
	props.insert( std::make_pair( "SFReco", Branch_Types() ) );	          props["SFReco"].f = 1.0;
	props.insert( std::make_pair( "SFTTVA", Branch_Types() ) );	          props["SFTTVA"].f = 1.0;
	props.insert( std::make_pair( "SFObjLoose", Branch_Types() ) );	          props["SFObjLoose"].f = 1.0;
	props.insert( std::make_pair( "SFObjTight", Branch_Types() ) );	          props["SFObjTight"].f = 1.0;

    };

    std::map< std::string, MiniNTupMaker::Branch_Types > props;

  };

  struct SorterEta {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       return  fabs(lep0.get()->props["Eta"].f) > fabs(lep1.get()->props["Eta"].f); /* sort in descending order of |eta| */
    }
  };

  struct SorterPt {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       return  lep0.get()->props["Pt"].f > lep1.get()->props["Pt"].f; /* sort in descending order of pT (get highest pT first) */
    }
  };

  struct SorterTrackIsoOverPt {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       if ( lep0.get()->props["TrackIsoOverPt"].f == lep1.get()->props["TrackIsoOverPt"].f ) { /* if they have same iso (aka 0 ), use pT as criterion */
         return lep0.get()->props["Pt"].f > lep1.get()->props["Pt"].f;
       }
       return  lep0.get()->props["TrackIsoOverPt"].f < lep1.get()->props["TrackIsoOverPt"].f; /* sort in ascending order of trackisooverpt (get more isolated first) */
    }
  };

  struct SorterDistanceClosestBJet {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       return  lep0.get()->props["deltaRClosestBJet"].f > lep1.get()->props["deltaRClosestBJet"].f; /* sort in descending order of DeltaR(lep, closest bjet) (get lep w/ maximal distance first) */
    }
  };

  struct SorterMassClosestBJet {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       return  lep0.get()->props["massClosestBJet"].f > lep1.get()->props["massClosestBJet"].f; /* sort in descending order of M(lep, closest bjet) (get lep w/ maximal mass w/ closest bjet first --> assume it's less likely to be fake) */
    }
  };

  struct Sorter3L {
    bool operator() ( const std::shared_ptr<leptonObj>& lep0, const std::shared_ptr<leptonObj>& lep1 ) const {
       return  lep0.get()->props["Rank3L"].i < lep1.get()->props["Rank3L"].i; /* sort in ascending order of 3L ranking (get lowest rank idx first) */
    }
  };

}


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

  /** Activate if want to define T&P leptons based on truth matching (NB: do this only on TTBar!)  */
  bool m_useTruthTP;

  /** Activate if want to define T&P leptons in the standard way (same criteria applied on data and MC). This is inspired by what done in the SUSY SS analysis */
  bool m_useNominalTP;

  /** Acivate if want to do T&P for 3L. This will ensure the T eletron definition does NOT include the chargeIDBDT selection */
  bool m_do3LTP;

  /** Choose which criterion to use to solve ambiguous case of both T (TM) leptons in T&P Fake SS CR

      Default is "OF" :
      -) for electrons, use only OF events where the muon is tagging the event (--> T & T.M.)
      -) for muons, use SF events, and solve the ambiguity (when necessary) by taking the lowest pT muon as probe
   */
  std::string m_ambiSolvingCrit;

  /** Choose which tight selection method to use in T&P algorithm: "CutBased" or "MVA". Default is "MVA" */
  std::string m_lepSelForTP;

  /** Activate if want to perform reco jet truth matching */
  bool m_jetTruthMatching;

private:

  /** Input TTree */

  TTree*          m_inputNTuple;

  /** Input TTree with the sum of weights */

  TTree*          m_sumWeightsTree;

  /** Histogram containing the total number of generated raw/weighted events */

  TH1F*           m_sumGenEventsHist;

  /** Output TTree (svc) */

  EL::NTupleSvc*  m_outputNTuple;

  /** Input TTree branches which need to be used by the algorithm */

  ULong64_t       m_EventNumber;
  UInt_t          m_RunNumber;
  Int_t           m_RunYear;
  Bool_t          m_passEventCleaning;
  UInt_t          m_mc_channel_number; /** for DATA, mc_channel_number=0 */

  Double_t        m_mcWeightOrg;
  Double_t	  m_pileupEventWeight_090;
  Double_t	  m_MV2c10_70_EventWeight;
  Double_t	  m_JVT_EventWeight;

  Int_t 	  m_dilep_type;
  Int_t 	  m_trilep_type;
  Int_t 	  m_quadlep_type;

  Int_t           m_nJets_OR;
  Int_t           m_nJets_OR_MV2c10_70;
  Int_t           m_nJets_OR_T;
  Int_t           m_nJets_OR_T_MV2c10_70;

  std::vector<std::string> m_LEP_INPUT_VARS  = {"ID/F","Pt/F","E/F","Eta/F","Phi/F","EtaBE2/F","sigd0PV/F","Z0SinTheta/F","isTightLH/B","isMediumLH/B","isLooseLH/B","isTight/B","isMedium/B","isLoose/B","isolationLooseTrackOnly/I","isolationLoose/I","isolationFixedCutTight/I","isolationFixedCutTightTrackOnly/I","isolationFixedCutLoose/I","topoEtcone20/F","ptVarcone20/F","ptVarcone30/F","promptLeptonIso_TagWeight/F","chargeIDBDTLoose/F","chargeIDBDTMedium/F","chargeIDBDTTight/F","isTrigMatch/B","isTrigMatchDLT/B","isPrompt/B","isBrems/B","isQMisID/B","truthType/I","truthOrigin/I","truthParentType/I","truthParentOrigin/I","truthParentPdgId/I","SFIDLoose/F","SFIDTight/F","SFTrigLoose/F","SFTrigTight/F","EffTrigLoose/F","EffTrigTight/F","SFIsoLoose/F","SFIsoTight/F","SFReco/F","SFTTVA/F","SFObjLoose/F","SFObjTight/F"};

  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep0_INPUT_branches;
  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep1_INPUT_branches;
  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep2_INPUT_branches;

  ULong64_t       m_totalEvents;
  Float_t         m_totalEventsWeighted;

  /** Trigger match decision per-lepton (for each chain) - BEFORE overlap removal! */

  // 2015

  std::vector<int> *m_electron_match_HLT_e24_lhmedium_L1EM20VH        = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e60_lhmedium                 = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e120_lhloose                 = nullptr; //!
  std::vector<int> *m_electron_match_HLT_2e12_lhloose_L12EM10VH       = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1 = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e7_medium_mu24               = nullptr; //!
  std::vector<int> *m_muon_match_HLT_mu20_iloose_L1MU15               = nullptr; //!
  std::vector<int> *m_muon_match_HLT_mu18_mu8noL1                     = nullptr; //!
  std::vector<int> *m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1     = nullptr; //!
  std::vector<int> *m_muon_match_HLT_e7_medium_mu24                   = nullptr; //!

  // 2016

  std::vector<int> *m_electron_match_HLT_e26_lhtight_nod0_ivarloose = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e60_lhmedium_nod0          = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e140_lhloose_nod0          = nullptr; //!
  std::vector<int> *m_electron_match_HLT_2e17_lhvloose_nod0         = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e17_lhloose_mu14           = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e17_lhloose_nod0_mu14      = nullptr; //!
  std::vector<int> *m_electron_match_HLT_e7_lhmedium_mu24           = nullptr; //!
  std::vector<int> *m_muon_match_HLT_mu26_ivarmedium                = nullptr; //!
  std::vector<int> *m_muon_match_HLT_mu22_mu8noL1                   = nullptr; //!
  std::vector<int> *m_muon_match_HLT_e17_lhloose_mu14               = nullptr; //!
  std::vector<int> *m_muon_match_HLT_e17_lhloose_nod0_mu14          = nullptr; //!

  // 2015 & 2016

  std::vector<int> *m_muon_match_HLT_mu50 = nullptr; //!

  /** Index of leptons which passed the overlap removal */

  std::vector<char> *m_electron_passOR = nullptr; //!
  std::vector<char> *m_muon_passOR     = nullptr; //!

  /** Some other lepton vector branches */

  std::vector<float> *m_electron_pt	      = nullptr; //!
  std::vector<float> *m_electron_eta	      = nullptr; //!
  std::vector<float> *m_electron_EtaBE2       = nullptr; //!
  std::vector<float> *m_electron_phi	      = nullptr; //!
  std::vector<float> *m_electron_E	      = nullptr; //!
  std::vector<int>   *m_electron_ID           = nullptr; //!
  std::vector<float> *m_electron_sigd0PV      = nullptr; //!
  std::vector<float> *m_electron_z0SinTheta   = nullptr; //!
  std::vector<float> *m_electron_topoetcone20 = nullptr; //!
  std::vector<float> *m_electron_ptvarcone20  = nullptr; //!
  std::vector<int>   *m_electron_truthType    = nullptr; //!
  std::vector<int>   *m_electron_truthOrigin  = nullptr; //!
  std::vector<float> *m_electron_PromptLeptonIso_TagWeight  = nullptr; //!
  std::vector<float> *m_electron_ChargeIDBDTLoose           = nullptr; //!
  std::vector<float> *m_electron_ChargeIDBDTMedium          = nullptr; //!
  std::vector<float> *m_electron_ChargeIDBDTTight           = nullptr; //!

  std::vector<float> *m_muon_pt	              = nullptr; //!
  std::vector<float> *m_muon_eta	      = nullptr; //!
  std::vector<float> *m_muon_phi	      = nullptr; //!
  std::vector<int>   *m_muon_ID               = nullptr; //!
  std::vector<float> *m_muon_sigd0PV	      = nullptr; //!
  std::vector<float> *m_muon_z0SinTheta       = nullptr; //!
  std::vector<float> *m_muon_ptvarcone30      = nullptr; //!
  std::vector<int>   *m_muon_truthType        = nullptr; //!
  std::vector<int>   *m_muon_truthOrigin      = nullptr; //!
  std::vector<float> *m_muon_PromptLeptonIso_TagWeight  = nullptr; //!

  /** Reco jets BEFORE overlap removal */

  std::vector<float>   *m_jet_pt = nullptr;  //!
  std::vector<float>   *m_jet_eta = nullptr; //!
  std::vector<float>   *m_jet_phi = nullptr; //!
  std::vector<float>   *m_jet_E   = nullptr;   //!
  std::vector<float>   *m_jet_flavor_weight_MV2c10 = nullptr;     //!
  std::vector<int>     *m_jet_flavor_truth_label = nullptr;       //!
  std::vector<int>     *m_jet_flavor_truth_label_ghost = nullptr; //!

  /** Indexes of jets that pass overlap removal */

  std::vector<short>   *m_selected_jets   = nullptr;  //!
  std::vector<short>   *m_selected_jets_T = nullptr;  //!

  /** Truth jets */

  std::vector<float>   *m_truth_jet_pt  = nullptr;  //!
  std::vector<float>   *m_truth_jet_eta = nullptr;  //!
  std::vector<float>   *m_truth_jet_phi = nullptr;  //!
  std::vector<float>   *m_truth_jet_e   = nullptr;  //!

  /** Extra branches to be stored in output TTree */

  float     m_weight_event;
  float     m_weight_event_trig_SLT;
  float     m_weight_event_lep;

  float     m_weight_lep_tag;
  float     m_weight_trig_tag;
  float     m_weight_lep_probe;
  float     m_weight_trig_probe;

  char	    m_isSS01;
  char	    m_isSS12;

  char      m_isQMisIDEvent_v2;

  char	    m_is_T_T;
  char	    m_is_T_AntiT;
  char	    m_is_AntiT_T;
  char	    m_is_AntiT_AntiT;
  char      m_is_Tel_AntiTmu;
  char      m_is_Tmu_AntiTel;
  char      m_is_AntiTel_Tmu;
  char      m_is_AntiTmu_Tel;

  char	    m_is_TMVA_TMVA;
  char	    m_is_TMVA_AntiTMVA;
  char	    m_is_AntiTMVA_TMVA;
  char	    m_is_AntiTMVA_AntiTMVA;
  char      m_is_TMVAel_AntiTMVAmu;
  char      m_is_TMVAmu_AntiTMVAel;
  char      m_is_AntiTMVAel_TMVAmu;
  char      m_is_AntiTMVAmu_TMVAel;

  int       m_nmuons;
  int       m_nelectrons;
  int       m_nleptons;

  /** Additional lepton flat branches which do not exist in GN1 */

  std::vector<std::string> m_LEP_OUTPUT_VARS  = {
      "isTightSelectedMVA/B",
      "isTightSelected/B",
      "deltaRClosestJet/F",
      "deltaRClosestBJet/F",
      "isPrompt_v2/B",
      "isBremsPrompt_v2/B",
      "isQMisID_v2/B",
      "isPromptQMisID_v2/B",
      "isBremsQMisID_v2/B",
      "isPhotonConv_v2/B",
      "isPromptPhotonConv_v2/B",
      "isRadiationPhotonConv_v2/B",
      "isOtherPhotonConv_v2/B",
      "isNonPrompt_v2/B",
      "isOtherFake_v2/B",
      "isUnknown_v2/B"
  };

  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep0_OUTPUT_branches;
  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep1_OUTPUT_branches;
  std::map< std::string, MiniNTupMaker::Branch_Types > m_lep2_OUTPUT_branches;

  char m_event_isTrigMatch_DLT;

  /** Some vector branches for leptons after OLR (pT-ordered) */

  std::vector<std::string> m_EL_VEC_VARS  = {
      // "ID/F",
      "Pt/F",
      // "Eta/F",
      "EtaBE2/F",
      "Phi/F",
      "isTightSelected/B",
      "isTightSelectedMVA/B",
      "promptLeptonIso_TagWeight/F",
      // "chargeIDBDTLoose/F",
      // "chargeIDBDTMedium/F",
      "chargeIDBDTTight/F",
      "sigd0PV/F",
      "Z0SinTheta/F",
      "deltaRClosestJet/F",
      "deltaRClosestBJet/F"
      // "ptVarcone20/F",
      // "topoEtcone20/F"
  };
  std::vector<std::string> m_MU_VEC_VARS  = {
      // "ID/F",
      "Pt/F",
      "Eta/F",
      "Phi/F",
      "isTightSelected/B",
      "isTightSelectedMVA/B",
      "promptLeptonIso_TagWeight/F",
      "sigd0PV/F",
      "Z0SinTheta/F",
      "deltaRClosestJet/F",
      "deltaRClosestBJet/F"
      // "ptVarcone30/F"
  };

  std::map< std::string, MiniNTupMaker::Branch_Types > m_electron_OR_branches;
  std::map< std::string, MiniNTupMaker::Branch_Types > m_muon_OR_branches;

  /** Tag & Probe flat,vector branches */

  std::vector<std::string> m_TPS      = { "Tag", "Probe" };
  std::vector<std::string> m_TRIGS    = { "SLT" /*, "DLT"*/ };
  std::vector<std::string> m_TP_VARS  = {
      "Pt/F",
      "Pt/VECF",
      "Eta/F",
      "Eta/VECF",
      "EtaBE2/F",
      "EtaBE2/VECF",
      // "ptVarcone20/F",
      // "ptVarcone30/F",
      // "topoEtcone20/F",
      "sigd0PV/F",
      "Z0SinTheta/F",
      "ID/F",
      "deltaRClosestJet/F",
      "deltaRClosestJet/VECF",
      "deltaRClosestBJet/F",
      "deltaRClosestBJet/VECF",
      "massClosestBJet/F",
      "isTrigMatch/B",
      "isTightSelected/B",
      "isTightSelectedMVA/B",
      "isolationLoose/I",
      "isTightLH/B",
      "promptLeptonIso_TagWeight/F",
      // "chargeIDBDTLoose/F",
      // "chargeIDBDTMedium/F",
      "chargeIDBDTTight/F",
      "truthType/I",
      "truthOrigin/I",
      "truthParentType/I",
      "truthParentOrigin/I",
      "truthParentPdgId/I",
      "isPrompt/B",
      "isBrems/B",
      "isQMisID/B",
      "isPrompt_v2/B",
      "isBremsPrompt_v2/B",
      "isQMisID_v2/B",
      "isPromptQMisID_v2/B",
      "isBremsQMisID_v2/B",
      "isPhotonConv_v2/B",
      "isPromptPhotonConv_v2/B",
      "isRadiationPhotonConv_v2/B",
      "isOtherPhotonConv_v2/B",
      "isNonPrompt_v2/B",
      "isOtherFake_v2/B",
      "isUnknown_v2/B"
  };

  char m_isBadTPEvent_SLT; /** No T&TM (SLT) leptons found */

  std::map< std::string, MiniNTupMaker::Branch_Types > m_TagProbe_branches;

  /** Jets AFTER overlap removal */

  std::vector<float> m_jet_OR_Pt;
  std::vector<float> m_jet_OR_Eta;
  std::vector<float> m_jet_OR_Phi;
  std::vector<float> m_jet_OR_E;
  std::vector<float> m_jet_OR_flavor_weight_MV2c10;
  std::vector<float> m_jet_OR_truthMatch_Pt;
  std::vector<float> m_jet_OR_truthMatch_Eta;
  std::vector<float> m_jet_OR_truthMatch_Phi;
  std::vector<float> m_jet_OR_truthMatch_E;
  std::vector<char>  m_jet_OR_truthMatch_isBJet;
  std::vector<char>  m_jet_OR_truthMatch_isCJet;
  std::vector<char>  m_jet_OR_truthMatch_isLFJet;
  std::vector<char>  m_jet_OR_truthMatch_isGluonJet;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

  /** Other private members */

  unsigned int m_effectiveTotEntries;  //!
  int          m_numEntry;             //!
  float        m_sumGenEvents;         //!
  float        m_sumGenEventsWeighted; //!

  std::shared_ptr<MiniNTupMaker::eventObj>                 m_event;   //!
  std::vector< std::shared_ptr<MiniNTupMaker::leptonObj> > m_leptons; //!
  std::vector< std::shared_ptr<MiniNTupMaker::jetObj> >    m_jets;    //!

  TRandom3* m_rand; //!

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
  EL::StatusCode defineTruthLepFlags ( std::shared_ptr<MiniNTupMaker::leptonObj> lep );
  EL::StatusCode checkIsTightLep ( std::shared_ptr<MiniNTupMaker::leptonObj> lep, const std::string& strategy = "CutBased" );
  EL::StatusCode decorateEvent ();
  EL::StatusCode decorateWeights ();

  EL::StatusCode getPostOLRIndex( int& idx, const unsigned int& pos, const std::string& lep_type );
  EL::StatusCode flatLepVars ();

  EL::StatusCode classify3L( std::vector< std::shared_ptr<MiniNTupMaker::leptonObj> >& leptons );

  EL::StatusCode triggerMatching ();
  EL::StatusCode findClosestJetLep ( const std::string& jetCollection = "jets" );

  /**
    * @brief  Set which lepton is tag and which is probe for the r/f efficiency measurement
    *
  */
  EL::StatusCode defineTagAndProbe ();

  EL::StatusCode fillTPFlatBranches ( std::shared_ptr<MiniNTupMaker::leptonObj> lep, const std::string& trig );
  EL::StatusCode storeLeptonBranches ();
  bool isQMisIDBDTLoose( std::shared_ptr<MiniNTupMaker::leptonObj> lepA, std::shared_ptr<MiniNTupMaker::leptonObj> lepB );
  EL::StatusCode setOutputBranches ();
  EL::StatusCode clearBranches ( const std::string& type );

  EL::StatusCode jetKinematics();
  EL::StatusCode jetTruthMatching();

};

#endif
