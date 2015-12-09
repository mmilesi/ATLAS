#ifndef HTopMultilepAnalysis_HTopMultilepTree_H
#define HTopMultilepAnalysis_HTopMultilepTree_H

// package include(s):
#include "xAODAnaHelpers/HelpTreeBase.h"

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODEgamma/Electron.h"
#include "xAODJet/Jet.h"
#include "xAODTau/TauJet.h"

// Infrastructure include(s):
#include "xAODRootAccess/TEvent.h"

// ROOT include(s):
#include "TTree.h"
#include "TFile.h"

class HTopMultilepTree : public HelpTreeBase
{

  private:

    /* event variables*/
    int           	m_is_mc;
    float         	m_ystar;
    unsigned int  	m_categoryFlag;
    unsigned int        m_nBjets_MV2c20_Fix70;
    int  	  	m_isSS01;
    int  	  	m_isSS12;
    std::vector<double> m_MMWeight;
    std::vector<double> m_FFWeight;
    float  		m_mll01;
    float  		m_mll02;
    float  		m_mll12;
    float  		m_mlll012;
    float               m_mT_lep0MET;
    float               m_mT_lep1MET;
    int  		m_isTT;
    int  		m_isTL;
    int  		m_isLT;
    int  		m_isLL;
    int  		m_isTM;
    int  		m_isMT;
    int  		m_isMM;
    int  		m_isTelLmu;
    int  		m_isLelTmu;
    int  		m_isTmuLel;
    int  		m_isLmuTel;
    int  		m_isTelMmu;
    int  		m_isMelTmu;
    int  		m_isTmuMel;
    int  		m_isMmuTel;
    int                 m_isNonTightEvent;
    int                 m_isProbeElEvent;
    int                 m_isProbeMuEvent;

    /* jet variables */
    std::vector<float> m_jet_m;

    /* muon variables */
    std::vector<int> m_muon_isTight;
    std::vector<int> m_muon_isMedium;
    std::vector<int> m_muon_isOS;
    std::vector<int> m_muon_isClosestSS;
    std::vector<int> m_muon_isTag;

    std::vector<int> m_muon_isTruthMatched;
    std::vector<int> m_muon_isChFlip;
    std::vector<int> m_muon_isBrem;
    std::vector<int> m_muon_truthType;
    std::vector<int> m_muon_truthPdgId;
    std::vector<int> m_muon_truthOrigin;
    std::vector<int> m_muon_truthStatus;
    std::vector<int> m_muon_ancestorTruthType;
    std::vector<int> m_muon_ancestorTruthPdgId;
    std::vector<int> m_muon_ancestorTruthOrigin;
    std::vector<int> m_muon_ancestorTruthStatus;

    /* muon TAG variables */
    std::vector<float> m_muon_tag_pt;
    std::vector<float> m_muon_tag_eta;
    std::vector<float> m_muon_tag_trkd0sig;
    std::vector<float> m_muon_tag_ptvarcone30;
    std::vector<float> m_muon_tag_topoetcone20;
    std::vector<int>   m_muon_tag_isIsolated_LooseTrackOnly;
    std::vector<int>   m_muon_tag_isIsolated_Loose;
    std::vector<int>   m_muon_tag_isIsolated_Tight;
    std::vector<int>   m_muon_tag_isIsolated_Gradient;
    std::vector<int>   m_muon_tag_isIsolated_GradientLoose;
    std::vector<int>   m_muon_tag_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_muon_tag_isIsolated_UserDefinedCut;
    std::vector<int>   m_muon_tag_isTight;
    std::vector<int>   m_muon_tag_isMedium;
    std::vector<int>   m_muon_tag_isTruthMatched;
    std::vector<int>   m_muon_tag_isChFlip;
    std::vector<int>   m_muon_tag_isBrem;
    std::vector<int>   m_muon_tag_truthType;
    std::vector<int>   m_muon_tag_truthPdgId;
    std::vector<int>   m_muon_tag_truthOrigin;
    std::vector<int>   m_muon_tag_truthStatus;
    std::vector<int>   m_muon_tag_ancestorTruthType;
    std::vector<int>   m_muon_tag_ancestorTruthPdgId;
    std::vector<int>   m_muon_tag_ancestorTruthOrigin;
    std::vector<int>   m_muon_tag_ancestorTruthStatus;
    /* muon PROBE variables */
    std::vector<float> m_muon_probe_pt;
    std::vector<float> m_muon_probe_eta;
    std::vector<float> m_muon_probe_trkd0sig;
    std::vector<float> m_muon_probe_ptvarcone30;
    std::vector<float> m_muon_probe_topoetcone20;
    std::vector<int>   m_muon_probe_isIsolated_LooseTrackOnly;
    std::vector<int>   m_muon_probe_isIsolated_Loose;
    std::vector<int>   m_muon_probe_isIsolated_Tight;
    std::vector<int>   m_muon_probe_isIsolated_Gradient;
    std::vector<int>   m_muon_probe_isIsolated_GradientLoose;
    std::vector<int>   m_muon_probe_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_muon_probe_isIsolated_UserDefinedCut;
    std::vector<int>   m_muon_probe_isTight;
    std::vector<int>   m_muon_probe_isMedium;
    std::vector<int>   m_muon_probe_isTruthMatched;
    std::vector<int>   m_muon_probe_isChFlip;
    std::vector<int>   m_muon_probe_isBrem;
    std::vector<int>   m_muon_probe_truthType;
    std::vector<int>   m_muon_probe_truthPdgId;
    std::vector<int>   m_muon_probe_truthOrigin;
    std::vector<int>   m_muon_probe_truthStatus;
    std::vector<int>   m_muon_probe_ancestorTruthType;
    std::vector<int>   m_muon_probe_ancestorTruthPdgId;
    std::vector<int>   m_muon_probe_ancestorTruthOrigin;
    std::vector<int>   m_muon_probe_ancestorTruthStatus;

    /* electron variables */
    std::vector<float> m_electron_calo_eta;
    std::vector<int>   m_electron_crack;
    std::vector<int>   m_electron_isTight;
    std::vector<int>   m_electron_isMedium;
    std::vector<int>   m_electron_isOS;
    std::vector<int>   m_electron_isClosestSS;
    std::vector<int>   m_electron_isTag;
    std::vector<int>   m_electron_isTruthMatched;
    std::vector<int>   m_electron_isChFlip;
    std::vector<int>   m_electron_isBrem;
    std::vector<int>   m_electron_truthType;
    std::vector<int>   m_electron_truthPdgId;
    std::vector<int>   m_electron_truthOrigin;
    std::vector<int>   m_electron_truthStatus;
    std::vector<int>   m_electron_ancestorTruthType;
    std::vector<int>   m_electron_ancestorTruthPdgId;
    std::vector<int>   m_electron_ancestorTruthOrigin;
    std::vector<int>   m_electron_ancestorTruthStatus;
    /* electron TAG variables */
    std::vector<float> m_electron_tag_pt;
    std::vector<float> m_electron_tag_eta;
    std::vector<int>   m_electron_tag_LHVeryLoose;
    std::vector<int>   m_electron_tag_LHLoose;
    std::vector<int>   m_electron_tag_LHMedium;
    std::vector<int>   m_electron_tag_LHTight;
    std::vector<int>   m_electron_tag_IsEMLoose;
    std::vector<int>   m_electron_tag_IsEMMedium;
    std::vector<int>   m_electron_tag_IsEMTight;
    std::vector<float> m_electron_tag_ptvarcone20;
    std::vector<float> m_electron_tag_topoetcone20;
    std::vector<int>   m_electron_tag_isIsolated_LooseTrackOnly;
    std::vector<int>   m_electron_tag_isIsolated_Loose;
    std::vector<int>   m_electron_tag_isIsolated_Tight;
    std::vector<int>   m_electron_tag_isIsolated_Gradient;
    std::vector<int>   m_electron_tag_isIsolated_GradientLoose;
    std::vector<int>   m_electron_tag_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_electron_tag_isIsolated_UserDefinedCut;
    std::vector<int>   m_electron_tag_isTight;
    std::vector<int>   m_electron_tag_isMedium;
    std::vector<int>   m_electron_tag_isTruthMatched;
    std::vector<int>   m_electron_tag_isChFlip;
    std::vector<int>   m_electron_tag_isBrem;
    std::vector<int>   m_electron_tag_truthType;
    std::vector<int>   m_electron_tag_truthPdgId;
    std::vector<int>   m_electron_tag_truthOrigin;
    std::vector<int>   m_electron_tag_truthStatus;
    std::vector<int>   m_electron_tag_ancestorTruthType;
    std::vector<int>   m_electron_tag_ancestorTruthPdgId;
    std::vector<int>   m_electron_tag_ancestorTruthOrigin;
    std::vector<int>   m_electron_tag_ancestorTruthStatus;
    /* electron PROBE variables */
    std::vector<float> m_electron_probe_pt;
    std::vector<float> m_electron_probe_eta;
    std::vector<int>   m_electron_probe_LHVeryLoose;
    std::vector<int>   m_electron_probe_LHLoose;
    std::vector<int>   m_electron_probe_LHMedium;
    std::vector<int>   m_electron_probe_LHTight;
    std::vector<int>   m_electron_probe_IsEMLoose;
    std::vector<int>   m_electron_probe_IsEMMedium;
    std::vector<int>   m_electron_probe_IsEMTight;
    std::vector<float> m_electron_probe_ptvarcone20;
    std::vector<float> m_electron_probe_topoetcone20;
    std::vector<int>   m_electron_probe_isIsolated_LooseTrackOnly;
    std::vector<int>   m_electron_probe_isIsolated_Loose;
    std::vector<int>   m_electron_probe_isIsolated_Tight;
    std::vector<int>   m_electron_probe_isIsolated_Gradient;
    std::vector<int>   m_electron_probe_isIsolated_GradientLoose;
    std::vector<int>   m_electron_probe_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_electron_probe_isIsolated_UserDefinedCut;
    std::vector<int>   m_electron_probe_isTight;
    std::vector<int>   m_electron_probe_isMedium;
    std::vector<int>   m_electron_probe_isTruthMatched;
    std::vector<int>   m_electron_probe_isChFlip;
    std::vector<int>   m_electron_probe_isBrem;
    std::vector<int>   m_electron_probe_truthType;
    std::vector<int>   m_electron_probe_truthPdgId;
    std::vector<int>   m_electron_probe_truthOrigin;
    std::vector<int>   m_electron_probe_truthStatus;
    std::vector<int>   m_electron_probe_ancestorTruthType;
    std::vector<int>   m_electron_probe_ancestorTruthPdgId;
    std::vector<int>   m_electron_probe_ancestorTruthOrigin;
    std::vector<int>   m_electron_probe_ancestorTruthStatus;

    /* lepton variables */
    int                m_nlep;
    std::vector<float> m_lepton_pt;
    std::vector<float> m_lepton_phi;
    std::vector<float> m_lepton_eta;
    std::vector<float> m_lepton_m;
    std::vector<float> m_lepton_charge;
    std::vector<int>   m_lepton_flavour;
    std::vector<int>   m_lepton_isTrigMatched;
    std::vector<int>   m_lepton_isIsolated_LooseTrackOnly;
    std::vector<int>   m_lepton_isIsolated_Loose;
    std::vector<int>   m_lepton_isIsolated_Tight;
    std::vector<int>   m_lepton_isIsolated_Gradient;
    std::vector<int>   m_lepton_isIsolated_GradientLoose;
    std::vector<int>   m_lepton_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_lepton_isIsolated_UserDefinedCut;
    std::vector<int>   m_lepton_isTight;
    std::vector<int>   m_lepton_isMedium;
    std::vector<int>   m_lepton_isOS;
    std::vector<int>   m_lepton_isClosestSS;
    std::vector<int>   m_lepton_isTag;
    std::vector<int>   m_lepton_isTruthMatched;
    std::vector<int>   m_lepton_isChFlip;
    std::vector<int>   m_lepton_isBrem;
    std::vector<int>   m_lepton_truthType;
    std::vector<int>   m_lepton_truthPdgId;
    std::vector<int>   m_lepton_truthOrigin;
    std::vector<int>   m_lepton_truthStatus;
    std::vector<int>   m_lepton_ancestorTruthType;
    std::vector<int>   m_lepton_ancestorTruthPdgId;
    std::vector<int>   m_lepton_ancestorTruthOrigin;
    std::vector<int>   m_lepton_ancestorTruthStatus;

    /* lepton TAG variables */
    std::vector<float> m_lepton_tag_pt;
    std::vector<float> m_lepton_tag_eta;
    std::vector<int>   m_lepton_tag_flavour;
    std::vector<int>   m_lepton_tag_charge;
    std::vector<int>   m_lepton_tag_isIsolated_LooseTrackOnly;
    std::vector<int>   m_lepton_tag_isIsolated_Loose;
    std::vector<int>   m_lepton_tag_isIsolated_Tight;
    std::vector<int>   m_lepton_tag_isIsolated_Gradient;
    std::vector<int>   m_lepton_tag_isIsolated_GradientLoose;
    std::vector<int>   m_lepton_tag_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_lepton_tag_isIsolated_UserDefinedCut;
    std::vector<int>   m_lepton_tag_isTrigMatched;
    std::vector<int>   m_lepton_tag_isTight;
    std::vector<int>   m_lepton_tag_isMedium;
    std::vector<int>   m_lepton_tag_isTruthMatched;
    std::vector<int>   m_lepton_tag_isChFlip;
    std::vector<int>   m_lepton_tag_isBrem;
    std::vector<int>   m_lepton_tag_truthType;
    std::vector<int>   m_lepton_tag_truthPdgId;
    std::vector<int>   m_lepton_tag_truthOrigin;
    std::vector<int>   m_lepton_tag_truthStatus;
    std::vector<int>   m_lepton_tag_ancestorTruthType;
    std::vector<int>   m_lepton_tag_ancestorTruthPdgId;
    std::vector<int>   m_lepton_tag_ancestorTruthOrigin;
    std::vector<int>   m_lepton_tag_ancestorTruthStatus;
    /* lepton PROBE variables */
    std::vector<float> m_lepton_probe_pt;
    std::vector<float> m_lepton_probe_eta;
    std::vector<int>   m_lepton_probe_flavour;
    std::vector<int>   m_lepton_probe_charge;
    std::vector<int>   m_lepton_probe_isIsolated_LooseTrackOnly;
    std::vector<int>   m_lepton_probe_isIsolated_Loose;
    std::vector<int>   m_lepton_probe_isIsolated_Tight;
    std::vector<int>   m_lepton_probe_isIsolated_Gradient;
    std::vector<int>   m_lepton_probe_isIsolated_GradientLoose;
    std::vector<int>   m_lepton_probe_isIsolated_UserDefinedFixEfficiency;
    std::vector<int>   m_lepton_probe_isIsolated_UserDefinedCut;
    std::vector<int>   m_lepton_probe_isTrigMatched;
    std::vector<int>   m_lepton_probe_isTight;
    std::vector<int>   m_lepton_probe_isMedium;
    std::vector<int>   m_lepton_probe_isTruthMatched;
    std::vector<int>   m_lepton_probe_isChFlip;
    std::vector<int>   m_lepton_probe_isBrem;
    std::vector<int>   m_lepton_probe_truthType;
    std::vector<int>   m_lepton_probe_truthPdgId;
    std::vector<int>   m_lepton_probe_truthOrigin;
    std::vector<int>   m_lepton_probe_truthStatus;
    std::vector<int>   m_lepton_probe_ancestorTruthType;
    std::vector<int>   m_lepton_probe_ancestorTruthPdgId;
    std::vector<int>   m_lepton_probe_ancestorTruthOrigin;
    std::vector<int>   m_lepton_probe_ancestorTruthStatus;

    /* tau variables */
    std::vector<int>   m_tau_isBDTTight;

  public:

    HTopMultilepTree( TTree* tree, TFile* file, xAOD::TEvent* event, xAOD::TStore* store, const float units = 1e3, bool debug = false, bool DC14 = false );
    ~HTopMultilepTree();

    void AddEventUser(const std::string detailStrUser = "");
    /*void AddTriggerUser(const std::string detailStrUser = "");*/
    void AddMuonsUser(const std::string detailStrUser = "");
    void AddElectronsUser(const std::string detailStrUser = "");
    void AddJetsUser(const std::string detailStrUser = "", const std::string = "jet" );
    void AddTausUser(const std::string detailStrUser = "");
    /*void AddMETUser(const std::string detailStrUser = "");*/

    void AddLeptons();

    void ClearEventUser();
    /*void ClearTriggerUser();*/
    void ClearMuonsUser();
    void ClearElectronsUser();
    void ClearJetsUser( const std::string = "jet" );
    void ClearTausUser();
    /*void ClearMETUser();*/
    void ClearLeptons();

    void FillEventUser( const xAOD::EventInfo* );
    /*void FillTriggerUser( const xAOD::EventInfo*  );*/
    void FillMuonsUser( const xAOD::Muon* );
    void FillElectronsUser( const xAOD::Electron*  );
    void FillJetsUser( const xAOD::Jet*, const std::string = "jet" );
    void FillTausUser( const xAOD::TauJet* );
    /*void FillMETUser( const xAOD::MissingETContainer*  );*/
    void FillLeptons( const xAOD::IParticleContainer*  );
};
#endif
