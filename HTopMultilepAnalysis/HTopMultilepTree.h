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

    float               m_mOSPair01;
    float               m_mOSPair02;
    int                 m_isOSPairSF01;
    int                 m_isOSPairSF02;

    float               m_mT_lep0MET;
    float               m_mT_lep1MET;

    int  		m_is_T_T;
    int  		m_is_T_AntiT;
    int  		m_is_AntiT_T;
    int  		m_is_AntiT_AntiT;

    int  		m_is_T_MAntiT;
    int  		m_is_MAntiT_T;
    int  		m_is_MAntiT_MAntiT;

    int  		m_is_M_M;
    int                 m_is_T_AntiM;
    int                 m_is_AntiM_T;
    int                 m_is_M_AntiM;
    int                 m_is_AntiM_M;
    int                 m_is_AntiM_AntiM;

    int  		m_is_Tel_AntiTmu;
    int  		m_is_AntiTel_Tmu;
    int  		m_is_Tmu_AntiTel;
    int  		m_is_AntiTmu_Tel;

    int  		m_is_Tel_MAntiTmu;
    int  		m_is_MAntiTel_Tmu;
    int  		m_is_Tmu_MAntiTel;
    int  		m_is_MAntiTmu_Tel;

    int  		m_is_Mel_AntiMmu;
    int  		m_is_AntiMel_Mmu;
    int  		m_is_Mmu_AntiMel;
    int  		m_is_AntiMmu_Mel;

    int  		m_is_Tel_AntiMmu;
    int  		m_is_AntiMel_Tmu;
    int  		m_is_Tmu_AntiMel;
    int  		m_is_AntiMmu_Tel;

    int                 m_isNonTightEvent;
    int                 m_isProbeElEvent;
    int                 m_isProbeMuEvent;

    // HTop scale factors per-event
    std::vector<float> m_weight_lepton_trig_HTop;
    std::vector<float> m_weight_lepton_reco_HTop;
    std::vector<float> m_weight_lepton_iso_HTop;
    std::vector<float> m_weight_lepton_ID_HTop;
    std::vector<float> m_weight_lepton_TTVA_HTop;

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
    std::vector<float> m_muon_tag_trkz0sintheta;
    std::vector<float> m_muon_tag_ptvarcone30;
    std::vector<float> m_muon_tag_topoetcone20;
    std::vector<int>   m_muon_tag_isIsolated_Loose;
    std::vector<int>   m_muon_tag_isIsolated_FixedCutTightTrackOnly;
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
    std::vector<float> m_muon_probe_trkz0sintheta;
    std::vector<float> m_muon_probe_ptvarcone30;
    std::vector<float> m_muon_probe_topoetcone20;
    std::vector<int>   m_muon_probe_isIsolated_Loose;
    std::vector<int>   m_muon_probe_isIsolated_FixedCutTightTrackOnly;
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
    std::vector<float> m_electron_tag_caloCluster_eta;
    std::vector<float> m_electron_tag_trkd0sig;
    std::vector<float> m_electron_tag_trkz0sintheta;
    std::vector<int>   m_electron_tag_LHLoose;
    std::vector<int>   m_electron_tag_LHMedium;
    std::vector<int>   m_electron_tag_LHTight;
    std::vector<int>   m_electron_tag_IsEMLoose;
    std::vector<int>   m_electron_tag_IsEMMedium;
    std::vector<int>   m_electron_tag_IsEMTight;
    std::vector<float> m_electron_tag_ptvarcone20;
    std::vector<float> m_electron_tag_topoetcone20;
    std::vector<int>   m_electron_tag_isIsolated_Loose;
    std::vector<int>   m_electron_tag_isIsolated_FixedCutTight;
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
    std::vector<float> m_electron_probe_caloCluster_eta;
    std::vector<float> m_electron_probe_trkd0sig;
    std::vector<float> m_electron_probe_trkz0sintheta;
    std::vector<int>   m_electron_probe_LHLoose;
    std::vector<int>   m_electron_probe_LHMedium;
    std::vector<int>   m_electron_probe_LHTight;
    std::vector<int>   m_electron_probe_IsEMLoose;
    std::vector<int>   m_electron_probe_IsEMMedium;
    std::vector<int>   m_electron_probe_IsEMTight;
    std::vector<float> m_electron_probe_ptvarcone20;
    std::vector<float> m_electron_probe_topoetcone20;
    std::vector<int>   m_electron_probe_isIsolated_Loose;
    std::vector<int>   m_electron_probe_isIsolated_FixedCutTight;
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
    std::vector<float> m_lepton_trkd0sig;
    std::vector<float> m_lepton_trkz0sintheta;
    std::vector<float> m_lepton_charge;
    std::vector<int>   m_lepton_flavour;
    std::vector<int>   m_lepton_isTrigMatched;
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

    /* Flat varables for leptons in the 3lep category */
    std::vector<float> m_lepton_3lepOS_pt;
    std::vector<float> m_lepton_3lepOS_phi;
    std::vector<float> m_lepton_3lepOS_eta;
    std::vector<float> m_lepton_3lepOS_m;
    std::vector<float> m_lepton_3lepOS_trkd0sig;
    std::vector<float> m_lepton_3lepOS_trkz0sintheta;
    std::vector<float> m_lepton_3lepOS_charge;
    std::vector<int>   m_lepton_3lepOS_flavour;
    std::vector<int>   m_lepton_3lepOS_isTrigMatched;
    std::vector<int>   m_lepton_3lepOS_isTight;
    std::vector<int>   m_lepton_3lepOS_isMedium;
    std::vector<float> m_lepton_3lepClosestSS_pt;
    std::vector<float> m_lepton_3lepClosestSS_phi;
    std::vector<float> m_lepton_3lepClosestSS_eta;
    std::vector<float> m_lepton_3lepClosestSS_m;
    std::vector<float> m_lepton_3lepClosestSS_trkd0sig;
    std::vector<float> m_lepton_3lepClosestSS_trkz0sintheta;
    std::vector<float> m_lepton_3lepClosestSS_charge;
    std::vector<int>   m_lepton_3lepClosestSS_flavour;
    std::vector<int>   m_lepton_3lepClosestSS_isTrigMatched;
    std::vector<int>   m_lepton_3lepClosestSS_isTight;
    std::vector<int>   m_lepton_3lepClosestSS_isMedium;
    std::vector<float> m_lepton_3lepOtherSS_pt;
    std::vector<float> m_lepton_3lepOtherSS_phi;
    std::vector<float> m_lepton_3lepOtherSS_eta;
    std::vector<float> m_lepton_3lepOtherSS_m;
    std::vector<float> m_lepton_3lepOtherSS_trkd0sig;
    std::vector<float> m_lepton_3lepOtherSS_trkz0sintheta;
    std::vector<float> m_lepton_3lepOtherSS_charge;
    std::vector<int>   m_lepton_3lepOtherSS_flavour;
    std::vector<int>   m_lepton_3lepOtherSS_isTrigMatched;
    std::vector<int>   m_lepton_3lepOtherSS_isTight;
    std::vector<int>   m_lepton_3lepOtherSS_isMedium;

    /* lepton TAG variables */
    std::vector<float> m_lepton_tag_pt;
    std::vector<float> m_lepton_tag_eta;
    std::vector<float> m_lepton_tag_trkd0sig;
    std::vector<float> m_lepton_tag_trkz0sintheta;
    std::vector<int>   m_lepton_tag_flavour;
    std::vector<int>   m_lepton_tag_charge;
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
    std::vector<float> m_lepton_probe_trkd0sig;
    std::vector<float> m_lepton_probe_trkz0sintheta;
    std::vector<int>   m_lepton_probe_flavour;
    std::vector<int>   m_lepton_probe_charge;
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
