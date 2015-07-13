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
    unsigned int  	m_nBjetsMedium;  
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
    int                 m_isNonTightEvent;
    int                 m_isProbeElEvent;
    int                 m_isProbeMuEvent;
     
    /* jet variables */
    std::vector<float> m_jet_m;    
    
    /* muon variables */
    std::vector<int> m_muon_isTight;
    std::vector<int> m_muon_isTruthMatched; 
    std::vector<int> m_muon_isChFlip; 
    std::vector<int> m_muon_isBrem;
    std::vector<int> m_muon_truthType;  
    std::vector<int> m_muon_truthPdgId;  
    std::vector<int> m_muon_truthOrigin; 
    std::vector<int> m_muon_truthStatus;  
    std::vector<int> m_muon_isOS;
    std::vector<int> m_muon_isClosestSS;
    std::vector<int> m_muon_isTag; 
    /* muon TAG variables */
    std::vector<float> m_muon_tag_pt;
    std::vector<float> m_muon_tag_eta;
    std::vector<int>   m_muon_tag_isIsolated; 
    std::vector<int>   m_muon_tag_isTrigMatched; 
    std::vector<int>   m_muon_tag_isTight; 
    std::vector<int>   m_muon_tag_isTruthMatched; 
    std::vector<int>   m_muon_tag_isChFlip; 
    std::vector<int>   m_muon_tag_isBrem;
    std::vector<int>   m_muon_tag_truthType; 
    std::vector<int>   m_muon_tag_truthPdgId; 
    std::vector<int>   m_muon_tag_truthOrigin; 
    std::vector<int>   m_muon_tag_truthStatus;  
    /* muon PROBE variables */
    std::vector<float> m_muon_probe_pt;
    std::vector<float> m_muon_probe_eta;
    std::vector<int>   m_muon_probe_isIsolated; 
    std::vector<int>   m_muon_probe_isTrigMatched; 
    std::vector<int>   m_muon_probe_isTight; 
    std::vector<int>   m_muon_probe_isTruthMatched; 
    std::vector<int>   m_muon_probe_isChFlip; 
    std::vector<int>   m_muon_probe_isBrem;
    std::vector<int>   m_muon_probe_truthType; 
    std::vector<int>   m_muon_probe_truthPdgId; 
    std::vector<int>   m_muon_probe_truthOrigin; 
    std::vector<int>   m_muon_probe_truthStatus;  
    
    /* electron variables */
    std::vector<float> m_electron_calo_eta;   
    std::vector<int>   m_electron_crack;   
    std::vector<int>   m_electron_isTight;
    std::vector<int>   m_electron_isTruthMatched; 
    std::vector<int>   m_electron_isChFlip; 
    std::vector<int>   m_electron_isBrem; 
    std::vector<int>   m_electron_truthType;
    std::vector<int>   m_electron_truthPdgId;
    std::vector<int>   m_electron_truthOrigin; 
    std::vector<int>   m_electron_truthStatus;
    std::vector<int>   m_electron_isOS;
    std::vector<int>   m_electron_isClosestSS;
    std::vector<int>   m_electron_isTag; 
    /* electron TAG variables */
    std::vector<float> m_electron_tag_pt;
    std::vector<float> m_electron_tag_eta;
    std::vector<int>   m_electron_tag_isIsolated; 
    std::vector<int>   m_electron_tag_isTrigMatched;     
    std::vector<int>   m_electron_tag_isTight; 
    std::vector<int>   m_electron_tag_isTruthMatched; 
    std::vector<int>   m_electron_tag_isChFlip; 
    std::vector<int>   m_electron_tag_isBrem;
    std::vector<int>   m_electron_tag_truthType;
    std::vector<int>   m_electron_tag_truthPdgId; 
    std::vector<int>   m_electron_tag_truthOrigin; 
    std::vector<int>   m_electron_tag_truthStatus;  
    /* electron PROBE variables */
    std::vector<float> m_electron_probe_pt;
    std::vector<float> m_electron_probe_eta;
    std::vector<int>   m_electron_probe_isIsolated; 
    std::vector<int>   m_electron_probe_isTrigMatched;     
    std::vector<int>   m_electron_probe_isTight; 
    std::vector<int>   m_electron_probe_isTruthMatched; 
    std::vector<int>   m_electron_probe_isChFlip; 
    std::vector<int>   m_electron_probe_isBrem;
    std::vector<int>   m_electron_probe_truthType; 
    std::vector<int>   m_electron_probe_truthPdgId; 
    std::vector<int>   m_electron_probe_truthOrigin; 
    std::vector<int>   m_electron_probe_truthStatus;         
       
    /* lepton variables */
    int                m_nlep;
    std::vector<float> m_lepton_pt;
    std::vector<float> m_lepton_phi;
    std::vector<float> m_lepton_eta;
    std::vector<float> m_lepton_m;
    std::vector<float> m_lepton_charge;
    std::vector<int>   m_lepton_flavour;  
    std::vector<int>   m_lepton_isTrigMatched; 
    std::vector<int>   m_lepton_isIsolated;
    std::vector<int>   m_lepton_isTight; 
    std::vector<int>   m_lepton_isTruthMatched; 
    std::vector<int>   m_lepton_isChFlip; 
    std::vector<int>   m_lepton_isBrem; 
    std::vector<int>   m_lepton_truthType; 
    std::vector<int>   m_lepton_truthPdgId;     
    std::vector<int>   m_lepton_truthOrigin; 
    std::vector<int>   m_lepton_truthStatus; 
    std::vector<int>   m_lepton_isOS;
    std::vector<int>   m_lepton_isClosestSS;
    std::vector<int>   m_lepton_isTag; 
    /* lepton TAG variables */
    std::vector<float> m_lepton_tag_pt;
    std::vector<float> m_lepton_tag_eta;
    std::vector<int>   m_lepton_tag_flavour;
    std::vector<int>   m_lepton_tag_charge;
    std::vector<int>   m_lepton_tag_isIsolated; 
    std::vector<int>   m_lepton_tag_isTrigMatched;     
    std::vector<int>   m_lepton_tag_isTight; 
    std::vector<int>   m_lepton_tag_isTruthMatched; 
    std::vector<int>   m_lepton_tag_isChFlip; 
    std::vector<int>   m_lepton_tag_isBrem;
    std::vector<int>   m_lepton_tag_truthType;
    std::vector<int>   m_lepton_tag_truthPdgId; 
    std::vector<int>   m_lepton_tag_truthOrigin; 
    std::vector<int>   m_lepton_tag_truthStatus;  
    /* lepton PROBE variables */
    std::vector<float> m_lepton_probe_pt;
    std::vector<float> m_lepton_probe_eta;
    std::vector<int>   m_lepton_probe_flavour;
    std::vector<int>   m_lepton_probe_charge; 
    std::vector<int>   m_lepton_probe_isIsolated; 
    std::vector<int>   m_lepton_probe_isTrigMatched;       
    std::vector<int>   m_lepton_probe_isTight; 
    std::vector<int>   m_lepton_probe_isTruthMatched; 
    std::vector<int>   m_lepton_probe_isChFlip; 
    std::vector<int>   m_lepton_probe_isBrem;
    std::vector<int>   m_lepton_probe_truthType;  
    std::vector<int>   m_lepton_probe_truthPdgId;
    std::vector<int>   m_lepton_probe_truthOrigin; 
    std::vector<int>   m_lepton_probe_truthStatus;         

    /* tau variables */
    std::vector<int>   m_tau_isBDTTight;    
    
  public:
    
    HTopMultilepTree( TTree* tree, TFile* file, xAOD::TEvent* event, xAOD::TStore* store, const float units = 1e3, bool debug = false, bool DC14 = false );
    ~HTopMultilepTree();

    void AddEventUser(const std::string detailStrUser = "");
    /*void AddTriggerUser(const std::string detailStrUser = "");*/
    void AddMuonsUser(const std::string detailStrUser = "");
    void AddElectronsUser(const std::string detailStrUser = "");
    void AddJetsUser(const std::string detailStrUser = "");
    void AddTausUser(const std::string detailStrUser = "");
    /*void AddMETUser(const std::string detailStrUser = "");*/    
    
    void AddLeptons();
    
    void ClearEventUser();
    /*void ClearTriggerUser();*/
    void ClearMuonsUser();   
    void ClearElectronsUser();  
    void ClearJetsUser();
    void ClearTausUser();
    /*void ClearMETUser();*/    
    void ClearLeptons();
    
    void FillEventUser( const xAOD::EventInfo* );
    /*void FillTriggerUser( const xAOD::EventInfo*  );*/
    void FillMuonsUser( const xAOD::Muon* );
    void FillElectronsUser( const xAOD::Electron*  );
    void FillJetsUser( const xAOD::Jet* );
    void FillTausUser( const xAOD::TauJet* );
    /*void FillMETUser( const xAOD::MissingETContainer*  );*/
    void FillLeptons( const xAOD::IParticleContainer*  ); 
};
#endif
