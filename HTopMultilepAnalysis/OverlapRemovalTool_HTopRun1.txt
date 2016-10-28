#ifndef ASSOCIATIONUTILS_OVERLAPREMOVALTOOL_HTOPRUN1_H
#define ASSOCIATIONUTILS_OVERLAPREMOVALTOOL_HTOPRUN1_H

// Framework includes
#include "AsgTools/AsgTool.h"

// EDM includes
#include "xAODBase/IParticle.h"

// Local includes
#include "AssociationUtils/OverlapRemovalTool.h"

class OverlapRemovalTool_HTopRun1 : public OverlapRemovalTool
{

  public:

    /* Constructor for standalone usage */
    OverlapRemovalTool_HTopRun1(const std::string& name);
    
    /* Destructor */
    ~OverlapRemovalTool_HTopRun1(){};
    
    /* Initialize the tool */
    StatusCode initialize();    
    
    StatusCode removeOverlaps(const xAOD::ElectronContainer* electrons,
                              const xAOD::MuonContainer* muons,
                              const xAOD::JetContainer* jets,
                              const xAOD::TauJetContainer* taus = 0,
                              const xAOD::PhotonContainer* photons = 0);
    
    StatusCode removeOverlaps(const xAOD::ElectronContainer* electrons,
                              const xAOD::MuonContainer* muons,
                              const xAOD::JetContainer* jets,
                              const xAOD::TauJetContainer* taus,
                              const xAOD::ElectronContainer* looseElectrons,
                              const xAOD::MuonContainer* looseMuons,
                              const xAOD::PhotonContainer* photons = 0);

    StatusCode removeEleMuonOverlap(const xAOD::ElectronContainer& electrons,
                                    const xAOD::MuonContainer& muons);

    StatusCode removeEleEleOverlap(const xAOD::ElectronContainer& electrons);

    StatusCode removeJetEleOverlap(const xAOD::JetContainer& jets,
    				   const xAOD::ElectronContainer& electrons);

    StatusCode removeMuonJetOverlap(const xAOD::MuonContainer& muons,
                                    const xAOD::JetContainer& jets);
    
    StatusCode removeTauEleOverlap(const xAOD::TauJetContainer& taus,
                                   const xAOD::ElectronContainer& electrons);
				   
    StatusCode removeTauMuonOverlap(const xAOD::TauJetContainer& taus,
                                    const xAOD::MuonContainer& muons);
   
    StatusCode removeJetTauOverlap(const xAOD::JetContainer& jets,
                                   const xAOD::TauJetContainer& taus);
	
    /* 
    re-implementing this since ASG AssociationUtils::OverlapRemoalTool is 
    using rapidity, not pseudorapidity 
    */				   
    inline double deltaRHTop(const xAOD::IParticle* p1, const xAOD::IParticle* p2)
    { 
       double dEta = p1->eta() - p2->eta();
       double dPhi = TVector2::Phi_mpi_pi(p1->phi() - p2->phi());
       return sqrt( dEta*dEta + dPhi*dPhi );
    };
				   
				   
  private:

    // electron-muon overlap cone (removes electron) 
    float m_electronMuonDR_Run1;
    
    // electron-electron overlap cone (removes electron)
    float m_electronElectronDR_Run1;
    
    /// jet-electron overlap cone (removes jet)
    float m_jetElectronDR_Run1;
    
    /// muon-jet overlap cone (removes muon)
    float m_muonJetDR_Run1; 
    
    /// jet-tau overlap cone (removes jet)
    float m_jetTauDR_Run1;
    
    /// tau-electron overlap cone (removes tau)
    float m_tauElectronDR_Run1; 
    
    /// tau-muon overlap cone (removes tau)
    float m_tauMuonDR_Run1;


}; 

#endif
