// EDM includes
#include "AthContainers/AuxElement.h"

// Local includes
#include "HTopMultilepAnalysis/OverlapRemovalTool_HTopRun1.h"

const float invGeV = 0.001;

//-----------------------------------------------------------------------------
// Standard constructor - will call the base class constructor!
//-----------------------------------------------------------------------------
OverlapRemovalTool_HTopRun1::OverlapRemovalTool_HTopRun1(const std::string& name) 
        : OverlapRemovalTool(name)
{
  // dR cones for defining overlap
  //
  declareProperty("ElectronMuonDRCone_Run1",     m_electronMuonDR_Run1	   = 0.1);
  declareProperty("ElectronElectronDRCone_Run1", m_electronElectronDR_Run1 = 0.3);
  declareProperty("JetElectronDRCone_Run1",      m_jetElectronDR_Run1      = 0.3);
  declareProperty("MuonJetDRCone_Run1",          m_muonJetDR_Run1          = 0.04);
  declareProperty("JetTauDRCone_Run1",           m_jetTauDR_Run1           = 0.3);
  declareProperty("TauElectronDRCone_Run1",      m_tauElectronDR_Run1      = 0.2);
  declareProperty("TauMuonDRCone_Run1",          m_tauMuonDR_Run1          = 0.2);
}

//-----------------------------------------------------------------------------
// Remove all overlapping objects according to the old Run1 HTop prescription
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::
removeOverlaps(const xAOD::ElectronContainer* electrons,
               const xAOD::MuonContainer* muons,
               const xAOD::JetContainer* jets,
               const xAOD::TauJetContainer* taus,
               const xAOD::PhotonContainer* photons)
{
  return removeOverlaps(electrons, muons, jets, taus,
                        electrons, muons, photons);
}
//-----------------------------------------------------------------------------
StatusCode OverlapRemovalTool_HTopRun1::
removeOverlaps(const xAOD::ElectronContainer* electrons,
               const xAOD::MuonContainer* muons,
               const xAOD::JetContainer* jets,
               const xAOD::TauJetContainer* taus,
               const xAOD::ElectronContainer* looseElectrons,
               const xAOD::MuonContainer* looseMuons,
               const xAOD::PhotonContainer* photons)
{
  
  // Check pointer validity. I can add more flexibility later,
  // but for now, if users don't want to use one of these object types,
  // they can use the lower-level methods instead.
  if ( !electrons || !muons || !jets ) {
    ATH_MSG_ERROR("Encountered NULL pointer in required object");
    return StatusCode::FAILURE;
  }
  if ( taus && !( looseElectrons && looseMuons )){
    ATH_MSG_ERROR("Taus provided but loose leptons are NULL!");
    return StatusCode::FAILURE;
  }

  // Reset all decorations to passing
  resetDecorations(*electrons);
  resetDecorations(*muons);
  resetDecorations(*jets);
  if ( photons ) resetDecorations(*photons);
  if ( taus )    resetDecorations(*taus);

  // remove overlaps in the following order
  //
  ATH_CHECK( removeEleMuonOverlap(*electrons, *muons) );
  ATH_CHECK( removeEleEleOverlap(*electrons) );
  ATH_CHECK( removeJetEleOverlap(*jets, *electrons) );
  ATH_CHECK( removeMuonJetOverlap(*muons, *jets) );
  if ( taus ) {
    ATH_CHECK( removeTauEleOverlap(*taus, *electrons) );
    ATH_CHECK( removeTauMuonOverlap(*taus, *muons) );
    ATH_CHECK( removeJetTauOverlap(*jets, *taus) );
  }

  return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------
// Remove electrons overlapping w/ muons
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeEleMuonOverlap(const xAOD::ElectronContainer& electrons,
                                                             const xAOD::MuonContainer& muons) 
{
  
  for ( const auto electron : electrons ) {
    if ( isSurvivingObject(electron) ) {
      if ( objectOverlaps(electron, muons, m_electronMuonDR_Run1 ) ) {
        ATH_MSG_DEBUG("  Found overlap electron w/ muon: " << electron->pt()*invGeV);
        setObjectFail(electron);
      }
    }
  }
  
  return StatusCode::SUCCESS;
}							     

//-----------------------------------------------------------------------------
// Remove electrons overlapping w/ electrons
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeEleEleOverlap(const xAOD::ElectronContainer& electrons)
{
 
  for ( const auto electron : electrons ) {
    if ( isSurvivingObject(electron) ) {
      if ( objectOverlaps(electron, electrons, m_electronMuonDR_Run1 ) ) {
        ATH_MSG_DEBUG("  Found overlap electron w/ electron: " << electron->pt()*invGeV);
        setObjectFail(electron);
      }
    }
  }

  return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------
// Remove jets overlapping w/ electrons
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeJetEleOverlap(const xAOD::JetContainer& jets,
    				                            const xAOD::ElectronContainer& electrons)
{
  
  for ( const auto jet : jets ) {
    if ( isSurvivingObject(jet) /*&& !isBJet(jet)*/ ) { // maybe we should remove the jet only if the jet is NOT btagged...
      if ( objectOverlaps(jet, electrons, m_jetElectronDR_Run1) ) {
        ATH_MSG_DEBUG("  Found overlap jet w/ electron: " << jet->pt()*invGeV);
        setObjectFail(jet);
      }
    }
  }
  
  return StatusCode::SUCCESS;
}
				   
//-----------------------------------------------------------------------------
// Remove muons overlapping w/ jets
//-----------------------------------------------------------------------------
				   

StatusCode OverlapRemovalTool_HTopRun1::removeMuonJetOverlap(const xAOD::MuonContainer& muons,
                                                             const xAOD::JetContainer& jets)
{
  
  for ( const auto muon : muons ) {
    if ( isSurvivingObject(muon) ) {
      for ( const auto jet : jets ) {  
        if ( isSurvivingObject(jet) /*&& isBJet(jet)*/ ) { // maybe we should remove the muon only if the jet is btagged...
          if ( deltaR(jet, muon) < ( m_muonJetDR_Run1 + 10.0 * ( invGeV / muon->pt() ) ) ) {
            ATH_MSG_DEBUG("  Found overlap muon w/ jet: " << muon->pt()*invGeV);
            setObjectFail(muon);
          }
        }
      }
    }
  }
  
  return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------
// Remove taus overlapping w/ electrons
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeTauEleOverlap(const xAOD::TauJetContainer& taus,
                                                            const xAOD::ElectronContainer& electrons)
{
  
  for ( const auto tau : taus ) {
    if ( isSurvivingObject(tau) ) { 
      if ( objectOverlaps(tau, electrons, m_tauElectronDR_Run1) ) {
        ATH_MSG_DEBUG("  Found overlap tau w/ electron: " << tau->pt()*invGeV);
        setObjectFail(tau);
      }
    }
  }
  
  return StatusCode::SUCCESS;
}							    

//-----------------------------------------------------------------------------
// Remove taus overlapping w/ muons
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeTauMuonOverlap(const xAOD::TauJetContainer& taus,
                                                             const xAOD::MuonContainer& muons)
{
  
  for ( const auto tau : taus ) {
    if ( isSurvivingObject(tau) ) { 
      if ( objectOverlaps(tau, muons, m_tauMuonDR_Run1) ) {
        ATH_MSG_DEBUG("  Found overlap tau w/ muon: " << tau->pt()*invGeV);
        setObjectFail(tau);
      }
    }
  }
    
  return StatusCode::SUCCESS;
}
//-----------------------------------------------------------------------------
// Remove jets overlapping w/ taus
//-----------------------------------------------------------------------------

StatusCode OverlapRemovalTool_HTopRun1::removeJetTauOverlap(const xAOD::JetContainer& jets,
                                                            const xAOD::TauJetContainer& taus)
{
  
  for ( const auto jet : jets ) {
    if ( isSurvivingObject(jet) ) { 
      if ( objectOverlaps(jet, taus, m_jetTauDR_Run1) ) {
        ATH_MSG_DEBUG("  Found overlap jet w/ tau: " << jet->pt()*invGeV);
        setObjectFail(jet);
      }
    }
  }

  return StatusCode::SUCCESS;
}

