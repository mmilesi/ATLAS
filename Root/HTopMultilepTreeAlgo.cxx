#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <xAODJet/JetContainer.h>
#include <xAODTracking/VertexContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <AthContainers/ConstDataVector.h>

#include <HTopMultilepAnalysis/HTopMultilepTree.h>
#include <HTopMultilepAnalysis/HTopMultilepTreeAlgo.h>

#include <xAODAnaHelpers/TreeAlgo.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>
#include <xAODAnaHelpers/tools/ReturnCheckConfig.h>

#include "TEnv.h"
#include "TSystem.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepTreeAlgo)

EL::StatusCode HTopMultilepTreeAlgo :: treeInitialize ()
{

  TTree * outTree = new TTree(m_name.c_str(),m_name.c_str());
  if ( !outTree ) {
    Error("treeInitialize()","Failed to instantiate output tree!");
    return EL::StatusCode::FAILURE;
  }

  // get the file we created already
  TFile* treeFile = wk()->getOutputFile ("tree");

  m_HTopTree = new HTopMultilepTree( outTree, treeFile, m_event, m_store, 1e0, m_debug, m_DC14 ); // 1e0 = MeV, default 1e3 = GeV

  // tell the tree to go into the file
  outTree->SetDirectory( treeFile );

  // choose if want to add tree to same directory as ouput histograms
  if ( m_outHistDir ) {
    wk()->addOutput( outTree );
  }

  m_HTopTree->AddEvent(m_evtDetailStr);
  if ( !m_trigDetailStr.empty() )        { m_HTopTree->AddTrigger    (m_trigDetailStr);    }
  if ( !m_muContainerName.empty() )      { m_HTopTree->AddMuons      (m_muDetailStr);	   }
  if ( !m_elContainerName.empty() )      { m_HTopTree->AddElectrons  (m_elDetailStr);	   }
  if ( !m_jetContainerName.empty() )     { m_HTopTree->AddJets       (m_jetDetailStr, "jet"); }
  if ( !m_trigJetContainerName.empty() ) { m_HTopTree->AddJets       (m_trigJetDetailStr, "trigJet"); }
  if ( !m_truthJetContainerName.empty() ){ m_HTopTree->AddJets       (m_truthJetDetailStr, "truthJet"); }
  if ( !m_fatJetContainerName.empty() )  { m_HTopTree->AddFatJets    (m_fatJetDetailStr);  }
  if ( !m_photonContainerName.empty() )  { m_HTopTree->AddPhotons    (m_photonDetailStr);  }
  if ( !m_tauContainerName.empty() )     { m_HTopTree->AddTaus       (m_tauDetailStr);     }
  if ( !m_METContainerName.empty() )     { m_HTopTree->AddMET	     (m_METDetailStr);     }
  if ( !m_lepContainerName.empty() )     { m_HTopTree->AddLeptons    ();		   }

  Info("treeInitialize()", "Successfully initialized output tree");

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepTreeAlgo :: execute ()
{

  // Get EventInfo and the PrimaryVertices
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");
  const xAOD::VertexContainer* vertices(nullptr);
  RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store, m_debug) , "");
  // get the hard-scatter primaryVertex
  const xAOD::Vertex* primaryVertex = HelperFunctions::getPrimaryVertex( vertices );

  m_HTopTree->FillEvent( eventInfo );

  // Fill trigger information
  if ( !m_trigDetailStr.empty() )    {
    m_HTopTree->FillTrigger( eventInfo );
  }

  // for the containers the were supplied, fill the appropiate vectors

  const xAOD::MuonContainer* inMuons(nullptr);
  const xAOD::ElectronContainer* inElectrons(nullptr);
  const xAOD::JetContainer* inJets(nullptr);
  const xAOD::TauJetContainer* inTaus(nullptr);
  const xAOD::MissingETContainer* inMETCont(nullptr);
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);

  if ( !m_muContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inMuons, m_muContainerName, m_event, m_store, m_debug) , "");

    // sort inMuons, and pass the reference to FillMuons()
    const xAOD::MuonContainer inMuonsSorted = HelperFunctions::sort_container_pt( inMuons );
    m_HTopTree->FillMuons( &inMuonsSorted, primaryVertex );

  }
  if ( !m_elContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inElectrons, m_elContainerName, m_event, m_store, m_debug) , "");

    // sort inElectrons, and pass the reference to FillElectrons()
    const xAOD::ElectronContainer inElectronsSorted = HelperFunctions::sort_container_pt( inElectrons );
    m_HTopTree->FillElectrons( &inElectronsSorted, primaryVertex );

  }
  if ( !m_jetContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inJets, m_jetContainerName, m_event, m_store, m_debug) , "");

    // sort inJets, and pass the reference to FillJets()
    const xAOD::JetContainer inJetsSorted = HelperFunctions::sort_container_pt( inJets );
    m_HTopTree->FillJets( &inJetsSorted, HelperFunctions::getPrimaryVertexLocation(vertices) );

    //m_HTopTree->FillJets( inJets, HelperFunctions::getPrimaryVertexLocation(vertices) );
  }
  if ( !m_tauContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inTaus, m_tauContainerName, m_event, m_store, m_debug) , "");

    // sort inTaus, and pass the reference to FillTaus()
    const xAOD::TauJetContainer inTausSorted = HelperFunctions::sort_container_pt( inTaus );
    m_HTopTree->FillTaus( &inTausSorted );

  }
  if ( !m_METContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inMETCont, m_METContainerName, m_event, m_store, m_debug) , "");
    m_HTopTree->FillMET( inMETCont );

  }
  if ( !m_lepContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, m_lepContainerName, m_event, m_store, m_debug) ,"");

    // sort inLeptons, and pass the reference to FillLeptons()
    const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );
    m_HTopTree->FillLeptons( &leptonsSorted );

  }

  // fill the tree
  m_HTopTree->Fill();

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepTreeAlgo :: finalize () {

  Info("finalize()", "Deleting tree instances...");

  if ( m_HTopTree ) { delete m_HTopTree; m_HTopTree = nullptr; }

  return EL::StatusCode::SUCCESS;
}

