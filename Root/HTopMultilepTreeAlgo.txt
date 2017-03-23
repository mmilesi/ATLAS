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

EL::StatusCode HTopMultilepTreeAlgo :: execute ()
{

  // what systematics do we need to process for this event?
  // handle the nominal case (merge all) on every event, always
  std::vector<std::string> event_systNames({""});
  std::vector<std::string> muSystNames;
  std::vector<std::string> elSystNames;
  std::vector<std::string> jetSystNames;
  std::vector<std::string> photonSystNames;

  // this is a temporary pointer that gets switched around to check each of the systematics
  std::vector<std::string>* systNames(nullptr);

  // note that the way we set this up, none of the below ##SystNames vectors contain the nominal case
  // TODO: do we really need to check for duplicates? Maybe, maybe not.
  if(!m_muSystsVec.empty()){
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(systNames, m_muSystsVec, 0, m_store, m_verbose) ,"");
    for(const auto& systName: *systNames){
      if (std::find(event_systNames.begin(), event_systNames.end(), systName) != event_systNames.end()) continue;
      event_systNames.push_back(systName);
      muSystNames.push_back(systName);
    }
  }

  if(!m_elSystsVec.empty()){
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(systNames, m_elSystsVec, 0, m_store, m_verbose) ,"");
    for(const auto& systName: *systNames){
      if (std::find(event_systNames.begin(), event_systNames.end(), systName) != event_systNames.end()) continue;
      event_systNames.push_back(systName);
      elSystNames.push_back(systName);
    }
  }

  if(!m_jetSystsVec.empty()){
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(systNames, m_jetSystsVec, 0, m_store, m_verbose) ,"");
    for(const auto& systName: *systNames){
      if (std::find(event_systNames.begin(), event_systNames.end(), systName) != event_systNames.end()) continue;
      event_systNames.push_back(systName);
      jetSystNames.push_back(systName);
    }
  }

  if(!m_photonSystsVec.empty()){
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(systNames, m_photonSystsVec, 0, m_store, m_verbose) ,"");
    for(const auto& systName: *systNames){
      if (std::find(event_systNames.begin(), event_systNames.end(), systName) != event_systNames.end()) continue;
      event_systNames.push_back(systName);
      photonSystNames.push_back(systName);
    }
  }

  TFile* treeFile = wk()->getOutputFile ("tree");

  // let's make the tdirectory and ttrees
  for(const auto& systName: event_systNames){
    // check if we have already created the tree
    if(m_trees.find(systName) != m_trees.end()) continue;
    std::string treeName = systName;
    if(systName.empty()) treeName = "nominal";

    Info("execute()", "Making tree %s/%s", m_name.c_str(), treeName.c_str());
    TTree * outTree = new TTree(treeName.c_str(),treeName.c_str());
    if ( !outTree ) {
      Error("execute()","Failed to instantiate output tree!");
      return EL::StatusCode::FAILURE;
    }

    //m_units = 1e0; // use MeV by default!

    m_trees[systName] = new HTopMultilepTree( outTree, treeFile, m_event, m_store, m_units, m_debug, m_DC14 );
    HTopMultilepTree* helpTree = dynamic_cast<HTopMultilepTree*>(m_trees[systName]);

    // tell the tree to go into the file
    outTree->SetDirectory( treeFile->GetDirectory(m_name.c_str()) );
    // choose if want to add tree to same directory as ouput histograms
    if ( m_outHistDir ) {
      if(m_trees.size() > 1) Warning("execute()", "You're running systematics! You may find issues in writing all of the output TTrees to the output histogram file... Set `m_outHistDir = false` if you run into issues!");
      wk()->addOutput( outTree );
    }

    // initialize all branch addresses since we just added this tree
    helpTree->AddEvent( m_evtDetailStr );
    if ( !m_trigDetailStr.empty() )       {   helpTree->AddTrigger    (m_trigDetailStr);    }
    if ( !m_muContainerName.empty() )     {   helpTree->AddMuons      (m_muDetailStr);      }
    if ( !m_elContainerName.empty() )     {   helpTree->AddElectrons  (m_elDetailStr);      }
    if ( !m_jetContainerName.empty() )    {   helpTree->AddJets       (m_jetDetailStr, "jet");     }
    if ( !m_trigJetContainerName.empty() ){   helpTree->AddJets       (m_trigJetDetailStr, "trigJet");     }
    if ( !m_truthJetContainerName.empty() ){  helpTree->AddJets       (m_truthJetDetailStr, "truthJet");     }
    if ( !m_fatJetContainerName.empty() ) {   helpTree->AddFatJets    (m_fatJetDetailStr);  }
    if ( !m_tauContainerName.empty() )    {   helpTree->AddTaus       (m_tauDetailStr);     }
    if ( !m_METContainerName.empty() )    {   helpTree->AddMET        (m_METDetailStr);     }
    if ( !m_photonContainerName.empty() ) {   helpTree->AddPhotons    (m_photonDetailStr);  }
    if ( !m_lepContainerName.empty() )    {   helpTree->AddLeptons    ();                   }
  }

  /* THIS IS WHERE WE START PROCESSING THE EVENT AND PLOTTING THINGS */

  // Get EventInfo and the PrimaryVertices
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");
  const xAOD::VertexContainer* vertices(nullptr);
  RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store, m_verbose) ,"");
  // get the primaryVertex
  const xAOD::Vertex* primaryVertex = HelperFunctions::getPrimaryVertex( vertices );

  for(const auto& systName: event_systNames){

    HTopMultilepTree* helpTree = dynamic_cast<HTopMultilepTree*>(m_trees[systName]);

    // assume the nominal container by default
    std::string muSuffix("");
    std::string elSuffix("");
    std::string jetSuffix("");
    std::string photonSuffix("");
    /*
       if we find the systematic in the corresponding vector, we will use that container's systematic version instead of nominal version
        NB: since none of these contain the "" (nominal) case because of how I filter it, we handle the merging.. why?
        - in each loop to make the ##systNames vectors, we check to see if the systName exists in event_systNames which is initialized
        -   to {""} - the nominal case. If the systName exists, we do not add it to the corresponding ##systNames vector, otherwise, we do.
        -   This precludes the nominal case in all of the ##systNames vectors, which means the default will always be to run nominal.
    */
    if (std::find(muSystNames.begin(), muSystNames.end(), systName) != muSystNames.end()) muSuffix = systName;
    if (std::find(elSystNames.begin(), elSystNames.end(), systName) != elSystNames.end()) elSuffix = systName;
    if (std::find(jetSystNames.begin(), jetSystNames.end(), systName) != jetSystNames.end()) jetSuffix = systName;
    if (std::find(photonSystNames.begin(), photonSystNames.end(), systName) != photonSystNames.end()) photonSuffix = systName;

    helpTree->FillEvent( eventInfo, m_event );

    // Fill trigger information
    if ( !m_trigDetailStr.empty() )    {
      helpTree->FillTrigger( eventInfo );
    }

    // Fill jet trigger information - this can be used if with layer/cleaning info we need to turn off some variables?
    /*if ( !m_trigJetDetailStr.empty() ) {
      helpTree->FillJetTrigger();
    }*/

    // for the containers the were supplied, fill the appropriate vectors
    if ( !m_muContainerName.empty() ) {
      const xAOD::MuonContainer* inMuon(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inMuon, m_muContainerName+muSuffix, m_event, m_store, m_verbose) ,"");
      // sort, and pass the reference to FillMuons()
      const xAOD::MuonContainer inMuonsSorted = HelperFunctions::sort_container_pt( inMuon );
      helpTree->FillMuons( &inMuonsSorted, primaryVertex );
    }

    if ( !m_elContainerName.empty() ) {
      const xAOD::ElectronContainer* inElec(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inElec, m_elContainerName+elSuffix, m_event, m_store, m_verbose) ,"");
      // sort, and pass the reference to FillElectrons()
      const xAOD::ElectronContainer inElectronsSorted = HelperFunctions::sort_container_pt( inElec );
      helpTree->FillElectrons( &inElectronsSorted, primaryVertex );

    }
    if ( !m_jetContainerName.empty() ) {
      const xAOD::JetContainer* inJets(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inJets, m_jetContainerName+jetSuffix, m_event, m_store, m_verbose) ,"");
      // sort, and pass the reference to FillJets()
      const xAOD::JetContainer inJetsSorted = HelperFunctions::sort_container_pt( inJets );
      helpTree->FillJets( &inJetsSorted );
    }
    if ( !m_trigJetContainerName.empty() ) {
      const xAOD::JetContainer* inTrigJets(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inTrigJets, m_trigJetContainerName, m_event, m_store, m_verbose) ,"");
      helpTree->FillJets( inTrigJets, HelperFunctions::getPrimaryVertexLocation(vertices), "trigJet" );
    }
    if ( !m_truthJetContainerName.empty() ) {
      const xAOD::JetContainer* inTruthJets(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inTruthJets, m_truthJetContainerName, m_event, m_store, m_verbose) ,"");
          helpTree->FillJets( inTruthJets, HelperFunctions::getPrimaryVertexLocation(vertices), "truthJet" );
    }
    if ( !m_fatJetContainerName.empty() ) {
      const xAOD::JetContainer* inFatJets(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inFatJets, m_fatJetContainerName, m_event, m_store, m_verbose) ,"");
      helpTree->FillFatJets( inFatJets );
    }
    if ( !m_tauContainerName.empty() ) {
      const xAOD::TauJetContainer* inTaus(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inTaus, m_tauContainerName, m_event, m_store, m_verbose) , "");
      // sort, and pass the reference to FillTaus()
      const xAOD::TauJetContainer inTausSorted = HelperFunctions::sort_container_pt( inTaus );
      helpTree->FillTaus( &inTausSorted );
    }
    if ( !m_METContainerName.empty() ) {
      const xAOD::MissingETContainer* inMETCont(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inMETCont, m_METContainerName, m_event, m_store, m_debug) , "");
      helpTree->FillMET( inMETCont );
    }
    if ( !m_photonContainerName.empty() ) {
      const xAOD::PhotonContainer* inPhotons(nullptr);
      RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inPhotons, m_photonContainerName+photonSuffix, m_event, m_store, m_verbose) ,"");
      helpTree->FillPhotons( inPhotons );
    }
    if ( !m_lepContainerName.empty() ) {
      ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);
      RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, m_lepContainerName, m_event, m_store, m_debug) ,"");
      // sort, and pass the reference to FillLeptons()
      const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );
      helpTree->FillLeptons( &leptonsSorted );
    }

    // fill the tree
    helpTree->Fill();

  }

  return EL::StatusCode::SUCCESS;

}
