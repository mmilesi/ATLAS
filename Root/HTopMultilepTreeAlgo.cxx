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

HTopMultilepTreeAlgo :: HTopMultilepTreeAlgo() {}

HTopMultilepTreeAlgo :: HTopMultilepTreeAlgo (const std::string name, const std::string configName) :
  TreeAlgo (name, configName), // probably no need to do this ...   
  m_name(name),
  m_configName(configName),
  m_helpTree(nullptr)
{}


EL::StatusCode HTopMultilepTreeAlgo :: initialize ()
{
  Info("initialize()", m_name.c_str());
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  this->treeInitialize();
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepTreeAlgo :: treeInitialize ()
{

  Info("treeInitialize()", "%s", m_name.c_str() );
  // needed here and not in initalize since this is called first
  Info("treeInitialize()", "Attempting to configure using: %s", m_configName.c_str());
  
  TTree * outTree = new TTree(m_name.c_str(),m_name.c_str());
  if ( !outTree ) {
    Error("treeInitialize()","Failed to instantiate output tree!");
    return EL::StatusCode::FAILURE;
  }

  // get the input from user which determines which branches are created (and more)!
  this->configure(); 
  
  // get the file we created already
  TFile *outTreeFile = wk()->getOutputFile ("tree");  
  m_helpTree =  new HTopMultilepTree( m_event, outTree, outTreeFile, 1e0, m_debug ); // 1e0 = MeV, default 1e3 = GeV
  // tell the tree to go into the file
  outTree->SetDirectory( outTreeFile );
  // uncomment if want to add to same file as ouput histograms
  wk()->addOutput( outTree ); 

  m_helpTree->AddEvent(m_evtDetailStr);

  if ( !m_muContainerName.empty() )     {   m_helpTree->AddMuons    (m_muDetailStr);      }
  if ( !m_elContainerName.empty() )     {   m_helpTree->AddElectrons(m_elDetailStr);      }
  if ( !m_jetContainerName.empty() )    {   m_helpTree->AddJets     (m_jetDetailStr);     }
  if ( !m_fatJetContainerName.empty() ) {   m_helpTree->AddFatJets  (m_fatJetDetailStr);  }
  if ( !m_tauContainerName.empty() )    {   m_helpTree->AddTaus     (m_tauDetailStr);     }
  if ( !m_lepContainerName.empty() )    {   m_helpTree->AddLeptons  ();                   }

  Info("treeInitialize()", "Successfully initialized output tree");
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepTreeAlgo :: configure ()
{
  
  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  RETURN_CHECK_CONFIG("HTopMultilepTreeAlgo::configure()", m_configName);

  // the file exists, use TEnv to read it off
  TEnv* config = new TEnv(m_configName.c_str());
  
  // read debug flag from .config file
  m_debug                   = config->GetValue("Debug" ,      false );
  
  m_evtDetailStr            = config->GetValue("EventDetailStr",       "");
  m_muDetailStr             = config->GetValue("MuonDetailStr",        "");
  m_elDetailStr             = config->GetValue("ElectronDetailStr",    "");
  m_jetDetailStr            = config->GetValue("JetDetailStr",         "");
  m_fatJetDetailStr         = config->GetValue("FatJetDetailStr",      "");
  m_tauDetailStr            = config->GetValue("TauDetailStr",         "kinematic");

  m_muContainerName         = config->GetValue("MuonContainerName",       "");
  m_elContainerName         = config->GetValue("ElectronContainerName",   "");
  m_jetContainerName        = config->GetValue("JetContainerName",        "");
  m_fatJetContainerName     = config->GetValue("FatJetContainerName",     "");
  m_tauContainerName        = config->GetValue("TauContainerName",        "");
  m_lepContainerName        = config->GetValue("LepContainerName",        "");

  Info("configure()", "Loaded in configuration values");

  // everything seems preliminarily ok, let's print config and say we were successful
  config->Print();
  
  delete config; config = nullptr;
  
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

  m_helpTree->FillEvent( eventInfo, m_event );

  // for the containers the were supplied, fill the appropiate vectors
  
  const xAOD::MuonContainer* inMuons(nullptr); 
  const xAOD::ElectronContainer* inElectrons(nullptr);
  const xAOD::JetContainer* inJets(nullptr);
  const xAOD::JetContainer* inFatJets(nullptr);
  const xAOD::TauJetContainer* inTaus(nullptr);
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);

  if ( !m_muContainerName.empty() ) {	
    
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inMuons, m_muContainerName, m_event, m_store, m_debug) , "");
    
    // sort inMuons, and pass the reference to FillMuons()
    const xAOD::MuonContainer inMuonsSorted = HelperFunctions::sort_container_pt( inMuons );
    m_helpTree->FillMuons( &inMuonsSorted, primaryVertex );

  }
  if ( !m_elContainerName.empty() ) { 	

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inElectrons, m_elContainerName, m_event, m_store, m_debug) , "");
    
    // sort inElectrons, and pass the reference to FillElectrons()
    const xAOD::ElectronContainer inElectronsSorted = HelperFunctions::sort_container_pt( inElectrons );
    m_helpTree->FillElectrons( &inElectronsSorted, primaryVertex );
    
  }
  if ( !m_jetContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inJets, m_jetContainerName, m_event, m_store, m_debug) , "");

    // sort inJets, and pass the reference to FillJets()
    const xAOD::JetContainer inJetsSorted = HelperFunctions::sort_container_pt( inJets );
    m_helpTree->FillJets( &inJetsSorted, HelperFunctions::getPrimaryVertexLocation(vertices) );

    //m_helpTree->FillJets( inJets, HelperFunctions::getPrimaryVertexLocation(vertices) );
  }
  if ( !m_fatJetContainerName.empty() ) {

    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inFatJets, m_fatJetContainerName, m_event, m_store, m_debug) , "");
    m_helpTree->FillFatJets( inFatJets );

  }
  if ( !m_tauContainerName.empty() ) {	
    
    RETURN_CHECK("HTopMultilepTreeAlgo::execute()", HelperFunctions::retrieve(inTaus, m_tauContainerName, m_event, m_store, m_debug) , "");
    
    // sort inTaus, and pass the reference to FillTaus()
    const xAOD::TauJetContainer inTausSorted = HelperFunctions::sort_container_pt( inTaus );
    m_helpTree->FillTaus( &inTausSorted );
    
  }
  if ( !m_lepContainerName.empty() ) {	
    
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, m_lepContainerName, m_event, m_store, m_debug) ,"");

    // sort inLeptons, and pass the reference to FillLeptons()  
    const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );
    m_helpTree->FillLeptons( &leptonsSorted );

  }    
    
  // fill the tree
  m_helpTree->Fill();
   
  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepTreeAlgo :: finalize () {  

  if ( m_helpTree ) { delete  m_helpTree; m_helpTree = nullptr; }
  
  return EL::StatusCode::SUCCESS;
}

