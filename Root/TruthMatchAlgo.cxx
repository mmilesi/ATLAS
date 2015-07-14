/****************************************************************************************
 *
 * An algorithm that performs truth matching and truth-match classification for leptons
 *
 * M. Milesi (marco.milesi@cern.ch)
 *
 ****************************************************************************************/

// c++ include(s):
#include <iostream>
#include <typeinfo>
#include <sstream>

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "xAODBase/ObjectType.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"

// package include(s):
#include "HTopMultilepAnalysis/TruthMatchAlgo.h"
#include "HTopMultilepAnalysis/MCTruthClassifierDefs.h" // until not made available from ASG release
#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <xAODAnaHelpers/tools/ReturnCheck.h>
#include <xAODAnaHelpers/tools/ReturnCheckConfig.h>

// external tools include(s):

// ROOT include(s):
#include "TEnv.h"
#include "TFile.h"
#include "TSystem.h"

// external tools include(s):

// this is needed to distribute the algorithm to the workers
ClassImp(TruthMatchAlgo)


TruthMatchAlgo :: TruthMatchAlgo () :
  m_cutflowHist(nullptr),
  m_cutflowHistW(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("TruthMatchAlgo()", "Calling constructor");
    
  m_inContainerName_Electrons   = "";     
  m_inContainerName_Muons       = "";    
  m_inContainerName_Leptons     = "";    
  
  m_doMuonTruthPartMatching     = false;
  m_doMuonTrackMatching         = true;
}

TruthMatchAlgo::~TruthMatchAlgo() {}

EL::StatusCode  TruthMatchAlgo :: configure ()
{

  if ( !getConfig().empty() ) {

    // read in user configuration from text file
    TEnv *config = new TEnv(getConfig(true).c_str());
    if ( !config ) {
      Error("TruthMatchAlgo()", "Failed to initialize reading of config file. Exiting." );
      return EL::StatusCode::FAILURE;
    }
    Info("configure()", "Configuing TruthMatchAlgo Interface. User configuration read from : %s", getConfig().c_str());
    
    // read debug flag from .config file
    m_debug			 = config->GetValue("Debug" ,	  m_debug );
    m_useCutFlow		 = config->GetValue("UseCutFlow", m_useCutFlow );
    
    m_inContainerName_Electrons  = config->GetValue("InputContainerElectrons", m_inContainerName_Electrons.c_str());
    m_inContainerName_Muons      = config->GetValue("InputContainerMuons",     m_inContainerName_Muons.c_str());
    m_inContainerName_Leptons	 = config->GetValue("InputContainerLeptons",   m_inContainerName_Leptons.c_str());
    
    m_doMuonTruthPartMatching	 = config->GetValue("DoMuonTruthPartMatching" , m_doMuonTruthPartMatching );
    m_doMuonTrackMatching	 = config->GetValue("DoMuonTrackMatching"  , m_doMuonTrackMatching );
    
    config->Print();
    
    Info("configure()", "TruthMatchAlgo Interface succesfully configured!");

    delete config; config = nullptr;
  }
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  Info("setupJob()", "Calling setupJob");

  job.useXAOD ();
  xAOD::Init( "TruthMatchAlgo" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthMatchAlgo :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  Info("histInitialize()", "Calling histInitialize");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthMatchAlgo :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  Info("fileExecute()", "Calling fileExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthMatchAlgo :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "Calling changeInput");

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info("initialize()", "Initializing TruthMatchAlgo Interface...");

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("TruthMatchAlgo::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");

  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );

  if ( m_useCutFlow ) {
    TFile *file = wk()->getOutputFile("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");
    m_cutflow_bin  = m_cutflowHist->GetXaxis()->FindBin(m_name.c_str());
    m_cutflowHistW->GetXaxis()->FindBin(m_name.c_str());
  }

  Info("initialize()", "Number of events: %lld ", m_event->getEntries() );

  if ( configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  m_numEvent            = 0;
  m_numObject           = 0;
  m_numEventPass        = 0;
  m_weightNumEventPass  = 0;
  m_numObjectPass       = 0;
    
  m_isTruthMatchedDecor = nullptr          ; m_isTruthMatchedDecor          = new SG::AuxElement::Decorator< char >("isTruthMatched");	        // has a lepton truth match
  m_truthPdgIdDecor = nullptr	           ; m_truthPdgIdDecor              = new SG::AuxElement::Decorator< int >("truthPdgId");		// pdgId of the match particle
  m_truthTypeDecor = nullptr	           ; m_truthTypeDecor               = new SG::AuxElement::Decorator< int >("truthType"); 	        // type of the parent particle (according to MCTruthClassifier) - this decorates only muons (info is originally available only for the track!)
  m_truthOriginDecor = nullptr	           ; m_truthOriginDecor             = new SG::AuxElement::Decorator< int >("truthOrigin"); 	        // origin of the parent particle - this decorates only muons (info is originally available only for the track!)
  m_truthStatusDecor = nullptr	           ; m_truthStatusDecor             = new SG::AuxElement::Decorator< int >("truthStatus"); 	        // status of the match particle
  m_isChFlipDecor = nullptr		   ; m_isChFlipDecor                = new SG::AuxElement::Decorator< char >("isChFlip");		// reco has opposite charge wrt to primitive truth match
  m_isBremDecor = nullptr 		   ; m_isBremDecor                  = new SG::AuxElement::Decorator< char >("isBrem");  		// reco is matched to a brem lepton
  
  m_mcEvtWeightAcc = nullptr		   ; m_mcEvtWeightAcc		    = new SG::AuxElement::Accessor< float >("mcEventWeight");
  m_isTruthMatchedAcc = nullptr		   ; m_isTruthMatchedAcc	    = new SG::AuxElement::Accessor< char >("isTruthMatched");		
  m_isChFlipAcc = nullptr 		   ; m_isChFlipAcc		    = new SG::AuxElement::Accessor< char >("isChFlip"); 			   
  m_isBremAcc = nullptr		           ; m_isBremAcc		    = new SG::AuxElement::Accessor< char >("isBrem");			    
  m_truthPLAcc = nullptr  		   ; m_truthPLAcc		    = new SG::AuxElement::Accessor< TruthLink_t >("truthParticleLink");
  m_truthTypeAcc = nullptr		   ; m_truthTypeAcc		    = new SG::AuxElement::ConstAccessor< int >("truthType");
  m_truthOriginAcc = nullptr		   ; m_truthOriginAcc		    = new SG::AuxElement::ConstAccessor< int >("truthOrigin");
  m_truthMatchProbabilityAcc = nullptr     ; m_truthMatchProbabilityAcc     = new SG::AuxElement::Accessor<float>("truthMatchProbability");
  
  Info("initialize()", "TruthMatchAlgo Interface succesfully initialized!" );

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if ( m_debug ) { Info("execute()", "Applying TruthMatchAlgo..."); }

  // retrieve event 
  //
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("TruthMatchAlgo::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");

  // MC event weight 
  //
  float mcEvtWeight(1.0);
  if ( ! (*m_mcEvtWeightAcc).isAvailable( *eventInfo ) ) {
    Error("execute()", "mcEventWeight is not available as decoration! Aborting" );
    return EL::StatusCode::FAILURE;
  }
  mcEvtWeight = (*m_mcEvtWeightAcc)( *eventInfo );
  
  m_numEvent++;

  // retrieve leptonsCDV from store
  //
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);  
  RETURN_CHECK("TruthMatchAlgo::execute()", HelperFunctions::retrieve(leptonsCDV, m_inContainerName_Leptons, m_event, m_store, m_debug) , "");
  
  if ( m_debug ) { Info("execute()"," number of leptons: %lu ", leptonsCDV->size() ); }

  // -------------------------------------
  // Truth matching for leptons
  // -------------------------------------
  
  if ( m_isMC ) { 
    for ( auto lep_itr : *(leptonsCDV) ) {

      if ( lep_itr->type() == xAOD::Type::Electron ) {

	if ( m_debug ) { Info("execute()"," truth matching reco electron, pT = %2f ", lep_itr->pt() / 1e3 ); }

	if ( this->applyTruthMatchingElectronMC15( lep_itr ) != EL::StatusCode::SUCCESS ) {
	  Error("execute()", "Problem with applyTruthMatchingElectronMC15()! Aborting" );
	  return EL::StatusCode::FAILURE;
	} 

      } else if ( lep_itr->type() == xAOD::Type::Muon ) {

	if ( m_debug ) { Info("execute()"," truth matching reco muon, pT = %2f ", lep_itr->pt() / 1e3 ); }
	if ( this->applyTruthMatchingMuonMC15( lep_itr ) != EL::StatusCode::SUCCESS ) {
	  Error("execute()", "Problem with applyTruthMatchingMuonMC15()! Aborting" );
	  return EL::StatusCode::FAILURE;
	} 

      } 

    } // end loop over leptons

  } // end check isMC
  
  m_numEventPass++;
  m_weightNumEventPass += mcEvtWeight;
   
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  if ( m_debug ) { Info("postExecute()", "Calling postExecute"); }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  Info("finalize()", "Deleting pointers...");

  delete m_isTruthMatchedDecor; m_isTruthMatchedDecor = nullptr;
  delete m_truthTypeDecor; m_truthTypeDecor = nullptr;
  delete m_truthPdgIdDecor; m_truthPdgIdDecor = nullptr;	 
  delete m_truthOriginDecor; m_truthOriginDecor = nullptr;	 
  delete m_truthStatusDecor; m_truthStatusDecor = nullptr;	 
  delete m_isChFlipDecor; m_isChFlipDecor = nullptr;		 
  delete m_isBremDecor; m_isBremDecor = nullptr;	 
  
  delete m_mcEvtWeightAcc; m_mcEvtWeightAcc = nullptr;		 
  delete m_isTruthMatchedAcc; m_isTruthMatchedAcc = nullptr; 
  delete m_truthTypeAcc; m_truthTypeAcc = nullptr;  
  delete m_truthOriginAcc; m_truthOriginAcc = nullptr;	    
  delete m_isChFlipAcc; m_isChFlipAcc = nullptr;	 
  delete m_isBremAcc; m_isBremAcc = nullptr;	 
  delete m_truthPLAcc; m_truthPLAcc = nullptr;   
  delete m_truthTypeAcc; m_truthTypeAcc = nullptr;  
  delete m_truthOriginAcc; m_truthOriginAcc = nullptr;	    
  delete m_truthMatchProbabilityAcc; m_truthMatchProbabilityAcc = nullptr;
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  Info("histFinalize()", "Calling histFinalize");
  
  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );
  }
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo ::  checkChargeFlipMC15 ( const xAOD::IParticle* recoPart, const xAOD::TruthParticle* matchTruth )
{
  
  // default decorations
  //
  (*m_isChFlipDecor)( *recoPart )  = 0;  
  (*m_isBremDecor)( *recoPart )    = 0;  

  float reco_charge(0.0);	
  if ( recoPart->type() == xAOD::Type::Electron ) {
    if ( m_debug ) { Info("checkChargeFlip()", "This reco lepton is an electron" ); }
    reco_charge = dynamic_cast<const xAOD::Electron*>(recoPart)->charge();
  } else if ( recoPart->type() == xAOD::Type::Muon ) {
    if ( m_debug ) { Info("checkChargeFlip()", "This reco lepton is a muon" ); }  
    reco_charge = dynamic_cast<const xAOD::Muon*>(recoPart)->charge();
  }
  if ( !reco_charge ) {
     Error("checkChargeFlip()", "Reco particle has zero charge. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;	    
  }

  xAOD::TruthParticle* primitiveTruth(nullptr);

  if ( ! (*m_truthTypeAcc).isAvailable( *recoPart ) ) {
     Error("checkChargeFlip()", "No accessor to truthType available for this reco lepton. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }

  // case 1: 
  //
  // lepton (in most cases, an electron) is matched to a truth lepton which is part of a bremmmstrahlung shower. 
  // In this case, we need to go back until we find the original lepton that radiated the photon.
  // The charge of this primitive lepton is the one to look at! 
  //
  // look at 'Background'-type el/mu (see MCTruthClassifier.h)
  //
  if ( (*m_truthTypeAcc)( *recoPart ) == 4 || (*m_truthTypeAcc)( *recoPart ) == 8 ) {

    if ( m_debug ) { Info("checkChargeFlip()", "This reco lepton (charge: %f ) is matched to a secondary truth lepton. Let's go back until we find the primitive", reco_charge ); }

    bool foundPrimitive(false), isBrem(false);
    primitiveTruth = const_cast<xAOD::TruthParticle*>( matchTruth->parent(0) );
    
    while ( !foundPrimitive ) {
     
      if ( primitiveTruth->prodVtx()->barcode() < -200000 ) { 
	
	if ( m_debug ) { Info("checkChargeFlip()", "Parent has pdgId: %i , prodVtx barcode: %i - Need to go backwards in the decay chain", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	
	primitiveTruth = const_cast<xAOD::TruthParticle*>( primitiveTruth->parent(0) );
	
	// do this only once
	//
	if ( !isBrem ) { 
	  isBrem = true;
	  (*m_isBremDecor)( *recoPart ) = 1; 
	}
	
      } else { 
      
	if ( m_debug ) { Info("checkChargeFlip()", "We found the primitive! pdgId: %i , prodVtx barcode: %i - Stop here", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	foundPrimitive = true; 
	
      }

    }
  
  }
  // case 2:
  //
  // lepton is matched to a truth lepton which is not produced in a secondary interaction. 
  //
  else 
  {
    primitiveTruth = const_cast<xAOD::TruthParticle*>( matchTruth );
  }

  float truth_charge	= primitiveTruth->charge();
  int truth_norm_charge = static_cast<int>( truth_charge / fabs(truth_charge) );   
  int reco_norm_charge  = static_cast<int>( reco_charge  / fabs(reco_charge)  ); 
  
  if ( ( reco_norm_charge * truth_norm_charge ) < 0 ) { 
    if ( m_debug ) { Info("checkChargeFlip()", "Reco norm charge: %i \n, Primitive truth charge: %f  norm charge: %i  pdgId: %i  prodVtxBarcode: %i \n It's charge flip!", reco_norm_charge, truth_charge, truth_norm_charge, primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
    (*m_isChFlipDecor)( *recoPart ) = 1; 
  }

  return StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo ::  applyTruthMatchingElectronMC15 ( const xAOD::IParticle* recoPart )
{
  
  // return immediately if input particle is not an electron 
  //
  if ( !( recoPart->type() == xAOD::Type::Electron ) ) {
      Warning("applyTruthMatchingElectronMC15()", "Not passing an electron! Won't try anything"); 
      return StatusCode::SUCCESS; 	  	  
  }

  // truth particle types are defined in HTopMultilepAnalysis/MCTruthClassifierDefs.h:
  //
  // further explaination can be found in: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/MCTruthClassifier
  // 

  // decorate reconstructed particle with default values
  //
  (*m_isTruthMatchedDecor)( *recoPart )	         = 0;
  (*m_truthPdgIdDecor)( *recoPart )              = 0;
  (*m_truthStatusDecor)( *recoPart )	         = -1;  

  // Now try to do the matching
  //
  // For electrons, the link to truth matching particle (and some useful info) is already saved in ElectronCollection
  //
  if ( ! (*m_truthPLAcc).isAvailable( *recoPart ) ) {
     Error("applyTruthMatchingElectronMC15()", "No link available to truth match for this reco electron. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }
  if ( ! (*m_truthPLAcc)( *recoPart ).isValid() ) {
     Error("applyTruthMatchingElectronMC15()", "Link to truth match for this reco electron is invalid. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }
  const xAOD::TruthParticle* matchTruthEl = *( (*m_truthPLAcc)(*recoPart) );

  // if there's no matching truth electron, abort.
  //
  if ( !matchTruthEl ) {
     Error( "applyTruthMatchingElectronMC15()", "This reco electron is not matched to a generic truth particle. This shouldn't happen. Aborting");
     return StatusCode::FAILURE;
  } 
  
  // decorate with true if the truth match is an electron
  //
  bool isTMElectron(false);
  if ( matchTruthEl->isElectron() ) {
      isTMElectron = true;
     (*m_isTruthMatchedDecor)( *recoPart ) = 1;
  }

  // store the pdgId of the match
  //
  static SG::AuxElement::Accessor< int > pdgIdAcc("pdgId");
  if ( pdgIdAcc.isAvailable( *matchTruthEl ) ) {
  
    if ( m_debug ) { Info( "applyTruthMatchingElectronMC15()", "decorating truthPdgId with value: %i", matchTruthEl->pdgId() ); }
    (*m_truthPdgIdDecor)( *recoPart ) = matchTruthEl->pdgId();
    
  }
  // store the status of the match
  //
  static SG::AuxElement::Accessor< int > statusAcc("status");
  if ( statusAcc.isAvailable( *matchTruthEl ) ) {
  
    if ( m_debug ) { Info( "applyTruthMatchingElectronMC15()", "decorating truthStatus with value: %i", matchTruthEl->status() ); }
    (*m_truthStatusDecor)( *recoPart ) = matchTruthEl->status();
  
  }
    
  // check if electron is charge flip  
  //
  if ( isTMElectron ) {
    if ( this->checkChargeFlipMC15( recoPart, matchTruthEl ) != EL::StatusCode::SUCCESS ) {
      Error("applyTruthMatchingElectronMC15()", "Problem with checkChargeFlipMC15(). Aborting"); 
      return EL::StatusCode::FAILURE;
    }
  }
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthMatchAlgo ::  applyTruthMatchingMuonMC15 ( const xAOD::IParticle* recoPart )
{
  
  // return immediately if input particle is not a muon 
  //
  if ( !( recoPart->type() == xAOD::Type::Muon ) ) {     
    Warning("applyTruthMatchingMuonMC15()", "Not passing a muon! Won't try anything"); 
    return EL::StatusCode::SUCCESS;
  }

  // Now try to do the matching
  //  
  // It can be either done by finding the link to MuonTruthParticles container, or by matching the muon track (default for |eta| < 2.5)
  // See the header file for more info.
  //
  const xAOD::TruthParticleContainer* muonTruthPartContainer(nullptr);
  RETURN_CHECK("TruthMatchAlgo::applyTruthMatchingMuonMC15()", HelperFunctions::retrieve(muonTruthPartContainer, "MuonTruthParticles", m_event, m_store, m_debug) , "");
  
  if ( m_doMuonTrackMatching  ) {
     
     if ( fabs( recoPart->eta() ) < 2.5 ) {
      
       if ( this->doMuonTrackMatching( recoPart ) != EL::StatusCode::SUCCESS ) {
         Error("applyTruthMatchingMuonMC15()", "Problem with doMuonTrackMatching() for this muon ( | eta | < 2.5 ). Aborting");
         return EL::StatusCode::FAILURE;
       }
       
     } else {
     
       if ( this->doMuonTruthPartMatching( recoPart ) != EL::StatusCode::SUCCESS ) {
         Error("applyTruthMatchingMuonMC15()", "Problem with doMuonTruthPartMatching() for this forward muon ( | eta| > 2.5 ). Aborting");
         return EL::StatusCode::FAILURE;
       }
       
     }  
     
  } else if ( m_doMuonTruthPartMatching && ( this->doMuonTruthPartMatching( recoPart ) != EL::StatusCode::SUCCESS ) ) {
     
     Error("applyTruthMatchingMuonMC15()", "Problem with doMuonTruthPartMatching(). Aborting");
     return EL::StatusCode::FAILURE;
  
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo :: doMuonTrackMatching( const xAOD::IParticle* recoPart )
{

   // decorate reconstructed particle with default values
   //
   (*m_isTruthMatchedDecor)( *recoPart )	  = 0;
   (*m_truthTypeDecor)( *recoPart )               = 0; // need it b/c for muons we need to pass from the track/muon truth container
   (*m_truthPdgIdDecor)( *recoPart )              = 0; 
   (*m_truthOriginDecor)( *recoPart )             = 0; // need it b/c for muons we need to pass from the track/muon truth container
   (*m_truthStatusDecor)( *recoPart )	          = -1;  
  
   // get the reco muon ID track particle 
   //
   const xAOD::TrackParticle* trk(nullptr); 
   ElementLink< xAOD::TrackParticleContainer > trkLink = dynamic_cast<const xAOD::Muon*>(recoPart)->inDetTrackParticleLink();
   
   if ( !trkLink.isValid() ) {
     Error("doMuonTrackMatching()", "Link to ID track particle for this reco muon is invalid. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;   
   }
   trk = *trkLink;

   // get the truth particle matching the ID track
   //
   if ( ! (*m_truthPLAcc).isAvailable( *trk ) ) {
      Error("doMuonTrackMatching()", "No link available to truth match for this reco muon's ID track. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;		 
   }
   if ( ! (*m_truthPLAcc)( *trk ).isValid() ) {
      Error("doMuonTrackMatching()", "Link to truth match for this reco muon's ID track is invalid. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;
   }
   const xAOD::TruthParticle* matchTruthMu = *( (*m_truthPLAcc)(*trk) );

   // if there is no matching truth particle for the ID track, abort 
   //
   if ( !matchTruthMu ) {
      if ( m_debug ) { Info("doMuonTrackMatching()", "No truth match for this reco muon's ID track. This shouldn't happen. Aborting"); }
      return StatusCode::SUCCESS;    
   } 

   // retrieve track truth MC probability (should be always available for ID tracks!) 
   //
   float trk_prob(-1.0); 
   if ( m_truthMatchProbabilityAcc->isAvailable( *trk ) ) { 
     trk_prob = (*m_truthMatchProbabilityAcc)( *trk ); 
   }

   // decorate with true if the truth match is a muon (NB: since we are looking at the track, this might not always be the case!!), 
   // and the track mc probability (when available) is high enough
   // 
   bool isTMMuon(false);
   if ( matchTruthMu->isMuon()  && ( trk_prob < 0.0 || trk_prob > 0.8 ) ) {
      isTMMuon = true;
     (*m_isTruthMatchedDecor)( *recoPart )     = 1;
   } 
 
   // store the type of the match: pass the track type info to the reco muon
   //
   if ( ! (*m_truthTypeAcc).isAvailable( *trk ) ) {
     Error("doMuonTrackMatching()", "No truth type info available for this muon's ID track matching truth particle. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthTrkMatchType = (*m_truthTypeAcc)(*trk);
   (*m_truthTypeDecor)( *recoPart ) = truthTrkMatchType;
  
   // store the pdgId of the match
   //
   static SG::AuxElement::Accessor< int > pdgIdAcc("pdgId");
   if ( pdgIdAcc.isAvailable( *matchTruthMu ) ) {
  
     if ( m_debug ) { Info( "doMuonTrackMatching()", "decorating truthPdgId with value: %i", matchTruthMu->pdgId() ); }
     (*m_truthPdgIdDecor)( *recoPart ) = matchTruthMu->pdgId();
   
   }
   // store the status of the match
   //
   static SG::AuxElement::Accessor< int > statusAcc("status");
   if ( statusAcc.isAvailable( *matchTruthMu ) ) {
  
     if ( m_debug ) { Info( "doMuonTrackMatching()", "decorating truthStatus with value: %i", matchTruthMu->status() ); }
     (*m_truthStatusDecor)( *recoPart ) = matchTruthMu->status();
  
   }

   // store the pdgId of the parent particle of the match: pass the track origin info to the reco muon
   //
   if ( ! (*m_truthOriginAcc).isAvailable( *trk ) ) {
     Error("doMuonTrackMatching()", "No truth origin info available for this muon's ID track matching truth particle. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthTrkMatchOrigin = (*m_truthOriginAcc)(*trk);
   (*m_truthOriginDecor)( *recoPart ) = truthTrkMatchOrigin;
   
   // check if muon is charge flip
   //  
   if ( isTMMuon ) {
     if ( this->checkChargeFlipMC15( recoPart, matchTruthMu ) != EL::StatusCode::SUCCESS ) {
       Error("doMuonTrackMatching()", "Problem with checkChargeFlipMC15(). Aborting"); 
       return EL::StatusCode::FAILURE;
     }
   }  
   
   return EL::StatusCode::SUCCESS;

}

EL::StatusCode TruthMatchAlgo :: doMuonTruthPartMatching ( const xAOD::IParticle* recoPart )
{

   // decorate reconstructed particle with default values
   //
   (*m_isTruthMatchedDecor)( *recoPart )	  = 0;
   (*m_truthTypeDecor)( *recoPart )               = 0; // need it b/c for muons we need to pass from the track/muon truth container
   (*m_truthPdgIdDecor)( *recoPart )              = 0; 
   (*m_truthOriginDecor)( *recoPart )             = 0; // need it b/c for muons we need to pass from the track/muon truth container
   (*m_truthStatusDecor)( *recoPart )	          = -1;  

   // get the truth muon matching the reco muon
   //
   const xAOD::TruthParticle* matchTruthMu(nullptr);

   if ( ! (*m_truthPLAcc).isAvailable( *recoPart ) ) {
      Error("doMuonTruthPartMatching()", "No link available to truth match for this reco muon. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;		 
   }
   if ( ! (*m_truthPLAcc)( *recoPart ).isValid() ) {
      Error("doMuonTruthPartMatching()", "Link to truth match for this reco muon is invalid. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;
   }
   matchTruthMu = *( (*m_truthPLAcc)(*recoPart) );

   // if there is no matching truth muon, abort
   //
   if ( !matchTruthMu ) {
      Error("doMuonTruthPartMatching()", "No truth matching for this reco muon. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;    
   }
   
   // decorate with this if the truth match is a muon ( should ALWAYS be the case, since the truth we are getting is in the "MuonTruthParticles" container.! )
   //
   bool isTMMuon(false);
   if ( matchTruthMu->isMuon() ) {
     isTMMuon = true;
     (*m_isTruthMatchedDecor)( *recoPart ) = 1;
   } 

   // store the type of the match: pass the truth muon type info to the reco muon
   //
   if ( ! (*m_truthTypeAcc).isAvailable( *matchTruthMu ) ) {
     Error("doMuonTruthPartMatching()", "No truth type info available for this muon's matching truth particle. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthMatchType = (*m_truthTypeAcc)(*matchTruthMu);
   (*m_truthTypeDecor)( *recoPart ) = truthMatchType;
  
   // store the pdgId of the match
   //
   static SG::AuxElement::Accessor< int > pdgIdAcc("pdgId");
   if ( pdgIdAcc.isAvailable( *matchTruthMu ) ) {
  
     if ( m_debug ) { Info( "doMuonTruthPartMatching()", "decorating truthPdgId with value: %i", matchTruthMu->pdgId() ); }
     (*m_truthPdgIdDecor)( *recoPart ) = matchTruthMu->pdgId();
   
   }
   // store the status of the match
   //
   static SG::AuxElement::Accessor< int > statusAcc("status");
   if ( statusAcc.isAvailable( *matchTruthMu ) ) {
  
     if ( m_debug ) { Info( "doMuonTruthPartMatching()", "decorating truthStatus with value: %i", matchTruthMu->status() ); }
     (*m_truthStatusDecor)( *recoPart ) = matchTruthMu->status();
  
   }

   // store the pdgId of the parent particle of the match: pass the truth muon origin info to the reco muon
   //
   if ( ! (*m_truthOriginAcc).isAvailable( *matchTruthMu ) ) {
     Error("doMuonTruthPartMatching()", "No truth origin info available for this muon's matching truth particle. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthMatchOrigin = (*m_truthOriginAcc)(*matchTruthMu);
   (*m_truthOriginDecor)( *recoPart ) = truthMatchOrigin;
   

   // check if muon is charge flip  
   //
   if ( isTMMuon ) {
     if ( this->checkChargeFlipMC15( recoPart, matchTruthMu ) != EL::StatusCode::SUCCESS ) {
       Error("doMuonTruthPartMatching()", "Problem with checkChargeFlipMC15(). Aborting"); 
       return EL::StatusCode::FAILURE;
     }
   }    
   
   return EL::StatusCode::SUCCESS;

}

