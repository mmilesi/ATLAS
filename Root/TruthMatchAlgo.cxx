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

  m_doDC14Matching              = false;
  m_doMC15Matching              = true;
  
  m_doMuonDeltaRMatching        = false;
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
    Info("configure()", "Configuing TruthMatchAlgo Interface. User configuration read from : %s \n", getConfig().c_str());
    
    // read debug flag from .config file
    m_debug			 = config->GetValue("Debug" ,	  m_debug );
    m_useCutFlow		 = config->GetValue("UseCutFlow", m_useCutFlow );
    
    m_inContainerName_Electrons  = config->GetValue("InputContainerElectrons", m_inContainerName_Electrons.c_str());
    m_inContainerName_Muons      = config->GetValue("InputContainerMuons",     m_inContainerName_Muons.c_str());
    m_inContainerName_Leptons	 = config->GetValue("InputContainerLeptons",   m_inContainerName_Leptons.c_str());
    
    m_doDC14Matching		 = config->GetValue("DoDC14Matching",  m_doDC14Matching );
    m_doMC15Matching		 = config->GetValue("DoMC15Matching",  m_doMC15Matching );
    
    m_doMuonDeltaRMatching	 = config->GetValue("DoMuonDeltaRMatching" , m_doMuonDeltaRMatching );
    m_doMuonTrackMatching	 = config->GetValue("DoMuonTrackMatching"  , m_doMuonTrackMatching );
    
    config->Print();
    
    Info("configure()", "TruthMatchAlgo Interface succesfully configured! \n");

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
  m_isTruthMatchedIsoDecor = nullptr	   ; m_isTruthMatchedIsoDecor       = new SG::AuxElement::Decorator< char >("isTruthMatchedIso");       // prompt leptons
  m_isTruthMatchedNonIsoDecor = nullptr	   ; m_isTruthMatchedNonIsoDecor    = new SG::AuxElement::Decorator< char >("isTruthMatchedNonIso");    // non-prompt leptons (from HF hadrons, or decays of hadrons in jets)
  m_isTruthMatchedSecondaryDecor = nullptr ; m_isTruthMatchedSecondaryDecor = new SG::AuxElement::Decorator< char >("isTruthMatchedSecondary"); // from secondary material interaction (e.g. conversion)
  m_isTruthMatchedNoProdVtxDecor = nullptr ; m_isTruthMatchedNoProdVtxDecor = new SG::AuxElement::Decorator< char >("isTruthMatchedNoProdVtx"); // matched to a lepton w/o production vertex
  m_isTruthMatchedUnknownDecor = nullptr   ; m_isTruthMatchedUnknownDecor   = new SG::AuxElement::Decorator< char >("isTruthMatchedUnknown");   // matched to an unknown truth particle
  m_isTruthMatchedOtherDecor = nullptr	   ; m_isTruthMatchedOtherDecor     = new SG::AuxElement::Decorator< char >("isTruthMatchedOther");     // matched to a non-lepton truth particle
  m_truthPdgIdDecor = nullptr	           ; m_truthPdgIdDecor              = new SG::AuxElement::Decorator< int >("truthPdgId");		// pdgId of the match particle
  m_truthTypeDecor = nullptr	           ; m_truthTypeDecor               = new SG::AuxElement::Decorator< int >("truthType"); 	        // type of the parent particle (according to MCTruthClassifier) - this decorates only muons (info is originally available only for the track!)
  m_truthOriginDecor = nullptr	           ; m_truthOriginDecor             = new SG::AuxElement::Decorator< int >("truthOrigin"); 	        // origin of the parent particle - this decorates only muons (info is originally available only for the track!)
  m_truthStatusDecor = nullptr	           ; m_truthStatusDecor             = new SG::AuxElement::Decorator< int >("truthStatus"); 	        // status of the match particle
  m_isChFlipDecor = nullptr		   ; m_isChFlipDecor                = new SG::AuxElement::Decorator< char >("isChFlip");		// reco has opposite charge wrt to primitive truth match
  m_isBremDecor = nullptr 		   ; m_isBremDecor                  = new SG::AuxElement::Decorator< char >("isBrem");  		// reco is matched to a brem lepton
  
  m_mcEvtWeightAcc = nullptr		   ; m_mcEvtWeightAcc		    = new SG::AuxElement::Accessor< float >("mcEventWeight");
  m_isTruthMatchedAcc = nullptr		   ; m_isTruthMatchedAcc	    = new SG::AuxElement::Accessor< char >("isTruthMatched");		
  m_isTruthMatchedIsoAcc = nullptr	   ; m_isTruthMatchedIsoAcc	    = new SG::AuxElement::Accessor< char >("isTruthMatchedIso");	
  m_isTruthMatchedNonIsoAcc = nullptr	   ; m_isTruthMatchedNonIsoAcc      = new SG::AuxElement::Accessor< char >("isTruthMatchedNonIso");	  
  m_isTruthMatchedSecondaryAcc = nullptr   ; m_isTruthMatchedSecondaryAcc   = new SG::AuxElement::Accessor< char >("isTruthMatchedSecondary");  
  m_isTruthMatchedNoProdVtxAcc = nullptr   ; m_isTruthMatchedNoProdVtxAcc   = new SG::AuxElement::Accessor< char >("isTruthMatchedNoProdVtx"); 
  m_isTruthMatchedUnknownAcc = nullptr	   ; m_isTruthMatchedUnknownAcc	    = new SG::AuxElement::Accessor< char >("isTruthMatchedUnknown"); 
  m_isTruthMatchedOtherAcc = nullptr	   ; m_isTruthMatchedOtherAcc	    = new SG::AuxElement::Accessor< char >("isTruthMatchedOther"); 
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
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("TruthMatchAlgo::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");

  bool isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );

  // MC event weight 
  float mcEvtWeight(1.0);
  if ( ! (*m_mcEvtWeightAcc).isAvailable( *eventInfo ) ) {
    Error("execute()", "mcEventWeight is not available as decoration! Aborting" );
    return EL::StatusCode::FAILURE;
  }
  mcEvtWeight = (*m_mcEvtWeightAcc)( *eventInfo );
  
  m_numEvent++;

  // retrieve leptonsCDV from store
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);  
  RETURN_CHECK("TruthMatchAlgo::execute()", HelperFunctions::retrieve(leptonsCDV, m_inContainerName_Leptons, m_event, m_store, m_debug) , "");
  
  if ( m_debug ) { Info("execute()"," number of leptons: %lu ", leptonsCDV->size() ); }

  // -------------------------------------
  // Truth matching for leptons
  // -------------------------------------
  
  if ( isMC ) { 
    for ( auto lep_itr : *(leptonsCDV) ) {

      if ( lep_itr->type() == xAOD::Type::Electron ) {

	if ( m_debug ) { Info("execute()"," truth matching reco electron, pT = %2f ", lep_itr->pt() / 1e3 ); }

	if ( m_doMC15Matching && ( this->applyTruthMatchingElectronMC15( lep_itr ) != EL::StatusCode::SUCCESS ) ) {
	  Error("execute()", "Problem with applyTruthMatchingElectronMC15()! Aborting" );
	  return EL::StatusCode::FAILURE;
	} 
	else if ( m_doDC14Matching && ( this->applyTruthMatchingDC14( lep_itr ) != EL::StatusCode::SUCCESS ) ) { 
	  Error("execute()", "Problem with applyTruthMatchingDC14()! Aborting" );
	  return EL::StatusCode::FAILURE;
	}
      } 
      else if ( lep_itr->type() == xAOD::Type::Muon ) {

	if ( m_debug ) { Info("execute()"," truth matching reco muon, pT = %2f ", lep_itr->pt() / 1e3 ); }
	if ( m_doMC15Matching && ( this->applyTruthMatchingMuonMC15( lep_itr ) != EL::StatusCode::SUCCESS ) ) {
	  Error("execute()", "Problem with applyTruthMatchingMuonMC15()! Aborting" );
	  return EL::StatusCode::FAILURE;
	} 
	else if ( m_doDC14Matching && ( this->applyTruthMatchingDC14( lep_itr ) != EL::StatusCode::SUCCESS ) ) { 
	  Error("execute()", "Problem with applyTruthMatchingDC14()! Aborting" );
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
  delete m_isTruthMatchedIsoDecor; m_isTruthMatchedIsoDecor = nullptr;	 
  delete m_isTruthMatchedNonIsoDecor; m_isTruthMatchedNonIsoDecor = nullptr;
  delete m_isTruthMatchedSecondaryDecor; m_isTruthMatchedSecondaryDecor = nullptr;
  delete m_isTruthMatchedUnknownDecor; m_isTruthMatchedUnknownDecor = nullptr;	 
  delete m_isTruthMatchedOtherDecor; m_isTruthMatchedOtherDecor = nullptr;
  delete m_truthTypeDecor; m_truthTypeDecor = nullptr;
  delete m_truthPdgIdDecor; m_truthPdgIdDecor = nullptr;	 
  delete m_truthOriginDecor; m_truthOriginDecor = nullptr;	 
  delete m_truthStatusDecor; m_truthStatusDecor = nullptr;	 
  delete m_isChFlipDecor; m_isChFlipDecor = nullptr;		 
  delete m_isBremDecor; m_isBremDecor = nullptr;	 
  
  delete m_mcEvtWeightAcc; m_mcEvtWeightAcc = nullptr;		 
  delete m_isTruthMatchedAcc; m_isTruthMatchedAcc = nullptr; 
  delete m_isTruthMatchedIsoAcc; m_isTruthMatchedIsoAcc = nullptr;
  delete m_isTruthMatchedNonIsoAcc; m_isTruthMatchedNonIsoAcc = nullptr;	 
  delete m_isTruthMatchedSecondaryAcc; m_isTruthMatchedSecondaryAcc = nullptr;
  delete m_isTruthMatchedUnknownAcc; m_isTruthMatchedUnknownAcc = nullptr;	
  delete m_isTruthMatchedOtherAcc; m_isTruthMatchedOtherAcc = nullptr;	
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

EL::StatusCode TruthMatchAlgo ::  applyTruthMatchingDC14 ( const xAOD::IParticle* recoPart )
{
    
  const xAOD::TruthEventContainer* truthEventContainer(nullptr);
  RETURN_CHECK("applyTruthMatchingDC14()", HelperFunctions::retrieve(truthEventContainer, "TruthEvent", m_event, m_store, m_debug) , "");

  // decorate reconstructed particle with default values
  (*m_isTruthMatchedDecor)( *recoPart )	         = 0;  
  (*m_isTruthMatchedIsoDecor)( *recoPart )       = 0;  
  (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 0;   
  (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 0;  
  (*m_isTruthMatchedNoProdVtxDecor)( *recoPart ) = 0;  
  (*m_isTruthMatchedOtherDecor)( *recoPart )	 = 0;  
  (*m_truthPdgIdDecor)( *recoPart )         = 0; 
  (*m_truthOriginDecor)( *recoPart )	 = 0;  
  (*m_truthStatusDecor)( *recoPart )	 = 0;  
  
  // this will be the best-matching truth particle (if ever found)
  const xAOD::TruthParticle* matchTruth(nullptr);

  xAOD::IParticle::FourMom_t recoPart4mom = recoPart->p4();	
     
  // loop over truth events
  for ( auto evt_it : *(truthEventContainer) ) {
    
    double minDR(0.2);
    bool checkOnlyElectrons(false); // this needed to make sure you pick the right electron match in case of brem electrons
    
    int nPart = evt_it->nTruthParticles();
    // loop over truth particles in this truth event
    for ( int iTPart = 0 ; iTPart < nPart ; ++iTPart ) {
    
      const xAOD::TruthParticle* truthPart = evt_it->truthParticle(iTPart);
      if ( !truthPart ) { continue; }
      
      // skip all the protons : they have null p4() (ROOT will complain when calculating pseudorapidity), and won't be ever interesting for us
      if ( truthPart->pdgId() == 2212 ) { continue; }

      // special treatment for electrons (see below)	
      if ( checkOnlyElectrons && !truthPart->isElectron() ) { continue; }

      //if ( iTPart < 10 ) {
      //  Info( "applyTruthMatchingDC14()", "candidate truth match:  "); 
      //  Info( "applyTruthMatchingDC14()", " \t type: %i - status %i ",         truthPart->pdgId(), truthPart->status() ); 
      //  Info( "applyTruthMatchingDC14()", " \t pT: %f eta: %f phi: %f m: %f ", truthPart->pt()/1e3 , truthPart->eta(), truthPart->phi(), truthPart->m()/1e3 );
      //  Info( "applyTruthMatchingDC14()", " \t cos(theta)*cos(theta): %f",     truthPart4mom.CosTheta()*truthPart4mom.CosTheta());   
      //}
       
      // set pt, eta, phi, m of truth 4 momentum (NB: FourMom_t it's a typedef of TLorentzVector)
      xAOD::IParticle::FourMom_t truthPart4mom = truthPart->p4();
      
      double thisDR = recoPart4mom.DeltaR( truthPart4mom );
      if ( thisDR < minDR ) { 
 
    	matchTruth = truthPart;

	// if the matching particle within thisDR is:
	//
	// 1.  a muon : no need to shrink the cone anymore! Keep first one as the match
	if ( matchTruth->isMuon() ) { break; } 
	// 2. an electron : in the next iterations, make sure you will only check matching against other electrons
	//                  ( you can have multiple truth electrons from brem, they can be close together, and you 
	//                    want to take the one with the smallest dR ) 
	if ( matchTruth->isElectron() ) { 
	  checkOnlyElectrons = true; 
	} 
     
     	minDR = thisDR; 
      } 	     
      
    } // close loop on truth particles

  } // close loop on truth events

  // if there is no matching truth lepton, just return
  if ( !matchTruth  ) {
     if ( m_debug ) { Info( "applyTruthMatchingDC14()", "No truth matching for this reco lepton"); }
     return EL::StatusCode::SUCCESS;    
  }

  // now try to find truth origin and type of the match

  if ( !( matchTruth->isElectron() || matchTruth->isMuon() ) ) {
    
    if ( m_debug ) { 
      Info( "applyTruthMatchingDC14()", "This reco lepton is truth matched neither to an electron, nor to a muon particle"); 
      Info( "applyTruthMatchingDC14()", "\t truth match type: %i - status %i ", matchTruth->pdgId(), matchTruth->status() ); 
    }
    (*m_isTruthMatchedOtherDecor)( *recoPart ) = 1;
  
  } else {
      
    // okay, this truth match is an electron/muon. Let's see where it comes from...
    
    // ...but first, let's decorate!
    (*m_isTruthMatchedDecor)( *recoPart ) = 1;  
      
    if ( m_debug ) { Info( "applyTruthMatchingDC14()", "\t match type: %i - status %i ", matchTruth->pdgId(), matchTruth->status() ); }
    
    // look at all its parents
    if ( m_debug ) { 
      for ( size_t iParent = 0; iParent < matchTruth->nParents() ; ++iParent ) {
    	 if ( m_debug ) { Info( "applyTruthMatchingDC14()", "\t parent idx %lu - type: %i - parent status %i ", iParent, matchTruth->parent(0)->pdgId(), matchTruth->parent(0)->status() ); }
      }
    }
    // from now on, we will assume each truth lepton has only one parent

    // check if it has a production vertex
    if ( matchTruth->hasProdVtx() ) {
       
      // check if it comes from a secondary (GEANT4) interaction in the material (e.g., a conversion)
      if ( matchTruth->prodVtx()->barcode() < -200000 ) {
    	 
    	 if ( m_debug ) { Info( "applyTruthMatchingDC14()", "This reco lepton is truth matched to a lepton from a secondary interaction in the detector"); }
    	 if ( m_debug ) { Info( "applyTruthMatchingDC14()", "\t parent type: %i - parent status %i ", matchTruth->parent(0)->pdgId(), matchTruth->parent(0)->status() ); }
    	 (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 1;

    	 if ( m_debug ) { Info( "applyTruthMatchingDC14()", "\t accessing isTruthMatchedSecondary %i ", (*m_isTruthMatchedSecondaryAcc)( *recoPart ) ) ; }
      
      } else {
    
    	// check the provenance
    	// keep on going backwards in the decay chain until not found :
    	// -)  a B/C (heavy flavour) hadron  ||
    	// -)  a generic light hadron (such as a pion in a jet) ||
    	// -)  the incoming parton ( status == 2 --> unstable for the generator ) or a photon from ISR/FSR  ( status == 1 --> on shell photons are stable for the generator ) : PrimaryInteraction (PI) 
    	bool foundHFHAD(false), foundLFHAD(false), foundPI(false);
	xAOD::TruthParticle* this_parent = const_cast<xAOD::TruthParticle*>(matchTruth->parent(0));
	unsigned int iGeneration(0);
	while ( !foundHFHAD && !foundLFHAD && !foundPI ) {

	   if ( !this_parent ) {
    	     if ( m_debug ) { Info( "applyTruthMatchingDC14()", "generation: %u - has no parent! Breaking", iGeneration); }
	     break;
	   }

    	   if ( m_debug ) { Info( "applyTruthMatchingDC14()", "generation: %u - truth match parent type: %i - status %i ", iGeneration, this_parent->pdgId(), this_parent->status() ); }
    	   if ( this_parent->isHeavyHadron() ) {
    	     foundHFHAD = true;
    	   }
    	   else if ( this_parent->isLightHadron() ) {
    	     foundLFHAD = true;
    	   }
    	   else if ( ( this_parent->isParton() && ( this_parent->status() == 2 || this_parent->status() == 3 ) ) || ( this_parent->isPhoton() && ( this_parent->status() == 1 ) ) ) {
    	     foundPI = true;
    	   }

	   ++iGeneration;
  	   
    	   // okay, if at the 7-th generation back in the chain we stil haven't found the intial parton (or a HF/LF hadron), let's break the loop
    	   if ( iGeneration > 6 ) {
	     if ( m_debug ) { Info( "applyTruthMatchingDC14()", "After %u generations back, we haven't reached the primary interaction vertex yet. We will flag this lepton as prompt.", iGeneration ); }
    	     foundPI = true;
	   }
	   
	   // re-assign the parent, going back of one generation in the decay chain
	   this_parent = const_cast<xAOD::TruthParticle*>(this_parent->parent(0));
	   
    	}

    	if ( foundHFHAD ) {
    	  if ( m_debug ) { Info( "applyTruthMatchingDC14()", "Truth lepton is non-prompt (HFHAD decay)"); }
    	  (*m_isTruthMatchedNonIsoDecor)( *recoPart )  = 1;
    	} else if ( foundLFHAD ) {
    	  if ( m_debug ) { Info( "applyTruthMatchingDC14()", "Truth lepton is non-prompt (generic LFHAD decay)"); }
    	  (*m_isTruthMatchedNonIsoDecor)( *recoPart )  = 1;
    	} else if ( foundPI ) {
    	  if ( m_debug ) { Info( "applyTruthMatchingDC14()", "Truth lepton is prompt"); }
    	  (*m_isTruthMatchedIsoDecor)( *recoPart )     = 1;
    	} else {
	  // should never reach this point, but let's decorate still...
          (*m_isTruthMatchedOtherDecor)( *recoPart )   = 1;
	}
    		      
      }
      
    } else {
      if ( m_debug ) { Info( "applyTruthMatchingDC14()", "Truth lepton has no production vertex... decorating with NoProdVtx "); }
      (*m_isTruthMatchedNoProdVtxDecor)( *recoPart ) = 1;
    }
    
  } // end of case: truth match is a lepton
  
  // store the pdgId and status of the match
  if ( matchTruth ) {
     if ( m_debug ) { Info( "applyTruthMatchingDC14()", "decorating truthPdgId with value : %i - truthStatus with value: %i", matchTruth->pdgId(), matchTruth->status() ); }
    (*m_truthPdgIdDecor)( *recoPart )  = matchTruth->pdgId();
    (*m_truthStatusDecor)( *recoPart ) = matchTruth->status();
    
     // store the pdgId of the parent particle of the match
     if ( matchTruth->parent(0) ) {
       if ( m_debug ) { Info( "applyTruthMatchingDC14()", "decorating truthOrigin with value: %i", matchTruth->parent(0)->pdgId() ); }
      (*m_truthOriginDecor)( *recoPart ) = matchTruth->parent(0)->pdgId();
     }
  }
     
  // check if lepton is charge flip. Do it only for leptons that are truth matched to leptons!  
  if ( ! (*m_isTruthMatchedOtherAcc).isAvailable( *recoPart ) ) {
     Error("checkChargeFlip()", "No accessor isTruthMatchedOther available for this reco lepton. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }  
  if ( ! (*m_isTruthMatchedOtherAcc)( *recoPart ) ) {
    if ( this->checkChargeFlip( recoPart, matchTruth ) != EL::StatusCode::SUCCESS ) {
      Error("applyTruthMatchingDC14()", "Problem with checkChargeFlip(). Aborting"); 
      return EL::StatusCode::FAILURE;
    }
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo ::  checkChargeFlip ( const xAOD::IParticle* recoPart, const xAOD::TruthParticle* matchTruth )
{
  
  // default decorations
  (*m_isChFlipDecor)( *recoPart )  = 0;  
  (*m_isBremDecor)( *recoPart )    = 0;  

  float reco_charge(0.0);	
  if ( recoPart->type() == xAOD::Type::Electron ) {
    reco_charge = dynamic_cast<const xAOD::Electron*>(recoPart)->charge();
  } else if ( recoPart->type() == xAOD::Type::Muon ) {
    reco_charge = dynamic_cast<const xAOD::Muon*>(recoPart)->charge();
  }
  if ( !reco_charge ) {
     Error("checkChargeFlip()", "Reco particle has zero charge. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;	    
  }

  xAOD::TruthParticle* primitiveTruth(nullptr);
  unsigned int iGeneration(0);

  if ( ! (*m_isTruthMatchedSecondaryAcc).isAvailable( *recoPart ) ) {
     Error("checkChargeFlip()", "No accessor isTruthMatchedSecondary available for this reco electron. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }

  // case 1: 
  // lepton (in most cases, an electron) is matched to a truth lepton which is part of a bremmmstrahlung shower. 
  // In this case, we need to go back until we find the original lepton that radiated the photon.
  // The charge of this primitive lepton is the one to look at! 
  if ( (*m_isTruthMatchedSecondaryAcc)( *recoPart ) ) {

    if ( m_debug ) { Info("checkChargeFlip()", "This reco lepton (charge: %f ) is matched to a secondary truth lepton (pdgId: %i , prodVtx barcode: %i). Let's go back until we find the primitive", reco_charge, matchTruth->pdgId(), matchTruth->prodVtx()->barcode() ); }

    bool foundPrimitive(false), isBrem(false);
    primitiveTruth = const_cast<xAOD::TruthParticle*>( matchTruth->parent(0) );
    
    while ( !foundPrimitive ) {
     
      if ( primitiveTruth->prodVtx()->barcode() < -200000 ) { 
	if ( m_debug ) { Info("checkChargeFlip()", "Parent has pdgId: %i , prodVtx barcode: %i - Need to go backwards in the decay chain", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	primitiveTruth = const_cast<xAOD::TruthParticle*>( primitiveTruth->parent(0) );
	// do this only once
	if ( !isBrem ) { 
	  isBrem = true;
	  (*m_isBremDecor)( *recoPart ) = 1; 
	}
      } else { 
	if ( m_debug ) { Info("checkChargeFlip()", "We found the primitive! pdgId: %i , prodVtx barcode: %i - Stop here", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	foundPrimitive = true; 
      }

      ++iGeneration;
  	   
      // okay, if at the 20-th generation back in the chain we stil haven't found the primitive lepton (i.e, the one that radiated the photon), let's break the loop
      if ( iGeneration > 19 ) {
	if ( m_debug ) { Info( "checkChargeFlip()", "After %u generations back, we haven't reached the primitive yet. Let's break the loop.", iGeneration ); }
	break;
      }

    }
  
  }
  // case 2:
  // lepton is matched to a truth lepton which is not produced in a secondary interaction (i.e., charge flip is due to charge mis-reconstruction). 
  else 
  {
    primitiveTruth = const_cast<xAOD::TruthParticle*>( matchTruth );
  }

  if ( primitiveTruth->isNeutral() ) {
    if ( m_debug ) { 
      Info("checkChargeFlip()", "primitive truth particle is neutral. PdgId : %i , Origin %i.  Returning", primitiveTruth->pdgId(), primitiveTruth->parent(0)->pdgId() ); 
    }
    return StatusCode::SUCCESS;
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

EL::StatusCode TruthMatchAlgo ::  checkChargeFlipMC15 ( const xAOD::IParticle* recoPart, const xAOD::TruthParticle* matchTruth )
{
  
  // default decorations
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
  unsigned int iGeneration(0);

  if ( ! (*m_isTruthMatchedSecondaryAcc).isAvailable( *recoPart ) ) {
     Error("checkChargeFlip()", "No accessor isTruthMatchedSecondary available for this reco electron. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }

  // case 1: 
  // lepton (in most cases, an electron) is matched to a truth lepton which is part of a bremmmstrahlung shower. 
  // In this case, we need to go back until we find the original lepton that radiated the photon.
  // The charge of this primitive lepton is the one to look at! 
  if ( (*m_isTruthMatchedSecondaryAcc)( *recoPart ) ) {

    if ( m_debug ) { Info("checkChargeFlip()", "This reco lepton (charge: %f ) is matched to a secondary truth lepton. Let's go back until we find the primitive", reco_charge ); }

    bool foundPrimitive(false), isBrem(false);
    primitiveTruth = const_cast<xAOD::TruthParticle*>( matchTruth->parent(0) );
    
    while ( !foundPrimitive ) {
     
      if ( primitiveTruth->prodVtx()->barcode() < -200000 ) { 
	if ( m_debug ) { Info("checkChargeFlip()", "Parent has pdgId: %i , prodVtx barcode: %i - Need to go backwards in the decay chain", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	primitiveTruth = const_cast<xAOD::TruthParticle*>( primitiveTruth->parent(0) );
	// do this only once
	if ( !isBrem ) { 
	  isBrem = true;
	  (*m_isBremDecor)( *recoPart ) = 1; 
	}
      } else { 
	if ( m_debug ) { Info("checkChargeFlip()", "We found the primitive! pdgId: %i , prodVtx barcode: %i - Stop here", primitiveTruth->pdgId(), primitiveTruth->prodVtx()->barcode() ); }
	foundPrimitive = true; 
      }

      ++iGeneration;
  	   
      // okay, if at the 20-th generation back in the chain we stil haven't found the primitive lepton (i.e, the one that radiated the photon), let's break the loop
      if ( iGeneration > 19 ) {
	if ( m_debug ) { Info( "checkChargeFlip()", "After %u generations back, we haven't reached the primitive yet. Let's break the loop.", iGeneration ); }
	break;
      }

    }
  
  }
  // case 2:
  // lepton is matched to a truth lepton which is not produced in a secondary interaction (i.e., charge flip is due to charge mis-reconstruction). 
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
  if ( !( recoPart->type() == xAOD::Type::Electron ) ) {
      Warning("applyTruthMatchingElectronMC15()", "Not passing an electron! Won't try anything"); 
      return StatusCode::SUCCESS; 	  	  
  }

  // truth particle types are defined in HTopMultilepAnalysis/MCTruthClassifierDefs.h:
  //
  // explaination can be found in: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/MCTruthClassifier
  // 

  // decorate reconstructed particle with default values
  //
  (*m_isTruthMatchedDecor)( *recoPart )	         = 0;
  (*m_isTruthMatchedIsoDecor)( *recoPart )       = 0;
  (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 0;
  (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 0;
  (*m_isTruthMatchedUnknownDecor)( *recoPart )   = 0;
  (*m_isTruthMatchedOtherDecor)( *recoPart )     = 0;
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
  
  // decorate with this if the truth match is an electron
  //
  if ( matchTruthEl->isElectron() ) {
     (*m_isTruthMatchedDecor)( *recoPart ) = 1;
  }
  
  if( ! (*m_truthTypeAcc).isAvailable( *recoPart ) ) {
     Error("applyTruthMatchingElectronMC15()", "No truth type info available for this electron muon. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }
  int truthType = (*m_truthTypeAcc)( *recoPart );

  if	  ( truthType == static_cast<int>(MCTruthPartClassifier::ParticleType::IsoElectron) )    {  (*m_isTruthMatchedIsoDecor)( *recoPart )	   = 1; }
  else if ( truthType == static_cast<int>(MCTruthPartClassifier::ParticleType::NonIsoElectron) ) {  (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 1; }
  else if ( truthType == static_cast<int>(MCTruthPartClassifier::ParticleType::BkgElectron) )    {  (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 1; }
  else if ( truthType == static_cast<int>(MCTruthPartClassifier::ParticleType::UnknownElectron) ){  (*m_isTruthMatchedUnknownDecor)( *recoPart )   = 1; }
  else   											 {  (*m_isTruthMatchedOtherDecor)( *recoPart )     = 1; }

  // store the pdgId of the match
  //
  if ( m_debug ) { Info( "applyTruthMatchingElectronMC15()", "decorating truthPdgId with value: %i", matchTruthEl->pdgId() ); }
  (*m_truthPdgIdDecor)( *recoPart ) = matchTruthEl->pdgId();
  // store the status of the match
  //
  if ( m_debug ) { Info( "applyTruthMatchingElectronMC15()", "decorating truthStatus with value: %i", matchTruthEl->status() ); }
  (*m_truthStatusDecor)( *recoPart ) = matchTruthEl->status();
    
  // check if lepton is charge flip  
  //
  if ( ! (*m_isTruthMatchedOtherAcc).isAvailable( *recoPart ) ) {
     Error("applyTruthMatchingElectronMC15()", "No accessor isTruthMatchedOther available for this reco lepton. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
  }  
  if ( !(*m_isTruthMatchedOtherAcc)( *recoPart ) ) {
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
  // It can be either done by DeltaR, or by track truth match (recommended)
  //
  // If using DeltaR matching, at the moment no origin and type of the matching truth muon is saved 
  //
  const xAOD::TruthParticleContainer* muonTruthPartContainer(nullptr);
  RETURN_CHECK("TruthMatchAlgo::applyTruthMatchingMuonMC15()", HelperFunctions::retrieve(muonTruthPartContainer, "MuonTruthParticles", m_event, m_store, m_debug) , "");
  
  if ( m_doMuonTrackMatching && ( this->doTrackProbMatching( recoPart ) != EL::StatusCode::SUCCESS ) ) {
     Error("applyTruthMatchingMuonMC15()", "Problem with doTrackProbMatching(). Aborting");
     return EL::StatusCode::FAILURE;
  } else if ( m_doMuonDeltaRMatching && ( this->doDeltaRMatching( muonTruthPartContainer, recoPart, 0.15) != EL::StatusCode::SUCCESS ) ) {
     Error("applyTruthMatchingMuonMC15()", "Problem with doDeltaRMatching(). Aborting");
     return EL::StatusCode::FAILURE;
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthMatchAlgo :: doTrackProbMatching( const xAOD::IParticle* recoPart )
{

   // decorate reconstructed particle with default values
   //
   (*m_isTruthMatchedDecor)( *recoPart )	  = 0;
   (*m_isTruthMatchedIsoDecor)( *recoPart )       = 0;
   (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 0;
   (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 0;
   (*m_isTruthMatchedUnknownDecor)( *recoPart )   = 0;
   (*m_isTruthMatchedOtherDecor)( *recoPart )     = 0;
   (*m_truthTypeDecor)( *recoPart )       = 0; // need it b/c for muons we need to pass from the track
   (*m_truthPdgIdDecor)( *recoPart )      = 0; 
   (*m_truthOriginDecor)( *recoPart )     = 0; // need it b/c for muons we need to pass from the track
   (*m_truthStatusDecor)( *recoPart )	  = -1;  
  
   // get the reco muon track particle 
   //
   const xAOD::TrackParticle* trk = dynamic_cast<const xAOD::Muon*>(recoPart)->primaryTrackParticle(); 
   if ( !trk ) {
     Error("doTrackProbMatching()", "No track particle available for this reco muon. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   } 

   if ( ! (*m_truthPLAcc).isAvailable( *trk ) ) {
      Error("doTrackProbMatching()", "No link available to truth match for this reco muon's track. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;		 
   }
   if ( ! (*m_truthPLAcc)( *trk ).isValid() ) {
      Error("doTrackProbMatching()", "Link to truth match for this reco muon's track is invalid. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;
   }
   const xAOD::TruthParticle* matchTruthMu = *( (*m_truthPLAcc)(*trk) );

   // if there is no matching truth track, abort 
   //
   if ( !matchTruthMu ) {
      if ( m_debug ) { Info("doTrackProbMatching()", "No truth match for this reco muon's track. This shouldn't happen. Aborting"); }
      return StatusCode::SUCCESS;    
   } 

   // retrieve track truth MC probability 
   // ( when available... this should be the case for SiliconAssociatedForwardMuon muons, whose track is a InDetTrackParticle )
   //
   float trk_prob(-1.0); 
   if ( m_truthMatchProbabilityAcc->isAvailable( *trk ) ) { 
     trk_prob = (*m_truthMatchProbabilityAcc)( *trk ); 
   }

   // decorate with this if the truth match is a muon, 
   // and the track mc probability (when available) is high enough
   // 
   if ( matchTruthMu->isMuon()  && ( trk_prob < 0.0 || trk_prob > 0.8 ) ) {
     (*m_isTruthMatchedDecor)( *recoPart )     = 1;
   } 
 
   // store the type of the parent particle of the match: pass the track type info to the reco muon
   //
   if ( ! (*m_truthTypeAcc).isAvailable( *trk ) ) {
     Error("doTrackProbMatching()", "No truth type info available for this muon's matching truth track. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthTrkMatchType = (*m_truthTypeAcc)(*trk);
   (*m_truthTypeDecor)( *recoPart ) = truthTrkMatchType;

   // store the pdgId of the match
   //
   if ( m_debug ) { Info( "doTrackProbMatching()", "decorating truthPdgId with value: %i", matchTruthMu->pdgId() ); }
   (*m_truthPdgIdDecor)( *recoPart ) = matchTruthMu->pdgId();
   
   // store the status of the match
   //
   if ( m_debug ) { Info( "doTrackProbMatching()", "decorating truthStatus with value: %i", matchTruthMu->status() ); }
   (*m_truthStatusDecor)( *recoPart ) = matchTruthMu->status();


   // store the pdgId of the parent particle of the match: pass the track origin info to the reco muon
   //
   if ( ! (*m_truthOriginAcc).isAvailable( *trk ) ) {
     Error("doTrackProbMatching()", "No truth origin info available for this muon's matching truth track. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }	    
   int truthTrkMatchOrigin = (*m_truthOriginAcc)(*trk);
   (*m_truthOriginDecor)( *recoPart ) = truthTrkMatchOrigin;

   if	   ( truthTrkMatchType == static_cast<int>(MCTruthPartClassifier::ParticleType::IsoMuon) )    {  (*m_isTruthMatchedIsoDecor)( *recoPart )	= 1; }
   else if ( truthTrkMatchType == static_cast<int>(MCTruthPartClassifier::ParticleType::NonIsoMuon) ) {  (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 1; }
   else if ( truthTrkMatchType == static_cast<int>(MCTruthPartClassifier::ParticleType::BkgMuon) )    {  (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 1; }
   else  											      {  (*m_isTruthMatchedOtherDecor)( *recoPart )     = 1; }
   
   // check if lepton is charge flip
   //  
   if ( ! (*m_isTruthMatchedOtherAcc).isAvailable( *recoPart ) ) {
     Error("doTrackProbMatching()", "No accessor isTruthMatchedOther available for this reco lepton. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }  
   if ( ! (*m_isTruthMatchedOtherAcc)(*recoPart) ) {
     if ( this->checkChargeFlipMC15( recoPart, matchTruthMu ) != EL::StatusCode::SUCCESS ) {
       Error("doTrackProbMatching()", "Problem with checkChargeFlipMC15(). Aborting"); 
       return EL::StatusCode::FAILURE;
     }
   }  
   
   return EL::StatusCode::SUCCESS;

}

EL::StatusCode TruthMatchAlgo :: doDeltaRMatching ( const xAOD::TruthParticleContainer* muonTruthPartContainer, const xAOD::IParticle* recoPart, double minDR )
{

   // decorate reconstructed particle with default values
   //
   (*m_isTruthMatchedDecor)( *recoPart )	  = 0;
   (*m_isTruthMatchedIsoDecor)( *recoPart )       = 0;
   (*m_isTruthMatchedNonIsoDecor)( *recoPart )    = 0;
   (*m_isTruthMatchedSecondaryDecor)( *recoPart ) = 0;
   (*m_isTruthMatchedOtherDecor)( *recoPart )     = 0;

   // this will be the best-matching truth muon (if ever found)
   //
   const xAOD::TruthParticle* matchTruthMu(nullptr);

   // set pt, eta, phi, m of reco 4 momentum (NB: FourMom_t it's a typedef of TLorentzVector)
   //
   xAOD::IParticle::FourMom_t recoPart4mom  = recoPart->p4();     
   ( recoPart4mom  ).SetPtEtaPhiM( recoPart->pt(), recoPart->eta(), recoPart->phi(), recoPart->m() );
   	
   for ( auto part_it : *(muonTruthPartContainer) ) {
     
     const xAOD::TruthParticle* truthMu = part_it;
     
     if( truthMu->pt() < 2e3 || fabs( truthMu->eta() ) > 2.5 ) { continue; }
             
     // set pt, eta, phi, m of truth 4 momentum (NB: FourMom_t it's a typedef of TLorentzVector)
     //
     xAOD::IParticle::FourMom_t truthMu4mom = truthMu->p4();

     double thisDR = recoPart4mom.DeltaR( truthMu4mom );
     if ( thisDR < minDR ) { 
       matchTruthMu = truthMu; 
       minDR = thisDR; 
     }	  
     
   } // close loop on truth muons

   // if there is no matching truth muon, abort
   //
   if ( !matchTruthMu ) {
      Error("doDeltaRMatching()", "No truth matching for this reco muon's track. This shouldn't happen. Aborting"); 
      return StatusCode::FAILURE;    
   }
   
   // decorate with this if the truth match is a muon
   //
   if ( matchTruthMu->isMuon() ) {
     (*m_isTruthMatchedDecor)( *recoPart ) = 1;
   } 

   // check if lepton is charge flip  
   //
   if ( ! (*m_isTruthMatchedOtherAcc).isAvailable( *recoPart ) ) {
     Error("doDeltaRMatching()", "No accessor isTruthMatchedOther available for this reco lepton. This shouldn't happen. Aborting"); 
     return StatusCode::FAILURE;		
   }  
   if ( ! (*m_isTruthMatchedOtherAcc)(*recoPart) ) {
     if ( this->checkChargeFlipMC15( recoPart, matchTruthMu ) != EL::StatusCode::SUCCESS ) {
       Error("doDeltaRMatching()", "Problem with checkChargeFlipMC15(). Aborting"); 
       return EL::StatusCode::FAILURE;
     }
   }    
   
   return EL::StatusCode::SUCCESS;

}

