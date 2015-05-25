// System include(s):
#include <memory>
#include <cstdlib>

// ROOT include(s):
#include <TFile.h>
#include <TError.h>
#include <TString.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/Electron.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"

#include <string>

#define CHECK( ARG )				      \
do {						      \
    const bool result = ARG;			      \
  if( ! result ) {				      \
      ::Error( APP_NAME, "Failed to execute: \"%s\"", \
#ARG ); 					      \
      return 1; 				      \
  }						      \
 } while( false )


int main( int argc, char* argv[] ) {

   // The application's name:
   const char* APP_NAME = argv[ 0 ];

   // Check if we received a file name:
   if( argc < 2 ) {
     Error( APP_NAME, "No file name received!" );
     Error( APP_NAME, "  Usage: %s [xAOD file name]", APP_NAME );
      return 1;
   }
   
   // Initialise the application:
   CHECK( xAOD::Init( APP_NAME ) );
   
   // Open the input file:
   const TString fileName = argv[ 1 ];
   Info( APP_NAME, "Opening file: %s", fileName.Data() );
   std::auto_ptr< TFile > ifile( TFile::Open( fileName, "READ" ) );
   CHECK( ifile.get() );

   // Create a TEvent object:
   //xAOD::TEvent event( xAOD::TEvent::kBranchAccess );
   xAOD::TEvent event( xAOD::TEvent::kClassAccess );
   CHECK( event.readFrom( ifile.get() ) );
   Info( APP_NAME, "Number of events in the file: %i",
	 static_cast< int >( event.getEntries() ) );

   bool debug = ( argv[2]  );

   // Decide how many events to run over:
   Long64_t entries = event.getEntries();
   if( argc > 3 ) {
      const Long64_t e = atoll( argv[ 3 ] );
      if( e < entries ) {
	 entries = e;
      }
   }
   

   // Loop over the events:
   for( Long64_t entry = 0; entry < entries; ++entry ) {
       
     // Tell the object which entry to look at:
     event.getEntry( entry );
     
     if( debug ){ Info(APP_NAME, "=================NEXT EVENT==========================" ); }
 
     //if ( (entry % 10000) == 0 ) {
	Info(APP_NAME, "Processed %lld events", entry);
     //}
     
     const xAOD::EventInfo* event_info(nullptr); 
     CHECK( event.retrieve( event_info, "EventInfo" ) ); 
     
     const xAOD::ElectronContainer* electrons(nullptr); 
     CHECK( event.retrieve(electrons, "ElectronCollection") );

     const xAOD::TruthEventContainer* truthEventContainer(nullptr);
     CHECK( event.retrieve(truthEventContainer, "TruthEvent") );


     // declare some particle decorations
     static SG::AuxElement::Decorator< char > isTruthMatchedDecor("isTruthMatched");
     static SG::AuxElement::Decorator< char > isTruthMatchedIsoDecor("isTruthMatchedIso");
     static SG::AuxElement::Decorator< char > isTruthMatchedNonIsoDecor("isTruthMatchedNonIso");
     static SG::AuxElement::Decorator< char > isTruthMatchedBkgDecor("isTruthMatchedBkg");
     static SG::AuxElement::Decorator< char > isTruthMatchedOtherDecor("isTruthMatchedOther");
     static SG::AuxElement::Decorator< char > isTruthMatchedChFlipDecor("isTruthMatchedChFlip");
     static SG::AuxElement::Decorator< int >  truthMatchedOriginDecor("truthMatchedOrigin");

     if ( debug ) { Info( APP_NAME, "************************* \n number of reco electrons: %lu \n *************************", electrons->size() ); }
     
     //Iterate over the electrons
     for ( auto el : *(electrons) ) { 
     	  
     	  // decorate reconstructed particle with default values
     	  isTruthMatchedDecor( *el )	    = 0;
     	  isTruthMatchedIsoDecor( *el )	    = 0;
     	  isTruthMatchedNonIsoDecor( *el )  = 0;
     	  isTruthMatchedBkgDecor( *el )	    = 0;
     	  isTruthMatchedOtherDecor( *el )   = 0;
     	  isTruthMatchedChFlipDecor( *el )  = 0;
     	  truthMatchedOriginDecor( *el )    = 0; // i.e., NonDefined  
	  
          // this will be the best-matching truth particle (if ever found)
          const xAOD::TruthParticle* matchTruth(nullptr);

          // set pt, eta, phi, m of reco 4 momentum (NB: FourMom_t it's a typedef of TLorentzVector)
	  if ( debug ) { Info( APP_NAME, "reco electron: " ); }
	  if ( debug ) { Info( APP_NAME, "pT: %f eta: %f phi: %f m: %f ", el->pt()/1e3 , el->eta(), el->phi(), el->m()/1e3 ); }
          xAOD::IParticle::FourMom_t el4mom(1., 1., 1., 1.);	
          ( el4mom  ).SetPtEtaPhiM( el->pt(), el->eta(), el->phi(), el->m() );
          
	  // loop over truth events
          for ( auto evt_it : *(truthEventContainer) ) {

	    double minDR(0.2);
	    
	    int nPart = evt_it->nTruthParticles();
            // loop over truth particles in this truth event
            for ( int iTPart = 0 ; iTPart < nPart ; ++iTPart ) {
          
              const xAOD::TruthParticle* truthPart = evt_it->truthParticle(iTPart);
            
	      if ( !truthPart ) { 
		continue;
	      }
	      
              // for the moment, let's do a simple deltaR match
               
              // set pt, eta, phi, m of truth 4 momentum (NB: FourMom_t it's a typedef of TLorentzVector)
	      
              xAOD::IParticle::FourMom_t truthPart4mom(1., 1., 1., 1.);
              ( truthPart4mom ).SetPtEtaPhiM( truthPart->pt(), truthPart->eta(), truthPart->phi(), truthPart->m() ); 
	      //if ( debug ) { Info( APP_NAME, "\t candidate truth match: "); }
	      //if ( debug ) { Info( APP_NAME, "\t pT: %f eta: %f phi: %f m: %f ", truthPart->pt()/1e3 , truthPart->eta(), truthPart->phi(), truthPart->m()/1e3 ); }
	      double thisDR = el4mom.DeltaR( truthPart4mom );
              if ( thisDR < minDR ) { 
                matchTruth = truthPart; 
                minDR = thisDR; 
              } 	 

            } // close loop on truth particles

          } // close loop on truth events

          // if there is no matching truth muon, just return
          if ( !matchTruth  ) {
             if ( debug ) { Info( APP_NAME, "No truth matching for this reco electron"); }
             continue;    
          }
	  
	  // /*
	  // now try to find truth origin and type of the match
	  
	  if ( ! matchTruth->isElectron() ) {
	    if ( debug ) { Info( APP_NAME, "This reco electron is truth matched to a NON-ELECTRON particle"); }
	    if ( debug ) { Info( APP_NAME, "\t match type: %i - status %i ", matchTruth->pdgId(), matchTruth->status() ); }
	    isTruthMatchedOtherDecor( *el ) = 1;
	  } else {
	    
	    // okay, this truth match is an electron. Let's see where it comes from...
	    
	    if ( debug ) { Info( APP_NAME, "\t match type: %i - status %i ", matchTruth->pdgId(), matchTruth->status() ); }
	    
	    // look at all its parents
	    if ( debug ) { 
	      for ( size_t iParent = 0; iParent < matchTruth->nParents() ; ++iParent ) {
	         if ( debug ) { Info( APP_NAME, "\t parent idx %lu - type: %i - parent status %i ", iParent, matchTruth->parent(0)->pdgId(), matchTruth->parent(0)->status() ); }
	      }
	    }
	    
	    // check if it has a production vertex
	    if ( matchTruth->hasProdVtx() ) {
	       
	      // check if it comes from a secondary (GEANT4) interaction
	      if ( matchTruth->prodVtx()->barcode() < -200000 ) {
	         
		 if ( debug ) { Info( APP_NAME, "This reco electron is truth matched to an electron from secondaries"); }
	         if ( debug ) { Info( APP_NAME, "\t parent type: %i - parent status %i ", matchTruth->parent(0)->pdgId(), matchTruth->parent(0)->status() ); }
	         isTruthMatchedBkgDecor( *el ) = 1;
	      
	      } else {
	        
		// check the provenance
		
		// keep on going backwards in the decay chain until not found :
		// -)  a B/C hadron  ||
		// -)  a generic hadron (like a pion in a jet) ||
		// -)  the incoming parton ( status == 2 --> unstable for the generator ) or a photon from ISR/FSR  ( status == 1 --> on shell photons are stable for the generator ) : PrimaryInteraction(PI) 
		bool foundHFHAD(false), foundHAD(false), foundPI(false);
		while ( !foundHFHAD && !foundHAD && !foundPI ) {
		
		   const xAOD::TruthParticle* parent = matchTruth->parent(0);
		   if ( debug ) { Info( APP_NAME, "truth match parent type: %i - status %i ", parent->pdgId(), parent->status() ); }
		   if ( parent->isHeavyHadron() ) {
		     foundHFHAD = true;
		     continue;
		   }
		   if ( parent->isHadron() ) {
		     foundHAD = true;
		     continue;
		   }
		   if ( ( parent->isParton() && ( parent->status() == 2 || parent->status() == 3 ) ) || ( parent->isPhoton() && ( parent->status() == 1 ) ) ) {
		     foundPI = true;
		     continue;
		   }
		
		   const xAOD::TruthParticle* grandparent = parent->parent(0);
		   if ( debug ) { Info( APP_NAME, "truth match grand parent type: %i - status %i ", grandparent->pdgId(), grandparent->status() ); }
		   if ( grandparent->isHeavyHadron() ) {
		     foundHFHAD = true;
		     continue;
		   }
		   if ( grandparent->isHadron() ) {
		     foundHAD = true;
		     continue;
		   }
		   if ( ( grandparent->isParton() && ( grandparent->status() == 2 || grandparent->status() == 3 ) ) || ( grandparent->isPhoton() && ( grandparent->status() == 1 ) ) ) {
		     foundPI = true;
		     continue;
		   }
		   
		   const xAOD::TruthParticle* grandgrandparent = grandparent->parent(0);
		   if ( debug ) { Info( APP_NAME, "truth match grand grand parent type: %i - status %i ", grandgrandparent->pdgId(), grandgrandparent->status() ); }
		   if ( grandgrandparent->isHeavyHadron() ) {
		     foundHFHAD = true;
		     continue;
		   }
		   if ( grandgrandparent->isHadron() ) {
		     foundHAD = true;
		     continue;
		   }
		   if ( ( grandgrandparent->isParton() && ( grandgrandparent->status() == 2 || grandgrandparent->status() == 3 ) ) || ( grandgrandparent->isPhoton() && ( grandgrandparent->status() == 1 ) ) ) {
		     foundPI = true;
		     continue;
		   }
		   
		   // okay, if at this stage we stil haven't found the intial parton, let's break the loop
		   if ( debug ) { Info( APP_NAME, "We haven't reached the primary interaction vertex yet. We will flag this lepton as prompt." ); }
		   foundPI = true;
		}
		
		if ( foundHFHAD ) {
		  if ( debug ) { Info( APP_NAME, "Truth electron is non-prompt (HFHAD decay)"); }
		  isTruthMatchedNonIsoDecor( *el )  = 1;
		} else if ( foundHAD ) {
		  if ( debug ) { Info( APP_NAME, "Truth electron is non-prompt (generic HAD decay)"); }
		  isTruthMatchedNonIsoDecor( *el )  = 1;
		} else if ( foundPI ) {
		  if ( debug ) { Info( APP_NAME, "Truth electron is prompt"); }
		  isTruthMatchedIsoDecor( *el )  = 1;
		}
			      
	      }
	      
	    } else {
	      if ( debug ) { Info( APP_NAME, "Truth electron has no production vertex..."); }
	    }
	    
	  }
	  
	  // */

     } // close loop over electrons

     if( (entry % 10000) == 0 ){
       Info( APP_NAME, "===>>>  done processing event #%lld ",entry);
     }

   } // close loop over events
   
   // Return gracefully:
   return 0;
}











