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
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/Egamma.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/Muon.h"
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
 
     if ( (entry % 10000) == 0 ) {
	Info(APP_NAME, "Processed %lld events", entry);
     }
     
     const xAOD::EventInfo* event_info = 0; 
     CHECK( event.retrieve( event_info, "EventInfo" ) ); 
     
     const xAOD::ElectronContainer* electrons = 0; 
     CHECK( event.retrieve(electrons, "ElectronCollection") );


     //Iterate over the electrons
     for ( auto el : *(electrons) ) { 

          // declare decorator
    	  static SG::AuxElement::Decorator< char > isTruthMatchedDecor("isTruthMatched");
    	  // decorate  with default value
    	  isTruthMatchedDecor( *el ) = 0;

    	  //
    	  // Now try to do the matching
    	  //
    	  // For electrons, the link to truth matching particle (and some useful info) is already saved in ElectronCollection
    	  //
    	  
    	  typedef ElementLink< xAOD::TruthParticleContainer > TruthLink_t; 
    	  static SG::AuxElement::Accessor< TruthLink_t > truthPLAcc("truthParticleLink");
    	  if ( ! truthPLAcc.isAvailable(*el) ) {
    	     Error(APP_NAME, "No link available to truth match for this reco electron. This shouldn't happen. Aborting"); 
    	     continue;		
    	  }
    	  if ( ! truthPLAcc(*el).isValid() ) {
    	     Error(APP_NAME, "Link to truth match for this reco electron is invalid. This shouldn't happen. Aborting"); 
    	     continue;		
    	  }
    	  const xAOD::TruthParticle* matchTruthEl = *( truthPLAcc(*el) );

    	  // if there is no matching truth electron, just return
    	  if ( !matchTruthEl ) {
    	     if ( debug ) { Info(APP_NAME, "No truth matching for this reco electron"); }
    	     continue;    
    	  } 
    	  
    	  // decorate reco object with matching  
    	  isTruthMatchedDecor( *el ) = 1;
	  
	  if ( debug ) { Info(APP_NAME, "reco electron has truth match!"); }
    	  
	  // accessors to truth properties
    	  static SG::AuxElement::ConstAccessor< int > truthTypeAcc("truthType");
    	  static SG::AuxElement::ConstAccessor< int > truthOriginAcc("truthOrigin");
    	   
    	  if( ! truthTypeAcc.isAvailable(*el) ) {
    	     Error(APP_NAME, "No truth type info available for this electron muon. This shouldn't happen. Aborting"); 
    	     continue;		
    	  }
    	  int truthMatchType = truthTypeAcc(*el);
    	  
    	  if( ! truthOriginAcc.isAvailable(*el) ) {
    	     Error(APP_NAME, "No truth origin info available for this truth electron. This shouldn't happen. Aborting"); 
    	     continue;		
    	  }	    
    	  int truthMatchOrigin = truthOriginAcc(*el);
	  
	  if ( debug ) { Info(APP_NAME, "truth match: \n type = %i \t origin %i !", truthMatchType, truthMatchOrigin); }

     }

     if( (entry % 10000) == 0 ){
       Info( APP_NAME, "===>>>  done processing event #%lld ",entry);
     }

   }
   
   // Return gracefully:
   return 0;
}
