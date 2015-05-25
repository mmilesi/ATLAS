//simple macro to test vertex container reading 

// System include(s):
#include <memory>
#include <cstdlib>

// ROOT include(s):
#include <TFile.h>
#include <TError.h>
#include <TString.h>

// Infrastructure include(s):
#ifdef ROOTCORE
 #include "xAODRootAccess/Init.h"
 #include "xAODRootAccess/TEvent.h"
 #include "xAODRootAccess/TStore.h"
#endif // ROOTCORE

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h" 
#include "xAODTracking/Vertex.h" 


//copied from CPAnalysisExamples/errorcheck.h
#define CHECK( ARG )                                  \
do {                                                  \
    const bool result = ARG;                          \
  if( ! result ) {                                    \
      ::Error( APP_NAME, "Failed to execute: \"%s\"", \
#ARG );                                               \
      return 1;                                       \
  }                                                   \
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
   xAOD::TEvent event( xAOD::TEvent::kClassAccess );
   CHECK( event.readFrom( ifile.get() ) );
   Info( APP_NAME, "Number of events in the file: %i",
         static_cast< int >( event.getEntries() ) );
   
   // Decide how many events to run over:
   Long64_t entries = event.getEntries();
   if( argc > 2 ) {
      const Long64_t e = atoll( argv[ 2 ] );
      if( e < entries ) {
         entries = e;
      }
   }

   // see later ...
   int Ntracks(3);
   
   // Loop over the events:
   for( Long64_t entry = 0; entry < entries; ++entry ) {
     
     // Tell the object which entry to look at:
     event.getEntry( entry );
     
     std::cout << "=================NEXT EVENT==========================" << std::endl;
     
     const xAOD::EventInfo* event_info = 0;  
     CHECK( event.retrieve( event_info, "EventInfo" ) ); 
      
     const xAOD::VertexContainer* vertexContainer;  
     CHECK( event.retrieve(vertexContainer, "PrimaryVertices") );
     
     // now do something with the vertex container!

     // get THE primary vertex
     const xAOD::Vertex* primaryVertex = 0;
     for( auto vtx_itr : *vertexContainer )
     {
	 if(vtx_itr->vertexType() != xAOD::VxType::VertexType::PriVtx) { continue; }
	 primaryVertex = vtx_itr;
     }

     // get number of tracks to PV
     if( primaryVertex ){
       std::cout << " Number of tracks associated to PV: " << static_cast<int>( primaryVertex->nTrackParticles() ) << std::endl; 
     } else {
       Error( APP_NAME, "could not retrieve primary vertex. Breaking event loop");
       break;
     }

     // count number of vertices w/ at least N tracks associated
     int nGoodVtx(0);
     // Loop over vertices in the container
     for( auto vtx_itr : *vertexContainer )
     {
        if( static_cast<int>( vtx_itr->nTrackParticles() ) < Ntracks ) { continue; }
        nGoodVtx++;
     }
     std::cout << " Number of vertices w/ at least " << Ntracks << " tracks: " << nGoodVtx << std::endl; 
     
     // find position of primary vertex in container
     if ( primaryVertex ){
       int location(0);
       for( auto vtx_itr : *vertexContainer )
       {
         if(vtx_itr->vertexType() == xAOD::VxType::VertexType::PriVtx) {
            break; 
         }
         location++;
       }
       std::cout << " Location of PV in container " << location << std::endl; 
     }

    
     Info( APP_NAME,
            "===>>>  done processing event #%i, "
            "run #%i %i events processed so far  <<<===",
            static_cast< int >( event_info->eventNumber() ),
            static_cast< int >( event_info->runNumber() ),
            static_cast< int >( entry + 1 ) );
     
   }

   // Return gracefully:
   return 0;
}
