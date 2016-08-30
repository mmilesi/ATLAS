// System include(s):
#include <memory>
#include <cstdlib>

// ROOT include(s):
#include <TFile.h>
#include <TError.h>
#include <TString.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"

// external tools include(s):
#include "AsgTools/AnaToolHandle.h"
#include "PileupReweighting/PileupReweightingTool.h"

#include <vector>
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
   
   // Initialise the application:
   CHECK( xAOD::Init( APP_NAME ) );
   
   asg::AnaToolHandle<CP::IPileupReweightingTool> pileup_tool_handle_1;
   asg::AnaToolHandle<CP::IPileupReweightingTool> pileup_tool_handle_2;

   std::vector<std::string> PRWFiles = {"/afs/cern.ch/user/m/mmilesi/public/merged_prw_mc15c.root"};
   std::vector<std::string> lumiCalcFiles = {"/afs/cern.ch/user/m/mmilesi/public/ilumicalc_histograms_None_276262-284484_OflLumi-13TeV-005.root","/afs/cern.ch/user/m/mmilesi/public/ilumicalc_histograms_None_297730-303560_OflLumi-13TeV-005.root"};
   
   CHECK( pileup_tool_handle_1.make("CP::PileupReweightingTool/myPRWTool") );
   if ( pileup_tool_handle_1.isUserConfigured() ) { Info("initialize()","tool handle has been already configured by the user"); } else { Info("initialize()","tool handle has NOT been configured by the user yet"); }
   if ( pileup_tool_handle_1.isConfigurable() )   { Info("initialize()","tool handle can be configured by ASG"); } else { Info("initialize()","tool handle can NOT be configured by ASG"); }
   CHECK( pileup_tool_handle_1.setProperty("ConfigFiles", PRWFiles) );
   CHECK( pileup_tool_handle_1.setProperty("LumiCalcFiles", lumiCalcFiles) );
   CHECK( pileup_tool_handle_1.setProperty("DefaultChannel", 363372 ) );
   CHECK( pileup_tool_handle_1.setProperty("DataScaleFactor", 1.0/1.16) );
   CHECK( pileup_tool_handle_1.setProperty("DataScaleFactorUP", 1.0) );
   CHECK( pileup_tool_handle_1.setProperty("DataScaleFactorDOWN", 1.0/1.23) );
   CHECK( pileup_tool_handle_1.initialize() );

   CHECK( pileup_tool_handle_2.make("CP::PileupReweightingTool/myPRWTool") );
   if ( pileup_tool_handle_2.isUserConfigured() ) { Info("initialize()","tool handle has been already configured by the user"); } else { Info("initialize()","tool handle has NOT been configured by the user yet"); }
   if ( pileup_tool_handle_2.isConfigurable() )   { Info("initialize()","tool handle can be configured by ASG"); } else { Info("initialize()","tool handle can NOT be configured by ASG"); }
   CHECK( pileup_tool_handle_2.setProperty("ConfigFiles", PRWFiles) );
   CHECK( pileup_tool_handle_2.setProperty("LumiCalcFiles", lumiCalcFiles) );
   CHECK( pileup_tool_handle_2.setProperty("DefaultChannel", 363372 ) );
   CHECK( pileup_tool_handle_2.setProperty("DataScaleFactor", 1.0/1.16) );
   CHECK( pileup_tool_handle_2.setProperty("DataScaleFactorUP", 1.0) );
   CHECK( pileup_tool_handle_2.setProperty("DataScaleFactorDOWN", 1.0/1.23) );
   CHECK( pileup_tool_handle_2.initialize() );   
   
   // Return gracefully:
   return 0;
}
