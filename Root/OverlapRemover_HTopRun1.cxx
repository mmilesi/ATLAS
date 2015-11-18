/*************************************************
 *
 * Interface to HTop overlap removal tool 
 *
 * This shares everything with the base xAH::OverlapRemover
 * interface, except for the initailize() method, which 
 * instantiate the HTop OLR tool 
 * (and copies it over to the base xAH tool instance)
 *
 * M. Milesi (marco.milesi@cern.ch)
 *
 ************************************************/

// c++ include(s):
#include <iostream>
#include <sstream>

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

// package include(s):

#include "HTopMultilepAnalysis/OverlapRemovalTool_HTopRun1.h"
#include "HTopMultilepAnalysis/OverlapRemover_HTopRun1.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include <xAODAnaHelpers/tools/ReturnCheck.h>

// ROOT include(s):
#include "TEnv.h"
#include "TSystem.h"

// this is needed to distribute the algorithm to the workers
ClassImp(OverlapRemover_HTopRun1)

EL::StatusCode OverlapRemover_HTopRun1 :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info("initialize()", "Initializing OverlapRemover_HTopRun1 Interface... ");

  if ( setCutFlowHist() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to setup cutflow histograms. Exiting." );
    return EL::StatusCode::FAILURE;
  }  

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  Info("initialize()", "Number of events in file: %lld ", m_event->getEntries() );

  if ( configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if ( setCounters() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly set event/object counters. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  // initialize HTop overlap removal tool
  //
  OverlapRemovalTool_HTopRun1* OLRTool = new OverlapRemovalTool_HTopRun1( "OverlapRemovalTool" );
  OLRTool->msg().setLevel( MSG::INFO ); // VERBOSE, INFO, DEBUG

  // set input object "selection" decoration
  //
  const std::string selected_label = ( m_useSelected ) ? "passSel" : "";  // set with decoration flag you use for selected objects if want to consider only selected objects in OR, otherwise it will perform OR on all objects
  RETURN_CHECK( "OverlapRemover_HTopRun1::initialize()", OLRTool->setProperty("InputLabel",  selected_label), "");
  RETURN_CHECK( "OverlapRemover_HTopRun1::initialize()", OLRTool->setProperty("OverlapLabel", "overlaps"), "Failed to set property OverlapLabel"); // tool will decorate objects with 'overlaps' boolean if they overlap 
  RETURN_CHECK( "OverlapRemover_HTopRun1::initialize()", OLRTool->initialize(), "Failed to properly initialize the OverlapRemovalTool.");

  // now copy HTop tool over to xAH member!
  //
  m_overlapRemovalTool = OLRTool;
  
  Info("initialize()", "OverlapRemover_HTopRun1 Interface succesfully initialized!" );

  return EL::StatusCode::SUCCESS;
}
