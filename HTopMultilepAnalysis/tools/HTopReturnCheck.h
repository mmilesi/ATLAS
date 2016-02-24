// Dear emacs, this is -*- c++ -*-
#ifndef  HTopMultilepAnalysis_Tools_HTopReturnCheck_H
#define  HTopMultilepAnalysis_Tools_HTopReturnCheck_H
 
// ROOT include(s):
#include <TError.h>

#include <xAODAnaHelpers/tools/Message.h>

/// Helper macro for checking return codes in a compact form in the code
///
/// This is pretty much a rip-off of the (ATH_)CHECK macros of the offline
/// code. It is used in the package in functions that return a TReturnCode,
/// and themselves call functions returning TReturnCode.
///
/// @param CONTEXT A context string to print an error message on failure
/// @param EXP The expression to execute in a checked manner
/// @param INFO Extra information to print with the error message to debug
///
#define HTOP_RETURN_CHECK( CONTEXT, EXP, INFO )                            \
   do {                                                                    \
      if( ! EXP ) {                                                        \
         ::Error( CONTEXT, XAOD_MESSAGE( "Failed to execute: %s\n\t%s\n" ),\
                  #EXP, INFO );                                            \
         return EL::StatusCode::FAILURE;                                   \
      }                                                                    \
   } while( false )
   
#endif // HTopMultilepAnalysis_Tools_HTopReturnCheck_H 
