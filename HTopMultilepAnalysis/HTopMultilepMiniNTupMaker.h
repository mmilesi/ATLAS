/**
 * @file   HTopMultilepMiniNTupMaker.h
 * @author Marco Milesi <marco.milesi@cern.ch>
 * @brief  EventLoop algorithm to skim/slim/augment HTop group ntuples
 *
 */

#ifndef HTopMultilepAnalysis_HTopMultilepMiniNTupMaker_H
#define HTopMultilepAnalysis_HTopMultilepMiniNTupMaker_H

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// EL include(s):
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>
#include <EventLoop/Worker.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Job.h>
#include <EventLoop/OutputStream.h>

// ROOT include(s):
#include "TTree.h"

class HTopMultilepMiniNTupMaker : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  
  std::string m_outputNTupStreamName;

private:
 
  /** Input TTree */

  TTree*          m_inputNTuple;
  
  /** Output TTree */
  
  EL::NTupleSvc*  m_outputNTuple;

  /** Input TTree branches */
  
  ULong64_t       m_EventNumber;
  
  Float_t	  m_lep_Pt_0;
  Float_t	  m_lep_E_0;
  Float_t	  m_lep_Eta_0;
  Float_t	  m_lep_Phi_0;
  Float_t	  m_lep_EtaBE2_0;
  Float_t	  m_lep_Pt_1;

  /** Extra branches for output TTree */

  Float_t         m_lep_Pt_0_Squared;
  Float_t         m_lep_Pt_01;


  unsigned int m_numEvent;  //!
  
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  
  // this is a standard constructor
  HTopMultilepMiniNTupMaker (std::string className = "HTopMultilepMiniNTupMaker"); 

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(HTopMultilepMiniNTupMaker, 1);
};

#endif
