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
  /** A comma-separated list of input branches to be activated */
  std::string m_inputBranches;

private:
 
  /** Input TTree */

  TTree*          m_inputNTuple;
  
  /** Output TTree */
  
  EL::NTupleSvc*  m_outputNTuple;

  /** Input TTree branches whcih need to be used by the algorithm */
  
  ULong64_t       m_EventNumber;
  UInt_t          m_RunNumber;

  Int_t 	  m_dilep_type;
  Int_t 	  m_trilep_type;
  
  Float_t	  m_lep_ID_0;
  Float_t	  m_lep_Pt_0;
  Float_t	  m_lep_E_0;
  Float_t	  m_lep_Eta_0;
  Float_t	  m_lep_Phi_0;
  Float_t	  m_lep_EtaBE2_0;
  Float_t	  m_lep_sigd0PV_0;
  Float_t	  m_lep_Z0SinTheta_0;
  Char_t	  m_lep_isTightLH_0;
  Char_t	  m_lep_isMediumLH_0;
  Char_t	  m_lep_isLooseLH_0;
  Char_t	  m_lep_isTight_0;
  Char_t	  m_lep_isMedium_0;
  Char_t	  m_lep_isLoose_0;
  Int_t 	  m_lep_isolationLooseTrackOnly_0;
  Int_t 	  m_lep_isolationLoose_0;
  Int_t 	  m_lep_isolationFixedCutTight_0;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_0;
  Int_t 	  m_lep_isolationFixedCutLoose_0;
  Char_t	  m_lep_isTrigMatch_0;
  Char_t	  m_lep_isPrompt_0;
  Char_t	  m_lep_isBremsElec_0;
  Char_t	  m_lep_isFakeLep_0;
  Float_t	  m_lep_SFIDLoose_0;
  Float_t	  m_lep_SFIDTight_0;
  Float_t	  m_lep_SFTrigLoose_0;
  Float_t	  m_lep_SFTrigTight_0;
  Float_t	  m_lep_SFIsoLoose_0;
  Float_t	  m_lep_SFIsoTight_0;
  Float_t	  m_lep_SFReco_0;
  Float_t	  m_lep_SFTTVA_0;
  Float_t	  m_lep_SFObjLoose_0;
  Float_t	  m_lep_SFObjTight_0;
  
  Float_t	  m_lep_ID_1;
  Float_t	  m_lep_Pt_1;
  Float_t	  m_lep_E_1;
  Float_t	  m_lep_Eta_1;
  Float_t	  m_lep_Phi_1;
  Float_t	  m_lep_EtaBE2_1;
  Float_t	  m_lep_sigd0PV_1;
  Float_t	  m_lep_Z0SinTheta_1;
  Char_t	  m_lep_isTightLH_1;
  Char_t	  m_lep_isMediumLH_1;
  Char_t	  m_lep_isLooseLH_1;
  Char_t	  m_lep_isTight_1;
  Char_t	  m_lep_isMedium_1;
  Char_t	  m_lep_isLoose_1;
  Int_t 	  m_lep_isolationLooseTrackOnly_1;
  Int_t 	  m_lep_isolationLoose_1;
  Int_t 	  m_lep_isolationFixedCutTight_1;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_1;
  Int_t 	  m_lep_isolationFixedCutLoose_1;
  Char_t	  m_lep_isTrigMatch_1;
  Char_t	  m_lep_isPrompt_1;
  Char_t	  m_lep_isBremsElec_1;
  Char_t	  m_lep_isFakeLep_1;
  Float_t	  m_lep_SFIDLoose_1;
  Float_t	  m_lep_SFIDTight_1;
  Float_t	  m_lep_SFTrigLoose_1;
  Float_t	  m_lep_SFTrigTight_1;
  Float_t	  m_lep_SFIsoLoose_1;
  Float_t	  m_lep_SFIsoTight_1;
  Float_t	  m_lep_SFReco_1;
  Float_t	  m_lep_SFTTVA_1;
  Float_t	  m_lep_SFObjLoose_1;
  Float_t	  m_lep_SFObjTight_1;

  Float_t	  m_lep_ID_2;
  Float_t	  m_lep_Pt_2;
  Float_t	  m_lep_E_2;
  Float_t	  m_lep_Eta_2;
  Float_t	  m_lep_Phi_2;
  Float_t	  m_lep_EtaBE2_2;
  Float_t	  m_lep_sigd0PV_2;
  Float_t	  m_lep_Z0SinTheta_2;
  Char_t	  m_lep_isTightLH_2;
  Char_t	  m_lep_isMediumLH_2;
  Char_t	  m_lep_isLooseLH_2;
  Char_t	  m_lep_isTight_2;
  Char_t	  m_lep_isMedium_2;
  Char_t	  m_lep_isLoose_2;
  Int_t 	  m_lep_isolationLooseTrackOnly_2;
  Int_t 	  m_lep_isolationLoose_2;
  Int_t 	  m_lep_isolationFixedCutTight_2;
  Int_t 	  m_lep_isolationFixedCutTightTrackOnly_2;
  Int_t 	  m_lep_isolationFixedCutLoose_2;
  Char_t	  m_lep_isTrigMatch_2;
  Char_t	  m_lep_isPrompt_2;
  Char_t	  m_lep_isBremsElec_2;
  Char_t	  m_lep_isFakeLep_2;
  Float_t	  m_lep_SFIDLoose_2;
  Float_t	  m_lep_SFIDTight_2;
  Float_t	  m_lep_SFTrigLoose_2;
  Float_t	  m_lep_SFTrigTight_2;
  Float_t	  m_lep_SFIsoLoose_2;
  Float_t	  m_lep_SFIsoTight_2;
  Float_t	  m_lep_SFReco_2;
  Float_t	  m_lep_SFTTVA_2;
  Float_t	  m_lep_SFObjLoose_2;
  Float_t	  m_lep_SFObjTight_2;  
  
  /** Extra branches to be stored in output TTree */
  
  char      isMC;

  char	    m_isSS01;
  char	    m_isSS12;
  
  char	    m_is_T_T;
  char	    m_is_T_AntiT;
  char	    m_is_AntiT_T;
  char	    m_is_AntiT_AntiT;
  
  int       m_nmuons;
  int       m_nelectrons;
  int       m_nleptons;
  
  char	    m_lep_isTightSelected_0;
  char	    m_lep_isTightSelected_1;
  char	    m_lep_isTightSelected_2;
  
  float     m_lep_Tag_Pt;
  float     m_lep_Tag_Eta;
  float     m_lep_Tag_EtaBE2;
  float     m_lep_Tag_sigd0PV;
  float     m_lep_Tag_Z0SinTheta;
  float	    m_lep_Tag_ID;
  char	    m_lep_Tag_isTrigMatch;
  char	    m_lep_Tag_isTightSelected;
    
  float     m_lep_Probe_Pt;
  float     m_lep_Probe_Eta;
  float     m_lep_Probe_EtaBE2;
  float     m_lep_Probe_sigd0PV;
  float     m_lep_Probe_Z0SinTheta;
  float	    m_lep_Probe_ID;
  char	    m_lep_Probe_isTrigMatch;
  char	    m_lep_Probe_isTightSelected;
    
  /** Other private members */
  
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
  
private:

  EL::StatusCode enableSelectedBranches ();
  
};

#endif
