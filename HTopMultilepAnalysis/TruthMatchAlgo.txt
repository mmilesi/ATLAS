#ifndef HTopMultilepAnalysis_TruthMatchAlgo_H
#define HTopMultilepAnalysis_TruthMatchAlgo_H

// EL include(s):
#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "AthContainers/AuxElement.h"
#include "AthLinks/ElementLink.h"

// EDM include(s):
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"

// external tools include(s):
#include "MCTruthClassifier/MCTruthClassifier.h"
#include "MCTruthClassifier/MCTruthClassifierDefs.h"

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// ROOT include(s):
#include "TH1D.h"

class TruthMatchAlgo : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  std::string    m_inContainerName_Electrons;
  std::string    m_inContainerName_Muons;
  std::string    m_inContainerName_Leptons;

  bool m_doMuonTruthPartMatching;
  bool m_doMuonTrackMatching;

private:

  int m_numEvent;           //!
  int m_numObject;          //!
  int m_numEventPass;       //!
  int m_weightNumEventPass; //!
  int m_numObjectPass;      //!

  bool m_isMC;              //!
  bool m_isDerivation;      //!

  // cutflow
  bool  m_useCutFlow;       //!
  TH1D* m_cutflowHist;      //!
  TH1D* m_cutflowHistW;     //!
  int   m_cutflow_bin;      //!

  /* Initialise decorators */
  SG::AuxElement::Decorator< char >* 	   m_isTruthMatchedDecor;	   //! /* has a lepton truth match */
  SG::AuxElement::Decorator< int >*  	   m_truthTypeDecor;		   //! /* type of the parent particle (according to MCTruthClassifier) - need it for muons since we have to retrieve this from the truth track */
  SG::AuxElement::Decorator< int >*  	   m_truthPdgIdDecor;		   //! /* pdgId of the match particle */
  SG::AuxElement::Decorator< int >*  	   m_truthOriginDecor;  	   //! /* origin of the parent particle - need it for muons since we have to retrieve this from the truth track */
  SG::AuxElement::Decorator< int >*  	   m_truthStatusDecor;  	   //! /* status of the match particle */
  SG::AuxElement::Decorator< char >* 	   m_isChFlipDecor;		   //! /* reco has opposite charge wrt to primitive truth ancestor */
  SG::AuxElement::Decorator< char >* 	   m_isBremDecor;		   //! /* reco is matched to a brem lepton (i.e, ancestor emitted a photon by brem) */
  SG::AuxElement::Decorator< int >*  	   m_ancestorTruthTypeDecor;	   //! /* type of the primitive ancestor (according to MCTruthClassifier) - need it for brem leptons*/
  SG::AuxElement::Decorator< int >*  	   m_ancestorTruthPdgIdDecor;	   //! /* pdgId of the primitive ancestor (according to MCTruthClassifier) - need it for brem leptons*/
  SG::AuxElement::Decorator< int >*  	   m_ancestorTruthOriginDecor;     //! /* origin of the primitive ancestor (according to MCTruthClassifier) - need it for brem leptons*/
  SG::AuxElement::Decorator< int >*  	   m_ancestorTruthStatusDecor;     //! /* status of the primitive ancestor (according to MCTruthClassifier) - need it for brem leptons*/

  /* Initialise accessors */
  SG::AuxElement::Accessor< float >*       m_mcEvtWeightAcc;		   //!
  SG::AuxElement::Accessor< char >*        m_isTruthMatchedAcc; 	   //!
  SG::AuxElement::Accessor< char >*        m_isChFlipAcc;		   //!
  SG::AuxElement::Accessor< char >*        m_isBremAcc; 		   //!
  typedef ElementLink< xAOD::TruthParticleContainer > TruthLink_t;
  SG::AuxElement::Accessor< TruthLink_t >* m_truthPLAcc; 	           //!
  SG::AuxElement::ConstAccessor< int >*    m_truthTypeAcc;		   //!  /* accessor to built-in xAOD attribute */
  SG::AuxElement::ConstAccessor< int >*    m_truthOriginAcc;		   //!  /* accessor to built-in xAOD attribute */
  SG::AuxElement::Accessor< float >*       m_truthMatchProbabilityAcc;     //!
  SG::AuxElement::Accessor< int >*         m_ancestorTruthTypeAcc;         //!
  SG::AuxElement::Accessor< int >*         m_ancestorTruthOriginAcc;       //!

  // MC Truth Classifier
  MCTruthClassifier *m_MCTClassifier; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  TruthMatchAlgo ();

  ~TruthMatchAlgo();

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

  // these are the functions not inherited from Algorithm

  /*
  / For muons, there are two options to perform the truth matching:
  /
  /  1) use the truth-link and get the info from the TruthParticle in the "MuonTruthParticles" container.
  /     Since this is filled only with truth muons one can only get the type and origin information for 'isolated' and 'not-isolated' truth muons.
  /
  /  2) retrieve the truth and origin information from the ID TrackParticle linked from the xAOD::Muon.
  /     This includes all types as defined by MCTruthClassifier for muons (i.e. also when the ID track is a pion).
  /     Clearly this is available only within ID coverage (|eta|<2.5).
  /
  */
  virtual EL::StatusCode applyTruthMatchingMuon ( const xAOD::IParticle* recoParticle );
  virtual EL::StatusCode doMuonTruthPartMatching ( const xAOD::IParticle* recoParticle );
  virtual EL::StatusCode doMuonTrackMatching ( const xAOD::IParticle* recoParticle );

  virtual EL::StatusCode applyTruthMatchingElectron ( const xAOD::IParticle* recoParticle );

  virtual EL::StatusCode checkChargeFlip ( const xAOD::IParticle* recoPart, const xAOD::TruthParticle* matchTruth );

  // this is needed to distribute the algorithm to the workers
  ClassDef(TruthMatchAlgo, 1);
};

#endif
