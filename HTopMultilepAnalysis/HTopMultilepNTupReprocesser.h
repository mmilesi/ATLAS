/**
 * @file   HTopMultilepNTupReprocesser.h
 * @author Marco Milesi <marco.milesi@cern.ch>
 * @brief  EventLoop algorithm to update information of an existing TTree. Used to augment trees w/ new weights (QMisID, Matrix Method)
 *
 */

#ifndef HTopMultilepAnalysis_HTopMultilepNTupReprocesser_H
#define HTopMultilepAnalysis_HTopMultilepNTupReprocesser_H

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// EL include(s):
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/Worker.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Job.h>
#include <EventLoop/OutputStream.h>

// ROOT include(s):
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"

namespace NTupReprocesser {

  class eventObj {

  public:
    eventObj():
    	  isMC(0),
    	  isSS01(0),
    	  dilep(0),
  	  TT(0),TAntiT(0),AntiTT(0),AntiTAntiT(0)
    { };

    char isMC;
    char isSS01;
    char dilep;
    char TT;
    char TAntiT;
    char AntiTT;
    char AntiTAntiT;
  
  };

  class leptonObj {

  public:
    leptonObj():
    	  pt(-1.0),
  	  eta(-999.0),
  	  etaBE2(-999.0),
  	  ID(0),flavour(0),
  	  charge(-999.0),
  	  tightselected(0),
	  trigmatched(0)
    { };

    float pt;
    float eta;
    float etaBE2;
    int ID;
    int flavour;
    float charge;
    char tightselected;
    char trigmatched;
  };

}

class HTopMultilepNTupReprocesser : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  /** The name of the output TTree */

  std::string m_outputNTupName;

  std::string m_outputNTupStreamName;

  /** A comma-separated list of input branches to be activated */

  std::string m_inputBranches;

  /** The weight to be computed and stored/updated in the ntuple. Can be one of {"MM", "QMisID"} */

  std::string m_weightToCalc;

  /** The path to the real/fake efficiency file directory */

  std::string m_REFF_dir;
  std::string m_FEFF_dir;
  std::string m_EFF_YES_TM_dir;
  std::string m_EFF_NO_TM_dir;
  std::string m_Efficiency_Filename;
  bool m_doMMClosure;
  bool m_useEtaParametrisation;

  /** Read different r/f rates depending on whether the lepton is trigger-matched or not */

  bool m_useTrigMatchingInfo;

  bool m_useTEfficiency;
  
  /** A list of systematics affecting the efficiency measurement, whoich will be eventually propagated to the final event weight.
      By default, it includes the statistical uncertainty on the efficiencies
  */
  
  std::string              m_systematics_list;
  std::vector<std::string> m_systematics;

  /** Use the QMisID-eff-scaled-real-efficiency as fake efficiency for electrons when running MM on DATA */

  bool m_useScaledFakeEfficiency;

  /** The path to the QMisID rates file directory */

  std::string m_QMisIDRates_dir;
  std::string m_QMisIDRates_Filename_T;
  std::string m_QMisIDRates_Filename_AntiT;
  bool m_useTAntiTRates;

private:

  /** Input TTree */

  TTree*          m_inputNTuple;

  /** Output TTree (svc) */

  EL::NTupleSvc*  m_outputNTuple;

  /** Input TTree branches which need to be used by the algorithm */

  ULong64_t       m_EventNumber;
  UInt_t          m_RunNumber;
  UInt_t          m_mc_channel_number; /** for DATA, mc_channel_number=0 */
  Int_t 	  m_dilep_type;
  Int_t           m_trilep_type;
  Char_t          m_isSS01;
  Char_t          m_is_T_T;
  Char_t          m_is_T_AntiT;
  Char_t          m_is_AntiT_T;
  Char_t          m_is_AntiT_AntiT;

  Float_t	  m_lep_ID_0;
  Float_t	  m_lep_Pt_0;
  Float_t	  m_lep_E_0;
  Float_t	  m_lep_Eta_0;
  Float_t	  m_lep_Phi_0;
  Float_t	  m_lep_EtaBE2_0;
  Char_t	  m_lep_isTightSelected_0;
  Char_t	  m_lep_isTrigMatch_0;

  Float_t	  m_lep_ID_1;
  Float_t	  m_lep_Pt_1;
  Float_t	  m_lep_E_1;
  Float_t	  m_lep_Eta_1;
  Float_t	  m_lep_Phi_1;
  Float_t	  m_lep_EtaBE2_1;
  Char_t	  m_lep_isTightSelected_1;
  Char_t	  m_lep_isTrigMatch_1;

  Float_t	  m_QMisIDWeight_NOMINAL_in;
  Float_t	  m_QMisIDWeight_UP_in;
  Float_t	  m_QMisIDWeight_DN_in;

  Float_t	  m_MMWeight_NOMINAL_in;
  
  /** Map containing input branches with the variations of MM weights for each systematic */
  
  std::map<std::string, std::vector<float> > m_MMWeight_in; 

  /** Value of these flags will be inferred from input tree content.
      If false (aka branch does not exist yet), will ADD new corresponding branch to output tree, otherwise (aka branch already exists) will UPDATE it
  */

  bool m_isQMisIDBranchIn;
  bool m_isMMBranchIn;

  /** Add/update weight or not */

  bool m_doQMisIDWeighting; //!
  bool m_doMMWeighting;     //!

  /** Extra branches to be stored in output TTree */

  Float_t	  m_QMisIDWeight_NOMINAL_out = 1.0;
  Float_t	  m_QMisIDWeight_UP_out = 1.0;
  Float_t	  m_QMisIDWeight_DN_out = 1.0;

  /** Output branch with nominal MM weight */

  Float_t	  m_MMWeight_NOMINAL_out = 1.0;

  /** Map containing output branches with the variations of MM weights for each systematic */

  std::map<std::string, std::vector<float> > m_MMWeight_out; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

  /** Other private members */

  unsigned int m_numEntry;   //!
  unsigned int m_count_inf;  //!
  
  /** This will be updated when looping over the input systematics for a given event, so all the methods know about it */
  
  std::string m_this_syst;   //!

  std::shared_ptr<NTupReprocesser::eventObj>                 m_event;   //!
  std::vector< std::shared_ptr<NTupReprocesser::leptonObj> > m_leptons; //!


  std::map< std::string, TH2D* > m_QMisID_hist_map; //!

  std::map< std::string, TH1D* > m_el_hist_map; //!
  std::map< std::string, TH1D* > m_mu_hist_map; //!

  std::map< std::string, TEfficiency* > m_el_teff_map; //!
  std::map< std::string, TEfficiency* > m_mu_teff_map; //!

  /** Normalisation factors for r/f efficiencies */

  double m_el_reff_tot = 1.0;
  double m_el_feff_tot = 1.0;
  double m_mu_reff_tot = 1.0;
  double m_mu_feff_tot = 1.0;

public:

  // this is a standard constructor
  HTopMultilepNTupReprocesser (std::string className = "HTopMultilepNTupReprocesser");

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
  ClassDef(HTopMultilepNTupReprocesser, 1);

  /** Template function to get a generic TObject from a TFile */

  template<typename T>
  T* get_object( TFile& file, const std::string& name ) {
    T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
    if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
    return obj;
  }

private:

  EL::StatusCode tokenize ( char separator, std::vector<std::string>& vec_tokens, const std::string& list ); 

  EL::StatusCode readRFEfficiencies ();
  EL::StatusCode getMMEfficiencyAndError ( std::shared_ptr<NTupReprocesser::leptonObj> lep,
  					   std::vector<float>& efficiency,
					   const std::string& type );
  EL::StatusCode getMMWeightAndError ( std::vector<float>& mm_weight,
  				       const std::vector<float>& r0, const std::vector<float>& r1,
				       const std::vector<float>& f0, const std::vector<float>& f1 );
  float matrix_equation ( const float& r0, const float& r1, const float& f0, const float& f1 );
  EL::StatusCode calculateMMWeights ();

  EL::StatusCode readQMisIDRates ();
  EL::StatusCode getQMisIDRatesAndError ( std::shared_ptr<NTupReprocesser::leptonObj> lep,
					  float& r, float& r_up, float& r_dn,
					  const std::string& selection );
  EL::StatusCode calculateQMisIDWeights ();

  EL::StatusCode enableSelectedBranches ();
  EL::StatusCode clearBranches ();
  EL::StatusCode resetDefaultWeights ();


};

#endif
