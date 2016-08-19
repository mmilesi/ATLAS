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
#include "TH1F.h"
#include "TH1D.h"

class eventObj {

public:
  eventObj():
        isSS01(0),
        dilep(0),
	TT(0),TL(0),LT(0),LL(0),
	weight_QMisID{1.0,0.0,0.0},
	weight_MM{1.0,0.0,0.0,0.0,0.0}
  { };

  char isSS01;
  char dilep;
  char TT;
  char TL;
  char LT;
  char LL;

  std::vector<float> weight_QMisID;
  std::vector<float> weight_MM;

};

class leptonObj {

public:
  leptonObj():
        pt(-1.0),
	eta(-999.0),
	etaBE2(-999.0),
	ID(0),flavour(0),
	charge(-999.0),
	tightselected(0)
  { };

  float pt;
  float eta;
  float etaBE2;
  int ID;
  int flavour;
  float charge;
  char tightselected;
};


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
  std::string m_RR_dir;
  std::string m_FR_dir;
  bool m_doMMClosure;
  bool m_useEtaParemetrisation;

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
  Char_t          m_isSS01;

  Float_t	  m_lep_ID_0;
  Float_t	  m_lep_Pt_0;
  Float_t	  m_lep_E_0;
  Float_t	  m_lep_Eta_0;
  Float_t	  m_lep_Phi_0;
  Float_t	  m_lep_EtaBE2_0;
  Char_t	  m_lep_isTightSelected_0;

  Float_t	  m_lep_ID_1;
  Float_t	  m_lep_Pt_1;
  Float_t	  m_lep_E_1;
  Float_t	  m_lep_Eta_1;
  Float_t	  m_lep_Phi_1;
  Float_t	  m_lep_EtaBE2_1;
  Char_t	  m_lep_isTightSelected_1;

  std::vector<float> *m_QMisIDWeight_in = nullptr; //!
  std::vector<float> *m_MMWeight_in     = nullptr; //!

  /** Add/update weight or not */

  bool m_doQMisIDWeighting; //!
  bool m_doMMWeighting;     //!

  /** Value of these flags will be inferred from input tree content.
      If false (aka branch does not exist yet), will ADD new corresponding branch to output tree, otherwise (aka branch already exists) will UPDATE it
  */
  bool m_isQMisIDBranchIn;
  bool m_isMMBranchIn;

  /** Extra branches to be stored in output TTree */

  std::vector<float> m_QMisIDWeight_out;
  std::vector<float> m_MMWeight_out;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

  /** Other private members */

  std::shared_ptr<eventObj>                 m_event;   //!
  std::vector< std::shared_ptr<leptonObj> > m_leptons; //!

  std::map< std::string, TH2D* > m_QMisID_hist_map; //!

  std::map< std::string, TH1D* > m_el_hist_map; //!
  std::map< std::string, TH1D* > m_mu_hist_map; //!

  /** Number of bins of the r/f efficiency input histograms */
  int  m_n_el_bins_pt_rr;
  int  m_n_el_bins_pt_fr;
  int  m_n_mu_bins_pt_rr;
  int  m_n_mu_bins_pt_fr;
  int  m_n_el_bins_eta;
  int  m_n_mu_bins_eta;

  /** Normalisation factors for r/f efficiencies */
  double m_el_rr_tot;
  double m_el_fr_tot;
  double m_mu_rr_tot;
  double m_mu_fr_tot;

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

  EL::StatusCode readRFEfficiencies ();
  EL::StatusCode getMMWeightAndError ();
  //std::vector<float>&  calculateMMWeights ();
  float calculateFinalMMWeight ( const std::string& region, float f1, float f2, float r1, float r2, const std::string& variation )

  EL::StatusCode calculateMMWeights ();

  EL::StatusCode readQMisIDRates ();
  EL::StatusCode getQMisIDRatesAndError ( std::shared_ptr<leptonObj> lep,
					  float& r, float& r_up, float& r_dn,
					  const std::string& selection );
  EL::StatusCode calculateQMisIDWeights ();


  EL::StatusCode enableSelectedBranches ();
  EL::StatusCode setOutputBranches ();
  EL::StatusCode clearBranches ( const std::string& type );

};

#endif
