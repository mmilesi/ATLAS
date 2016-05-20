/**
 * @file   RealFakeEffFitter.cc
 * @author Marco Milesi <marco.milesi@cern.ch>
 * @date   10 March 2016
 * @brief  ROOT macro to measure real/fake efficiencies via a binned maximum likelihood fit.
 *
 * Use ROOT::TMinuit class
 * @see https://root.cern.ch/doc/master/classTMinuit.html
 *
 */

#include <TROOT.h>
#include "TSystem.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "AtlasStyle.C"
#include "TLine.h"

// -------------------------------------------------------------------------------------

// GLOBAL VARIABLES, FUNCTIONS (TMinuit fiited function must be global, unfortunately...)

bool g_verbose(false);

int g_nPtBins_Squared;

/**
  Estimated parameters' indexes
  These vectors contain the corresponding index in the total par[] vector of the TMinuit object
  --> can use these quantities in the likelihood function!
*/
std::vector<int> g_r_idxs;
std::vector<int> g_f_idxs;
std::vector<int> g_RR_idxs;
std::vector<int> g_FR_idxs;
std::vector<int> g_RF_idxs;
//std::vector<int> g_FF_idxs;
std::vector<int> g_r1_idxs;  // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
std::vector<int> g_r2_idxs;  // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin
std::vector<int> g_f1_idxs;  // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
std::vector<int> g_f2_idxs;  // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin};

std::vector<double> g_fitted_r1_el; // fitted r value for leading e in every 2D bin
std::vector<double> g_fitted_r2_el; // fitted r value for subleading e in every 2D bin
std::vector<double> g_fitted_r1_mu; // fitted r value for leading mu in every 2D bin
std::vector<double> g_fitted_r2_mu; // fitted r value for subleading mu in every 2D bin

std::vector<double> g_TT_resized;
std::vector<double> g_TL_resized;
std::vector<double> g_LT_resized;
std::vector<double> g_LL_resized;

bool g_elec(false);
bool g_muon(false);

bool g_reff(false);
bool g_feff(false);
bool g_RR(false);
bool g_RF(false);
bool g_FR(false);
//bool g_FF(false);

/**
  Define global likelihood function to be minimised
*/
void myLikelihood( int& nDim, double* gout, double& result, double par[], int flg ) {

   if ( g_verbose ) { std::cout << "" << std::endl; }

   double likelihood_TT(0.0), likelihood_TL(0.0), likelihood_LT(0.0), likelihood_LL(0.0);
   double likelihood(0.0);
   double obs_TT(-1.0), obs_TL(-1.0), obs_LT(-1.0), obs_LL(-1.0);
   double exp_TT(-1.0), exp_TL(-1.0), exp_LT(-1.0), exp_LL(-1.0);

   double r1(0.0), r2(0.0), f1(0.0), f2(0.0);
   double RR(0.0), RF(0.0), FR(0.0), FF(0.0);

   // Compute the likelihood by summing TL,LT.. sub-blocks in every bin of the (linearised) 2D pT histogram
   //
   for ( int ibin = 0; ibin < g_nPtBins_Squared; ++ibin  ) {

     // Read the likelihood parameters from the TMinuit par[] vector
     //

     double fitted_r1(0.0), fitted_r2(0.0);
     if ( g_elec && g_feff ) { fitted_r1 = g_fitted_r1_el.at(ibin); fitted_r2 = g_fitted_r2_el.at(ibin); }
     if ( g_muon && g_feff ) { fitted_r1 = g_fitted_r1_mu.at(ibin); fitted_r2 = g_fitted_r2_mu.at(ibin); }

     // This assumes real efficiencies are fitted first
     // When fitting fake efficiencies, the previously fitted values for r are used

     r1 = ( g_reff )  ? par[ g_r1_idxs.at(ibin) ] : fitted_r1;
     r2 = ( g_reff )  ? par[ g_r2_idxs.at(ibin) ] : fitted_r2;
     f1 = ( g_feff )  ? par[ g_f1_idxs.at(ibin) ] : 0.0;
     f2 = ( g_feff )  ? par[ g_f2_idxs.at(ibin) ] : 0.0;
     RR = (  g_RR  )  ? par[ g_RR_idxs.at(ibin) ] : 0.0;
     RF = (  g_RF  )  ? par[ g_RF_idxs.at(ibin) ] : 0.0;
     FR = (  g_FR  )  ? par[ g_FR_idxs.at(ibin) ] : 0.0;
     //FF = (  g_FF  )  ? par[ g_FF_idxs.at(ibin) ] : 0.0; // can we set this to 0 by brute force?
     FF = 0.0;

     if ( g_verbose ) {
       std::cout << "" << std::endl;
       Info("myLikelihood()","Ingredients of likelihood - global bin (%i):", ibin );
       std::cout << "\tr1 = " << r1 << std::endl;
       std::cout << "\tr2 = " << r2 << std::endl;
       std::cout << "\tf1 = " << f1 << std::endl;
       std::cout << "\tf2 = " << f2 << std::endl;
       std::cout << "\tRR = " << RR << std::endl;
       std::cout << "\tRF = " << RF << std::endl;
       std::cout << "\tFR = " << FR << std::endl;
       std::cout << "\tFF = " << FF << std::endl;
       std::cout << "" << std::endl;
     }

     // TT block
     //
     obs_TT  = g_TT_resized.at(ibin);
     exp_TT  = r1 * r2 * RR + r1 * f2 * RF + r2 * f1 * FR + f1 * f2 * FF;
     if ( g_verbose ) {
       std::cout << "" << std::endl;
       Info("myLikelihood()","Adding term in likelihood for observed selection: TT in global bin (%i)", ibin );
       std::cout << "\tobs_TT = " << obs_TT << " - exp_TT = " << exp_TT << " - likelihood block = " << ( obs_TT * log( exp_TT ) - ( exp_TT ) ) << std::endl;
       std::cout << "" << std::endl;
     }
     likelihood_TT = ( obs_TT * log( exp_TT ) - ( exp_TT ) );

     // TL block
     //
     obs_TL  = g_TL_resized.at(ibin);
     exp_TL  = r1 * ( 1 - r2 ) * RR + r1 * ( 1 - f2 ) * RF + f1 * ( 1 - r2 ) * FR + f1 * ( 1 - f2 ) * FF;
     if ( g_verbose ) {
       std::cout << "" << std::endl;
       Info("myLikelihood()","Adding term in likelihood for observed selection: TL in global bin (%i)", ibin );
       std::cout << "\tobs_TL = " << obs_TL << " - exp_TL = " << exp_TL << " - likelihood block = " << ( obs_TL * log( exp_TL ) - ( exp_TL ) ) << std::endl;
       std::cout << "" << std::endl;
     }
     likelihood_TL = ( obs_TL * log( exp_TL ) - ( exp_TL ) );

     // LT block
     //
     obs_LT  = g_LT_resized.at(ibin);
     exp_LT  = r2 * ( 1 - r1 ) * RR + f2 * ( 1 - r1 ) * RF + r2 * ( 1 - f1 ) * FR + f2 * ( 1 - f1 ) * FF;
     if ( g_verbose ) {
       std::cout << "" << std::endl;
       Info("myLikelihood()","Adding term in likelihood for observed selection: LT in global bin (%i)", ibin );
       std::cout << "\tobs_LT = " << obs_LT << " - exp_LT = " << exp_LT << " - likelihood block = " << ( obs_LT * log( exp_LT ) - ( exp_LT ) ) << std::endl;
       std::cout << "" << std::endl;
     }
     likelihood_LT = ( obs_LT * log( exp_LT ) - ( exp_LT ) );

     // LL block
     //
     obs_LL  = g_LL_resized.at(ibin);
     exp_LL  = ( 1 - r1 ) * ( 1 - r2 ) * RR + ( 1 - r1 ) * ( 1 - f2 ) * RF + ( 1 - f1 ) * ( 1 - r2 ) * FR + ( 1 - f1 ) * ( 1 - f2 ) * FF;
     if ( g_verbose ) {
       std::cout << "" << std::endl;
       Info("myLikelihood()","Adding term in likelihood for observed selection: LL in global bin (%i)", ibin );
       std::cout << "\tobs_LL = " << obs_LL << " - exp_LL = " << exp_LL << " - likelihood block = " << ( obs_LL * log( exp_LL ) - ( exp_LL ) ) << std::endl;
       std::cout << "" << std::endl;
     }
     likelihood_LL = ( obs_LL * log( exp_LL ) - ( exp_LL ) );

     if ( g_verbose ) { Info("myLikelihood()","Total likelihood block L(%i) = %2f.", ibin, ( likelihood_TT + likelihood_TL + likelihood_LT + likelihood_LL ) ); }

     likelihood += ( likelihood_TT + likelihood_TL + likelihood_LT + likelihood_LL );
   }

   if ( g_verbose ) { Info("myLikelihood()","===> Final likelihood - 2 * log(L) = %2f.", - 2.0 * likelihood ); }

   result = - 2.0 * likelihood;
}

/**
  Function to get histograms from file
*/
template<typename T>
T* get_object( TFile& file, const std::string& name ) {
    T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
    if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
    return obj;
}

/**
  Function to print content of a container
*/
template<typename T>
void printContainer( T& container, const std::string& message ) {
  std::cout << "" << std::endl;
  Info("printContainer()","%s",message.c_str() );
  unsigned int idx(0);
  for ( const auto& itr : container ) {
     std::cout << "\t[" << idx << "] = " << itr << std::endl;
     ++idx;
  }
  std::cout << "" << std::endl;
}

/**
  Function to print content of a container of errors
*/
void printErrorContainer( std::vector<std::tuple<double,double,double> >& container, const std::string& message ) {
  std::cout << "" << std::endl;
  Info("printErrorContainer()","%s",message.c_str() );
  unsigned int idx(0);
  for ( const auto& itr : container ) {
     std::cout << "\t[" << idx << "] - err UP = " << std::get<0>(itr) << " err DN = " << std::get<1>(itr) << " err CENT = " << std::get<2>(itr) << std::endl;
     ++idx;
  }
  std::cout << "" << std::endl;
}

// -------------------------------------------------------------------------------------

/**
  The main class for the fit
*/
class LHFitter {

  public :

    bool m_doSubtraction;
    bool m_doRebinning;

    enum kFlavour {
      ELEC  = 0,
      MUON  = 1,
    };

    enum kEfficiency {
      REAL  = 0,
      FAKE  = 1,
    };

    enum kVerbosity {
      NONE    = 0,
      DEBUG   = 1,
      VERBOSE = 2,
    };

    /**
      Class constructor:
      Set the type of efficiency to be fit, and the lepton flavour
    */
    LHFitter( kFlavour FLAVOUR, kEfficiency EFFICIENCY );

    /**
      Class destructor:
    */
    ~LHFitter();

    /**
      -) Read the yields for TT, TL, LT, LL, depending on the flavour and efficiency to be fit
      -) Reserve memory for parameters' vectors
      -) Set initial guesses for parameters r, f, RR, RF, FR, FF
      -) Create the TMinuit object and set NPAR, depending on initialisation
      -) Set the array of parameters for the fit (NB: every time this method gets called, the global parameters will be reset!)
    */
    void initialise();

    /**
      -) Calls MIGRAD and MINOS to perform the fit
      -) Save the fitted parameters into vectors
    */
    void fit();

    /**
      Set path for tag and probe files
    */
    inline void setTagAndProbePath( const std::string& path ) { m_input_path_TP = path; };

    /**
      Set path for input files with TL, LT... yields
    */
    inline void setInputHistPath( const std::string& path ) { m_input_path_yields = path; };

    inline void setVerbosity( kVerbosity VERBOSITY ) {
	if ( VERBOSITY == kVerbosity::DEBUG )   { m_debug = true; }
	if ( VERBOSITY == kVerbosity::VERBOSE ) { m_debug = true; m_verbose = true; }
    };

    /**
      Function to set bin grouping when rebinning w/ fixed bin size
    */
    inline void setBinGrouping( const int& groupsize ) {  m_nBinGroup = groupsize; };

    /**
      Function to set variable bins when rebinning w/ variable bin size
    */
    inline void setVariableBins( double *varbins, const int& nvarbins ) {
      m_useVariableBins = true;
      m_newNBins = nvarbins;
      m_newBins  = new double[m_newNBins+1];
      for ( int ibin(0); ibin < m_newNBins+1; ++ibin ) {
	  m_newBins[ibin] = varbins[ibin];
      }
    };

    static void useMC() {
      Info("useMC()","Fitting efficiencies on MC...");
      s_useMC = true;
    };

  private :

    /**
      Set this flag when doing closure on MC
    */
    static bool s_useMC;

    /**
      Upper/lower limits for fit parameters and step size
      Declared static so they can be initialised once and for all
    */
    static double s_step_eff;
    static double s_up_eff;
    static double s_dn_eff;
    static double s_step_yields;
    static double s_up_yields;
    static double s_dn_yields;

    kFlavour    m_flavour;
    kEfficiency m_efficiency;

    std::string m_flavour_str;
    std::string m_efficiency_str;

    int     m_nBinGroup; /** Set this if rebinning w/ fixed bin size */
    bool    m_useVariableBins;
    int     m_newNBins;  /** Number of bins to set if rebinning with variable bin size */
    double* m_newBins;   /** Array with variable bin size    */

    bool m_debug;
    bool m_verbose;

    int m_nPtBins_Linear;   /** including underflow and overflow */
    int m_nPtBins_Squared;  /** including underflow and overflow */

    std::vector<std::string> m_obs_selection;

    /**
      The TMinuit object
    */
    TMinuit* m_myFitter;

    /**
      The number of parameters of the fit
    */
    int m_NPAR;

    /**
      Container for the input histograms
    */
    std::vector<TH2D*> m_histograms;

    /**
      Observables
    */
    std::vector<double> m_TT;
    std::vector<double> m_TL;
    std::vector<double> m_LT;
    std::vector<double> m_LL;

    /**
      Initial values for parameters
    */
    std::vector<double> m_r_init;
    double		m_r_init_avg;
    std::vector<double> m_f_init;
    double		m_f_init_avg;
    std::vector<double> m_RR_init;
    std::vector<double> m_RF_init;
    std::vector<double> m_FR_init;
    //std::vector<double> m_FF_init;

    /**
      Final parameter values
    */
    std::vector<double> m_r_vals;
    std::vector<double> m_f_vals;
    std::vector<double> m_RR_vals;
    std::vector<double> m_RF_vals;
    std::vector<double> m_FR_vals;
    //std::vector<double> m_FF_vals;

    /**
      Final parameter asymmetric errors

      Values in the tuple correspond to:

      [0] --> eplus ( MINOS error UP )
      [1] --> eminus ( MINOS error DN )
      [2] --> eparab ( 'parabolic' error (from error matrix) )

    */
    std::vector<std::tuple<double,double,double> > m_r_errs;
    std::vector<std::tuple<double,double,double> > m_f_errs;
    std::vector<std::tuple<double,double,double> > m_RR_errs;
    std::vector<std::tuple<double,double,double> > m_RF_errs;
    std::vector<std::tuple<double,double,double> > m_FR_errs;
    //std::vector<std::tuple<double,double,double> > m_FF_errs;

    std::string m_input_path_TP;
    std::string m_input_path_yields;

    /**
      Function to get number of bins below (and including) diagonal in N X N 2D hist
    */
    inline int areaGrid( const int& n );

    /**
      Get tag and probe efficiencies from input file
    */
    void readTagAndProbeEff();

    /**
      Get input histograms for TT, TL, LT, LL from ROOT file
      Read the TT, TL, LT, LL events bin-by-bin from input histograms
    */
    void getHists();

    /**
      Find initial guesses for parameters. Use inputs from tag-and-probe measurement
    */
    void getEducatedGuess();

    /**
      Reset global flags
    */
    void resetGlobFlags();

    /**
      Update the content of error containers w/ info from MINOS
    */
    void getParametersAndErrors();

    /**
      Save efficiencies and their errors in a ROOT/text file
    */
    void  saveEfficiencies();

};

bool   LHFitter::s_useMC = false;

double LHFitter::s_step_eff    = 1e-3;
double LHFitter::s_up_eff      = 1.0;
double LHFitter::s_dn_eff      = 0.0;
double LHFitter::s_step_yields = 1e-3;
double LHFitter::s_up_yields   = 1e6;
double LHFitter::s_dn_yields   = 0.0;

// ----------------------------------------------------------------------------------------------------------------------

LHFitter :: LHFitter( kFlavour FLAVOUR, kEfficiency EFFICIENCY ) :
  m_debug(false),
  m_verbose(false),
  m_doSubtraction(false),
  m_doRebinning(false),
  m_nBinGroup(1),
  m_useVariableBins(false),
  m_newBins(nullptr),
  m_myFitter(nullptr),
  m_obs_selection({"TT","TL","LT","LL"}),
  m_r_init_avg(0.0),
  m_f_init_avg(0.0)
{

  m_flavour = FLAVOUR;
  m_efficiency = EFFICIENCY;

  if ( m_efficiency == kEfficiency::REAL ) { m_efficiency_str = "REAL"; }
  if ( m_efficiency == kEfficiency::FAKE ) { m_efficiency_str = "FAKE"; }
  if ( m_flavour == kFlavour::ELEC ) { m_flavour_str = "ELECTRON"; }
  if ( m_flavour == kFlavour::MUON ) { m_flavour_str = "MUON"; }

  Info("LHFitter()","Creating class instance to fit %s efficiency for %s... \n", m_efficiency_str.c_str(), m_flavour_str.c_str() );

  // Every time an instance is created, the global variables must be reset, no matter what
  //
  this->resetGlobFlags();
}

// ----------------------------------------------------------------------------------------------------------------------

LHFitter::~LHFitter() {
  if ( m_newBins != nullptr )  delete[] m_newBins;
}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: initialise() {

  Info("initialise()","Setting up fit...");

  // Get input histograms and store them
  //
  this->getHists();

  // Reserve memory for (initial) parameter containers
  //
  m_r_init.reserve(m_nPtBins_Linear);
  m_f_init.reserve(m_nPtBins_Linear);
  m_RR_init.reserve(m_nPtBins_Squared);
  m_RF_init.reserve(m_nPtBins_Squared);
  m_FR_init.reserve(m_nPtBins_Squared);
  //m_FF_init.reserve(m_nPtBins_Squared);

  if ( m_verbose ) {
     printContainer( m_TT, "Printing content of TT:" );
     printContainer( m_TL, "Printing content of TL:" );
     printContainer( m_LT, "Printing content of LT:" );
     printContainer( m_LL, "Printing content of LL:" );
  }

  std::string eff_type("");
  std::string flavour("");

  // Set the number of parameteres of the fit
  //
  int NFLAV(1);
  int NCOMP(0);
  if ( m_efficiency == kEfficiency::REAL ) {
    NCOMP = 1; // RR only
  }
  if ( m_efficiency == kEfficiency::FAKE ) {
   //NCOMP = 3; // RF, FR, FF
    NCOMP = 2; // RF, FR
  }

  std::cout << "" << std::endl;
  std::cout << " Number of 1D pT bins (including O/Flow) = "<<  m_nPtBins_Linear << std::endl;
  std::cout << " Number of effective 2D pT bins (including O/Flow) = "<<  m_nPtBins_Squared << std::endl;
  std::cout << " ----------------------------"<<  std::endl;
  std::cout << " Number of flavour bins = "<<  NFLAV << std::endl;
  std::cout << " Number of RR, RF, FR, FF bins = "<<  NCOMP << std::endl;
  std::cout << "" << std::endl;

  // Total number of parameters to be estimated in the fit
  //
  // 1st factor in sum: YIELDS
  // 2nd factor in sum: EFFICIENCIES
  //
  const int NPAR_YIELDS = m_nPtBins_Squared * NFLAV * NCOMP;
  const int NPAR_EFF	= m_nPtBins_Linear;

  m_NPAR = NPAR_YIELDS + NPAR_EFF;

  const int m_NOBS = m_obs_selection.size() * m_nPtBins_Squared;

  Info("initialise()","Number of observables in fit ===> %i", m_NOBS);
  Info("initialise()","Number of free parameters in fit ===> %i (efficiencies) + %i (yields) = %i\n", NPAR_EFF, NPAR_YIELDS, m_NPAR);

  int ierflg(0);

  // Create the TMinuit object
  //
  m_myFitter = new TMinuit(m_NPAR);

  // Set the fitting function
  //
  m_myFitter->SetFCN(myLikelihood);

  // Get an educated guess for the initial parameters, from the tag-and-probe measurement:
  //
  // r(i), f(i)  ( in 1D bins of pT )
  //
  // RR(j,k), RF(j,k), FR(j,k), FF(j,k) ( in 2D bins of lead-sublead pT )
  //
  this->getEducatedGuess();

  if ( m_debug ) {
     printContainer( m_r_init, "Printing content of r_init:" );
     printContainer( m_f_init, "Printing content of f_init:" );
     if ( m_efficiency == kEfficiency::REAL ) printContainer( m_RR_init,"Printing content of RR_init:" );
     if ( m_efficiency == kEfficiency::FAKE ) printContainer( m_RF_init,"Printing content of RF_init:" );
     if ( m_efficiency == kEfficiency::FAKE ) printContainer( m_FR_init,"Printing content of FR_init:" );
     //if ( m_efficiency == kEfficiency::FAKE ) printContainer( m_FF_init,"Printing content of FF_init:" );
     printContainer( g_TT_resized, "Printing content of TT (resized):" );
     printContainer( g_TL_resized, "Printing content of TL (resized):" );
     printContainer( g_LT_resized, "Printing content of LT (resized):" );
     printContainer( g_LL_resized, "Printing content of LL (resized):" );
  }

  // Set the parameters of the fit into the toatl par[] array of the TMinuit object
  // The order will be:
  // -) r  ( x nPtBins_Linear )
  // -) f  ( x nPtBins_Linear )
  // -) RR ( x nPtBins_Squared )
  // -) RF ( x nPtBins_Squared )
  // -) FR ( x nPtBins_Squared )
  // -) FF ( x nPtBins_Squared )

  // Set parameters for r,f efficiencies. This is just a 1D loop over pT bins
  //
  std::string param_name("");
  int offset(0); // needed for setting parameter offsets (TMinuit can accept only one huge array of parameters...)
  int param_idx(0);

  // --------------------
  // Set parameters for r
  // --------------------

  if ( m_efficiency == kEfficiency::REAL ) {

    g_reff = true;

    for ( auto ibin(0); ibin < m_nPtBins_Linear; ++ibin ) {
      param_name = "real efficiency - bin [" + std::to_string(ibin) + "]";
      param_idx  = offset + ibin;
      m_myFitter->mnparm( param_idx, param_name.c_str(), m_r_init.at(ibin), s_step_eff, s_dn_eff, s_up_eff, ierflg );
      g_r_idxs.push_back( param_idx );
    }

    // A trick: save the index of the bin for the real efficiency for every (lead, sublead) lepton pair
    // It will contain the same indexes of r_idxs, but each one repeated nPtBins_Linear times
    //
    for ( int ibiny(0); ibiny < m_nPtBins_Linear; ++ibiny ) {
      for ( int ibinx(0); ibinx < m_nPtBins_Linear; ++ibinx ) {
        if ( ibiny > ibinx ) { continue; }
        g_r1_idxs.push_back( ibinx );
        g_r2_idxs.push_back( ibiny );
      }
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Linear;
  }

  // --------------------
  // Set parameters for f
  // --------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    g_feff = true;

    for ( auto ibin(0); ibin < m_nPtBins_Linear; ++ibin ) {
      param_name = "fake efficiency - bin [" + std::to_string(ibin) + "]";
      param_idx  = offset + ibin;
      m_myFitter->mnparm( param_idx, param_name.c_str(), m_f_init.at(ibin), s_step_eff, s_dn_eff, s_up_eff, ierflg );
      g_f_idxs.push_back( param_idx );
    }

    // A trick: save the index of the bin for the real efficiency for every (lead, sublead) lepton pair
    // It will contain the same indexes of f_idxs, but each one repeated nPtBins_Linear times
    //
    for ( int ibiny(0); ibiny < m_nPtBins_Linear; ++ibiny ) {
      for ( int ibinx(0); ibinx < m_nPtBins_Linear; ++ibinx ) {
        if ( ibiny > ibinx ) { continue; }
        g_f1_idxs.push_back( ibinx );
        g_f2_idxs.push_back( ibiny );
      }
    }
     // Set the offset for total parameter index
    //
    offset += m_nPtBins_Linear;
  }

  // ---------------------
  // Set parameters for RR
  // ---------------------

  if ( m_efficiency == kEfficiency::REAL ) {

    g_RR = true;

    for ( auto ibin(0); ibin < m_nPtBins_Squared; ++ibin ) {
      param_name = "RR - bin [" + std::to_string(ibin) + "]";
      param_idx = offset + ibin; // can use any of the input histograms: take the first by default
      m_myFitter->mnparm( param_idx, param_name.c_str(), m_RR_init.at(ibin), s_step_yields, s_dn_yields, s_up_yields, ierflg );
      g_RR_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }
  // ---------------------
  // Set parameters for RF
  // ---------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    g_RF = true;

    for ( auto ibin(0); ibin < m_nPtBins_Squared; ++ibin ) {
    	param_name = "RF - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	m_myFitter->mnparm( param_idx, param_name.c_str(), m_RF_init.at(ibin), s_step_yields, s_dn_yields, s_up_yields, ierflg );
    	g_RF_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }
  // ---------------------
  // Set parameters for FR
  // ---------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    g_FR = true;

    for ( auto ibin(0); ibin < m_nPtBins_Squared; ++ibin ) {
    	param_name = "FR - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	m_myFitter->mnparm( param_idx, param_name.c_str(), m_FR_init.at(ibin), s_step_yields, s_dn_yields, s_up_yields, ierflg );
    	g_FR_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }

  // ---------------------
  // Set parameters for FF
  // ---------------------

  /*
  if ( m_efficiency == kEfficiency::FAKE ) {

    g_FF = true;

    for ( auto ibin(0); ibin < m_nPtBins_Squared; ++ibin ) {
    	param_name = "FF - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	m_myFitter->mnparm( param_idx, param_name.c_str(), m_FF_init.at(ibin), s_step_yields, s_dn_yields, s_up_yields, ierflg );
    	g_FF_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }
  */

  // Set default values for final parameter/error vectors
  //
  int idx1(0);
  while ( idx1 < m_nPtBins_Linear ) {
    m_r_vals.push_back(0.0);
    m_f_vals.push_back(0.0);
    m_r_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    m_f_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    ++idx1;
  }
  int idx2(0);
  while ( idx2 < m_nPtBins_Squared ) {
    m_RR_vals.push_back(0.0);
    m_RF_vals.push_back(0.0);
    m_FR_vals.push_back(0.0);
    //m_FF_vals.push_back(0.0);
    m_RR_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    m_RF_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    m_FR_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    //m_FF_errs.push_back(std::make_tuple( 0.0, 0.0, 0.0 ));
    ++idx2;
  }

  if ( m_verbose ) {
    if ( m_efficiency == kEfficiency::REAL )  printContainer( g_r_idxs, "Printing content of r_idxs:" );
    if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_f_idxs, "Printing content of f_idxs:" );
    if ( m_efficiency == kEfficiency::REAL )  printContainer( g_r1_idxs, "Printing content of r1_idxs:" );
    if ( m_efficiency == kEfficiency::REAL )  printContainer( g_r2_idxs, "Printing content of r2_idxs:" );
    if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_f1_idxs, "Printing content of f1_idxs:" );
    if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_f2_idxs, "Printing content of f2_idxs:" );
    if ( m_efficiency == kEfficiency::REAL )  printContainer( g_RR_idxs, "Printing content of RR_idxs:" );
    if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_RF_idxs, "Printing content of RF_idxs:" );
    if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_FR_idxs, "Printing content of FR_idxs:" );
    //if ( m_efficiency == kEfficiency::FAKE )  printContainer( g_FF_idxs, "Printing content of FF_idxs:" );
  }

  Info("initialise()","Initialisation done!");

}
// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: fit() {

  std::cout << "" << std::endl;
  Info("fit()","Performing the fit...\n" );

  double arglist[m_NPAR];
  int ierflg(0);

  arglist[0] = 1000; // maximum number of iterations
  m_myFitter->mnexcm("MIGRAD",arglist,2/*0*/,ierflg);

  std::cout << "\n\n" << std::endl;

  if ( !m_myFitter->fCstatu.Contains("CONVERGED") ) {
    Error("fit()","No convergence at fitting! Minuit return string: %s", m_myFitter->fCstatu.Data() );
    //exit(-1);
  }

  arglist[0] = 0.0;
  arglist[1] = 1.0;
  m_myFitter->mnexcm("MINOS",arglist,1,ierflg); // calculate MINOS errors for all the parameters

  if ( m_verbose ) {
    if ( m_efficiency == kEfficiency::REAL ) {
      printContainer( m_r_vals, "Printing content of r_vals (DEFAULT):" );
      printErrorContainer( m_r_errs, "Printing content of r_errs (DEFAULT):" );
      printContainer( m_RR_vals, "Printing content of RR_vals (DEFAULT):" );
      printErrorContainer( m_RR_errs, "Printing content of RR_errs (DEFAULT):" );
    }
    if ( m_efficiency == kEfficiency::FAKE ) {
      printContainer( m_f_vals, "Printing content of f_vals (DEFAULT):" );
      printErrorContainer( m_f_errs, "Printing content of f_errs (DEFAULT):" );
      printContainer( m_RF_vals, "Printing content of RF_vals (DEFAULT):" );
      printErrorContainer( m_RF_errs, "Printing content of RF_errs (DEFAULT):" );
      printContainer( m_FR_vals, "Printing content of FR_vals (DEFAULT):" );
      printErrorContainer( m_FR_errs, "Printing content of FR_errs (DEFAULT):" );
    }
  }

  // Get final parameters, parabolic error and the asymmetric errors from MINOS
  //
  this->getParametersAndErrors();

  if ( m_debug ) {
    if ( m_efficiency == kEfficiency::REAL ) {
      printContainer( m_r_vals, "Printing content of r_vals (POST-FIT):" );
      printErrorContainer( m_r_errs, "Printing content of r_errs (POST-FIT):" );
      printContainer( m_RR_vals, "Printing content of RR_vals (POST-FIT):" );
      printErrorContainer( m_RR_errs, "Printing content of RR_errs (POST-FIT):" );
      if ( m_flavour == kFlavour::ELEC ) {
        printContainer( g_fitted_r1_el, "Printing content of g_fitted_r1_el (POST-FIT):" );
        printContainer( g_fitted_r2_el, "Printing content of g_fitted_r2_el (POST-FIT):" );
      }
      if ( m_flavour == kFlavour::MUON ) {
        printContainer( g_fitted_r1_mu, "Printing content of g_fitted_r1_mu (POST-FIT):" );
        printContainer( g_fitted_r2_mu, "Printing content of g_fitted_r2_mu (POST-FIT):" );
      }
    }
    if ( m_efficiency == kEfficiency::FAKE ) {
      printContainer( m_f_vals, "Printing content of f_vals (POST-FIT):" );
      printErrorContainer( m_f_errs, "Printing content of f_errs (POST-FIT):" );
      printContainer( m_RF_vals, "Printing content of RF_vals (POST-FIT):" );
      printErrorContainer( m_RF_errs, "Printing content of RF_errs (POST-FIT):" );
      printContainer( m_FR_vals, "Printing content of FR_vals (POST-FIT):" );
      printErrorContainer( m_FR_errs, "Printing content of FR_errs (POST-FIT):" );
    }
  }

  // Save results to output
  //
  this->saveEfficiencies();

  // Print fit statistic
  //
  double best_min;  // the best function value found so far
  double est_vdist; // the estimated vertical distance remaining to minimum
  double err_def;   // the value of UP defining parameter uncertainties
  int nvpar;	    // the number of currently variable parameters
  int nparx;	    // the highest (external) parameter number defined by user
  int icstat;	    // a status integer indicating how good is the covariance matrix:
  		    //  0= not calculated at all
        	    //  1= approximation only, not accurate
        	    //  2= full matrix, but forced positive-definite
        	    //  3= full accurate covariance matrix

  m_myFitter->mnstat(best_min,est_vdist,err_def,nvpar,nparx,icstat);

  std::cout << "" << std::endl;
  Info("fit()","************************************************" );
  std::cout << "" << std::endl;

  switch (icstat) {
   case 0 : Error("fit()","No covariance matrix was calculated! Exiting..." );
  	    exit(-1);
  	    break;
   case 1 : Warning("fit()","An approximated covariance matrix was calculated! Not accurate..." );
  	    break;
   case 2 : Warning("fit()","Full covariance matrix was calculated, but forced to be positive-definite..." );
  	    break;
   case 3 : Info("fit()","Full covariance matrix was calculated :)" );
  	    break;
  }

  Info("fit()","Minimum of likelihood function = %f", best_min );
  Info("fit()","Estimated vert. distance to min. = %f", est_vdist );
  Info("fit()","Value of UP defining parameter unc. = %f", err_def );
  Info("fit()","Number of variable parameters = %i", nvpar );
  Info("fit()","Highest number of parameters defined by user = %i", nparx );
  std::cout << "" << std::endl;
  Info("fit()","************************************************\n" );

}

// ----------------------------------------------------------------------------------------------------------------------

int LHFitter :: areaGrid( const int& n ) {
  int area(0);
  for ( int i = 0; i < n; ++i ) {
    area += n - i;
  }
  return area;
}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: resetGlobFlags(){

  Info("resetGlobFlags()","Resetting all global variables..." );

  g_nPtBins_Squared = 0;

  g_r_idxs.clear();
  g_f_idxs.clear();
  g_RR_idxs.clear();
  g_FR_idxs.clear();
  g_RF_idxs.clear();
  //g_FF_idxs.clear();
  g_r1_idxs.clear();
  g_r2_idxs.clear();
  g_f1_idxs.clear();
  g_f2_idxs.clear();

  // these should not be reset...
  //g_fitted_r1_el.clear();
  //g_fitted_r2_el.clear();
  //g_fitted_r1_mu.clear();
  //g_fitted_r2_mu.clear();

  g_TT_resized.clear();
  g_TL_resized.clear();
  g_LT_resized.clear();
  g_LL_resized.clear();

  g_elec = ( m_flavour == kFlavour::ELEC );
  g_muon = ( m_flavour == kFlavour::MUON );

  g_reff = false;
  g_feff = false;
  g_RR   = false;
  g_RF   = false;
  g_FR   = false;
  //g_FF   = false;

}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: readTagAndProbeEff() {

  Info("readTagAndProbeEff()", "Reading histograms with tag and probe efficiencies from input...");

  std::string path     = m_input_path_TP + "Rates.root";
  std::string path_avg = m_input_path_TP + "AvgRates.root";

  TFile *file = TFile::Open(path.c_str());
  if ( !file->IsOpen() ) {
    SysError("readTagAndProbeEff()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
    exit(-1);
  }

  TFile *file_avg = TFile::Open(path_avg.c_str());
  if ( !file_avg->IsOpen() ) {
    SysError("readTagAndProbeEff()", "Failed to open ROOT file from path: %s . Aborting", path_avg.c_str() );
    exit(-1);
  }

  std::string filename_r(""), filename_r_avg("");
  std::string filename_f(""), filename_f_avg("");

  std::string eff_type = ( !s_useMC ) ? "observed" : "expected";

  if ( m_flavour == kFlavour::ELEC ) {
    filename_r = filename_r_avg = "El_ProbePt_Real_Efficiency_" + eff_type;
    filename_f = filename_f_avg = "El_ProbePt_Fake_Efficiency_" + eff_type;
  }
  if ( m_flavour == kFlavour::MUON ) {
    filename_r = filename_r_avg = "Mu_ProbePt_Real_Efficiency_" + eff_type;
    filename_f = filename_f_avg = "Mu_ProbePt_Fake_Efficiency_" + eff_type;
  }

  TH1D *reff = get_object<TH1D>( *file, filename_r );
  TH1D *feff = get_object<TH1D>( *file, filename_f );

  // Do not read underflow, but do read overflow!
  //
  for ( int ibin(1); ibin <= reff->GetNbinsX()+1; ++ibin ) {
    m_r_init.push_back( reff->GetBinContent(ibin) );
  }
  for ( int ibin(1); ibin <= feff->GetNbinsX()+1; ++ibin ) {
    m_f_init.push_back( feff->GetBinContent(ibin) );
  }

  // Save average efficiencies

  TH1D *reff_avg = get_object<TH1D>( *file_avg, filename_r_avg );
  TH1D *feff_avg = get_object<TH1D>( *file_avg, filename_f_avg );

  m_r_init_avg = reff_avg->GetBinContent(1);
  m_f_init_avg = feff_avg->GetBinContent(1);

}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: getHists() {

  Info("getHists()","Extracting histograms for TT, TL, LT, LL yields from input ROOT files...");

  std::string path("");
  std::string charge("");
  std::string flav("");

  if      ( m_efficiency == kEfficiency::REAL ) charge = "OS";
  else if ( m_efficiency == kEfficiency::FAKE ) charge = "SS";

  if      ( m_flavour == kFlavour::ELEC ) flav = "ElEl";
  else if ( m_flavour == kFlavour::MUON ) flav = "MuMu";

  for ( const auto& sl : m_obs_selection ) {

    path = m_input_path_yields + charge + "_" + flav + "_" + sl + "/" + charge + "_" + flav + "_" + sl + "_Lep0Pt_VS_Lep1Pt.root";

    if ( m_debug ) { Info("getHists()","Reading histogram from:\t %s", path.c_str() ); }

    TFile *file = TFile::Open(path.c_str());
    if ( !file->IsOpen() ) {
      SysError("getHists()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
      exit(-1);
    }

    std::string eff_type = ( !s_useMC ) ? "observed" : "expected";
    TH2D *hist = get_object<TH2D>( *file, eff_type.c_str() );

    // -) Do ( !prompt & charge flip ) subtraction in OS
    // -) Do ( prompt & charge flip ) subtraction in SS
    //
    if  ( !s_useMC && m_doSubtraction ) {

      Info("getHists()", "Subtracting bkgs to data...");

      TH2D *hist_to_sub = get_object<TH2D>( *file, "expected" );
      hist->Add(hist_to_sub, -1.0);

      // Set bin content to 0 if subtraction gives negative yield
      //
      for ( int ibiny(0); ibiny < hist->GetNbinsX()+2; ++ibiny ) {
        for ( int ibinx(0); ibinx < hist->GetNbinsX()+2; ++ibinx ) {
  	    if ( hist->GetBinContent( hist->GetBin(ibinx,ibiny)  ) < 0 ) { hist->SetBinContent( hist->GetBin(ibinx,ibiny), 0.0); }
  	}
      }
    }

    std::string new_name = charge + "_" + flav + "_" + sl;

    if ( m_debug ) { Info("getHists()","Storing histogram w/ name:\t %s for later use", new_name.c_str() ); }

    hist->SetName(new_name.c_str());

    hist->SetDirectory(0);

    if ( m_doRebinning ) {

      Info("getHists()","Rebinning histogram...");

      TH2D* htemp = dynamic_cast<TH2D*>( hist->Clone() );

      if ( m_nBinGroup > 1 ) {
	  Info("getHists()","Using fixed bin size, grouping %i bins", m_nBinGroup );
	  hist = dynamic_cast<TH2D*>( htemp->Rebin2D( m_nBinGroup, m_nBinGroup, htemp->GetName() ) );
      } else {
	  Info("getHists()","Using variable bin size:");
	  std::cout << "new nbins: " << m_newNBins << " - bin low edges: " << std::endl;
	  for ( int ibin(0); ibin < m_newNBins+1; ++ibin ) {
	      std::cout << m_newBins[ibin] << " ";
	  }
	  std::cout << std::endl;

	  // Need this trick to create 2D rebinned histogram w/ variable bin size
	  TH2D *hrebinned = new TH2D("rebinned",hist->GetTitle(), m_newNBins, m_newBins, m_newNBins, m_newBins );
	  hrebinned->SetName(hist->GetName());
	  hrebinned->SetDirectory(0);
	  TAxis *xaxis = hist->GetXaxis();
	  TAxis *yaxis = hist->GetYaxis();
	  for ( int ibiny(0); ibiny < hist->GetNbinsY()+2; ++ibiny ) {
            for ( int ibinx(0); ibinx < hist->GetNbinsX()+2; ++ibinx ) {
	      hrebinned->Fill( xaxis->GetBinCenter(ibinx), yaxis->GetBinCenter(ibiny), hist->GetBinContent(ibinx,ibiny) );
	    }
	  }
	  m_histograms.push_back(hrebinned);
	  continue;
      }

    }

    m_histograms.push_back(hist);

  }

  m_nPtBins_Linear  = ( m_histograms.at(0) )->GetNbinsX() + 1;
  m_nPtBins_Squared = areaGrid( m_nPtBins_Linear );

  // set it as a global variable
  //
  g_nPtBins_Squared = m_nPtBins_Squared;

  if ( m_debug ) { Info("getHists()","Reading observed histograms for TT, TL, LT, LL yields:"); }

  // UGLYYYYYYYYYY
  int idx(0);
  while ( idx < (m_nPtBins_Linear + 1) * (m_nPtBins_Linear + 1) ) {
    m_TT.push_back(0.0);
    m_TL.push_back(0.0);
    m_LT.push_back(0.0);
    m_LL.push_back(0.0);
    ++idx;
  }

  for ( auto hist : m_histograms ) {

    // To get overflow as well (but no underflow)
    //
    int firstxbin = 1;
    int firstybin = 1;
    int lastxbin  = m_nPtBins_Linear;
    int lastybin  = m_nPtBins_Linear;

    if ( m_debug ) { std::cout << "\t" << hist->GetName() << " - Integral: " << hist->Integral(firstxbin,lastxbin,firstybin,lastybin) << std::endl; }

    std::string charge("");
    if      ( m_efficiency == kEfficiency::REAL ) charge = "OS";
    else if ( m_efficiency == kEfficiency::FAKE ) charge = "SS";

    std::string flav("");
    if      ( m_flavour == kFlavour::ELEC ) flav = "ElEl";
    else if ( m_flavour == kFlavour::MUON ) flav = "MuMu";

    std::string TT_name = charge + "_" + flav + "_TT";
    std::string TL_name = charge + "_" + flav + "_TL";
    std::string LT_name = charge + "_" + flav + "_LT";
    std::string LL_name = charge + "_" + flav + "_LL";

    if ( strcmp( hist->GetName(), TT_name.c_str() ) == 0 ) {
      for ( int ibiny(0); ibiny <= m_nPtBins_Linear; ++ibiny ) {
        for ( int ibinx(0); ibinx <= m_nPtBins_Linear; ++ibinx ) {
          if ( m_verbose ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
          m_TT.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
        }
      }
    } else if ( strcmp( hist->GetName(), TL_name.c_str() ) == 0 ) {
      for ( int ibiny(0); ibiny <= m_nPtBins_Linear; ++ibiny ) {
        for ( int ibinx(0); ibinx <= m_nPtBins_Linear; ++ibinx ) {
          if ( m_verbose ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
          m_TL.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
        }
      }
    } else if ( strcmp( hist->GetName(), LT_name.c_str() ) == 0 ) {
      for ( int ibiny(0); ibiny <= m_nPtBins_Linear; ++ibiny ) {
        for ( int ibinx(0); ibinx <= m_nPtBins_Linear; ++ibinx ) {
          if ( m_verbose ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
          m_LT.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
        }
      }
    } else if ( strcmp( hist->GetName(), LL_name.c_str() ) == 0 ) {
      for ( int ibiny(0); ibiny <= m_nPtBins_Linear; ++ibiny ) {
        for ( int ibinx(0); ibinx <= m_nPtBins_Linear; ++ibinx ) {
          if ( m_verbose ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
          m_LL.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
        }
      }
    }

  }

}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: getEducatedGuess( ) {

  Info("getEducatedGuess()","Getting initial values for fit parameters...");

  // TEMP - hardcode efficiency
  //double r =  0.65;
  //double f = 0.172;

  //for ( int idx(0); idx < nPtBins_Linear-1; ++idx ) {
  //  m_r_init.push_back(r);
  //  m_f_init.push_back(f);
  //}

  // Read tag-and-probe efficiencies from input file
  //
  this->readTagAndProbeEff();

  // loop over 2D pT histogram, and depending on the pT value of the leading/subleading lepton,
  // read the corresponding tag-and-probe efficiency, read the TT, TL...yields per bin (set previously in the vectors),
  // and finally compute nRR, nRF...in every bin

  // Skip the (0,j) and (i,0) bins --> they are underflow
  // Skip the (i,j) bins where i < j (empty by construction)

  double r1(-1.0), r2(-1.0);
  double f1(-1.0), f2(-1.0);

  double nTT(0.0), nTL(0.0), nLT(0.0), nLL(0.0);

  double RR_default(1e-3), RF_default(1e-3), FR_default(1e-3), FF_default(1e-3);

  // Drop tail of m_r_init, m_f_init vectors if their size is greater than the actual number of bins
  while ( m_r_init.size() > m_nPtBins_Linear ) { m_r_init.pop_back(); }
  while ( m_f_init.size() > m_nPtBins_Linear ) { m_f_init.pop_back(); }

  // Add elements to tail of vector ( pushing back the average efficiency ) if size is smaller than the actual number of bins
  while ( m_r_init.size() < m_nPtBins_Linear ) { m_r_init.push_back( m_r_init_avg ); }
  while ( m_f_init.size() < m_nPtBins_Linear ) { m_f_init.push_back( m_f_init_avg ); }

  if ( m_flavour == kFlavour::MUON || m_efficiency == kEfficiency::FAKE || m_doRebinning ) {
    if ( m_doRebinning ) {
      std::cout << "" << std::endl;
      Info("getEducatedGuess()","===> REBINNING ACTIVATED! Will be using average r,f to determine initial RR, RF, FR, FF yields..." );
      Info("getEducatedGuess()","<r> = %.2f - <f> = %.2f", m_r_init_avg, m_f_init_avg );

      std::cout << "m_r_init.size() = " << m_r_init.size() << " - m_nPtBins_Linear = " << m_nPtBins_Linear << std::endl;
      std::cout << "m_f_init.size() = " << m_f_init.size() << " - m_nPtBins_Linear = " << m_nPtBins_Linear << std::endl;

      for ( int i(0); i < (int)m_r_init.size(); ++i ) {
	  std::cout << "idx [" << i << "] - r init = " << m_r_init.at(i) << std::endl;
      }
      for ( int i(0); i < (int)m_f_init.size(); ++i ) {
	  std::cout << "idx [" << i << "] - r init = " << m_f_init.at(i) << std::endl;
      }

      for ( int ibin(0); ibin < m_nPtBins_Linear; ++ibin ) {
        m_r_init.at(ibin) = m_r_init_avg;
	std::cout << "idx [" << ibin << "] - r init = " << m_r_init.at(ibin) << " - f init = " << m_f_init.at(ibin) << std::endl;
        m_f_init.at(ibin) = m_f_init_avg;
      }
    } else {
      std::cout << "" << std::endl;
      Info("getEducatedGuess()","===> Will be using average f to determine initial RR, RF, FR, FF yields..." );
      Info("getEducatedGuess()","<f> = %.2f", m_f_init_avg );
    }
  }

  int glob_bin_idx(-1);
  for ( int ibiny(1); ibiny <= m_nPtBins_Linear; ++ibiny ) {

    r2 = m_r_init.at(ibiny-1);
    f2 = m_f_init.at(ibiny-1);
    if ( m_flavour == kFlavour::MUON || m_efficiency == kEfficiency::FAKE  ) f2 = m_f_init_avg;

    for ( int ibinx(1); ibinx <= m_nPtBins_Linear; ++ibinx ) {

      r1 = m_r_init.at(ibinx-1);
      f1 = m_f_init.at(ibinx-1);
      if ( m_flavour == kFlavour::MUON || m_efficiency == kEfficiency::FAKE ) f1 = m_f_init_avg;

      // Skip the above-diagonal elements
      //
      if ( ibiny > ibinx ) { continue; }

      glob_bin_idx = m_histograms.at(0)->GetBin(ibinx,ibiny); // can use any of the input histograms: take the first by default

      nTT = m_TT.at(glob_bin_idx);
      nTL = m_TL.at(glob_bin_idx);
      nLT = m_LT.at(glob_bin_idx);
      nLL = m_LL.at(glob_bin_idx);

      if ( m_verbose ) {
        std::cout << "global bin: " << glob_bin_idx << std::endl;
        std::cout << "nTT: " << nTT << std::endl;
        std::cout << "nTL: " << nTL << std::endl;
        std::cout << "nLT: " << nLT << std::endl;
        std::cout << "nLL: " << nLL << std::endl;
        std::cout << "r1: " << r1 << std::endl;
        std::cout << "r2: " << r2 << std::endl;
        std::cout << "f1: " << f1 << std::endl;
        std::cout << "f2: " << f2 << std::endl;
        std::cout << "" << std::endl;
      }

      double alpha = ( r1 - f1 ) * ( r2 - f2 );

      double RR = ( 1.0 / alpha ) * ( ( 1 - f1 ) * ( 1 - f2 ) * nTT + ( f1 - 1 ) * f2 * nTL + ( f2 - 1 ) * f1 * nLT + f1 * f2 * nLL );
      double RF = ( 1.0 / alpha ) * ( ( f1 - 1 ) * ( 1 - r2 ) * nTT + ( 1 - f1 ) * r2 * nTL + ( 1 - r2 ) * f1 * nLT - f1 * r2 * nLL );
      double FR = ( 1.0 / alpha ) * ( ( r1 - 1 ) * ( 1 - f2 ) * nTT + ( 1 - r1 ) * f2 * nTL + ( 1 - f2 ) * r1 * nLT - r1 * f2 * nLL );
      //double FF = ( 1.0 / alpha ) * ( ( 1 - r1 ) * ( 1 - r2 ) * nTT + ( r1 - 1 ) * r2 * nTL + ( r2 - 1 ) * r1 * nLT + r1 * r2 * nLL );

      // Reset yields to zero if unphisically negative
      // (Actually to 1 to avoid warnings b/c of lower physical limit set on yields...)
      //
      if ( RR >= 0.0 ) { m_RR_init.push_back( RR ); } else { m_RR_init.push_back( 1.0 ); }
      if ( RF >= 0.0 ) { m_RF_init.push_back( RF ); } else { m_RF_init.push_back( 1.0 ); }
      if ( FR >= 0.0 ) { m_FR_init.push_back( FR ); } else { m_FR_init.push_back( 1.0 ); }
      //if ( FF >= 0.0 ) { m_FF_init.push_back( FF ); } else { m_FF_init.push_back( 1.0 ); }

      // Copy only relevant bins for observed yields into these containers
      // (will be the ones used in the likelihood)
      //
      g_TT_resized.push_back( nTT );
      g_TL_resized.push_back( nTL );
      g_LT_resized.push_back( nLT );
      g_LL_resized.push_back( nLL );
    }
  }

}

// ----------------------------------------------------------------------------------------------------------------------

void LHFitter :: getParametersAndErrors() {

  Info("getParametersAndErrors()","Retrieveing fitted parameters and their errors...");

  int offset(0); // needed for dealing with parameter offsets (TMinuit can accept only one huge array of parameters...)
  int param_idx(0);

  double globcc;

  // --------------------
  // Get errors for r
  // --------------------

  if ( m_efficiency == kEfficiency::REAL ) {

    for ( auto idx(0); idx < m_r_errs.size(); ++idx ) {
       param_idx = offset + idx;
       m_myFitter->GetParameter( param_idx, m_r_vals.at(idx), std::get<2>(m_r_errs.at(idx)) );
       m_myFitter->mnerrs( param_idx, std::get<0>(m_r_errs.at(idx)), std::get<1>(m_r_errs.at(idx)), std::get<2>(m_r_errs.at(idx)), globcc );
    }

    for ( int ibiny(1); ibiny <= m_nPtBins_Linear; ++ibiny ) {
      for ( int ibinx(1); ibinx <= m_nPtBins_Linear; ++ibinx ) {
        if ( ibiny > ibinx ) { continue; }
	if ( m_flavour == kFlavour::ELEC ) {
          g_fitted_r1_el.push_back( m_r_vals.at(ibinx-1) );
          g_fitted_r2_el.push_back( m_r_vals.at(ibiny-1) );
	}
	if ( m_flavour == kFlavour::MUON ) {
          g_fitted_r1_mu.push_back( m_r_vals.at(ibinx-1) );
          g_fitted_r2_mu.push_back( m_r_vals.at(ibiny-1) );
	}
      }
    }

    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Linear;
  }

  // --------------------
  // Get errors for f
  // --------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    for ( auto idx(0); idx < m_f_errs.size(); ++idx ) {
       param_idx = offset + idx;
       m_myFitter->GetParameter( param_idx, m_f_vals.at(idx), std::get<2>(m_f_errs.at(idx)) );
       m_myFitter->mnerrs( param_idx, std::get<0>(m_f_errs.at(idx)), std::get<1>(m_f_errs.at(idx)), std::get<2>(m_f_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Linear;
  }

  // ---------------------
  // Get errors for RR
  // ---------------------

  if ( m_efficiency == kEfficiency::REAL ) {

    for ( auto idx(0); idx < m_RR_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        m_myFitter->GetParameter( param_idx, m_RR_vals.at(idx), std::get<2>(m_RR_errs.at(idx)) );
        m_myFitter->mnerrs( param_idx, std::get<0>(m_RR_errs.at(idx)), std::get<1>(m_RR_errs.at(idx)), std::get<2>(m_RR_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }

  // ---------------------
  // Get errors for RF
  // ---------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    for ( auto idx(0); idx < m_RF_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        m_myFitter->GetParameter( param_idx, m_RF_vals.at(idx), std::get<2>(m_RF_errs.at(idx)) );
        m_myFitter->mnerrs( param_idx, std::get<0>(m_RF_errs.at(idx)), std::get<1>(m_RF_errs.at(idx)), std::get<2>(m_RF_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }

  // ---------------------
  // Get errors for FR
  // ---------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    for ( auto idx(0); idx < m_FR_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        m_myFitter->GetParameter( param_idx, m_FR_vals.at(idx), std::get<2>(m_FR_errs.at(idx)) );
        m_myFitter->mnerrs( param_idx, std::get<0>(m_FR_errs.at(idx)), std::get<1>(m_FR_errs.at(idx)), std::get<2>(m_FR_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }

  // ---------------------
  // Get errors for FF
  // ---------------------

  /*
  if ( m_efficiency == kEfficiency::FAKE ) {

    for ( auto idx(0); idx < m_FF_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        m_myFitter->GetParameter( param_idx, m_FF_vals.at(idx), std::get<2>(m_FF_errs.at(idx)) );
        m_myFitter->mnerrs( param_idx, std::get<0>(m_FF_errs.at(idx)), std::get<1>(m_FF_errs.at(idx)), std::get<2>(m_FF_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += m_nPtBins_Squared;
  }
  */
}

// ----------------------------------------------------------------------------------------------------------------------

void  LHFitter :: saveEfficiencies() {

  Info("saveEfficiencies()","Saving fitted efficiency to output...");

  //TH1::ResetBit(TH1::kCanRebin);

  std::string outfilename("LH_efficiencies");

  if      ( m_efficiency == kEfficiency::REAL ) { outfilename += "_real"; }
  else if ( m_efficiency == kEfficiency::FAKE ) { outfilename += "_fake"; }

  if      ( m_flavour == kFlavour::MUON ) { outfilename += "_mu"; }
  else if ( m_flavour == kFlavour::ELEC ) { outfilename += "_el"; }

  std::string rootfilename = outfilename + ".root";
  TFile outfile(rootfilename.c_str(),"RECREATE");

  std::string txtfilename = outfilename + ".txt";
  std::ofstream outtextfile(txtfilename.c_str(), std::ios::trunc);
  outtextfile << "Efficiencies for FF amd Matrix Method from LH fit \n";

  // --------------------
  // real efficiency
  // --------------------

  if ( m_efficiency == kEfficiency::REAL ) {

    TH1D *r_hist(nullptr);
    if ( m_doRebinning && m_useVariableBins ) {
	r_hist = new TH1D( "r_hist", "real efficiency", m_newNBins, m_newBins );
    } else {
      double rmin = ( m_histograms.at(0) )->GetXaxis()->GetBinLowEdge(1);
      double rmax = ( m_histograms.at(0) )->GetXaxis()->GetBinLowEdge( ( m_histograms.at(0) )->GetNbinsX()+1 );
      r_hist = new TH1D( "r_hist", "real efficiency", m_nPtBins_Linear-1, rmin, rmax );
    }
    r_hist->SetCanExtend(TH1::kXaxis);

    r_hist->GetYaxis()->SetTitle("Real efficiency");
    r_hist->GetYaxis()->SetRangeUser(0.0,1.0);

    std::string xtitle("");
    if ( m_flavour == kFlavour::MUON ) { xtitle += "muon pT [GeV]"; }
    if ( m_flavour == kFlavour::ELEC ) { xtitle += "electron pT [GeV]"; }

    r_hist->GetXaxis()->SetTitle(xtitle.c_str());

    outtextfile << "Real efficiency - " << xtitle << "\n";
    for ( auto idx(0); idx < m_nPtBins_Linear; ++idx ) {
       if ( m_debug ) { std:: cout << "Bin idx " << idx+1 << " central value: " << ( m_histograms.at(0) )->GetXaxis()->GetBinCenter(idx + 1) << " - value to fill in: " << m_r_vals.at(idx) << std::endl; }
       r_hist->Fill( ( m_histograms.at(0) )->GetXaxis()->GetBinCenter(idx + 1), m_r_vals.at(idx) );
       r_hist->SetBinError( idx + 1, std::get<2>(m_r_errs.at(idx)) );
       outtextfile << "{ Bin nr: " << idx << ", efficiency = " <<  m_r_vals.at(idx) << " + " << std::get<0>(m_r_errs.at(idx)) << " - " << fabs( std::get<1>(m_r_errs.at(idx)) ) << " ( +- " << std::get<2>(m_r_errs.at(idx)) << " ) }\n";
    }

    r_hist->Write();
  }

  // --------------------
  // fake efficiency
  // --------------------

  if ( m_efficiency == kEfficiency::FAKE ) {

    TH1D *f_hist(nullptr);
    if ( m_doRebinning && m_useVariableBins ) {
	f_hist = new TH1D( "f_hist", "fake efficiency", m_newNBins, m_newBins );
    } else {
      double rmin = ( m_histograms.at(0) )->GetXaxis()->GetBinLowEdge(1);
      double rmax = ( m_histograms.at(0) )->GetXaxis()->GetBinLowEdge( ( m_histograms.at(0) )->GetNbinsX()+1 );
      f_hist = new TH1D( "f_hist", "fake efficiency", m_nPtBins_Linear-1, rmin, rmax );
    }
    f_hist->SetCanExtend(TH1::kXaxis);

    f_hist->GetYaxis()->SetTitle("Fake efficiency");
    f_hist->GetYaxis()->SetRangeUser(0.0,1.0);

    std::string xtitle("");
    if ( m_flavour == kFlavour::MUON ) { xtitle += "muon pT [GeV]"; }
    if ( m_flavour == kFlavour::ELEC ) { xtitle += "electron pT [GeV]"; }

    f_hist->GetXaxis()->SetTitle(xtitle.c_str());

    outtextfile << "Fake efficiency - " << xtitle << "\n";
    for ( auto idx(0); idx < m_nPtBins_Linear; ++idx ) {
       if ( m_debug ) { std:: cout << "Bin idx " << idx+1 << " central value: " << ( m_histograms.at(0) )->GetXaxis()->GetBinCenter(idx + 1) << " - value to fill in: " << m_f_vals.at(idx) << std::endl; }
       f_hist->Fill( ( m_histograms.at(0) )->GetXaxis()->GetBinCenter(idx + 1), m_f_vals.at(idx) );
       f_hist->SetBinError( idx + 1, std::get<2>(m_f_errs.at(idx)) );
       outtextfile << "{ Bin nr: " << idx << ", efficiency = " <<  m_f_vals.at(idx) << " + " << std::get<0>(m_f_errs.at(idx)) << " - " << fabs( std::get<1>(m_f_errs.at(idx)) ) << " ( +- " << std::get<2>(m_f_errs.at(idx)) << " ) }\n";
    }

    f_hist->Write();
  }

  outfile.Close();
  outtextfile.close();

}

// ----------------------------------------------------------------------------------------------------------------------

// ------------------------------------------------------
// Call this function from command line
// ------------------------------------------------------

int main( int argc, char **argv ) {

    std::cout << "Starting up..." << std::endl;

    TH1::SetDefaultSumw2(kTRUE);

    //gSystem->Load( "libCintex.so" );
    //gROOT->ProcessLine("Cintex::Cintex::Enable();");

    ///if ( argc<2 ) {
    //   std::cout << "Missing input parameters: " << std::endl;
    //    std::cout << "[1] input path" << std::endl;
    //    return 0;
    //}

    //const char* input_name = argv[1];
    //std::string input(input_name);

    // ----------------------------------------------------

    // DO THE FIT ON DATA
    //const std::string tp_path("../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/");
    //const std::string input_path("../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/");

    // DO THE FIT ON TTBAR MC
    const std::string tp_path("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/");
    const std::string input_path("../OutputPlots_MMClosureRates_LHFit_25ns_v7_FinalSelection_NominalBinning/");
    LHFitter::useMC();

    // real efficiency - e

    LHFitter real_e( LHFitter::kFlavour::ELEC, LHFitter::kEfficiency::REAL );
    real_e.setVerbosity(LHFitter::kVerbosity::DEBUG);
    real_e.setTagAndProbePath(tp_path);
    real_e.setInputHistPath(input_path);
    real_e.m_doSubtraction = true;
    // REBINNING
    real_e.m_doRebinning = true;
    //real_e.setBinGrouping(2);
    double real_e_new_bins[8] = {10.0,15.0,20.0,25.0,30.0,40.0,60.0,200.0};
    real_e.setVariableBins( real_e_new_bins, 7 );
    real_e.initialise();
    real_e.fit();
///*
    // real efficiency - mu

    LHFitter real_mu( LHFitter::kFlavour::MUON, LHFitter::kEfficiency::REAL );
    //real_mu.setVerbosity(LHFitter::kVerbosity::DEBUG);
    real_mu.setTagAndProbePath(tp_path);
    real_mu.setInputHistPath(input_path);
    real_mu.m_doSubtraction = true;
    // REBINNING
    real_mu.m_doRebinning = true;
    //real_mu.setBinGrouping(2);
    double real_mu_new_bins[8] = {10.0,15.0,20.0,25.0,30.0,40.0,60.0,200.0};
    real_mu.setVariableBins( real_mu_new_bins, 7 );
    real_mu.initialise();
    real_mu.fit();

    // fake efficiency - e

    LHFitter fake_e( LHFitter::kFlavour::ELEC, LHFitter::kEfficiency::FAKE );
    fake_e.setVerbosity(LHFitter::kVerbosity::DEBUG);
    fake_e.setTagAndProbePath(tp_path);
    fake_e.setInputHistPath(input_path);
    fake_e.m_doSubtraction = true;
    // REBINNING
    fake_e.m_doRebinning = true;
    //fake_e.setBinGrouping(2);
    double fake_e_new_bins[6] = {10.0,15.0,20.0,25.0,40.0,200.0};
    fake_e.setVariableBins( fake_e_new_bins, 5 );
    fake_e.initialise();
    fake_e.fit();

    // fake efficiency - mu

    LHFitter fake_mu( LHFitter::kFlavour::MUON, LHFitter::kEfficiency::FAKE );
    fake_mu.setVerbosity(LHFitter::kVerbosity::DEBUG);
    fake_mu.setTagAndProbePath(tp_path);
    fake_mu.setInputHistPath(input_path);
    fake_mu.m_doSubtraction = true;
    // REBINNING
    fake_mu.m_doRebinning = true;
    fake_mu.m_doRebinning = true;
    //fake_mu.setBinGrouping(2);
    double fake_mu_new_bins[6] = {10.0,15.0,20.0,25.0,35.0,200.0};
    fake_mu.setVariableBins( fake_mu_new_bins, 5 );
    fake_mu.initialise();
    fake_mu.fit();

//*/
    std::cout << "End of program!" << std::endl;
    return 0;

}
