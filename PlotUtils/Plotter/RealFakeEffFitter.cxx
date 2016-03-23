/**
 * @file   RealFakeEffFitter.cxx
 * @author Marco Milesi <marco.milesi@cern.ch>
 * @date   10 March 2016
 * @brief  ROOT macro to measure real/fake efficiencies via a binned maximum likelihood fit.
 *
 * Use ROOT::TMinuit class
 * @see https://root.cern.ch/doc/master/classTMinuit.html
 *
 */

// TO DO:
/*

Use clever choice for initial parameters of the fit:

-) for r/f efficiencies: read the ones we already have

-) for RR,RF yields: use the observed TT,TL,LL events, and use the inverted matrix w/ the r/f eff. we already have

*/

#include <iostream>
#include <string>
#include <vector>
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
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


bool g_debug(true);
bool g_verbose(false);

bool g_doRebinning(true);
int g_nBINS = 1;//6;
//double g_BINS[7] = {10.0,15.0,20.0,25.0,35.0,50.0,150.0};
double g_BINS[2] = {10.0,150.0};

bool g_useUOflow(false);

std::vector<TH1D*> g_histograms_pt;

std::vector<std::string> g_charge = {"OS","SS"};
std::vector<std::string> g_flavour = {"ElEl","MuMu"};
std::vector<std::string> g_obs_selection = {"TT","TL","LL"};
std::vector<std::string> g_true_selection = {"RR","RF"/*,"FF"*/};
std::vector<std::string> g_efficiencies = {"r_eff_el","r_eff_mu","f_eff_el","f_eff_mu"};

std::vector<std::string> g_param_names;

// ------------------------------------
// Function to get histograms from file
// ------------------------------------

template<typename T>
T* get_object( TFile& file, const std::string& name ) {
  T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
  if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
  return obj;
}

void getHists(const std::string& input_path, int& nbins ) {

  std::string path("");

  for ( const auto& ch : g_charge ) {

    for ( const auto& fl : g_flavour ) {

      for ( const auto& sl : g_obs_selection ) {

	path = input_path + ch + "_" + fl + "_" + sl + "/" + ch + "_" + fl + "_" + sl + "_LepPt.root";

	if ( g_debug ) { Info("getHists()","Reading histogram from:\t %s", path.c_str() ); }

	TFile *file = TFile::Open(path.c_str());
    	if ( !file->IsOpen() ) {
    	  SysError("getHists()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
    	  exit(-1);
    	}

        TH1D *hist = get_object<TH1D>( *file, "observed" );

	// Do prompt and charge flip subtraction in SS
	//
        if  ( ch == "SS" ) {
          TH1D *hist_to_sub = get_object<TH1D>( *file, "expected" );
	  hist->Add(hist_to_sub, -1.0);
          for ( int ibin(0); ibin < nbins; ++ibin ) {
            if ( hist->GetBinContent(ibin) < 0 ) { hist->SetBinContent(ibin,0.0); }
	  }
	}

	std::string new_name = ch + "_" + fl + "_" + sl;

	if ( g_debug ) { Info("getHists()","Saving histogram w/ name:\t %s", new_name.c_str() ); }

	hist->SetName(new_name.c_str());
	  
	hist->SetDirectory(0);
	
	if ( g_doRebinning ) {
	  if ( g_debug ) { Info("getHists()","Rebinning histograms..."); }
	  TH1D* htemp = dynamic_cast<TH1D*>( hist->Clone() );
	  hist = dynamic_cast<TH1D*>( htemp->Rebin( g_nBINS, new_name.c_str(), g_BINS) );
	}
	
	if ( nbins == -1 ) { nbins = ( g_useUOflow ) ? hist->GetNbinsX()+2 : hist->GetNbinsX(); }
	
	g_histograms_pt.push_back(hist); 
	
      }

    }

  }

}

// -------------------------------------------------
// Define expected events with the matrix equation
// -------------------------------------------------

double getExpected( const std::string& obs_selection, const double& r,  const double& f, const double& RR, const double& RF, const double& FF = 0.0 ) {

  double exp(-1.0);

  if ( obs_selection == "TT" ) {
    
    if ( g_verbose ) { Info("getExpected()","Adding term in likelihood for observed selection:\t %s", obs_selection.c_str() ); }
    
    exp = r * r * RR + 2.0 * r * f * RF + f * f * FF;

  } else if ( obs_selection == "TL" ) {
    
    if ( g_verbose ) { Info("getExpected()","Adding term in likelihood for observed selection:\t %s", obs_selection.c_str() ); }

    exp = 2.0 * r * ( 1 - r ) * RR + 2.0 * ( r * ( 1 - f) + ( 1 - r ) * f ) * RF + 2.0 * f * ( 1 - f ) * FF;

  } else if  ( obs_selection == "LL" ) {
    
    if ( g_verbose ) { Info("getExpected()","Adding term in likelihood for observed selection:\t %s", obs_selection.c_str() ); }

    exp = ( 1 -r ) * ( 1 -r ) * RR + 2.0 * ( 1 -r ) * ( 1 - f ) * RF + ( 1 - f ) * ( 1 - f ) * FF;

  }

  return exp;
}


void getParamIndex( const std::string& input_string, const int& ibin,  int& idx, const std::string& param_id = "" ) {

  std::string search_string("");
  
  if ( param_id == "r_eff" || param_id == "f_eff" ) {

    if ( input_string.find("ElEl") != std::string::npos ) {
      search_string = param_id + "_el_pT";
    } else if ( input_string.find("MuMu") != std::string::npos ) {
      search_string = param_id + "_mu_pT";
    }

  } else {
    
    search_string = input_string;
    
    // Replace the last two characters of the input string
    //
    search_string.replace( search_string.end()-2, search_string.end(), param_id);
  }

  // Add bin index to name
  //
  search_string += ( "_" + std::to_string(ibin) );

  // Get corresponding index in list of parameters
  //
  auto it = std::find( g_param_names.begin(), g_param_names.end(), search_string );
  idx = std::distance(g_param_names.begin(),it);

  if ( g_verbose ) { Info("getParamIndex()","Index of parameter with name:\t %s --> %i", search_string.c_str(), idx ); }

}

// -------------------------------------------------
// Define global likelihood function to be minimised
// -------------------------------------------------

void myLikelihood( int& nDim, double* gout, double& result, double par[], int flg ) {

   double likelihood(-999.0);
   double obs(-1.0);

   int idx_RR(-1);
   int idx_RF(-1);
   int idx_FR(-1);
   //int idx_FF(-1);
   int idx_reff(-1);
   int idx_feff(-1);

   double exp(-1);

   for ( auto hist : g_histograms_pt ) {

     if ( g_verbose ) { Info("myLikelihood()","Getting observed yield from histogram:\t %s", hist->GetName() ); }

     std::string this_hist_name = hist->GetName();

     size_t pos = std::distance(this_hist_name.begin(),this_hist_name.end()-2);
     std::string obs_selection = this_hist_name.substr(pos);
     if ( g_verbose ) { Info("myLikelihood()","Observed selection:\t %s", obs_selection.c_str() ); }

     static int nbins = ( g_useUOflow ) ? hist->GetNbinsX()+2 : hist->GetNbinsX()+1;
     int ibin  = ( g_useUOflow ) ? 0 : 1;
     for ( ; ibin < nbins; ++ibin ) {

       obs = hist->GetBinContent(ibin);
       
       int bin_idx = ( g_useUOflow ) ? ibin : ibin - 1;
              
       if ( g_verbose ) { Info("myLikelihood()","\t bin %i (pT = [%.2f,%.2f]) - yield = %f.", bin_idx, hist->GetBinLowEdge(ibin), hist->GetBinLowEdge(ibin+1), obs ); }

       // Need to find correspondent set of parameters
       //
       getParamIndex( this_hist_name, bin_idx, idx_RR, "RR" );
       getParamIndex( this_hist_name, bin_idx, idx_RF, "RF" );
       //getParamIndex( this_hist_name, bin_idx, idx_FF, "FF" );

       getParamIndex( this_hist_name, bin_idx, idx_reff, "r_eff" );
       getParamIndex( this_hist_name, bin_idx, idx_feff, "f_eff" );

       exp = getExpected(obs_selection, par[idx_reff], par[idx_feff], par[idx_RR], par[idx_RF] /*, par[idx_FF] */ );

       likelihood += ( obs * log( exp ) - ( exp ) );
     }

   }

   result = likelihood;
}

// ------------------------------------------------------
// Set the free parameteres to the TMinuit object,
// specifying physical ranges and initial value/step size
// ------------------------------------------------------

void setParameters( TMinuit *myFitter, const int& nbins ) {

   int idx_r_el_pt(0);
   int idx_r_mu_pt(1);
   int idx_f_el_pt(2);
   int idx_f_mu_pt(3);

   double start(-1.0);
   double step(1e-6);
   double up(-1.0);
   double dn(-1.0);
   int ierflg(0);

   for ( int ibin(0); ibin < nbins; ++ibin ) {

     std::string r_eff_el_pt = "r_eff_el_pT_" + std::to_string(ibin);
     std::string r_eff_mu_pt = "r_eff_mu_pT_" + std::to_string(ibin);
     std::string f_eff_el_pt = "f_eff_el_pT_" + std::to_string(ibin);
     std::string f_eff_mu_pt = "f_eff_mu_pT_" + std::to_string(ibin);

     g_param_names.push_back(r_eff_el_pt);
     g_param_names.push_back(r_eff_mu_pt);
     g_param_names.push_back(f_eff_el_pt);
     g_param_names.push_back(f_eff_mu_pt);

   }

   for ( const auto& ch : g_charge ) {
     for ( const auto& fl : g_flavour ) {
       for ( const auto& sl : g_true_selection ) {
	 for ( int ibin = 0; ibin < nbins; ++ibin ) {
	   std::string param_name = ch + "_" + fl + "_" + sl + "_" + std::to_string(ibin);
	   g_param_names.push_back(param_name);
         }
       }
     }
   }

   if ( g_debug ) { Info("setParameters()","Adding parameters to fit function..."); }
   int ipar(0);
   for ( auto& par : g_param_names ) {
   
     bool isEff = ( par.find("r_eff_") != std::string::npos || par.find("f_eff_") != std::string::npos );
     
     start = isEff ? 0.001 : 1.001;
     up    = isEff ? 0.999 : 1e4;
     dn    = 0.0;
     
     myFitter->mnparm( ipar, par.c_str(), start, step, dn, up, ierflg );
     ++ipar;
   }

}

void fitEff( std::string input_path = "" ) {

   TH1::SetDefaultSumw2(kTRUE);

   // Read the input histograms, and set the number of bins
   //
   int nbins(-1);   
   getHists( input_path, nbins );

   // Set the number of parameteres of the fit
   //
   int NFLAV   = g_flavour.size();
   int NCHARGE = g_charge.size();
   int NCOMP   = g_true_selection.size();
   int NEFF    = g_efficiencies.size();
   
   // Total number of parameters to be estimated in the fit
   //
   const int NPAR = (nbins * NFLAV * NCHARGE * NCOMP) + (nbins * NEFF);

   Info("fitEff()","\n\n Number of free parameters in fit ===> %i\n\n", NPAR);
  
   double arglist[NPAR];
   int ierflg = 0;

   // Create the TMinuit object
   //
   TMinuit *myFitter = new TMinuit(NPAR);

   // Set the fitting function
   //
   myFitter->SetFCN(myLikelihood);

   // Set the parameters of the fit
   //
   setParameters(myFitter, nbins);

   //arglist[0] = 1;
   //myFitter->mnexcm("SET ERR",arglist,1,ierflg);

   arglist[0] = 100; // maximum number of iterations
   myFitter->mnexcm("MIGRAD",arglist,0,ierflg);

   if ( !myFitter->fCstatu.Contains("CONVERGED") ) {
     Error("fitEff()","No convergence at fitting! Minuit return string: %s", myFitter->fCstatu.Data() );
     exit(-1);;
   }

   // Print fit statistic
   //
   double best_min;  // the best function value found so far
   double est_vdist; // the estimated vertical distance remaining to minimum
   double err_def;   // the value of UP defining parameter uncertainties
   int nvpar;        // the number of currently variable parameters
   int nparx;        // the highest (external) parameter number defined by user
   int icstat;       // a status integer indicating how good is the covariance matrix:
                     //  0= not calculated at all
  		     //  1= approximation only, not accurate
  		     //  2= full matrix, but forced positive-definite
  		     //  3= full accurate covariance matrix

   myFitter->mnstat(best_min,est_vdist,err_def,nvpar,nparx,icstat);

   Info("fitEff()","************************************************" );
   Info("fitEff()","" );
   
   switch (icstat) {
    case 0 : Error("fitEff()","No covariance matrix was calculated!Exiting..." );
             exit(-1);
             break;       // and exits the switch
    case 1 : Warning("fitEff()","An approximated covariance matrix was calculated!Not accurate..." );
             break; 
    case 2 : Warning("fitEff()","Full covariance matrix was calculated, but forced to be positive-definite..." );
             break; 
    case 3 : Info("fitEff()","Full covariance matrix was calculated :)" );
             break; 
   }

   Info("fitEff()","Minimum of function = %f", best_min );
   Info("fitEff()","Estimated vert. distance to min. = %f", est_vdist ); 
   Info("fitEff()","Value of UP defining parameter unc. = %f", err_def ); 
   Info("fitEff()","Number of variable parameters = %i", nvpar );
   Info("fitEff()","Highest number of parameters defined by user = %i", nparx );
   Info("fitEff()","" );
   Info("fitEff()","************************************************" );

/*
   TH1D* rate =new TH1D("Likelihood",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);

   double pvalue, perror, plow,phigh;
   TString pname;

   int lu=0;

   if( amin < 99999) {
     for (int i=0;i<max;i++)
       {
         lu=i;
         if ((i>=Crack-1)) lu=i+1;

         pMinuit->mnpout(i,pname,pvalue,perror,plow,phigh,ierflg);
         pRate->SetBinContent(lu+1,pvalue);
         pRate->SetBinError(lu+1,perror);

       }
   }
*/

   delete myFitter; myFitter = nullptr;

}
