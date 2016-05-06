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

bool g_debug(true);
bool g_verbose(false);

std::vector<std::string> g_charge = {"OS"/*,"SS"*/};
std::vector<std::string> g_flavour = {"MuMu"};//{"ElEl","MuMu"};
std::vector<std::string> g_obs_selection = {"TT","TL","LT","LL"};
std::vector<std::string> g_true_selection = {"RR",/*"RF","FR","FF"*/};
std::vector<std::string> g_efficiencies = {"r",/*"f"*/};

bool doMuonEff(false);
bool doElectronEff(false);

bool found_par( const std::string& name ) {
  return ( std::find( g_true_selection.begin(), g_true_selection.end(), name ) != g_true_selection.end() ||
           std::find( g_efficiencies.begin(), g_efficiencies.end(), name )     != g_efficiencies.end()
	 );
}

int nPtBins_Linear(6); // including underflow and overflow
int nPtBins_Squared = nPtBins_Linear * nPtBins_Linear; // including underflow and overflow

// Do MC subtraction?
//
bool doSubtraction(false);

bool g_doRebinning(false);
int  g_nBinGroup(4);//6;

// observed quantities
//
std::vector<double> TT(nPtBins_Squared,0.0);
std::vector<double> TL(nPtBins_Squared,0.0);
std::vector<double> LT(nPtBins_Squared,0.0);
std::vector<double> LL(nPtBins_Squared,0.0);

std::vector<double> TT_resized;
std::vector<double> TL_resized;
std::vector<double> LT_resized;
std::vector<double> LL_resized;

// Estimated parameters' indexes
// These vectors contain the corresponding index in the total par[] vector of the TMinuit object
// --> can use these quantities in the likelihood function!
//
std::vector<int> r_idxs;
std::vector<int> f_idxs;
std::vector<int> RR_idxs;
std::vector<int> FR_idxs;
std::vector<int> RF_idxs;
std::vector<int> FF_idxs;
std::vector<int> r1_idxs;  // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
std::vector<int> r2_idxs;  // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin
std::vector<int> f1_idxs;  // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
std::vector<int> f2_idxs;  // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin

// Initial values for parameters
//
std::vector<double> r_init; double r_init_avg(0.0);
std::vector<double> f_init; double f_init_avg(0.0);
std::vector<double> RR_init;
std::vector<double> RF_init;
std::vector<double> FR_init;
std::vector<double> FF_init;

// -----------------------------------------
// Function to get number of bins below
// (and including) diagonal in N X N 2D hist
// -----------------------------------------

int areaGrid( const int& n ) {
  int area(0);
  for ( int i = 0; i < n; ++i ) {
    area += n - i;
  }
  return area;
}

// Final parameter values
//
std::vector<double> r_vals( nPtBins_Linear - 1, 0.0 );
std::vector<double> f_vals( nPtBins_Linear - 1, 0.0 );
std::vector<double> RR_vals( areaGrid(nPtBins_Linear-1), 0.0 );
std::vector<double> RF_vals( areaGrid(nPtBins_Linear-1), 0.0 );
std::vector<double> FR_vals( areaGrid(nPtBins_Linear-1), 0.0 );
std::vector<double> FF_vals( areaGrid(nPtBins_Linear-1), 0.0 );

// Final parameter asymmetric errors
//
// Values in the tuple correspond to:
//
// 0) eplus ( MINOS error UP )
// 1) eminus ( MINOS error DN )
// 2) eparab ( 'parabolic' error (from error matrix) )
//
std::vector<std::tuple<double,double,double> > r_errs( nPtBins_Linear - 1, std::make_tuple( 0.0, 0.0, 0.0 ) );
std::vector<std::tuple<double,double,double> > f_errs( nPtBins_Linear - 1, std::make_tuple( 0.0, 0.0, 0.0 ) );
std::vector<std::tuple<double,double,double> > RR_errs( areaGrid(nPtBins_Linear-1), std::make_tuple( 0.0, 0.0, 0.0 ) );
std::vector<std::tuple<double,double,double> > RF_errs( areaGrid(nPtBins_Linear-1), std::make_tuple( 0.0, 0.0, 0.0 ) );
std::vector<std::tuple<double,double,double> > FR_errs( areaGrid(nPtBins_Linear-1), std::make_tuple( 0.0, 0.0, 0.0 ) );
std::vector<std::tuple<double,double,double> > FF_errs( areaGrid(nPtBins_Linear-1), std::make_tuple( 0.0, 0.0, 0.0 ) );

// ------------------------------------
// Function to get histograms from file
// ------------------------------------

template<typename T>
T* get_object( TFile& file, const std::string& name ) {
    T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
    if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
    return obj;
}

// ----------------------------------------
// Function to print content of a container
// ----------------------------------------

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

// --------------------------------------------------
// Function to print content of a container of errors
// --------------------------------------------------

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

// ----------------------------------------------
// Get tag and probe efficiencies from input file
// ----------------------------------------------

const std::string g_tag_and_probe_path("../OutputPlots_MMRates_25ns_v7_FinalSelection_LHComparison/Rates_NoSub_LHInput/");

void readTagAndProbeEff( const std::string& input_path ) {

    std::string path = input_path + "Rates.root";
    std::string path_avg = input_path + "AvgRates.root";

    TFile *file = TFile::Open(path.c_str());
    if ( !file->IsOpen() ) {
      SysError("readTagAndProbeEff()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
      exit(-1);
    }
    ///*
    TFile *file_avg = TFile::Open(path_avg.c_str());
    if ( !file_avg->IsOpen() ) {
      SysError("readTagAndProbeEff()", "Failed to open ROOT file from path: %s . Aborting", path_avg.c_str() );
      exit(-1);
    }
    //*/
    std::string filename_r(""), filename_r_avg("");
    std::string filename_f(""), filename_f_avg("");

    if ( doElectronEff ) {
      filename_r = filename_r_avg = "El_ProbePt_Real_Efficiency_observed";
      filename_f = filename_f_avg = "El_ProbePt_Fake_Efficiency_observed";
    }
    if ( doMuonEff ) {
      filename_r = filename_r_avg = "Mu_ProbePt_Real_Efficiency_observed";
      filename_f = filename_f_avg = "Mu_ProbePt_Fake_Efficiency_observed";
    }

    TH1D *reff = get_object<TH1D>( *file, filename_r.c_str() );
    TH1D *feff = get_object<TH1D>( *file, filename_f.c_str() );

    // Do not read underflow, but do read overflow!
    //
    for ( int ibin(1); ibin <= reff->GetNbinsX()+1; ++ibin ) {
      r_init.push_back( reff->GetBinContent(ibin) );
    }
    for ( int ibin(1); ibin <= feff->GetNbinsX()+1; ++ibin ) {
      f_init.push_back( feff->GetBinContent(ibin) );
    }

    // Save average efficiencies
    ///*
    TH1D *reff_avg = get_object<TH1D>( *file_avg, filename_r_avg.c_str() );
    TH1D *feff_avg = get_object<TH1D>( *file_avg, filename_f_avg.c_str() );

    r_init_avg = reff_avg->GetBinContent(1);
    f_init_avg = feff_avg->GetBinContent(1);
    //*/
}

// -----------------------------------------
// Get input histograms from ROOT file
// -----------------------------------------

std::vector<TH2D*> g_histograms;

void getHists(const std::string& input_path ) {

  std::string path("");

  for ( const auto& ch : g_charge ) {

    for ( const auto& fl : g_flavour ) {

      for ( const auto& sl : g_obs_selection ) {

        path = input_path + ch + "_" + fl + "_" + sl + "/" + ch + "_" + fl + "_" + sl + "_Lep0Pt_VS_Lep1Pt.root";

        if ( g_debug ) { Info("getHists()","Reading histogram from:\t %s", path.c_str() ); }

        TFile *file = TFile::Open(path.c_str());
        if ( !file->IsOpen() ) {
          SysError("getHists()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
          exit(-1);
        }

        TH2D *hist = get_object<TH2D>( *file, "observed" );

        // -) Do ( !prompt & charge flip ) subtraction in OS
        // -) Do ( prompt & charge flip ) subtraction in SS
        //
        if  ( doSubtraction ) {

	  TH2D *hist_to_sub = get_object<TH2D>( *file, "expected" );
          hist->Add(hist_to_sub, -1.0);

          // Set bin content to 0 if subtraction gives negative yield
	  //
          for ( int ibinx(0); ibinx < nPtBins_Linear; ++ibinx ) {
            for ( int ibiny(0); ibiny < nPtBins_Linear; ++ibiny ) {
                if ( hist->GetBinContent( hist->GetBin(ibinx,ibiny)  ) < 0 ) { hist->SetBinContent( hist->GetBin(ibinx,ibiny), 0.0); }
            }
          }
        }

        std::string new_name = ch + "_" + fl + "_" + sl;

        if ( g_debug ) { Info("getHists()","Storing histogram w/ name:\t %s for later use", new_name.c_str() ); }

        hist->SetName(new_name.c_str());

        hist->SetDirectory(0);

        if ( g_doRebinning ) {

	  if ( g_debug ) { Info("getHists()","Rebinning histogram..."); }

	  TH2D* htemp = dynamic_cast<TH2D*>( hist->Clone() );
          hist = dynamic_cast<TH2D*>( htemp->Rebin2D( g_nBinGroup, g_nBinGroup, new_name.c_str() ) );

	  if ( g_debug ) { Info("getHists()","Re-setting total number of bins..."); }

          // TO DO
	}

        g_histograms.push_back(hist);

      }

    }

  }

}

// -----------------------------------------
// Read the TT, TL, LT, LL events bin-by-bin
// from input histograms
// -----------------------------------------

void readObsYields() {

    if ( g_debug ) { Info("readObsYields()","Reading observed histograms:"); }

    for ( auto hist : g_histograms ) {

	// To get overflow as well
	//
	int firstxbin = 1;
	int firstybin = 1;
	int lastxbin  = nPtBins_Linear - 1;
	int lastybin  = nPtBins_Linear - 1;

	if ( g_debug ) { std::cout << "\t" << hist->GetName() << " - Integral: " << hist->Integral(firstxbin,lastxbin,firstybin,lastybin) << std::endl; }

        std::string TT_name("");
        std::string TL_name("");
        std::string LT_name("");
        std::string LL_name("");
        if ( std::find( g_flavour.begin(), g_flavour.end(), "ElEl" ) != g_flavour.end() ) {
	  TT_name = "OS_ElEl_TT";
	  TL_name = "OS_ElEl_TL";
	  LT_name = "OS_ElEl_LT";
	  LL_name = "OS_ElEl_LL";
	}
        if ( std::find( g_flavour.begin(), g_flavour.end(), "MuMu" ) != g_flavour.end() ) {
	  TT_name = "OS_MuMu_TT";
	  TL_name = "OS_MuMu_TL";
	  LT_name = "OS_MuMu_LT";
	  LL_name = "OS_MuMu_LL";
	}

	if ( strcmp( hist->GetName(), TT_name.c_str() ) == 0 ) {
	    for ( int ibiny(0); ibiny <= nPtBins_Linear - 1; ++ibiny ) {
		for ( int ibinx(0); ibinx <= nPtBins_Linear - 1; ++ibinx ) {
		    if ( g_debug ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
		    TT.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
		}
	    }
	} else if ( strcmp( hist->GetName(), TL_name.c_str() ) == 0 ) {
	    for ( int ibiny(0); ibiny <= nPtBins_Linear - 1; ++ibiny ) {
		for ( int ibinx(0); ibinx <= nPtBins_Linear - 1; ++ibinx ) {
		    if ( g_debug ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
		    TL.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
		}
	    }
	} else if ( strcmp( hist->GetName(), LT_name.c_str() ) == 0 ) {
	    for ( int ibiny(0); ibiny <= nPtBins_Linear - 1; ++ibiny ) {
		for ( int ibinx(0); ibinx <= nPtBins_Linear - 1; ++ibinx ) {
		    if ( g_debug ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
		    LT.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
		}
	    }
	} else if ( strcmp( hist->GetName(), LL_name.c_str() ) == 0 ) {
	    for ( int ibiny(0); ibiny <= nPtBins_Linear - 1; ++ibiny ) {
		for ( int ibinx(0); ibinx <= nPtBins_Linear - 1; ++ibinx ) {
		    if ( g_debug ) { std::cout << "(" << ibinx << "," << ibiny << ") - global bin: " << hist->GetBin(ibinx,ibiny) << std::endl; }
		    LL.at( hist->GetBin(ibinx,ibiny) )  = hist->GetBinContent( hist->GetBin(ibinx,ibiny) );
		}
	    }
	}

    }

}

// -----------------------------------------
// Find initial guesses for parameters
// Use inputs from tag-and-probe measurement
// -----------------------------------------

void getEducatedGuess( ) {

    // TEMP - hardcode efficiency
    //double r =  0.65;
    //double f = 0.172;

    //for ( int idx(0); idx < nPtBins_Linear-1; ++idx ) {
    //  r_init.push_back(r);
    //  f_init.push_back(f);
    //}

    // Read tag-and-probe efficiencies from input file
    //
    readTagAndProbeEff(g_tag_and_probe_path);

    // loop over 2D pT histogram, and depending on the pT value of the leading/subleading lepton,
    // read the corresponding tag-and-probe efficiency, read the TT, TL...yields per bin (set previously in the vectors),
    // and finally compute nRR, nRF...in every bin

    // Skip the (0,j) and (i,0) bins --> they are underflow
    // Skip the (i,j) bins where i < j (empty by construction)

    double r1(-1.0), r2(-1.0);
    double f1(-1.0), f2(-1.0);

    double nTT(0.0), nTL(0.0), nLT(0.0), nLL(0.0);

    double RR_default(1e-3), RF_default(1e-3), FR_default(1e-3), FF_default(1e-3);

    if ( doMuonEff ) {
	std::cout << "" << std::endl;
	Info("getEducatedGuess()","Measuring MUON efficiency! ===> using average f as input..." );
	Info("getEducatedGuess()","f_init_avg = %.2f", f_init_avg );
    }

    int glob_bin_idx(-1);
    for ( int ibiny(1); ibiny <= nPtBins_Linear - 1; ++ibiny ) {

	r2 = r_init.at(ibiny-1);
	f2 = ( !doMuonEff ) ? f_init.at(ibiny-1) : f_init_avg;

	for ( int ibinx(1); ibinx <= nPtBins_Linear - 1; ++ibinx ) {

	    r1 = r_init.at(ibinx-1);
	    f1 = ( !doMuonEff ) ? f_init.at(ibinx-1) : f_init_avg;

            // Skip the above-diagonal elements
	    //
            if ( ibiny > ibinx ) { continue; }

	    glob_bin_idx = g_histograms.at(0)->GetBin(ibinx,ibiny); // can use any of the input histograms: take the first by default

	    nTT = TT.at(glob_bin_idx);
	    nTL = TL.at(glob_bin_idx);
	    nLT = LT.at(glob_bin_idx);
	    nLL = LL.at(glob_bin_idx);

            /*
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
            */

	    double alpha = ( r1 - f1 ) * ( r2 - f2 );

	    // TEMP: how to deal with r = f case???
	    //
	    //if ( alpha == 0.0 ) alpha = 0.01;

	    RR_init.push_back( ( 1.0 / alpha ) * ( ( 1 - f1 ) * ( 1 - f2 ) * nTT + ( f1 - 1 ) * f2 * nTL + ( f2 - 1 ) * f1 * nLT + f1 * f2 * nLL ) );
	    RF_init.push_back( ( 1.0 / alpha ) * ( ( f1 - 1 ) * ( 1 - r2 ) * nTT + ( 1 - f1 ) * r2 * nTL + ( 1 - r2 ) * f1 * nLT - f1 * r2 * nLL ) );
	    FR_init.push_back( ( 1.0 / alpha ) * ( ( r1 - 1 ) * ( 1 - f2 ) * nTT + ( 1 - r1 ) * f2 * nTL + ( 1 - f2 ) * r1 * nLT - r1 * f2 * nLL ) );
	    FF_init.push_back( ( 1.0 / alpha ) * ( ( 1 - r1 ) * ( 1 - r2 ) * nTT + ( r1 - 1 ) * r2 * nTL + ( r2 - 1 ) * r1 * nLT + r1 * r2 * nLL ) );

	    // Copy only relevant bins for observed yields into these containers
	    // (will be the ones used in the likelihood)
	    //
	    TT_resized.push_back( nTT );
	    TL_resized.push_back( nTL );
	    LT_resized.push_back( nLT );
	    LL_resized.push_back( nLL );
	}
    }

}

// ------------------------------------------------------
// Set the free parameteres to the TMinuit object,
// specifying physical ranges and initial value/step size
// ------------------------------------------------------

void setParameters( TMinuit *myFitter ) {

  int ierflg(0);

  // Set parameters for r,f efficiencies. This is just a 1D loop over pT bins
  //
  std::string param_name("");
  int offset(0); // needed for setting parameter offsets (TMinuit can accept only one huge array of parameters...)
  int param_idx(0);

  double step_eff(1e-3);
  double up_eff(0.999);
  double dn_eff(0.0);

  // --------------------
  // Set parameters for r
  // --------------------

  if ( found_par("r") ) {

    for ( auto ibin(0); ibin < r_init.size(); ++ibin ) {
      param_name = "real efficiency - bin [" + std::to_string(ibin) + "]";
      param_idx  = offset + ibin;
      myFitter->mnparm( param_idx, param_name.c_str(), r_init.at(ibin), step_eff, dn_eff, up_eff, ierflg );
      r_idxs.push_back( param_idx );
    }

    // A trick: save the index of the bin for the real efficiency for every (lead, sublead) lepton pair
    // It will contain the same indexes of r_idxs, but each one repeated nPtBins_Linear times
    //
    for ( int ibiny = 0; ibiny < nPtBins_Linear - 1; ++ibiny ) {
	for ( int ibinx = 0; ibinx < nPtBins_Linear - 1; ++ibinx ) {
            if ( ibiny > ibinx ) { continue; }
    	    r1_idxs.push_back( ibinx );
    	    r2_idxs.push_back( ibiny );
	}
    }
     // Set the offset for total parameter index
    //
    offset += nPtBins_Linear - 1;
  }

  // --------------------
  // Set parameters for f
  // --------------------

  if ( found_par("f") ) {

    for ( auto ibin(0); ibin < f_init.size(); ++ibin ) {
      param_name = "fake efficiency - bin [" + std::to_string(ibin) + "]";
      param_idx  = offset + ibin;
      myFitter->mnparm( param_idx, param_name.c_str(), f_init.at(ibin), step_eff, dn_eff, up_eff, ierflg );
      f_idxs.push_back( param_idx );
    }

    // A trick: save the index of the bin for the real efficiency for every (lead, sublead) lepton pair
    // It will contain the same indexes of f_idxs, but each one repeated nPtBins_Linear times
    //
    for ( int ibiny = 0; ibiny < nPtBins_Linear - 1; ++ibiny ) {
	for ( int ibinx = 0; ibinx < nPtBins_Linear - 1; ++ibinx ) {
            if ( ibiny > ibinx ) { continue; }
    	    f1_idxs.push_back( ibinx );
    	    f2_idxs.push_back( ibiny );
	}
    }
     // Set the offset for total parameter index
    //
    offset += nPtBins_Linear - 1;
  }

  double step_yields(1e-3);
  double up_yields(1e6);
  double dn_yields(0.0);

  // ---------------------
  // Set parameters for RR
  // ---------------------

  if ( found_par("RR") ) {

    for ( auto ibin(0); ibin < RR_init.size(); ++ibin ) {
    	param_name = "RR - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	myFitter->mnparm( param_idx, param_name.c_str(), RR_init.at(ibin), step_yields, dn_yields, up_yields, ierflg );
    	RR_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }
  // ---------------------
  // Set parameters for RF
  // ---------------------

  if ( found_par("RF") ) {

    for ( auto ibin(0); ibin < RF_init.size(); ++ibin ) {
    	param_name = "RF - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	myFitter->mnparm( param_idx, param_name.c_str(), RF_init.at(ibin), step_yields, dn_yields, up_yields, ierflg );
    	RF_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }
  // ---------------------
  // Set parameters for FR
  // ---------------------

  if ( found_par("FR") ) {

    for ( auto ibin(0); ibin < FR_init.size(); ++ibin ) {
    	param_name = "FR - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	myFitter->mnparm( param_idx, param_name.c_str(), FR_init.at(ibin), step_yields, dn_yields, up_yields, ierflg );
    	FR_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }

  // ---------------------
  // Set parameters for FF
  // ---------------------

  if ( found_par("FF") ) {

    for ( auto ibin(0); ibin < FF_init.size(); ++ibin ) {
    	param_name = "FF - bin [" + std::to_string(ibin) + "]";
    	param_idx = offset + ibin; // can use any of the input histograms: take the first by default
    	myFitter->mnparm( param_idx, param_name.c_str(), FF_init.at(ibin), step_yields, dn_yields, up_yields, ierflg );
    	FF_idxs.push_back( param_idx );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }
}

// -------------------------------------------------
// Define global likelihood function to be minimised
// -------------------------------------------------

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
   for ( int ibin = 0; ibin < areaGrid( nPtBins_Linear - 1 ); ++ibin  ) {

     //. Read the likelihood parameters from the TMinuit par[] vector
     //
     r1 = ( found_par("r") )  ? par[ r1_idxs.at(ibin) ] : 0.0;
     r2 = ( found_par("r") )  ? par[ r2_idxs.at(ibin) ] : 0.0;
     f1 = ( found_par("f") )  ? par[ f1_idxs.at(ibin) ] : 0.0;
     f2 = ( found_par("f") )  ? par[ f2_idxs.at(ibin) ] : 0.0;
     RR = ( found_par("RR") ) ? par[ RR_idxs.at(ibin) ] : 0.0;
     RF = ( found_par("RF") ) ? par[ RF_idxs.at(ibin) ] : 0.0;
     FR = ( found_par("FR") ) ? par[ FR_idxs.at(ibin) ] : 0.0;
     FF = ( found_par("FF") ) ? par[ FF_idxs.at(ibin) ] : 0.0; // can we set this to 0 by brute force?

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
     obs_TT  = TT_resized.at(ibin);
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
     obs_TL  = TL_resized.at(ibin);
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
     obs_LT  = LT_resized.at(ibin);
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
     obs_LL  = LL_resized.at(ibin);
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


// ---------------------------------------------------------
// Update the content of error containers w/ info from MINOS
// ---------------------------------------------------------

void getParametersAndErrors( TMinuit *myFitter ) {

  int offset(0); // needed for dealing with parameter offsets (TMinuit can accept only one huge array of parameters...)
  int param_idx(0);

  double globcc;

  // --------------------
  // Get errors for r
  // --------------------

  if ( found_par("r") ) {

    for ( auto idx(0); idx < r_errs.size(); ++idx ) {
       param_idx = offset + idx;
       myFitter->GetParameter( param_idx, r_vals.at(idx), std::get<2>(r_errs.at(idx)) );
       myFitter->mnerrs( param_idx, std::get<0>(r_errs.at(idx)), std::get<1>(r_errs.at(idx)), std::get<2>(r_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += nPtBins_Linear - 1;
  }

  // --------------------
  // Get errors for f
  // --------------------

  if ( found_par("f") ) {

    for ( auto idx(0); idx < f_errs.size(); ++idx ) {
       param_idx = offset + idx;
       myFitter->GetParameter( param_idx, f_vals.at(idx), std::get<2>(f_errs.at(idx)) );
       myFitter->mnerrs( param_idx, std::get<0>(f_errs.at(idx)), std::get<1>(f_errs.at(idx)), std::get<2>(f_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += nPtBins_Linear - 1;
  }

  // ---------------------
  // Get errors for RR
  // ---------------------

  if ( found_par("RR") ) {

    for ( auto idx(0); idx < RR_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        myFitter->GetParameter( param_idx, RR_vals.at(idx), std::get<2>(RR_errs.at(idx)) );
        myFitter->mnerrs( param_idx, std::get<0>(RR_errs.at(idx)), std::get<1>(RR_errs.at(idx)), std::get<2>(RR_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }

  // ---------------------
  // Get errors for RF
  // ---------------------

  if ( found_par("RF") ) {

    for ( auto idx(0); idx < RF_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        myFitter->GetParameter( param_idx, RF_vals.at(idx), std::get<2>(RF_errs.at(idx)) );
        myFitter->mnerrs( param_idx, std::get<0>(RF_errs.at(idx)), std::get<1>(RF_errs.at(idx)), std::get<2>(RF_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }

  // ---------------------
  // Get errors for FR
  // ---------------------

  if ( found_par("FR") ) {

    for ( auto idx(0); idx < FR_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        myFitter->GetParameter( param_idx, FR_vals.at(idx), std::get<2>(FR_errs.at(idx)) );
        myFitter->mnerrs( param_idx, std::get<0>(FR_errs.at(idx)), std::get<1>(FR_errs.at(idx)), std::get<2>(FR_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }

  // ---------------------
  // Get errors for FF
  // ---------------------

  if ( found_par("FF") ) {

    for ( auto idx(0); idx < FF_errs.size(); ++idx ) {
    	param_idx = offset + idx;
        myFitter->GetParameter( param_idx, FF_vals.at(idx), std::get<2>(FF_errs.at(idx)) );
        myFitter->mnerrs( param_idx, std::get<0>(FF_errs.at(idx)), std::get<1>(FF_errs.at(idx)), std::get<2>(FF_errs.at(idx)), globcc );
    }
    // Set the offset for total parameter index
    //
    offset += areaGrid(nPtBins_Linear-1);
  }

}

// ------------------------------------------------------
// Save efficiencies and their errors in a ROOT file
// ------------------------------------------------------

void  saveEfficiencies() {

  //TH1::ResetBit(TH1::kCanRebin);

  std::string outfilename("LH_efficiencies");
  if ( doMuonEff )     { outfilename += "_mu"; }
  if ( doElectronEff ) { outfilename += "_el"; }

  std::string rootfilename = outfilename + ".root";
  TFile outfile(rootfilename.c_str(),"RECREATE");

  std::string txtfilename = outfilename + ".txt";
  std::ofstream outtextfile(txtfilename.c_str(), std::ios::trunc);
  outtextfile << "Efficiencies for FF amd Matrix Method - LH fit \n";

  // --------------------
  // real efficiency
  // --------------------

  if ( found_par("r") ) {

    TH1D *r_hist = new TH1D( "r_hist", "real efficiency", nPtBins_Linear - 1, 10, 135 );
    r_hist->GetYaxis()->SetTitle("Real efficiency");
    r_hist->GetYaxis()->SetRangeUser(0.0,1.0);

    std::string xtitle("");
    if ( doMuonEff )     { xtitle += "muon pT [GeV]"; }
    if ( doElectronEff ) { xtitle += "electron pT [GeV]"; }

    r_hist->GetXaxis()->SetTitle(xtitle.c_str());

    outtextfile << "Real efficiency - " << xtitle << "\n";
    for ( auto idx(0); idx < r_vals.size(); ++idx ) {
       r_hist->SetBinContent( idx + 1, r_vals.at(idx) );
       r_hist->SetBinError( idx + 1, std::get<2>(r_errs.at(idx)) );
       outtextfile << "{ Bin nr: " << idx << ", efficiency = " <<  r_vals.at(idx) << " + " << std::get<0>(r_errs.at(idx)) << " - " << std::get<1>(r_errs.at(idx)) << " }\n";
    }
    r_hist->Write();
  }

  // --------------------
  // fake efficiency
  // --------------------

  if ( found_par("f") ) {

    TH1D *f_hist = new TH1D( "f_hist", "fake efficiency", nPtBins_Linear - 1, 10, 135 );
    f_hist->GetYaxis()->SetTitle("Fake efficiency");
    f_hist->GetYaxis()->SetRangeUser(0.0,1.0);

    std::string xtitle("");
    if ( doMuonEff )     { xtitle += "electron pT [GeV]"; }
    if ( doElectronEff ) { xtitle += "muon pT [GeV]"; }

    f_hist->GetXaxis()->SetTitle(xtitle.c_str());

    for ( auto idx(0); idx < f_vals.size(); ++idx ) {
       f_hist->SetBinContent( idx + 1, f_vals.at(idx) );
       f_hist->SetBinError( idx + 1, std::get<2>(f_errs.at(idx)) );
    }
    f_hist->Write();
  }

  outfile.Close();
  outtextfile.close();

}

// ------------------------------------------------------
// The main function
// ------------------------------------------------------

void fitEff( std::string input_path = "" ) {

  TH1::SetDefaultSumw2(kTRUE);

  doMuonEff     = ( std::find( g_flavour.begin(), g_flavour.end(), "MuMu" ) != g_flavour.end() );
  doElectronEff = ( std::find( g_flavour.begin(), g_flavour.end(), "ElEl" ) != g_flavour.end() );

  // Reserve memory for input containers
  //
  TT_resized.reserve(areaGrid(nPtBins_Linear-1));
  TL_resized.reserve(areaGrid(nPtBins_Linear-1));
  LT_resized.reserve(areaGrid(nPtBins_Linear-1));
  LL_resized.reserve(areaGrid(nPtBins_Linear-1));

  // Reserve memory for (initial) parameter containers
  //
  r_init.reserve(nPtBins_Linear-1);
  f_init.reserve(nPtBins_Linear-1);
  RR_init.reserve(areaGrid(nPtBins_Linear-1));
  RF_init.reserve(areaGrid(nPtBins_Linear-1));
  FR_init.reserve(areaGrid(nPtBins_Linear-1));
  FF_init.reserve(areaGrid(nPtBins_Linear-1));

  // Reserve memory for (final) parameter error containers
  //

  // Reserve memory for parameters' index containers
  //
  r_idxs.reserve(nPtBins_Linear-1);
  f_idxs.reserve(nPtBins_Linear-1);
  RR_idxs.reserve(areaGrid(nPtBins_Linear-1));
  FR_idxs.reserve(areaGrid(nPtBins_Linear-1));
  RF_idxs.reserve(areaGrid(nPtBins_Linear-1));
  FF_idxs.reserve(areaGrid(nPtBins_Linear-1));
  r1_idxs.reserve(areaGrid(nPtBins_Linear-1)); // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
  r2_idxs.reserve(areaGrid(nPtBins_Linear-1)); // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin
  f1_idxs.reserve(areaGrid(nPtBins_Linear-1)); // x-coordinate  (--> leading pT axis) of the leading lepton in the (x_i,y_j) bin
  f2_idxs.reserve(areaGrid(nPtBins_Linear-1)); // y-coordinate  (--> subleading pT axis) of the subleading lepton in the (x_i,y_j) bin

  // Get input histograms and store them
  //
  getHists( input_path );

  // Store the observed yields for TT, TL... histograms
  //
  readObsYields();

  if ( g_debug ) {
     printContainer( TT, "Printing content of TT:" );
     printContainer( TL, "Printing content of TL:" );
     printContainer( LT, "Printing content of LT:" );
     printContainer( LL, "Printing content of LL:" );
  }

  // Set the number of parameteres of the fit
  //
  int NFLAV   = g_flavour.size();
  int NCHARGE = g_charge.size();
  int NCOMP   = g_true_selection.size();
  int NEFF    = g_efficiencies.size();

  Info("fitEff()","");
  std::cout << "Number of pT bins = "<<  nPtBins_Linear << std::endl;
  std::cout << "----------------------------"<<  std::endl;
  std::cout << "Number of flavour bins = "<<  NFLAV << std::endl;
  std::cout << "Number of charge bins = "<<  NCHARGE << std::endl;
  std::cout << "Number of RR, RF... bins = "<<  NCOMP << std::endl;
  std::cout << "Number of efficiency bins = "<<  NEFF << std::endl;

  // Total number of parameters to be estimated in the fit
  //
  // 1st factor in sum: YIELDS
  // 2nd factor in sum: EFFICIENCIES
  //
  const int NPAR_YIELDS = areaGrid( nPtBins_Linear - 1 ) * NFLAV * NCHARGE * NCOMP; // get effective "area" of 2D hist (excluding underflows)
  const int NPAR_EFF	= ( nPtBins_Linear - 1 ) * NEFF; // do not take underflow into account

  const int NPAR = NPAR_YIELDS + NPAR_EFF;

  Info("fitEff()","\n\n Number of free parameters in fit ===> %i (efficiencies) + %i (yields) = %i\n\n", NPAR_EFF, NPAR_YIELDS, NPAR);

  double arglist[NPAR];
  int ierflg = 0;

  // Create the TMinuit object
  //
  TMinuit *myFitter = new TMinuit(NPAR);

  // Set the fitting function
  //
  myFitter->SetFCN(myLikelihood);

  // Get an educated guess for the initial parameters, from the tag-and-probe measurement:
  //
  // r(i), f(i)  ( in 1D bins of pT )
  //
  // RR(j,k), RF(j,k), FR(j,k), FF(j,k) ( in 2D bins of lead-sublead pT )
  //
  getEducatedGuess();

  if ( g_debug ) {
     printContainer( r_init, "Printing content of r_init:" );
     printContainer( f_init, "Printing content of f_init:" );
     if ( found_par("RR") ) printContainer( RR_init,"Printing content of RR_init:" );
     if ( found_par("RF") ) printContainer( RF_init,"Printing content of RF_init:" );
     if ( found_par("FR") ) printContainer( FR_init,"Printing content of FR_init:" );
     if ( found_par("FF") ) printContainer( FF_init,"Printing content of FF_init:" );
     printContainer( TT_resized, "Printing content of TT (resized):" );
     printContainer( TL_resized, "Printing content of TL (resized):" );
     printContainer( LT_resized, "Printing content of LT (resized):" );
     printContainer( LL_resized, "Printing content of LL (resized):" );
  }

  // Set the parameters of the fit into the toatl par[] array of the TMinuit object
  // The order will be:
  // -) r  ( x nPtBins_Linear )
  // -) f  ( x nPtBins_Linear )
  // -) RR ( x nPtBins_Squared )
  // -) RF ( x nPtBins_Squared )
  // -) FR ( x nPtBins_Squared )
  // -) FF ( x nPtBins_Squared )
  //
  setParameters(myFitter);

  if ( g_debug ) {
     if ( found_par("r") )  printContainer( r_idxs, "Printing content of r_idxs:" );
     if ( found_par("f") )  printContainer( f_idxs, "Printing content of f_idxs:" );
     if ( found_par("r") )  printContainer( r1_idxs, "Printing content of r1_idxs:" );
     if ( found_par("r") )  printContainer( r2_idxs, "Printing content of r2_idxs:" );
     if ( found_par("f") )  printContainer( f1_idxs, "Printing content of f1_idxs:" );
     if ( found_par("f") )  printContainer( f2_idxs, "Printing content of f2_idxs:" );
     if ( found_par("RR") ) printContainer( RR_idxs, "Printing content of RR_idxs:" );
     if ( found_par("RF") ) printContainer( RF_idxs, "Printing content of RF_idxs:" );
     if ( found_par("FR") ) printContainer( FR_idxs, "Printing content of FR_idxs:" );
     if ( found_par("FF") ) printContainer( FF_idxs, "Printing content of FF_idxs:" );
  }

  // ----------------------------------------------------------------------------------

  arglist[0] = 1000; // maximum number of iterations
  myFitter->mnexcm("MIGRAD",arglist,2/*0*/,ierflg);

  std::cout << "\n\n" << std::endl;

  if ( !myFitter->fCstatu.Contains("CONVERGED") ) {
    Error("fitEff()","No convergence at fitting! Minuit return string: %s", myFitter->fCstatu.Data() );
    exit(-1);;
  }

  // ----------------------------------------------------------------------------------

  arglist[0] = 0.0;
  arglist[1] = 1.0;
  myFitter->mnexcm("MINOS",arglist,1,ierflg); // calculate MINOS errors for all the parameters

  if ( g_debug ) {
     if ( found_par("r") )  printContainer( r_vals, "Printing content of r_vals (DEFAULT):" );
     if ( found_par("r") )  printErrorContainer( r_errs, "Printing content of r_errs (DEFAULT):" );
     if ( found_par("RR") ) printContainer( RR_vals, "Printing content of RR_vals (DEFAULT):" );
     if ( found_par("RR") ) printErrorContainer( RR_errs, "Printing content of RR_errs (DEFAULT):" );
  }

  // Get final parameters, parabolic error and the asymmetric errors from MINOS
  //
  getParametersAndErrors(myFitter);

  if ( g_debug ) {
     if ( found_par("r") )  printContainer( r_vals, "Printing content of r_vals (POST-FIT):" );
     if ( found_par("r") )  printErrorContainer( r_errs, "Printing content of r_errs (POST-FIT):" );
     if ( found_par("RR") ) printContainer( RR_vals, "Printing content of RR_vals (POST-FIT):" );
     if ( found_par("RR") ) printErrorContainer( RR_errs, "Printing content of RR_errs (POST-FIT):" );
  }

  // ----------------------------------------------------------------------------------

  saveEfficiencies();

  // ----------------------------------------------------------------------------------

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

  myFitter->mnstat(best_min,est_vdist,err_def,nvpar,nparx,icstat);

  Info("fitEff()","************************************************" );
  Info("fitEff()","" );

  switch (icstat) {
   case 0 : Error("fitEff()","No covariance matrix was calculated!Exiting..." );
  	    exit(-1);
  	    break;	 // and exits the switch
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

}

// ------------------------------------------------------
// Call this function from command line
// ------------------------------------------------------

int main(int argc, char **argv){

    //gSystem->Load( "libCintex.so" );
    //gROOT->ProcessLine("Cintex::Cintex::Enable();");

    if ( argc<2 ) {
        std::cout << "Missing input parameters: " << std::endl;
        std::cout << "[1] input path" << std::endl;
        return 0;
    }

    const char* input_name = argv[1];
    std::string input(input_name);

    fitEff(input);

    std::cout << "end" << std::endl;
    return 0;

}
