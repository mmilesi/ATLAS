#include <vector>
#include <math.h>
#include <iostream>
#include <sstream>

#include "TError.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TTree.h"

/* ***************
/
/ Global variables
/
/
*************** */

bool g_debug(true);
bool g_verbose(true);

std::map< std::string, TH1D* > g_el_hist_map;
std::map< std::string, TH1D* > g_mu_hist_map;

int g_n_el_bins_pt_rr(0);
int g_n_el_bins_pt_fr(0);
int g_n_mu_bins_pt_rr(0);
int g_n_mu_bins_pt_fr(0);

bool g_ratesEta;
int g_n_el_bins_eta(0);
int g_n_mu_bins_eta(0);

double g_el_rr_tot(1.0);
double g_el_fr_tot(1.0);
double g_mu_rr_tot(1.0);
double g_mu_fr_tot(1.0);

bool g_doClosure;

/* ***********************************
/
/ Function to get histograms from file
/
*********************************** */

TH1D* get_hist( TFile& file, const std::string& name ) {
  TH1D* hist = (TH1D*)( file.Get(name.c_str()) );
  if ( !hist ) {
    SysError("get_hist()", "Histogram %s not found. Aborting", name.c_str() );
    exit(-1);
  }
  return hist;
}

/* **********************************************
/
/ Function to read rates
/ (NB: rates can be read from 2 different files!)
/
********************************************** */

void read_rates(const std::string rr_dir, const std::string fr_dir )
{

  //std::string glob_path("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/");
  std::string glob_path("");

  // ******************************************************
  //
  // 1. 'REAL' efficiency

  Info("read_rates()", "REAL efficiency from directory: %s ", rr_dir.c_str() );

  std::string path_R_el = glob_path + rr_dir + "/Rates.root";// "/AvgRates.root";
  TFile *file_R_el = TFile::Open(path_R_el.c_str());
  if ( !file_R_el->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with R efficiency from path: %s . Aborting", path_R_el.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "ELECTRON REAL efficiency: %s ", path_R_el.c_str() );
  }
  std::string path_R_mu = glob_path + rr_dir + "/Rates.root"; //"/AvgRates.root";
  TFile *file_R_mu = TFile::Open(path_R_mu.c_str());
  if ( !file_R_mu->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with R efficiency from path: %s . Aborting", path_R_mu.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "MUON REAL efficiency: %s ", path_R_mu.c_str() );
  }

  // special case: use charge flip (i.e. "fake" rate) as real rate (only for electrons)!
  bool CFRateAsRR = ( rr_dir.find("ChFlip") != std::string::npos ) ? true : false;
  if ( CFRateAsRR ) { Info("read_rates()", "---> For ELECTRONS, will be treating charge flip efficiency as REAL efficiency!!"); }

  std::string rate_type = ( !g_doClosure ) ? "observed" : "expected";

  // ELECTRONS
  //
  std::string histname_el_pt_rr    = ( !CFRateAsRR ) ?  "El_ProbePt_Real_Efficiency_" + rate_type : "El_ProbePt_Fake_Efficiency_" + rate_type;
  std::string histname_el_eta_rr   = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_Efficiency_" + rate_type : "El_ProbeEta_Fake_Efficiency_" + rate_type;

  std::string histname_el_pt_r_T   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_T_" + rate_type : "El_ProbePt_Fake_T_" + rate_type;
  std::string histname_el_pt_r_L   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_L_" + rate_type : "El_ProbePt_Fake_L_" + rate_type;

  // get real efficiency histograms
  //
  TH1D *hist_el_pt_rr    = get_hist( *file_R_el, histname_el_pt_rr );
  TH1D *hist_el_eta_rr   = ( g_ratesEta ) ? get_hist( *file_R_el, histname_el_eta_rr ) : nullptr;

  TH1D *hist_el_pt_r_T   = get_hist( *file_R_el, histname_el_pt_r_T );
  TH1D *hist_el_pt_r_L   = get_hist( *file_R_el, histname_el_pt_r_L );

  // MUONS
  //
  std::string histname_mu_pt_rr    = "Mu_ProbePt_Real_Efficiency_" + rate_type;
  std::string histname_mu_eta_rr   = "Mu_ProbeEta_Real_Efficiency_" + rate_type;

  std::string histname_mu_pt_r_T   = "Mu_ProbePt_Real_T_" + rate_type;
  std::string histname_mu_pt_r_L   = "Mu_ProbePt_Real_L_" + rate_type;

  // get real efficiency histograms
  //
  TH1D *hist_mu_pt_rr    = get_hist( *file_R_mu, histname_mu_pt_rr );
  TH1D *hist_mu_eta_rr   = ( g_ratesEta ) ? get_hist( *file_R_mu, histname_mu_eta_rr ) : nullptr;

  TH1D *hist_mu_pt_r_T   = get_hist( *file_R_mu, histname_mu_pt_r_T );
  TH1D *hist_mu_pt_r_L   = get_hist( *file_R_mu, histname_mu_pt_r_L );

  // 2. FAKE efficiency

  std::string fake_dir(rr_dir);
  if ( !fr_dir.empty() ) {
     Warning("read_rates()", "FAKE efficiency is going to be read from %s. Check whether it's really what you want...", fr_dir.c_str());
     fake_dir = fr_dir;
  } else {
     Info("read_rates()", "FAKE efficiency from same directory" );
  }

  std::string path_F_el = glob_path + fake_dir + "/Rates.root"; // "/AvgRates.root";
  TFile *file_F_el = TFile::Open(path_F_el.c_str());
  if ( !file_F_el->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with F efficiency from path: %s . Aborting", path_F_el.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "ELECTRON FAKE efficiency: %s ", path_F_el.c_str() );
  }

  std::string path_F_mu = glob_path + fake_dir + "/Rates.root"; // "/AvgRates.root";
  TFile *file_F_mu = TFile::Open(path_F_mu.c_str());
  if ( !file_F_mu->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with F efficiency from path: %s . Aborting", path_F_mu.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "MUON FAKE efficiency: %s ", path_F_mu.c_str() );
  }

  // ELECTRONS
  //
  std::string histname_el_pt_fr(""), histname_el_eta_fr(""), histname_el_pt_f_T(""), histname_el_pt_f_L("");
  if ( g_doClosure ) {
    histname_el_pt_fr    = "El_ProbePt_Fake_Efficiency_" + rate_type;
    histname_el_eta_fr   = "El_ProbeEta_Fake_Efficiency_" + rate_type;
    histname_el_pt_f_T   = "El_ProbePt_Fake_T_" + rate_type;
    histname_el_pt_f_L   = "El_ProbePt_Fake_L_" + rate_type;
  } else {
    //histname_el_pt_fr    = "El_ProbePt_ScaledFake_Efficiency_" + rate_type; // Use the QMisID-eff-scaled fake efficiency for electrons when running MM on DATA
    //histname_el_eta_fr   = "El_ProbeEta_ScaledFake_Efficiency_" + rate_type;
    //histname_el_pt_f_T   = "El_ProbePt_ScaledFake_T_" + rate_type;
    //histname_el_pt_f_L   = "El_ProbePt_ScaledFake_L_" + rate_type;
    histname_el_pt_fr    = "El_ProbePt_Fake_Efficiency_" + rate_type; // Use the QMisID-eff-scaled fake efficiency for electrons when running MM on DATA
    histname_el_eta_fr   = "El_ProbeEta_Fake_Efficiency_" + rate_type;
    histname_el_pt_f_T   = "El_ProbePt_Fake_T_" + rate_type;
    histname_el_pt_f_L   = "El_ProbePt_Fake_L_" + rate_type;
  }

  // get fake efficiency histograms
  //
  TH1D *hist_el_pt_fr    = get_hist( *file_F_el, histname_el_pt_fr );
  TH1D *hist_el_eta_fr   = ( g_ratesEta ) ? get_hist( *file_F_el, histname_el_eta_fr ) : nullptr;

  TH1D *hist_el_pt_f_T   = get_hist( *file_F_el, histname_el_pt_f_T );
  TH1D *hist_el_pt_f_L   = get_hist( *file_F_el, histname_el_pt_f_L );

  // MUONS
  //
  std::string histname_mu_pt_fr    = "Mu_ProbePt_Fake_Efficiency_" + rate_type;
  std::string histname_mu_eta_fr   = "Mu_ProbeEta_Fake_Efficiency_" + rate_type;

  std::string histname_mu_pt_f_T   = "Mu_ProbePt_Fake_T_" + rate_type;
  std::string histname_mu_pt_f_L   = "Mu_ProbePt_Fake_L_" + rate_type;

  // get fake efficiency histograms
  //
  TH1D *hist_mu_pt_fr    = get_hist( *file_F_mu, histname_mu_pt_fr );
  TH1D *hist_mu_eta_fr   = ( g_ratesEta ) ? get_hist( *file_F_mu, histname_mu_eta_fr ) : nullptr;

  TH1D *hist_mu_pt_f_T   = get_hist( *file_F_mu, histname_mu_pt_f_T );
  TH1D *hist_mu_pt_f_L   = get_hist( *file_F_mu, histname_mu_pt_f_L );

  // ***********************************************************************

  // fill a map for later usage
  //
  g_el_hist_map["pt_rr"]   = hist_el_pt_rr;
  g_mu_hist_map["pt_rr"]   = hist_mu_pt_rr;
  g_el_hist_map["pt_fr"]   = hist_el_pt_fr;
  g_mu_hist_map["pt_fr"]   = hist_mu_pt_fr;

  g_el_hist_map["eta_rr"]  = hist_el_eta_rr;
  g_mu_hist_map["eta_rr"]  = hist_mu_eta_rr;
  g_el_hist_map["eta_fr"]  = hist_el_eta_fr;
  g_mu_hist_map["eta_fr"]  = hist_mu_eta_fr;

  // pt hist has two different binning for r/f
  //
  g_n_el_bins_pt_rr =  hist_el_pt_rr->GetNbinsX()+1;
  g_n_el_bins_pt_fr =  hist_el_pt_fr->GetNbinsX()+1;
  g_n_mu_bins_pt_rr =  hist_mu_pt_rr->GetNbinsX()+1;
  g_n_mu_bins_pt_fr =  hist_mu_pt_fr->GetNbinsX()+1;

  // eta hist has same binning for r/f
  //
  if ( g_ratesEta ) {
    g_n_el_bins_eta   =  hist_el_eta_rr->GetNbinsX()+1;
    g_n_mu_bins_eta   =  hist_mu_eta_rr->GetNbinsX()+1;
  }

  // Calculate normalisation factor for (pT X eta) 1D efficiencies case.
  //
  // This factor is the same for eta and pT r/f histograms (it's just Integral(N) / Integral(D) for the efficiency definition ): use pT
  //
  g_el_rr_tot = ( hist_el_pt_r_T->Integral(1,hist_el_pt_r_T->GetNbinsX()+1) ) / ( hist_el_pt_r_L->Integral(1,hist_el_pt_r_L->GetNbinsX()+1) );
  g_el_fr_tot = ( hist_el_pt_f_T->Integral(1,hist_el_pt_f_T->GetNbinsX()+1) ) / ( hist_el_pt_f_L->Integral(1,hist_el_pt_f_L->GetNbinsX()+1) );
  g_mu_rr_tot = ( hist_mu_pt_r_T->Integral(1,hist_mu_pt_r_T->GetNbinsX()+1) ) / ( hist_mu_pt_r_L->Integral(1,hist_mu_pt_r_L->GetNbinsX()+1) );
  g_mu_fr_tot = ( hist_mu_pt_f_T->Integral(1,hist_mu_pt_f_T->GetNbinsX()+1) ) / ( hist_mu_pt_f_L->Integral(1,hist_mu_pt_f_L->GetNbinsX()+1) );

  std::cout << "\n" << std::endl;
  Info("read_rates()", "MUON REAL efficiency - pT histogram name: %s ", histname_mu_pt_rr.c_str() );
  Info("read_rates()", "MUON FAKE efficiency - pT histogram name: %s ", histname_mu_pt_fr.c_str() );
  if ( g_ratesEta ) {
    Info("read_rates()", "MUON REAL efficiency - eta histogram name: %s ", histname_mu_eta_rr.c_str() );
    Info("read_rates()", "MUON FAKE efficiency - eta histogram name: %s ", histname_mu_eta_fr.c_str() );
  }
  std::cout << "            --------------------------------------------" << std::endl;
  Info("read_rates()", "ELECTRON REAL efficiency - pT histogram name: %s ", histname_el_pt_rr.c_str() );
  Info("read_rates()", "ELECTRON FAKE efficiency - pT histogram name: %s ", histname_el_pt_fr.c_str() );
  if ( g_ratesEta ) {
    Info("read_rates()", "ELECTRON REAL efficiency - eta histogram name: %s ", histname_el_eta_rr.c_str() );
    Info("read_rates()", "ELECTRON FAKE efficiency - eta histogram name: %s ", histname_el_eta_fr.c_str() );
  }
  std::cout << "\n" << std::endl;

}

/* **************************
/
/ Just algebraic manipulation
/
************************** */

double  scaleRateToEfficiency( double rate )
{
  if ( rate < 0 ) { rate = 0.0; }

  double eff = ( rate / (rate+1.0) );

  return eff;
}

/* ********************************************************
/
/ Function to calculate r/f weights and their unceratinties
/
******************************************************** */
///*
std::vector<double>  calc_weights( std::map< std::string, TH1D* >& histograms,
				   float pt,
				   float eta,
				   bool isFakeLep,
				   int n_bins_pt_fr,
				   int n_bins_pt_rr,
				   int n_bins_eta,
				   float fr_tot,
				   float rr_tot
                                  )
{
  // As a first thing, convert pT in GeV!
  //
  pt = pt/1e3;

  std::vector<double> weights(3,0.0); //initialized with zeroes

  weights.at(0) = 1.0;
  double error(0.0);

  float this_low_edge_pt(-1.0);
  float this_up_edge_pt(-1.0);

  float this_low_edge_eta(-999.0);
  float this_up_edge_eta(-999.0);

  // case 1) : lepton is fake: choose correct pt histogram
  //
  if ( isFakeLep ) {

    if ( g_verbose ) {
       Info("calc_weights()", "\tReading fake efficiency...");
    }

    // loop over number of pt bins
    // do not consider underflow, i.e. 0th bin
    //
    for ( int p = 1; p <= n_bins_pt_fr; p++ ) {

      this_low_edge_pt = (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p);
      this_up_edge_pt  = (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p+1);

      if ( g_verbose ) { Info("calc_weights()","\t\tpT bin %i : [%.0f,%.0f] GeV", p, this_low_edge_pt, this_up_edge_pt ); }

      if ( pt >= this_low_edge_pt && pt < this_up_edge_pt ) {

    	double fr_pt     = (histograms.find("pt_fr")->second)->GetBinContent(p);
    	double fr_pt_err = (histograms.find("pt_fr")->second)->GetBinError(p);

  	if ( g_verbose ) { Info("calc_weights()", "\t\tLepton pT = %.3f GeV ==> Reading fake efficiency in pT bin [%.0f,%.0f] GeV: fr_pt = %.3f", pt, this_low_edge_pt, this_up_edge_pt, fr_pt ); }

	double fr_eta(1.0);
	double fr_eta_err(1.0);
        if ( g_ratesEta ) {
	    // loop over number of eta bins
	    // do not consider underflow, i.e. 0th bin
	    //
	    for ( int e = 1; e <= n_bins_eta; e++ ) {
		this_low_edge_eta = (histograms.find("eta_fr")->second)->GetXaxis()->GetBinLowEdge(e);
		this_up_edge_eta  = (histograms.find("eta_fr")->second)->GetXaxis()->GetBinLowEdge(e+1);

		if ( g_verbose ) { Info("calc_weights()","\t\t|eta| bin %i : [%.3f,%.3f]", e, this_low_edge_eta, this_up_edge_eta ); }

		if ( fabs(eta) >= this_low_edge_eta && fabs(eta) < this_up_edge_eta ) {

		    fr_eta     = (histograms.find("eta_fr")->second)->GetBinContent(e);
		    fr_eta_err = (histograms.find("eta_fr")->second)->GetBinError(e);

		    if ( g_verbose ) {
			Info("calc_weights()", "\t\tLepton |eta| = %.3f ==> Reading fake efficiency in |eta| bin [%.3f,%.3f]: fr_eta = %.3f, norm factor = %.3f", fabs(eta), this_low_edge_eta, this_up_edge_eta, fr_eta, fr_tot );
		    }

		    break;
		}
	    }
	}

    	// nominal
    	//
    	weights.at(0) = fr_pt;
    	error	      = fr_pt_err;

    	// up syst
    	//
    	weights.at(1) = ( fr_pt + error );

    	// down syst
    	//
    	if ( fr_pt - error > 0 ) { weights.at(2) = ( fr_pt - error );}
    	else		         { weights.at(2) = 0.0; }

	if ( g_ratesEta ) {
	    weights.at(0) = ( fr_pt * fr_eta ) / fr_tot;
	    // (assuming  fr_pt,fr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back later in the def of weight...)
	    //
	    error  = sqrt( (fr_eta*fr_pt_err)*(fr_eta*fr_pt_err) + (fr_pt*fr_eta_err)*(fr_pt*fr_eta_err) );
	    weights.at(1) = ( (fr_pt * fr_eta) + error ) / fr_tot;
	    if ( (fr_pt * fr_eta) - error > 0 ) { weights.at(2) = ( (fr_pt * fr_eta) - error ) / fr_tot;}
	    else                                { weights.at(2) = 0.0; }
	}

        break;
      }

    } // close loop on pT bins: fake lepton

  // lepton is real: choose correct pt histogram
  //
  } else {

    if ( g_verbose ) { Info("calc_weights()", "\tReading real efficiency..."); }

    // loop over number of pt bins
    // do not consider underflow, i.e. 0th bin
    //
    for ( int p = 1; p <= n_bins_pt_rr; p++ ) {

      this_low_edge_pt = (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p);
      this_up_edge_pt  = (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p+1);

      if ( g_verbose ) { Info("calc_weights()","\t\tpT bin %i : [%.0f,%.0f] GeV", p, this_low_edge_pt, this_up_edge_pt ); }

      if ( pt >= this_low_edge_pt && pt < this_up_edge_pt ) {

    	double rr_pt      = (histograms.find("pt_rr")->second)->GetBinContent(p);
    	double rr_pt_err  = (histograms.find("pt_rr")->second)->GetBinError(p);

  	if ( g_verbose ) { Info("calc_weights()", "\t\tLepton pT = %.3f GeV ==> Reading real efficiency in pT bin [%.0f,%.0f] GeV: fr_pt = %.3f", pt, this_low_edge_pt, this_up_edge_pt, rr_pt ); }

	double rr_eta(1.0);
	double rr_eta_err(1.0);
        if ( g_ratesEta ) {
	    // loop over number of eta bins
	    // do not consider underflow, i.e. 0th bin
	    //
	    for ( int e = 1; e <= n_bins_eta; e++ ) {
		this_low_edge_eta = (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e);
		this_up_edge_eta  = (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e+1);

		if ( g_verbose ) { Info("calc_weights()","\t\t|eta| bin %i : [%.3f,%.3f]", e, this_low_edge_eta, this_up_edge_eta ); }

		if ( fabs(eta) >= this_low_edge_eta && fabs(eta) < this_up_edge_eta ) {

		    rr_eta     = (histograms.find("eta_rr")->second)->GetBinContent(e);
		    rr_eta_err = (histograms.find("eta_rr")->second)->GetBinError(e);

		    if ( g_verbose ) {
			Info("calc_weights()", "\t\tLepton |eta| = %.3f ==> Reading real efficiency in |eta| bin [%.3f,%.3f]: fr_eta = %.3f, norm factor = %.3f", fabs(eta), this_low_edge_eta, this_up_edge_eta, rr_eta, rr_tot );
		    }

		    break;
		}
	    }
	}

    	// nominal
    	//
    	weights.at(0) = rr_pt;
    	error	      = rr_pt_err;

    	// up syst
    	//
    	weights.at(1) = ( rr_pt + error );

    	// down syst
    	//
    	if ( rr_pt - error > 0 ) { weights.at(2) = ( rr_pt - error ); }
    	else		         { weights.at(2) = 0.0; }

	if ( g_ratesEta ) {
	    weights.at(0) = ( rr_pt * rr_eta ) / rr_tot;
	    // (assuming  rr_pt,rr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back later in the def of weight...)
	    //
	    error  = sqrt( (rr_eta*rr_pt_err)*(rr_eta*rr_pt_err) + (rr_pt*rr_eta_err)*(rr_pt*rr_eta_err) );
	    weights.at(1) = ( (rr_pt * rr_eta) + error ) / rr_tot;
	    if ( (rr_pt * rr_eta) - error > 0 ) { weights.at(2) = ( (rr_pt * rr_eta) - error ) / rr_tot;}
	    else                                { weights.at(2) = 0.0; }
	}

        break;

      }
    } // close loop on pT bins: real lepton

  } // close check isFakeLep

  // Now converting rates to the efficiencies for the MM/FF
  //
  // --> Not needed anymore: now reading efficiency directly!
  //
  // if ( g_verbose ) {
  //    if ( !isFakeLep ) { Info("calc_weights()", "MM/FF REAL rate ==> r = %.3f ( r_up = %.3f , r_dn = %.3f )", weights.at(0), weights.at(1), weights.at(2) ); }
  //    if ( isFakeLep )  { Info("calc_weights()", "MM/FF FAKE rate ==> f = %.3f ( f_up = %.3f , f_dn = %.3f )", weights.at(0), weights.at(1), weights.at(2) ); }
  // }
  //
  //weights.at(0) = scaleRateToEfficiency(weights.at(0));
  //weights.at(1) = scaleRateToEfficiency(weights.at(1));
  //weights.at(2) = scaleRateToEfficiency(weights.at(2));

  if ( g_verbose ) {
      if ( !isFakeLep ) { Info("calc_weights()", "\t\tEffective REAL efficiency ==> r = %.3f ( r_up = %.3f , r_dn = %.3f )", weights.at(0), weights.at(1), weights.at(2) ); }
      if ( isFakeLep )  { Info("calc_weights()", "\t\tEffective FAKE efficiency ==> f = %.3f ( f_up = %.3f , f_dn = %.3f )", weights.at(0), weights.at(1), weights.at(2) ); }
  }

  return weights;
}

/* ***************************************************************************************
/
/ This function saves the final MM/FF event weight, depending on the event type (TT,TL...)
/
*************************************************************************************** */

double calc_final_event_weight( std::string region, double f1, double f2, double r1, double r2, const std::string variation )
{

   double weight = 1.0;
   double alpha  = 1.0 / ( (r1-f1) * (r2-f2) );

   // TEMP!!
   //if ( r1 > 0.955 ) r1 = 0.955;
   //if ( r2 > 0.955 ) r2 = 0.955;

   if      ( region=="TT" ) { weight = 1 - ( r1 * r2 * (1-f1) * (1-f2) * alpha ); }
   else if ( region=="TL" ) { weight = r1 * r2 * f2 * (1-f1) * alpha;  }
   else if ( region=="LT" ) { weight = r1 * r2 * f1 * (1-f2) * alpha;  }
   else if ( region=="LL" ) { weight = -1 * r1 * r2 * f1 * f2 * alpha; }

   // The above formulas are equivalent to the following:
   //
   //double weight2 = 1.0;
   //if	   ( region=="TT" )   { weight2 = alpha * ( r1 * f2 * ( (f1 - 1) * (1 - r2) ) + r2 * f1 * ( (r1 - 1) * (1 - f2) ) + f1 * f2 * ( (1 - r1) * (1 - r2) ) ); }
   //else if ( region=="TL" ) { weight2 = alpha * ( r1 * f2 * ( (1 - f1) * r2 ) + r2 * f1 * ( (1 - r1) * f2 ) + f1 * f2 * ( (r1 - 1) * r2 ) ); }
   //else if ( region=="LT" ) { weight2 = alpha * ( r1 * f2 * ( (1 - r2) * f1 ) + r2 * f1 * ( (1 - f2) * r1 ) + f1 * f2 * ( (r2 - 1) * r1 ) ); }
   //else if ( region=="LL" ) { weight2 = alpha * ( r1 * f2 * ( -1.0 * f1 * r2 ) + r2 * f1 * ( -1.0 * r1 * f2 ) + f1 * f2 * ( r1 * r2 ) ); }

   if ( g_debug ) { Info("calc_final_event_weight()", "In region %s : weight (%s) = %.3f , alpha = %.3f ", region.c_str(), variation.c_str(), weight, alpha); }

   return weight;
}

/* ********************************************
/
/ Update MMWeight in TTree by reading new rates
/
******************************************** */

void recomputeMMW( std::vector<double>* MMW_out,  /* pass it by pointer, as you are going to modify it! */
     		   Int_t	        isTT,
     		   Int_t	        isTL,
     		   Int_t	        isLT,
     		   Int_t	        isLL,
     		   std::vector<float>&  lep_pt,
     		   std::vector<float>&  lep_eta,
     		   std::vector<int>&    lep_flavour
     		 )
{

  // ********************************************************

  if ( g_debug ) {
    Info("recomputeMMW()", "Lepton 1 - flavour: %i , pT = %.3f GeV , eta = %.3f", lep_flavour.at(0), lep_pt.at(0)/1e3, lep_eta.at(0) );
    Info("recomputeMMW()", "Lepton 2 - flavour: %i , pT = %.3f GeV , eta = %.3f", lep_flavour.at(1), lep_pt.at(1)/1e3, lep_eta.at(1) );
  }

  // real and fake rates w/ syst variations
  //
  std::vector<double> r1, r2, f1, f2;
  double r1up,r1dn, r2up, r2dn, f1up, f1dn, f2up, f2dn;

  bool isFakeLep = true;

  // NB: input lep_* vectors are pT-ordered.
  //
  if ( lep_flavour.at(0) == 11 ) {
      r1 = calc_weights( g_el_hist_map, lep_pt.at(0), lep_eta.at(0), !isFakeLep, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_n_el_bins_eta, g_el_fr_tot, g_el_rr_tot );
      f1 = calc_weights( g_el_hist_map, lep_pt.at(0), lep_eta.at(0), isFakeLep,  g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_n_el_bins_eta, g_el_fr_tot, g_el_rr_tot );
  } else if ( lep_flavour.at(0) == 13 ) {
      r1 = calc_weights( g_mu_hist_map, lep_pt.at(0), lep_eta.at(0), !isFakeLep, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_n_mu_bins_eta, g_mu_fr_tot, g_mu_rr_tot );
      f1 = calc_weights( g_mu_hist_map, lep_pt.at(0), lep_eta.at(0), isFakeLep,  g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_n_mu_bins_eta, g_mu_fr_tot, g_mu_rr_tot );
  }
  if ( lep_flavour.at(1) == 11 ) {
      r2 = calc_weights( g_el_hist_map, lep_pt.at(1), lep_eta.at(1), !isFakeLep, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_n_el_bins_eta, g_el_fr_tot, g_el_rr_tot );
      f2 = calc_weights( g_el_hist_map, lep_pt.at(1), lep_eta.at(1), isFakeLep,  g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_n_el_bins_eta, g_el_fr_tot, g_el_rr_tot );
  } else if ( lep_flavour.at(1) == 13 ) {
      r2 = calc_weights( g_mu_hist_map, lep_pt.at(1), lep_eta.at(1), !isFakeLep, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_n_mu_bins_eta, g_mu_fr_tot, g_mu_rr_tot );
      f2 = calc_weights( g_mu_hist_map, lep_pt.at(1), lep_eta.at(1), isFakeLep,  g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_n_mu_bins_eta, g_mu_fr_tot, g_mu_rr_tot );
  }

  if ( g_debug ) {
    std::cout << "" << std::endl;
    Info("recomputeMMW()", "Lepton 1 - effective real and fake eff. (nominal): r1 = %.3f , f1 = %.3f ", r1.at(0), f1.at(0) );
    Info("recomputeMMW()", "Lepton 2 - effective real and fake eff. (nominal): r2 = %.3f , f2 = %.3f ", r2.at(0), f2.at(0) );
    std::cout << "" << std::endl;
  }

  // ***************************************************************************************************************

   std::string region("");
   if      ( isTT ) { region = "TT"; }
   else if ( isTL ) { region = "TL"; }
   else if ( isLT ) { region = "LT"; }
   else if ( isLL ) { region = "LL"; }

  double mm_weight(0.0);

  if ( (r1.at(0) == 0) || (r2.at(0) == 0) || (r1.at(0) <= f1.at(0)) || (r2.at(0) <= f2.at(0)) ) {
      // event will get null weight - will basically cancel out this event at plotting
      //
      if ( g_debug ) {
        Warning("recomputeMMW()", "Warning! The Matrix Method cannot be applied because : \nr1 = %.3f , r2 = %.3f, \nf1 = %.3f , f2 = %.3f \ngiven that \npt1 = %.3f , eta1 = %.3f,\npt2 = %.3f , eta2 = %.3f", r1.at(0), r2.at(0),  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
        Warning("recomputeMMW()", "applying MM weight = 0 ...");
      }
  } else {
      // calculate nominal MM weight
      //
      mm_weight      = calc_final_event_weight( region, f1.at(0), f2.at(0), r1.at(0), r2.at(0), "nominal" );
  }

  // update branch for nominal MM weight
  //
  MMW_out->at(0) = mm_weight;

  if ( mm_weight != 0.0 ) {

    // calculate MM weight with systematics
    //
    r1up = ( r1.at(1) > 1.0 ) ? 1.0 :  r1.at(1) ;
    r2up = ( r2.at(1) > 1.0 ) ? 1.0 :  r2.at(1) ;
    r1dn = r1.at(2);
    r2dn = r2.at(2);

    f1up = f1.at(1);
    f2up = f2.at(1);
    f1dn = ( f1.at(2) < 0.0 ) ? 0.0 :  f1.at(2) ;
    f2dn = ( f2.at(2) < 0.0 ) ? 0.0 :  f2.at(2) ;

    // rup syst (save relative vweight wrt. nominal)
    //
    MMW_out->at(1) = ( calc_final_event_weight( region, f1.at(0), f2.at(0), r1up, r2up, "r1, r2 UP" ) / mm_weight );
    // fdn syst (save relative vweight wrt. nominal)
    //
    MMW_out->at(4) = ( calc_final_event_weight( region, f1dn, f2dn, r1.at(0), r2.at(0), "f1, f2 DN" ) / mm_weight );

    if ( (r1dn > f1.at(0)) && (r2dn > f2.at(0)) ) {
      // rdn syst (save relative vweight wrt. nominal)
      //
	MMW_out->at(2) = ( calc_final_event_weight( region, f1.at(0), f2.at(0), r1dn, r2dn, "r1, r2 DN" ) / mm_weight );
    } else {
      if ( g_debug ) {
         Warning("recomputeMMW()", "Warning! Systematic MMWeight_rdn cannot be calculated because : \nr1dn = %.3f , r2dn = %.3f, \nf1 = %.3f , f2 = %.3f \ngiven that \npt1 = %.3f , eta1 = %.3f, \npt2 = %.3f , eta2 = %.3f ", r1dn, r2dn,  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
      }
    }

    if ( (r1.at(0) > f1up) && (r2.at(0) > f2up) ) {
      // fup syst (save relative vweight wrt. nominal)
      //
	MMW_out->at(3) = ( calc_final_event_weight( region, f1up, f2up, r1.at(0),  r2.at(0), "f1, f2 UP" ) / mm_weight );
    } else {
      if ( g_debug ) {
         Warning("recomputeMMW()", "Warning! Systematic MMWeight_fup cannot be calculated because : \nr1dn = %.3f , r2dn = %.3f, \nf1 = %.3f , f2 = %.3f \ngiven that \npt1 = %.3f , eta1 = %.3f, \npt2 = %.3f , eta2 = %.3f ", r1dn, r2dn,  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
      }
    }
  }

}

/* ******************
/
/ The 'main' function
/
****************** */

void modifyttree_MM(std::string filename, std::string outfilename, std::string addWeight, std::string RR_dir, std::string doClosure, std::string ratesEta = "NO",
		    std::string  NENTRIES = "ALL", std::string treename = "physics", std::string FR_dir = "")
{

  Info("modifytree_MM()","Starting off...");

  //Get in file, in tree and set top branch address
  //
  TFile *infile = new TFile(filename.c_str());
  TTree *intree = (TTree*)infile->Get(treename.c_str());

  Long64_t nentries;

  if ( NENTRIES == "ALL" ) {
      nentries = intree->GetEntries();
  } else {
      std::stringstream ss; ss << NENTRIES;
      int n_e;              ss >> n_e;
      nentries = n_e;
  }

  // Set global variable for reading MC rates, if needed
  //
  g_doClosure = ( doClosure.compare("YES") == 0 );

  // Read efficiencies also parametrised in eta
  //
  g_ratesEta = ( ratesEta.compare("YES") == 0 );

  if ( addWeight.compare("NO") == 0 ) {

    // TO BE MODIFIED ACCORDINGLY TO YOUR NEEDS (name and type of the variables)
    //
    std::string in_eventNumber_name("EventNumber");
    std::string in_nlep_name("nleptons");
    std::string in_isSS01_name("isSS01");
    std::string in_isTT_name("is_T_T");
    std::string in_isTL_name("is_T_AntiT");
    std::string in_isLT_name("is_AntiT_T");
    std::string in_isLL_name("is_AntiT_AntiT");
    std::string in_lep_pt_0_name("lep_Pt_0");
    std::string in_lep_pt_1_name("lep_Pt_1");
    std::string in_lep_eta_0_name("lep_Eta_0");
    std::string in_lep_eta_1_name("lep_Eta_1");
    std::string in_lep_etaBE2_0_name("lep_EtaBE2_0");
    std::string in_lep_etaBE2_1_name("lep_EtaBE2_1");
    std::string in_lep_flavour_0_name("lep_ID_0");
    std::string in_lep_flavour_1_name("lep_ID_1");
    std::string in_lep_isT_0_name("lep_isTightSelected_0");
    std::string in_lep_isT_1_name("lep_isTightSelected_1");
    std::string in_MMWeight_name("MMWeight");

    ULong64_t  eventNumber_in; eventNumber_in = -1;
    Int_t      nlep_in;        nlep_in = -1;
    Char_t     isSS01_in;      isSS01_in = -1;
    Char_t     isTT_in;        isTT_in = -1;
    Char_t     isTL_in;        isTL_in = -1;
    Char_t     isLT_in;        isLT_in = -1;
    Char_t     isLL_in;        isLL_in = -1;
    Float_t    lep_pt_0_in;      lep_pt_0_in = 0;
    Float_t    lep_pt_1_in;      lep_pt_1_in = 0;
    Float_t    lep_eta_0_in;	 lep_eta_0_in = 0;
    Float_t    lep_eta_1_in;	 lep_eta_1_in = 0;
    Float_t    lep_etaBE2_0_in;  lep_etaBE2_0_in = 0;
    Float_t    lep_etaBE2_1_in;  lep_etaBE2_1_in = 0;
    Float_t    lep_flavour_0_in; lep_flavour_0_in = 0;
    Float_t    lep_flavour_1_in; lep_flavour_1_in = 0;
    Char_t     lep_isT_0_in;     lep_isT_0_in = 0;
    Char_t     lep_isT_1_in;     lep_isT_1_in = 0;
    std::vector<double>*   MMWeight_in;    MMWeight_in = 0;

    // List of input branches
    //
    TBranch      *b_EventNumber = 0;      //!
    TBranch	 *b_nleptons  = 0;        //!
    TBranch      *b_isSS01      = 0;      //!
    TBranch	 *b_is_T_T = 0;    	  //!
    TBranch	 *b_is_T_AntiT = 0;    	  //!
    TBranch	 *b_is_AntiT_T = 0;    	  //!
    TBranch	 *b_is_AntiT_AntiT = 0;   //!
    TBranch	 *b_lep_Pt_0 = 0;         //!
    TBranch	 *b_lep_Pt_1 = 0;         //!
    TBranch	 *b_lep_Eta_0 = 0;        //!
    TBranch	 *b_lep_Eta_1 = 0;        //!
    TBranch	 *b_lep_EtaBE2_0 = 0;     //!
    TBranch	 *b_lep_EtaBE2_1 = 0;     //!
    TBranch	 *b_lep_ID_0 = 0;         //!
    TBranch	 *b_lep_ID_1 = 0;         //!
    TBranch      *b_lep_isTightSelected_0 = 0;  //!
    TBranch      *b_lep_isTightSelected_1 = 0;  //!
    TBranch      *b_MMWeight = 0;         //!

    // Before cloning input TTree, tell ROOT to process all the in branches,
    // except for the one you want to change
    //
    intree->SetBranchStatus("*",1);
    intree->SetBranchStatus(in_MMWeight_name.c_str(),0);

    std::string out_MMWeight_name("MMWeight");
    std::vector<double>*   MMWeight_out;  MMWeight_out = 0;

    // Create a new file + a clone of in tree in new file
    //
    TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");
    TTree *outtree = intree->CloneTree(0); //clone only the structure

    // MUST re-activate the branch(es) previously deactivated before calling SetBranchAddress()!!
    //
    intree->SetBranchStatus(in_MMWeight_name.c_str(),1);

    // Get branches from input TTree
    //
    intree->SetBranchAddress(in_eventNumber_name.c_str(), &eventNumber_in, &b_EventNumber);
    intree->SetBranchAddress(in_nlep_name.c_str(), &nlep_in, &b_nleptons);
    intree->SetBranchAddress(in_isSS01_name.c_str(), &isSS01_in, &b_isSS01);
    intree->SetBranchAddress(in_isTT_name.c_str(), &isTT_in, &b_is_T_T);
    intree->SetBranchAddress(in_isTL_name.c_str(), &isTL_in, &b_is_T_AntiT);
    intree->SetBranchAddress(in_isLT_name.c_str(), &isLT_in, &b_is_AntiT_T);
    intree->SetBranchAddress(in_isLL_name.c_str(), &isLL_in, &b_is_AntiT_AntiT);
    intree->SetBranchAddress(in_lep_pt_0_name.c_str(), &lep_pt_0_in, &b_lep_Pt_0);
    intree->SetBranchAddress(in_lep_pt_1_name.c_str(), &lep_pt_1_in, &b_lep_Pt_1);
    intree->SetBranchAddress(in_lep_eta_0_name.c_str(), &lep_eta_0_in, &b_lep_Eta_0);
    intree->SetBranchAddress(in_lep_eta_1_name.c_str(), &lep_eta_1_in, &b_lep_Eta_1);
    intree->SetBranchAddress(in_lep_etaBE2_0_name.c_str(), &lep_etaBE2_0_in, &b_lep_EtaBE2_0);
    intree->SetBranchAddress(in_lep_etaBE2_1_name.c_str(), &lep_etaBE2_1_in, &b_lep_EtaBE2_1);
    intree->SetBranchAddress(in_lep_flavour_0_name.c_str(), &lep_flavour_0_in, &b_lep_ID_0);
    intree->SetBranchAddress(in_lep_flavour_1_name.c_str(), &lep_flavour_1_in, &b_lep_ID_1);
    intree->SetBranchAddress(in_lep_isT_0_name.c_str(), &lep_isT_0_in, &b_lep_isTightSelected_0);
    intree->SetBranchAddress(in_lep_isT_1_name.c_str(), &lep_isT_1_in, &b_lep_isTightSelected_1);
    intree->SetBranchAddress(in_MMWeight_name.c_str(), &MMWeight_in, &b_MMWeight);

    //Info("modifytree_MM()", "--> lep_pt_0 before SetBranchAddress() %p\n", intree->GetBranch(in_lep_pt_0_name.c_str())->GetAddress());
    //Info("modifytree_MM()", "--> lep_pt_0 after SetBranchAddress() %p\n", intree->GetBranch(in_lep_pt_0_name.c_str())->GetAddress());

    // Set the "new" branches in the output TTree
    //
    outtree->Branch(out_MMWeight_name.c_str(), &MMWeight_out);

    // read r/f rates from ROOT histograms
    //
    Info("modifytree_MM()","Reading r/f rates from ROOT file(s)..");
    read_rates(RR_dir,FR_dir);

    // Loop over entries in TTree
    //
    Info("modifytree_MM()","Begin loop on input tree entries...\n");
    int count_bad(0);
    Long64_t i = 0;
    for ( ; i < nentries; i++ ) {

      // Print out every N events to see where we are
      //
      if ( i > 0 && ( static_cast<int>(i) % 20000 == 0 ) ) { Info("modifytree_MM()","\t Processed %lld entries",i); }

      // Now, in the input tree, reset to 1 the status of the branches you
      // deactivated before cloning
      //
      if ( !intree->GetBranchStatus(in_MMWeight_name.c_str()) ) {
        intree->SetBranchStatus(in_MMWeight_name.c_str(),1);
      }

      intree->GetEntry(i);

      if ( g_debug ) { Info("modifytree_MM()","\t Processing entry: %lld - eventNumber: %lli \n",i, eventNumber_in); }

      // A security check...
      //
      if ( !MMWeight_in ) {
        Info("modifytree_MM()","\t --> MMWeight_in is NULL!! Skipping event...  \n");
        ++count_bad;
        continue;
      }

      // To start off, copy the in branch into the new
      // (then it will be overridden , if necessary)
      //
      *MMWeight_out = *MMWeight_in;

      // and now, recompute the MM weights!
      // Do it only for 2lep events
      //
      if ( nlep_in == 2 ) {

        if ( g_debug ) {
          int idx_in(0);
          for ( const auto& itr : *MMWeight_in ) {
	      if ( idx_in == 0 ) Info("modifytree_MM()","IN MMWeight[%i] = %.3f", idx_in, itr );
	      else                Info("modifytree_MM()","IN MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx_in, idx_in, itr, ( itr * *(MMWeight_out->begin()) )  );
	      ++idx_in;
          }
	  std::cout << "" << std::endl;
        }

        int TT =  ( isTT_in );
        int TL =  ( isTL_in );
        int LT =  ( isLT_in );
        int LL =  ( isLL_in );

        // Create these vectors to avoid changing the signature of the function
        //
        std::vector<float> lep_pt_vec      = { lep_pt_0_in, lep_pt_1_in };
        std::vector<int> lep_flavour_vec   = { abs(static_cast<int>(lep_flavour_0_in)), abs(static_cast<int>(lep_flavour_1_in)) };
        std::vector<float> lep_eta_vec;
	float eta0 = ( fabs(lep_flavour_0_in) == 13 ) ? lep_eta_0_in : lep_etaBE2_0_in; lep_eta_vec.push_back(eta0);
	float eta1 = ( fabs(lep_flavour_1_in) == 13 ) ? lep_eta_1_in : lep_etaBE2_1_in; lep_eta_vec.push_back(eta1);

        recomputeMMW( MMWeight_out, TT, TL, LT, LL, lep_pt_vec, lep_eta_vec, lep_flavour_vec );

        if ( g_debug ) {
	  std::cout << "" << std::endl;
          int idx_out(0);
          for ( const auto& itr : *MMWeight_out ) {
	      if ( idx_out == 0 ) Info("modifytree_MM()","OUT MMWeight[%i] = %.3f", idx_out, itr );
	      else                Info("modifytree_MM()","OUT MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx_out, idx_out, itr, ( itr * *(MMWeight_out->begin()) )  );
	      ++idx_out;
          }
	  std::cout << "" << std::endl;
        }

        if ( g_debug ) { Info("modifytree_MM()","\n\n"); }

      }

      // to avoid overriding new branch (in has same name) ?
      intree->SetBranchStatus(in_MMWeight_name.c_str(),0);

      outtree->Fill();

    }

    Info("modifytree_MM()","End of loop!\n ---> total number of processed events: %lld \n ---> number of skipped events: %i \n", i, count_bad );

    outfile->Write();
    outfile->Close();

    // Since we passed the address of a local variable we need
    // to remove it.
    intree->ResetBranchAddresses();

  } else if ( addWeight.compare("YES") == 0 ) {

    // TO BE MODIFIED ACCORDINGLY TO YOUR NEEDS (name and type of the variables)
    //
    std::string in_eventNumber_name("EventNumber");
    std::string in_nlep_name("nleptons");
    std::string in_isSS01_name("isSS01");
    std::string in_isTT_name("is_T_T");
    std::string in_isTL_name("is_T_AntiT");
    std::string in_isLT_name("is_AntiT_T");
    std::string in_isLL_name("is_AntiT_AntiT");
    std::string in_lep_pt_0_name("lep_Pt_0");
    std::string in_lep_pt_1_name("lep_Pt_1");
    std::string in_lep_eta_0_name("lep_Eta_0");
    std::string in_lep_eta_1_name("lep_Eta_1");
    std::string in_lep_etaBE2_0_name("lep_EtaBE2_0");
    std::string in_lep_etaBE2_1_name("lep_EtaBE2_1");
    std::string in_lep_flavour_0_name("lep_ID_0");
    std::string in_lep_flavour_1_name("lep_ID_1");
    std::string in_lep_isT_0_name("lep_isTightSelected_0");
    std::string in_lep_isT_1_name("lep_isTightSelected_1");

    ULong64_t  eventNumber_in; eventNumber_in = -1;
    Int_t      nlep_in;        nlep_in = -1;
    Char_t     isSS01_in;      isSS01_in = -1;
    Char_t     isTT_in;        isTT_in = -1;
    Char_t     isTL_in;        isTL_in = -1;
    Char_t     isLT_in;        isLT_in = -1;
    Char_t     isLL_in;        isLL_in = -1;
    Float_t    lep_pt_0_in;      lep_pt_0_in = 0;
    Float_t    lep_pt_1_in;      lep_pt_1_in = 0;
    Float_t    lep_eta_0_in;     lep_eta_0_in = 0;
    Float_t    lep_eta_1_in;     lep_eta_1_in = 0;
    Float_t    lep_etaBE2_0_in;  lep_etaBE2_0_in = 0;
    Float_t    lep_etaBE2_1_in;  lep_etaBE2_1_in = 0;
    Float_t    lep_flavour_0_in; lep_flavour_0_in = 0;
    Float_t    lep_flavour_1_in; lep_flavour_1_in = 0;
    Char_t     lep_isT_0_in;     lep_isT_0_in = 0;
    Char_t     lep_isT_1_in;     lep_isT_1_in = 0;

    // List of input branches
    //
    TBranch      *b_EventNumber = 0;      //!
    TBranch	 *b_nleptons  = 0;        //!
    TBranch      *b_isSS01      = 0;      //!
    TBranch	 *b_is_T_T = 0;    	  //!
    TBranch	 *b_is_T_AntiT = 0;    	  //!
    TBranch	 *b_is_AntiT_T = 0;    	  //!
    TBranch	 *b_is_AntiT_AntiT = 0;   //!
    TBranch	 *b_lep_Pt_0 = 0;         //!
    TBranch	 *b_lep_Pt_1 = 0;         //!
    TBranch	 *b_lep_Eta_0 = 0;        //!
    TBranch	 *b_lep_Eta_1 = 0;        //!
    TBranch	 *b_lep_EtaBE2_0 = 0;     //!
    TBranch	 *b_lep_EtaBE2_1 = 0;     //!
    TBranch	 *b_lep_ID_0 = 0;         //!
    TBranch	 *b_lep_ID_1 = 0;         //!
    TBranch      *b_lep_isTightSelected_0 = 0;  //!
    TBranch      *b_lep_isTightSelected_1 = 0;  //!

    // Before cloning input TTree, tell ROOT to process all the input branches
    //
    intree->SetBranchStatus("*",1);

    // Create a new file + a clone of input tree in new file
    //
    TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");
    TTree *outtree = intree->CloneTree(0); //clone only the structure

    // Get branches from input TTree
    //
    intree->SetBranchAddress(in_eventNumber_name.c_str(), &eventNumber_in, &b_EventNumber);
    intree->SetBranchAddress(in_nlep_name.c_str(), &nlep_in, &b_nleptons);
    intree->SetBranchAddress(in_isSS01_name.c_str(), &isSS01_in, &b_isSS01);
    intree->SetBranchAddress(in_isTT_name.c_str(), &isTT_in, &b_is_T_T);
    intree->SetBranchAddress(in_isTL_name.c_str(), &isTL_in, &b_is_T_AntiT);
    intree->SetBranchAddress(in_isLT_name.c_str(), &isLT_in, &b_is_AntiT_T);
    intree->SetBranchAddress(in_isLL_name.c_str(), &isLL_in, &b_is_AntiT_AntiT);
    intree->SetBranchAddress(in_lep_pt_0_name.c_str(), &lep_pt_0_in, &b_lep_Pt_0);
    intree->SetBranchAddress(in_lep_pt_1_name.c_str(), &lep_pt_1_in, &b_lep_Pt_1);
    intree->SetBranchAddress(in_lep_eta_0_name.c_str(), &lep_eta_0_in, &b_lep_Eta_0);
    intree->SetBranchAddress(in_lep_eta_1_name.c_str(), &lep_eta_1_in, &b_lep_Eta_1);
    intree->SetBranchAddress(in_lep_etaBE2_0_name.c_str(), &lep_etaBE2_0_in, &b_lep_EtaBE2_0);
    intree->SetBranchAddress(in_lep_etaBE2_1_name.c_str(), &lep_etaBE2_1_in, &b_lep_EtaBE2_1);
    intree->SetBranchAddress(in_lep_flavour_0_name.c_str(), &lep_flavour_0_in, &b_lep_ID_0);
    intree->SetBranchAddress(in_lep_flavour_1_name.c_str(), &lep_flavour_1_in, &b_lep_ID_1);
    intree->SetBranchAddress(in_lep_isT_0_name.c_str(), &lep_isT_0_in, &b_lep_isTightSelected_0);
    intree->SetBranchAddress(in_lep_isT_1_name.c_str(), &lep_isT_1_in, &b_lep_isTightSelected_1);

    // Set the "new" branches in the output TTree
    //
    std::vector<double>  MMWeight_out;
    outtree->Branch("MMWeight", &MMWeight_out);

    // read r/f rates from ROOT histograms
    //
    Info("modifytree_MM()","Reading r/f rates from ROOT file(s)..");
    read_rates(RR_dir,FR_dir);

    // Loop over entries in TTree
    //
    Info("modifytree_MM()","Begin loop on input tree entries...\n");
    int count_bad(0);
    Long64_t i = 0;
    for ( ; i < nentries; i++ ) {

      // Print out every N events to see where we are
      //
      if ( i > 0 && ( static_cast<int>(i) % 20000 == 0 ) ) { Info("modifytree_MM()","\t Processed %lld entries",i); }

      intree->GetEntry(i);

      if ( g_debug ) { Info("modifytree_MM()","\t Processing entry: %lld - eventNumber: %lli \n",i, eventNumber_in); }

      MMWeight_out.assign(5,1.0);

      if ( g_debug ) {
        int idx(0);
        for ( const auto& itr : MMWeight_out ) {
	   if ( idx == 0 ) Info("modifytree_MM()","Default MMWeight[%i] = %.3f", idx, itr );
	   else                Info("modifytree_MM()","Default MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(MMWeight_out.begin()) )  );
	   ++idx;
        }
	std::cout << "" << std::endl;
      }

      // and now, recompute the MM weights!
      // Do it only for 2lep events
      //
      if ( nlep_in == 2 ) {

        int TT =  ( isTT_in );
        int TL =  ( isTL_in );
        int LT =  ( isLT_in );
        int LL =  ( isLL_in );

        // Create these vectors to avoid changing the signature of the function
        //
        std::vector<float> lep_pt_vec      = { lep_pt_0_in, lep_pt_1_in };
        std::vector<int> lep_flavour_vec   = { abs(static_cast<int>(lep_flavour_0_in)), abs(static_cast<int>(lep_flavour_1_in)) };
        std::vector<float> lep_eta_vec;
	float eta0 = ( fabs(lep_flavour_0_in) == 13 ) ? lep_eta_0_in : lep_etaBE2_0_in; lep_eta_vec.push_back(eta0);
	float eta1 = ( fabs(lep_flavour_1_in) == 13 ) ? lep_eta_1_in : lep_etaBE2_1_in; lep_eta_vec.push_back(eta1);

        recomputeMMW( &MMWeight_out, TT, TL, LT, LL, lep_pt_vec, lep_eta_vec, lep_flavour_vec );

        if ( g_debug ) {
	  std::cout << "" << std::endl;
          int idx_out(0);
          for ( const auto& itr : MMWeight_out ) {
	      if ( idx_out == 0 ) Info("modifytree_MM()","OUT MMWeight[%i] = %.3f", idx_out, itr );
	      else                Info("modifytree_MM()","OUT MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx_out, idx_out, itr, ( itr * *(MMWeight_out.begin()) )  );
	      ++idx_out;
          }
	  std::cout << "" << std::endl;
        }

        if ( g_debug ) { Info("modifytree_MM()","\n\n"); }

      }

      outtree->Fill();

    }

    Info("modifytree_MM()","End of loop!\n ---> total number of processed events: %lld \n ---> number of skipped events: %i \n", i, count_bad );

    outfile->Write();
    outfile->Close();

    // Since we passed the address of a local variable we need
    // to remove it.
    intree->ResetBranchAddresses();

  }

  // //delete a branch:
  // TFile f("myfile.root","update");
  // TTree *T = (TTree*)f.Get("treename");
  // TBranch *b = T->GetBranch("name of branch to delete");
  // T->GetListOfBranches()->Remove(b);
  // TLeaf* l= T->GetLeaf("name of branch to delete");
  // T->GetListOfLeaves()->Remove(l)
  // T->Write();

}
