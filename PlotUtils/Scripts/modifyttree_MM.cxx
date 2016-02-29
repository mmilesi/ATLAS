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

bool g_debug(false);
bool g_verbose(false);

std::map< std::string, TH1D* > g_el_hist_map;
std::map< std::string, TH1D* > g_mu_hist_map;

int g_n_el_bins_eta(0);
int g_n_mu_bins_eta(0);
int g_n_el_bins_pt_rr(0);
int g_n_el_bins_pt_fr(0);
int g_n_mu_bins_pt_rr(0);
int g_n_mu_bins_pt_fr(0);

double g_el_rr_tot(9e5);
double g_el_fr_tot(9e5);
double g_mu_rr_tot(9e5);
double g_mu_fr_tot(9e5);

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

void read_rates(const std::string rr_dir, const std::string fr_dir = "")
{

  std::string glob_path("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/");

  // ******************************************************
  //
  // 1. 'REAL' rate

  Info("read_rates()", "REAL rate from directory: %s ", rr_dir.c_str() );

  std::string path_R_el = glob_path + rr_dir + "/Rates.root";// "/AvgRates.root";
  TFile *file_R_el = TFile::Open(path_R_el.c_str());
  if ( !file_R_el->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with R rate from path: %s . Aborting", path_R_el.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "ELECTRON REAL rate: %s ", path_R_el.c_str() );
  }
  std::string path_R_mu = glob_path + rr_dir + "/Rates.root"; //"/AvgRates.root";
  TFile *file_R_mu = TFile::Open(path_R_mu.c_str());
  if ( !file_R_mu->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with R rate from path: %s . Aborting", path_R_mu.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "MUON REAL rate: %s ", path_R_mu.c_str() );
  }

  // special case: use charge flip (i.e. "fake" rate) as real rate (only for electrons)!
  bool CFRateAsRR = ( rr_dir.find("ChFlip") != std::string::npos ) ? true : false;
  if ( CFRateAsRR ) { Info("read_rates()", "---> For ELECTRONS, will be treating charge flip rate as REAL rate!!"); }

  // ELECTRONS
  //
  std::string histname_el_eta_rr   = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_Rate_observed" : "El_ProbeEta_Fake_Rate_observed";
  std::string histname_el_pt_rr    = ( !CFRateAsRR ) ?  "El_ProbePt_Real_Rate_observed" : "El_ProbePt_Fake_Rate_observed";
  std::string histname_el_eta_r_T  = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_T_observed" : "El_ProbeEta_Fake_Rate_observed";
  std::string histname_el_eta_r_L  = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_L_observed" : "El_ProbeEta_Fake_Rate_observed";
  std::string histname_el_pt_r_T   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_T_observed" : "El_ProbePt_Fake_Rate_observed";
  std::string histname_el_pt_r_L   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_L_observed" : "El_ProbePt_Fake_Rate_observed";

  // get real rate histograms
  //
  TH1D *hist_el_eta_rr   = get_hist( *file_R_el, histname_el_eta_rr );
  TH1D *hist_el_pt_rr    = get_hist( *file_R_el, histname_el_pt_rr );
  TH1D *hist_el_eta_r_T  = get_hist( *file_R_el, histname_el_eta_r_T );
  TH1D *hist_el_eta_r_L  = get_hist( *file_R_el, histname_el_eta_r_L );
  TH1D *hist_el_pt_r_T   = get_hist( *file_R_el, histname_el_pt_r_T );
  TH1D *hist_el_pt_r_L   = get_hist( *file_R_el, histname_el_pt_r_L );

  // MUONS
  //
  std::string histname_mu_eta_rr   = "Mu_ProbeEta_Real_Rate_observed";
  std::string histname_mu_pt_rr    = "Mu_ProbePt_Real_Rate_observed";
  std::string histname_mu_eta_r_T  = "Mu_ProbeEta_Real_T_observed";
  std::string histname_mu_eta_r_L  = "Mu_ProbeEta_Real_L_observed";
  std::string histname_mu_pt_r_T   = "Mu_ProbePt_Real_T_observed";
  std::string histname_mu_pt_r_L   = "Mu_ProbePt_Real_L_observed";

  // get real rate histograms
  //
  TH1D *hist_mu_eta_rr   = get_hist( *file_R_mu, histname_mu_eta_rr );
  TH1D *hist_mu_pt_rr    = get_hist( *file_R_mu, histname_mu_pt_rr );
  TH1D *hist_mu_eta_r_T  = get_hist( *file_R_mu, histname_mu_eta_r_T );
  TH1D *hist_mu_eta_r_L  = get_hist( *file_R_mu, histname_mu_eta_r_L );
  TH1D *hist_mu_pt_r_T   = get_hist( *file_R_mu, histname_mu_pt_r_T );
  TH1D *hist_mu_pt_r_L   = get_hist( *file_R_mu, histname_mu_pt_r_L );

  // 2. FAKE rate

  std::string fake_dir(rr_dir);
  if ( !fr_dir.empty() ) {
     Warning("read_rates()", "FAKE rate is going to be read from %s. Check whether it's really what you want...", fr_dir.c_str());
     fake_dir = fr_dir;
  } else {
     Info("read_rates()", "FAKE rate from same directory" );
  }

  std::string path_F_el = glob_path + fake_dir + "/Rates.root"; // "/AvgRates.root";
  TFile *file_F_el = TFile::Open(path_F_el.c_str());
  if ( !file_F_el->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with F rate from path: %s . Aborting", path_F_el.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "ELECTRON FAKE rate: %s ", path_F_el.c_str() );
  }

  std::string path_F_mu = glob_path + fake_dir + "/Rates.root"; // "/AvgRates.root";
  TFile *file_F_mu = TFile::Open(path_F_mu.c_str());
  if ( !file_F_mu->IsOpen() ) {
    SysError("read_rates()", "Failed to open ROOT file with F rate from path: %s . Aborting", path_F_mu.c_str() );
    exit(-1);
  } else {
    Info("read_rates()", "MUON FAKE rate: %s ", path_F_mu.c_str() );
  }

  // ELECTRONS
  //
  std::string histname_el_eta_fr   = "El_ProbeEta_Fake_Rate_observed";
  std::string histname_el_pt_fr    = "El_ProbePt_Fake_Rate_observed";
  std::string histname_el_eta_f_T  = "El_ProbeEta_Fake_T_observed";
  std::string histname_el_eta_f_L  = "El_ProbeEta_Fake_L_observed";
  std::string histname_el_pt_f_T   = "El_ProbePt_Fake_T_observed";
  std::string histname_el_pt_f_L   = "El_ProbePt_Fake_L_observed";

  // get fake rate histograms
  //
  TH1D *hist_el_eta_fr   = get_hist( *file_F_el, histname_el_eta_fr );
  TH1D *hist_el_pt_fr    = get_hist( *file_F_el, histname_el_pt_fr );
  TH1D *hist_el_eta_f_T  = get_hist( *file_F_el, histname_el_eta_f_T );
  TH1D *hist_el_eta_f_L  = get_hist( *file_F_el, histname_el_eta_f_L );
  TH1D *hist_el_pt_f_T   = get_hist( *file_F_el, histname_el_pt_f_T );
  TH1D *hist_el_pt_f_L   = get_hist( *file_F_el, histname_el_pt_f_L );

  // MUONS
  //
  std::string histname_mu_eta_fr   = "Mu_ProbeEta_Fake_Rate_observed";
  std::string histname_mu_pt_fr    = "Mu_ProbePt_Fake_Rate_observed";
  std::string histname_mu_eta_f_T  = "Mu_ProbeEta_Fake_T_observed";
  std::string histname_mu_eta_f_L  = "Mu_ProbeEta_Fake_L_observed";
  std::string histname_mu_pt_f_T   = "Mu_ProbePt_Fake_T_observed";
  std::string histname_mu_pt_f_L   = "Mu_ProbePt_Fake_L_observed";

  // get fake rate histograms
  //
  TH1D *hist_mu_eta_fr   = get_hist( *file_F_mu, histname_mu_eta_fr );
  TH1D *hist_mu_pt_fr    = get_hist( *file_F_mu, histname_mu_pt_fr );
  TH1D *hist_mu_eta_f_T  = get_hist( *file_F_mu, histname_mu_eta_f_T );
  TH1D *hist_mu_eta_f_L  = get_hist( *file_F_mu, histname_mu_eta_f_L );
  TH1D *hist_mu_pt_f_T   = get_hist( *file_F_mu, histname_mu_pt_f_T );
  TH1D *hist_mu_pt_f_L   = get_hist( *file_F_mu, histname_mu_pt_f_L );

  // ***********************************************************************

  // fill a map for later usage
  //
  g_el_hist_map["eta_rr"]  = hist_el_eta_rr;
  g_el_hist_map["pt_rr"]   = hist_el_pt_rr;
  g_el_hist_map["eta_r_T"] = hist_el_eta_r_T;
  g_el_hist_map["eta_r_L"] = hist_el_eta_r_L;
  g_el_hist_map["pt_r_T"]  = hist_el_pt_r_T;
  g_el_hist_map["pt_r_L"]  = hist_el_pt_r_L;

  g_mu_hist_map["eta_rr"]  = hist_mu_eta_rr;
  g_mu_hist_map["pt_rr"]   = hist_mu_pt_rr;
  g_mu_hist_map["eta_r_T"] = hist_mu_eta_r_T;
  g_mu_hist_map["eta_r_L"] = hist_mu_eta_r_L;
  g_mu_hist_map["pt_r_T"]  = hist_mu_pt_r_T;
  g_mu_hist_map["pt_r_L"]  = hist_mu_pt_r_L;

  g_el_hist_map["eta_fr"]  = hist_el_eta_fr;
  g_el_hist_map["pt_fr"]   = hist_el_pt_fr;
  g_el_hist_map["eta_f_T"] = hist_el_eta_f_T;
  g_el_hist_map["eta_f_L"] = hist_el_eta_f_L;
  g_el_hist_map["pt_f_T"]  = hist_el_pt_f_T;
  g_el_hist_map["pt_f_L"]  = hist_el_pt_f_L;

  g_mu_hist_map["eta_fr"]  = hist_mu_eta_fr;
  g_mu_hist_map["pt_fr"]   = hist_mu_pt_fr;
  g_mu_hist_map["eta_f_T"] = hist_mu_eta_f_T;
  g_mu_hist_map["eta_f_L"] = hist_mu_eta_f_L;
  g_mu_hist_map["pt_f_T"]  = hist_mu_pt_f_T;
  g_mu_hist_map["pt_f_L"]  = hist_mu_pt_f_L;

  // eta hist has same binning for r/f
  //
  g_n_el_bins_eta   =  hist_el_eta_rr->GetNbinsX();
  g_n_mu_bins_eta   =  hist_mu_eta_rr->GetNbinsX();

  // pt hist has two different binning for r/f
  //
  g_n_el_bins_pt_rr =  hist_el_pt_rr->GetNbinsX();
  g_n_el_bins_pt_fr =  hist_el_pt_fr->GetNbinsX();
  g_n_mu_bins_pt_rr =  hist_mu_pt_rr->GetNbinsX();
  g_n_mu_bins_pt_fr =  hist_mu_pt_fr->GetNbinsX();

  // normalistaion factor is the same for eta and pt r/f histograms: use eta
  //
  g_el_rr_tot = ( hist_el_eta_r_T->Integral() ) / ( hist_el_eta_r_L->Integral() );
  g_el_fr_tot = ( hist_el_eta_f_T->Integral() ) / ( hist_el_eta_f_L->Integral() );
  g_mu_rr_tot = ( hist_mu_eta_r_T->Integral() ) / ( hist_mu_eta_r_L->Integral() );
  g_mu_fr_tot = ( hist_mu_eta_f_T->Integral() ) / ( hist_mu_eta_f_L->Integral() );

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
				   int n_bins_eta,
				   int n_bins_pt_fr,
				   int n_bins_pt_rr,
				   double fr_tot,
				   double rr_tot
				  )
{
  // Read the real/fake rates from input histograms
  //
  // Will eventually convert these to real/fake FACTORS

  // As a first thing, convert pT in GeV!
  //
  pt = pt/1e3;

  std::vector<double> weights(3,0.0); //initialized with zeroes

  weights.at(0) = 1.0;
  double error(0.0);

  // loop over number of eta bins
  // do not consider underflow, i.e. 0th bin
  //
  for ( int e = 1; e <= n_bins_eta; e++ ) {

    // check whether the eta under question is in *this* eta range
    //
    if ( g_verbose ) {
      Info("calc_weights()", "\tlepton fabs(eta) = %f", fabs(eta) );
      Info("calc_weights()", "\tbin %i : lower edge = %f, upper edge = %f", e, (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e),(histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e+1) );
    }

    if ( ( fabs(eta) >= (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e) ) && ( fabs(eta) < (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e+1) ) ) {

      // case 1) : lepton is fake: choose correct pt histogram
      //
      if ( isFakeLep ) {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
        for ( int p = 1; p <= n_bins_pt_fr; p++ ) {

          if ( g_verbose ) {
            Info("calc_weights()", "\t\tlepton pT = %f", pt );
            Info("calc_weights()", "\t\tbin %i : lower edge = %f, upper edge = %f", p,(histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p), (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p+1) );
          }

     	  if ( ( pt >= (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p) ) && ( pt < (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p+1) ) ) {

	    // combine eta and pt rates
	    //
	    double fr_pt  = (histograms.find("pt_fr")->second)->GetBinContent(p);
	    double fr_eta = (histograms.find("eta_fr")->second)->GetBinContent(e);

	    double fr_pt_err  = (histograms.find("pt_fr")->second)->GetBinError(p);
	    double fr_eta_err = (histograms.find("eta_fr")->second)->GetBinError(e);

            if ( g_verbose ) {
	       Info("calc_weights()", "\t\tFake lepton"  );
	       Info("calc_weights()", "\t\tfr_pt = %f", fr_pt );
	       Info("calc_weights()", "\t\tfr_eta = %f", fr_eta );
	    }

	    // nominal
	    //
     	    weights.at(0) = ( fr_pt * fr_eta ) / fr_tot;

	    // (assuming  fr_pt,fr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back later in the def of weight...
	    //
	    error	  = sqrt( (fr_eta*fr_pt_err)*(fr_eta*fr_pt_err) + (fr_pt*fr_eta_err)*(fr_pt*fr_eta_err) );

	    // up syst
	    //
     	    weights.at(1) = ( (fr_pt * fr_eta) + error ) / fr_tot;

     	    // down syst
     	    //
	    if ( (fr_pt * fr_eta) - error > 0 ) { weights.at(2) = ( (fr_pt * fr_eta) - error ) / fr_tot;}
	    else                                { weights.at(2) = 0.0; }

     	  }

	} // close loop on pT bins: fake lepton

      // lepton is real: choose correct pt histogram
      //
      }	else {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
        for ( int p = 1; p <= n_bins_pt_rr; p++ ) {

     	  if ( ( pt >= (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p) ) && ( pt < (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p+1) ) ) {


	    // combine eta and pt rates
	    //
	    double rr_pt  = (histograms.find("pt_rr")->second)->GetBinContent(p);
	    double rr_eta = (histograms.find("eta_rr")->second)->GetBinContent(e);

	    double rr_pt_err  = (histograms.find("pt_rr")->second)->GetBinError(p);
	    double rr_eta_err = (histograms.find("eta_rr")->second)->GetBinError(e);

            if ( g_verbose ) {
	       Info("calc_weights()", "Real lepton"  );
	       Info("calc_weights()", "rr_pt = %f", rr_pt );
	       Info("calc_weights()", "rr_eta = %f", rr_eta );
	    }

	    // nominal
	    //
	    weights.at(0) = ( rr_pt * rr_eta ) / rr_tot;

	    // (assuming  rr_pt,rr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back in the def of weight...
	    //
	    error	  = sqrt( (rr_eta*rr_pt_err)*(rr_eta*rr_pt_err) + (rr_pt*rr_eta_err)*(rr_pt*rr_eta_err) );

     	    // up syst
	    //
     	    weights.at(1) = ( (rr_pt * rr_eta) + error ) / rr_tot;

	    // down syst
	    //
     	    if ( (rr_pt * rr_eta) - error > 0 ) { weights.at(2) = ( (rr_pt * rr_eta) - error ) / rr_tot; }
	    else                                { weights.at(2) = 0.0; }

     	  }
        } // close loop on pT bins: real lepton

      } // close check isFakeLep

    } // close check on eta bin

  } // close loop on eta bins

  // Now converting rates to the factors for the MM/FF
  //
  if ( g_verbose ) { Info("calc_weights()", "Rates = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  weights.at(0) = scaleRateToEfficiency(weights.at(0));
  weights.at(1) = scaleRateToEfficiency(weights.at(1));
  weights.at(2) = scaleRateToEfficiency(weights.at(2));

  if ( g_verbose ) { Info("calc_weights()", "MM/FF efficiency = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  return weights;
}
//*/
/* ********************************************************
/
/ Function to calculate r/f weights and their unceratinties
/ (DIFFERENT TREATMENT FOR OVERFLOW BIN)
/
******************************************************** */
/*
std::vector<double>  calc_weights( std::map< std::string, TH1D* >& histograms,
				   float pt,
				   float eta,
				   bool isFakeLep,
				   int n_bins_eta,
				   int n_bins_pt_fr,
				   int n_bins_pt_rr,
				   double fr_tot,
				   double rr_tot
				  )
{

  // Read the real/fake rates from input histograms
  //
  // Will eventually convert these to real/fake FACTORS

  // As a first thing, convert pT in GeV!
  //
  pt = pt/1e3;

  std::vector<double> weights(3,0.0); //initialized with zeroes

  weights.at(0) = 1.0;
  double error(0.0);

  TH1D* h_eta_rr = histograms.find("eta_rr")->second;
  TH1D* h_eta_fr = histograms.find("eta_fr")->second;
  TH1D* h_pt_rr  = histograms.find("pt_rr")->second;
  TH1D* h_pt_fr  = histograms.find("pt_fr")->second;

  // loop over number of eta bins
  // do not consider underflow, i.e. 0th bin
  //
  // NB: this works assuming the binning for h_eta_rr and h_eta_fr is the same
  //
  bool isEtaPastBinOverFlow(false);
  for ( int e = 1; e <= n_bins_eta; e++ ) {

    // if next bin is the overflow bin, and eta under question
    // is >= lower edge of o.f. bin, switch on a flag
    //
    if ( h_eta_rr->IsBinOverflow(e+1) && fabs(eta) >= h_eta_rr->GetXaxis()->GetBinLowEdge(e+1) ) {
      isEtaPastBinOverFlow = true;
    }

    // check whether the eta under question is in *this* eta range
    //
    if ( ( ( fabs(eta) >= h_eta_rr->GetXaxis()->GetBinLowEdge(e) ) && ( fabs(eta) < h_eta_rr->GetXaxis()->GetBinLowEdge(e+1) ) ) || isEtaPastBinOverFlow ) {

      // case 1) : lepton is fake: choose correct pt histogram
      //
      if ( isFakeLep ) {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
	bool isPtPastBinOverFlow(false);
        for ( int p = 1; p <= n_bins_pt_fr; p++ ) {

	  // if next bin is the overflow bin, and pt under question
	  // is >= lower edge of o.f. bin, switch on a flag
	  //
	  if ( h_pt_fr->IsBinOverflow(p+1) && pt >= h_pt_fr->GetXaxis()->GetBinLowEdge(p+1) ) {
	    isPtPastBinOverFlow = true;
	  }

	  if ( ( ( pt >= h_pt_fr->GetXaxis()->GetBinLowEdge(p) ) && ( pt < h_pt_fr->GetXaxis()->GetBinLowEdge(p+1) ) ) || isPtPastBinOverFlow ) {

	    // combine eta and pt rates
	    // (NB: if eta/pt under question are >= o.f. bin, apply the rate of the last bin before o.f.)
	    //
	    double fr_pt  = ( !isPtPastBinOverFlow )  ? h_pt_fr->GetBinContent(p)  : h_pt_fr->GetBinContent(n_bins_pt_fr);
	    double fr_eta = ( !isEtaPastBinOverFlow ) ? h_eta_fr->GetBinContent(e) : h_eta_fr->GetBinContent(n_bins_eta);

	    double fr_pt_err  = ( !isPtPastBinOverFlow )  ? h_pt_fr->GetBinError(p)  : h_pt_fr->GetBinError(n_bins_pt_fr);
	    double fr_eta_err = ( !isEtaPastBinOverFlow ) ? h_eta_fr->GetBinError(e) : h_eta_fr->GetBinError(n_bins_eta);

	    // nominal
	    //
	    weights.at(0) = ( fr_pt * fr_eta ) / fr_tot;

	    // (assuming  fr_pt,fr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back later in the def of weight...
	    //
	    error  = sqrt( (fr_eta*fr_pt_err)*(fr_eta*fr_pt_err) + (fr_pt*fr_eta_err)*(fr_pt*fr_eta_err) );

	    // up syst
	    //
	    weights.at(1) = ( (fr_pt * fr_eta) + error ) / fr_tot;

	    // down syst
	    //
	    if ( (fr_pt * fr_eta) - error > 0 ) { weights.at(2) = ( (fr_pt * fr_eta) - error ) / fr_tot;}
	    else                                { weights.at(2) = 0.0; }

	  }

	} // close loop on pT bins: fake lepton

	// lepton is real: choose correct pt histogram
	//
      } else {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
	bool isPtPastBinOverFlow(false);
        for ( int p = 1; p <= n_bins_pt_rr; p++ ) {

	  // if next bin is the overflow bin, and pt under question
	  // is >= lower edge of o.f. bin, switch on a flag
	  //
	  if ( h_pt_rr->IsBinOverflow(p+1) && pt >= h_pt_rr->GetXaxis()->GetBinLowEdge(p+1) ) {
	    isPtPastBinOverFlow = true;
	  }

	  if ( ( ( pt >= h_pt_rr->GetXaxis()->GetBinLowEdge(p) ) && ( pt < h_pt_rr->GetXaxis()->GetBinLowEdge(p+1) ) ) || isPtPastBinOverFlow ) {

	    // combine eta and pt rates
	    // (NB: if eta/pt under question are >= o.f. bin, apply the rate of the last bin before o.f.)
	    //
	    double rr_pt  = ( !isPtPastBinOverFlow )  ? h_pt_rr->GetBinContent(p)  : h_pt_rr->GetBinContent(n_bins_pt_rr);
	    double rr_eta = ( !isEtaPastBinOverFlow ) ? h_eta_rr->GetBinContent(e) : h_eta_rr->GetBinContent(n_bins_eta);

	    double rr_pt_err  = ( !isPtPastBinOverFlow )  ? h_pt_rr->GetBinError(p)  : h_pt_rr->GetBinError(n_bins_pt_rr);
	    double rr_eta_err = ( !isEtaPastBinOverFlow ) ? h_eta_rr->GetBinError(e) : h_eta_rr->GetBinError(n_bins_eta);

	    // nominal
	    //
	    weights.at(0) = ( rr_pt * rr_eta ) / rr_tot;

	    // (assuming  rr_pt,rr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back in the def of weight...
	    //
	    error  = sqrt( (rr_eta*rr_pt_err)*(rr_eta*rr_pt_err) + (rr_pt*rr_eta_err)*(rr_pt*rr_eta_err) );

	    // up syst
	    //
	    weights.at(1) = ( (rr_pt * rr_eta) + error ) / rr_tot;

	    // down syst
	    //
	    if ( (rr_pt * rr_eta) - error > 0 ) { weights.at(2) = ( (rr_pt * rr_eta) - error ) / rr_tot; }
	    else                                { weights.at(2) = 0.0; }

	  }
        } // close loop on pT bins: real lepton

      } // close check isFakeLep

    } // close check on eta bin

  } // close loop on eta bins

  // Now converting rates to the factors for the MM/FF
  //
  if ( g_debug ) { Info("calc_weights()", "Rates = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  weights.at(0) = scaleRateToEfficiency(weights.at(0));
  weights.at(1) = scaleRateToEfficiency(weights.at(1));
  weights.at(2) = scaleRateToEfficiency(weights.at(2));

  if ( g_debug ) { Info("calc_weights()", "MM/FF factor = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  return weights;
}
*/

/* ***************************************************************************************
/
/ This function saves the final MM/FF event weight, depending on the event type (TT,TL...)
/
*************************************************************************************** */

double calc_final_event_weight( std::string region, double f1, double f2, double r1, double r2 )
{

   double weight = 1.0;
   double alpha  = 1.0 / ( (r1-f1) * (r2-f2) );

   if      ( region=="TT" ) { weight = 1 - ( r1 * r2 * (1-f1) * (1-f2) * alpha ); }
   else if ( region=="TL" ) { weight = r1 * r2 * f2 * (1-f1) * alpha;  }
   else if ( region=="LT" ) { weight = r1 * r2 * f1 * (1-f2) * alpha;  }
   else if ( region=="LL" ) { weight = -1 * r1 * r2 * f1 * f2 * alpha; }

   if ( g_debug ) { Info("calc_final_event_weight()", "In region %s : \n weight = %.15f , alpha = %.15f ", region.c_str(), weight, alpha); }

   return weight;
}

/* ********************************************
/
/ Update MMWeight in TTree by reading new rates
/
******************************************** */

void recomputeMMW( std::vector<double>* MMW_new,  /* pass it by pointer, as you are going to modify it! */
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

  // real and fake rates w/ syst variations
  //
  std::vector<double> r1, r2, f1, f2;
  double r1up,r1dn, r2up, r2dn, f1up, f1dn, f2up, f2dn;

  bool isFakeLep = true;

  if ( g_verbose ) {
    Info("recomputeMMW()", "\n Lepton 1 \n flavour: %i \n pT = %.2f [GeV] \n eta = %.2f \n", lep_flavour.at(0), lep_pt.at(0)/1e3, lep_eta.at(0) );
    Info("recomputeMMW()", "\n Lepton 2 \n flavour: %i \n pT = %.2f [GeV] \n eta = %.2f \n", lep_flavour.at(1), lep_pt.at(1)/1e3, lep_eta.at(1) );
  }

  // NB: input lep_* vectors are pT-ordered.
  //
  if ( lep_flavour.at(0) == 11 ) {
     r1 = calc_weights( g_el_hist_map, lep_pt.at(0), lep_eta.at(0), !isFakeLep, g_n_el_bins_eta, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_el_fr_tot, g_el_rr_tot );
     f1 = calc_weights( g_el_hist_map, lep_pt.at(0), lep_eta.at(0), isFakeLep,  g_n_el_bins_eta, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_el_fr_tot, g_el_rr_tot );
  } else if ( lep_flavour.at(0) == 13 ) {
     r1 = calc_weights( g_mu_hist_map, lep_pt.at(0), lep_eta.at(0), !isFakeLep, g_n_mu_bins_eta, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_mu_fr_tot, g_mu_rr_tot );
     f1 = calc_weights( g_mu_hist_map, lep_pt.at(0), lep_eta.at(0), isFakeLep,  g_n_mu_bins_eta, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_mu_fr_tot, g_mu_rr_tot );
  }

  if ( lep_flavour.at(1) == 11 ) {
     r2 = calc_weights( g_el_hist_map, lep_pt.at(1), lep_eta.at(1), !isFakeLep, g_n_el_bins_eta, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_el_fr_tot, g_el_rr_tot );
     f2 = calc_weights( g_el_hist_map, lep_pt.at(1), lep_eta.at(1), isFakeLep,  g_n_el_bins_eta, g_n_el_bins_pt_fr, g_n_el_bins_pt_rr, g_el_fr_tot, g_el_rr_tot );
  } else if ( lep_flavour.at(1) == 13 ) {
     r2 = calc_weights( g_mu_hist_map, lep_pt.at(1), lep_eta.at(1), !isFakeLep, g_n_mu_bins_eta, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_mu_fr_tot, g_mu_rr_tot );
     f2 = calc_weights( g_mu_hist_map, lep_pt.at(1), lep_eta.at(1), isFakeLep,  g_n_mu_bins_eta, g_n_mu_bins_pt_fr, g_n_mu_bins_pt_rr, g_mu_fr_tot, g_mu_rr_tot );
  }

  if ( g_debug ) {
    Info("recomputeMMW()", "\n Lepton 1 \n flavour: %i \n pT = %.2f [GeV] \n eta = %.2f \n ****** \n Nominal real and fake rates: \n r1 = %f , f1 = %f ", lep_flavour.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), r1.at(0), f1.at(0) );
    Info("recomputeMMW()", "\n Lepton 2 \n flavour: %i \n pT = %.2f [GeV] \n eta = %.2f \n ****** \n Nominal real and fake rates: \n r2 = %f , f2 = %f ", lep_flavour.at(1), lep_pt.at(1)/1e3, lep_eta.at(1), r2.at(0), f2.at(0) );
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
        Warning("recomputeMMW()", "Warning! The Matrix Method cannot be applied because : \n r1 = %f , r2 = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", r1.at(0), r2.at(0),  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
        Warning("recomputeMMW()", "applying MM weight = 0 ...");
      }
  } else {
      // calculate nominal MM weight
      //
      mm_weight      = calc_final_event_weight( region, f1.at(0), f2.at(0), r1.at(0), r2.at(0) );
  }

  // update branch for nominal MM weight
  //
  MMW_new->at(0) = mm_weight;

  if ( g_debug ) { Info("recomputeMMW()", "MM final weight = %f ", mm_weight ); }

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

    // rup syst
    //
    MMW_new->at(1) = ( calc_final_event_weight( region, f1.at(0), f2.at(0), r1up, r2up ) / mm_weight );
    // fdn syst
    //
    MMW_new->at(4) = ( calc_final_event_weight( region, f1dn, f2dn, r1.at(0), r2.at(0) ) / mm_weight );

    if ( (r1dn > f1.at(0)) && (r2dn > f2.at(0)) ) {
      // rdn syst
      //
      MMW_new->at(2) = ( calc_final_event_weight(region, f1.at(0), f2.at(0), r1dn, r2dn) / mm_weight );
    } else {
      if ( g_debug ) {
         Warning("recomputeMMW()", "Warning! Systematic MMWeight_rdn cannot be calculated because : \n r1dn = %f , r2dn = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", r1dn, r2dn,  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
      }
    }

    if ( (r1.at(0) > f1up) && (r2.at(0) > f2up) ) {
      // fup syst
      //
      MMW_new->at(3) = ( calc_final_event_weight(region, f1up, f2up, r1.at(0),  r2.at(0)) / mm_weight );
    } else {
      if ( g_debug ) {
         Warning("recomputeMMW()", "Warning! Systematic MMWeight_fup cannot be calculated because : \n r1dn = %f , r2dn = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", r1dn, r2dn,  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
      }
    }
  }

}

/* ******************
/
/ The 'main' function
/
****************** */

void modifyttree_MM(std::string filename = "input.root", std::string  NENTRIES = "ALL", std::string treename= "physics", std::string newfilename = "output.root")
{
  // This script loads a tree, clones it, removes a branch and substitutes it with another.
  // The branch can also have the same name and in this way you can change for example the type of the variable or the content.

  Info("modifytree()","Starting off...");

  //Get old file, old tree and set top branch address
  //
  TFile *oldfile = new TFile(filename.c_str());
  TTree *oldtree = (TTree*)oldfile->Get(treename.c_str());

  Long64_t nentries;

  if ( NENTRIES == "ALL" ) {
     nentries = oldtree->GetEntries();
  } else {
     std::stringstream ss; ss << NENTRIES;
     int n_e;              ss >> n_e;
     nentries = n_e;
  }

  // TO BE MODIFIED ACCORDINGLY TO YOUR NEEDS (name and type of the variables)
  //
  std::string old_eventNumber_name("eventNumber");
  std::string old_nlep_name("nlep");
  std::string old_isSS01_name("isSS01");
  std::string old_isTT_name("isTT");
  std::string old_isTL_name("isTL");
  std::string old_isLT_name("isLT");
  std::string old_isLL_name("isLL");
  std::string old_MMWeight_name("MMWeight");
  std::string old_lep_pt_name("lep_pt");
  std::string old_lep_eta_name("lep_eta");
  std::string old_lep_flavour_name("lep_flavour");
  std::string old_lep_isT_name("lep_isTightSelected");

  Long64_t               eventNumber_old; eventNumber_old = -1;
  Int_t                  nlep_old;        nlep_old = -1;
  Int_t                  isSS01_old;      isSS01_old = -1;
  Int_t 	         isTT_old;        isTT_old = -1;
  Int_t 	         isTL_old;        isTL_old = -1;
  Int_t 	         isLT_old;        isLT_old = -1;
  Int_t 	         isLL_old;        isLL_old = -1;
  std::vector<double>*   MMWeight_old;    MMWeight_old = 0;
  std::vector<float>*    lep_pt_old;      lep_pt_old = 0;
  std::vector<float>*    lep_eta_old;     lep_eta_old = 0;
  std::vector<int>*      lep_flavour_old; lep_flavour_old = 0;
  std::vector<int>*      lep_isT_old;     lep_isT_old = 0;

  // List of old branches
  //
  TBranch        *b_eventNumber = 0;      //!
  TBranch	 *b_nlep_old    = 0;      //!
  TBranch        *b_isSS01_old = 0;       //!
  TBranch	 *b_isTT = 0;    	  //!
  TBranch	 *b_isTL = 0;    	  //!
  TBranch	 *b_isLT = 0;    	  //!
  TBranch	 *b_isLL = 0;    	  //!
  TBranch        *b_MMWeight_old = 0;     //!
  TBranch	 *b_lep_pt_old = 0;       //!
  TBranch	 *b_lep_eta_old = 0;      //!
  TBranch	 *b_lep_flavour_old = 0;  //!
  TBranch        *b_lep_isT_old = 0;      //!

  // Before cloning input TTree, tell ROOT to process all the old branches,
  // except for the one you want to change
  //
  oldtree->SetBranchStatus("*",1);
  oldtree->SetBranchStatus(old_MMWeight_name.c_str(),0);

  std::string new_MMWeight_name("MMWeight");
  std::vector<double>*   MMWeight_new;  MMWeight_new = 0;

  // Create a new file + a clone of old tree in new file
  //
  TFile *newfile = new TFile(newfilename.c_str(),"RECREATE");
  TTree *newtree = oldtree->CloneTree(0); //clone only the structure

  // MUST re-activate the branch(es) previously deactivated before calling SetBranchAddress()!!
  //
  oldtree->SetBranchStatus(old_MMWeight_name.c_str(),1);

  // Get branches from input TTree
  //
  oldtree->SetBranchAddress(old_eventNumber_name.c_str(), &eventNumber_old, &b_eventNumber);
  oldtree->SetBranchAddress(old_nlep_name.c_str(), &nlep_old, &b_nlep_old);
  oldtree->SetBranchAddress(old_isSS01_name.c_str(), &isSS01_old, &b_isSS01_old);
  oldtree->SetBranchAddress(old_isTT_name.c_str(), &isTT_old, &b_isTT);
  oldtree->SetBranchAddress(old_isTL_name.c_str(), &isTL_old, &b_isTL);
  oldtree->SetBranchAddress(old_isLT_name.c_str(), &isLT_old, &b_isLT);
  oldtree->SetBranchAddress(old_isLL_name.c_str(), &isLL_old, &b_isLL);
  oldtree->SetBranchAddress(old_MMWeight_name.c_str(), &MMWeight_old, &b_MMWeight_old);
  oldtree->SetBranchAddress(old_lep_pt_name.c_str(), &lep_pt_old, &b_lep_pt_old);
  oldtree->SetBranchAddress(old_lep_eta_name.c_str(), &lep_eta_old, &b_lep_eta_old);
  oldtree->SetBranchAddress(old_lep_flavour_name.c_str(), &lep_flavour_old, &b_lep_flavour_old);
  oldtree->SetBranchAddress(old_lep_isT_name.c_str(), &lep_isT_old, &b_lep_isT_old);

  //Info("modifytree()", "--> lep_pt before SetBranchAddress() %p\n", oldtree->GetBranch(old_lep_pt_name.c_str())->GetAddress());
  //Info("modifytree()", "--> lep_pt after SetBranchAddress() %p\n", oldtree->GetBranch(old_lep_pt_name.c_str())->GetAddress());

  // Set the "new" branches in the output TTree
  //
  newtree->Branch(new_MMWeight_name.c_str(), &MMWeight_new);

  // read r/f rates from ROOT histograms
  //
  //std::string RR_dir("GOOD_STUFF/OutputPlots_MMRates_v021_Madgraph_Observed");
  //std::string RR_dir("OutputPlots_MMRates_v025");
  //std::string RR_dir("OutputPlots_MMRates_v028_FINAL");
  std::string RR_dir("OutputPlots_MMRates_v027_20GeVpT");
 
  // when using ch-flip rate as RR (for electrons)
  //std::string RR_dir("PLOTS/PLOTS_013/TEST_13F_2/OutputPlots_ChFlipBkgRates_13F");
  //std::string FR_dir("PLOTS/PLOTS_013/TEST_13F_2/OutputPlots_NonPromptBkgRates_13F");

  Info("modifytree()","Reading r/f rates from ROOT file(s)..");
  read_rates(RR_dir);
  //read_rates(RR_dir,FR_dir); // pass also FR_dir if different than RR_dir

  // Loop over entries in TTree
  //
  Info("modifytree()","Begin loop on input tree entries...\n");
  int count_bad(0);
  Long64_t i = 0;
  for ( ; i < nentries; i++ ) {

    // Print out every N events to see where we are
    //
    if ( i > 0 && ( static_cast<int>(i) % 20000 == 0 ) ) { Info("modifytree()","\t Processed %lld entries",i); }

    // Now, in the input tree, reset to 1 the status of the branches you
    // deactivated before cloning
    //
    if ( !oldtree->GetBranchStatus(old_MMWeight_name.c_str()) ) {
      oldtree->SetBranchStatus(old_MMWeight_name.c_str(),1);
    }

    oldtree->GetEntry(i);

    if ( g_debug ) { Info("modifytree()","\t Processing entry: %lld - eventNumber: %lli \n",i, eventNumber_old); }

    // A security check...
    //
    if ( !MMWeight_old ) {
      Info("modifytree()","\t --> MMWeight_old is NULL!! Skipping event...  \n");
      ++count_bad;
      continue;
    }

    // To start off, copy the old branch into the new
    // (then it will be overridden , if necessary)
    //
    *MMWeight_new = *MMWeight_old;

    // and now, recompute the MM weights!
    // Do it only for 2lepSS events
    //
    if ( nlep_old == 2 && isSS01_old == 1 ) {

      if ( g_debug ) {
        int idx_old(0);
        for ( const auto& itr : *MMWeight_old ) {
    	  Info("modifytree()","\t\t OLD MMWeight[%i] = %f", idx_old, itr );
    	  ++idx_old;
        }
      }

      // Just recompute MMWeight w/ new rates
      //
      //recomputeMMW( MMWeight_new, isTT_old, isTL_old, isLT_old, isLL_old, *lep_pt_old, *lep_eta_old, *lep_flavour_old );

      bool TT =  ( lep_isT_old->at(0) == 1 && lep_isT_old->at(1) == 1 );
      bool TL =  ( lep_isT_old->at(0) == 1 && lep_isT_old->at(1) == 0 );
      bool LT =  ( lep_isT_old->at(0) == 0 && lep_isT_old->at(1) == 1 );
      bool LL =  ( lep_isT_old->at(0) == 0 && lep_isT_old->at(1) == 0 );

      recomputeMMW( MMWeight_new, TT, TL, LT, LL, *lep_pt_old, *lep_eta_old, *lep_flavour_old );

      if ( g_debug ) {
        int idx_new(0);
        for ( const auto& itr : *MMWeight_new ) {
    	  Info("modifytree()","\t\t NEW MMWeight[%i] = %f", idx_new, itr );
    	  ++idx_new;
        }
      }

      if ( g_debug ) { Info("modifytree()","\n\n"); }

    }

    // to avoid overriding new branch (old has same name) ?
    oldtree->SetBranchStatus(old_MMWeight_name.c_str(),0);

    newtree->Fill();

  }

  Info("modifytree()","End of loop!\n ---> total number of processed events: %lld \n ---> number of skipped events: %i \n", i, count_bad );

  newfile->Write();
  newfile->Close();

  // Since we passed the address of a local variable we need
  // to remove it.
  oldtree->ResetBranchAddresses();

  delete oldfile;
  delete newfile;

  // //delete a branch:
  // TFile f("myfile.root","update");
  // TTree *T = (TTree*)f.Get("treename");
  // TBranch *b = T->GetBranch("name of branch to delete");
  // T->GetListOfBranches()->Remove(b);
  // TLeaf* l= T->GetLeaf("name of branch to delete");
  // T->GetListOfLeaves()->Remove(l)
  // T->Write();

}
