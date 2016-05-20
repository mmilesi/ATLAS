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

  std::string rate_type = ( !g_doClosure ) ? "observed" : "expected";

  // ELECTRONS
  //
  std::string histname_el_eta_rr   = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_Rate_" + rate_type : "El_ProbeEta_Fake_Rate_" + rate_type;
  std::string histname_el_pt_rr    = ( !CFRateAsRR ) ?  "El_ProbePt_Real_Rate_" + rate_type : "El_ProbePt_Fake_Rate_" + rate_type;
  std::string histname_el_eta_r_T  = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_T_" + rate_type : "El_ProbeEta_Fake_Rate_" + rate_type;
  std::string histname_el_eta_r_L  = ( !CFRateAsRR ) ?  "El_ProbeEta_Real_L_" + rate_type : "El_ProbeEta_Fake_Rate_" + rate_type;
  std::string histname_el_pt_r_T   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_T_" + rate_type : "El_ProbePt_Fake_Rate_" + rate_type;
  std::string histname_el_pt_r_L   = ( !CFRateAsRR ) ?  "El_ProbePt_Real_L_" + rate_type : "El_ProbePt_Fake_Rate_" + rate_type;

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
  std::string histname_mu_eta_rr   = "Mu_ProbeEta_Real_Rate_" + rate_type;
  std::string histname_mu_pt_rr    = "Mu_ProbePt_Real_Rate_" + rate_type;
  std::string histname_mu_eta_r_T  = "Mu_ProbeEta_Real_T_" + rate_type;
  std::string histname_mu_eta_r_L  = "Mu_ProbeEta_Real_L_" + rate_type;
  std::string histname_mu_pt_r_T   = "Mu_ProbePt_Real_T_" + rate_type;
  std::string histname_mu_pt_r_L   = "Mu_ProbePt_Real_L_" + rate_type;

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
  std::string histname_el_eta_fr   = "El_ProbeEta_Fake_Rate_" + rate_type;
  std::string histname_el_pt_fr    = "El_ProbePt_Fake_Rate_" + rate_type;
  std::string histname_el_eta_f_T  = "El_ProbeEta_Fake_T_" + rate_type;
  std::string histname_el_eta_f_L  = "El_ProbeEta_Fake_L_" + rate_type;
  std::string histname_el_pt_f_T   = "El_ProbePt_Fake_T_" + rate_type;
  std::string histname_el_pt_f_L   = "El_ProbePt_Fake_L_" + rate_type;

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
  std::string histname_mu_eta_fr   = "Mu_ProbeEta_Fake_Rate_" + rate_type;
  std::string histname_mu_pt_fr    = "Mu_ProbePt_Fake_Rate_" + rate_type;
  std::string histname_mu_eta_f_T  = "Mu_ProbeEta_Fake_T_" + rate_type;
  std::string histname_mu_eta_f_L  = "Mu_ProbeEta_Fake_L_" + rate_type;
  std::string histname_mu_pt_f_T   = "Mu_ProbePt_Fake_T_" + rate_type;
  std::string histname_mu_pt_f_L   = "Mu_ProbePt_Fake_L_" + rate_type;

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
  g_n_el_bins_eta   =  hist_el_eta_rr->GetNbinsX()+1;
  g_n_mu_bins_eta   =  hist_mu_eta_rr->GetNbinsX()+1;

  // pt hist has two different binning for r/f
  //
  g_n_el_bins_pt_rr =  hist_el_pt_rr->GetNbinsX()+1;
  g_n_el_bins_pt_fr =  hist_el_pt_fr->GetNbinsX()+1;
  g_n_mu_bins_pt_rr =  hist_mu_pt_rr->GetNbinsX()+1;
  g_n_mu_bins_pt_fr =  hist_mu_pt_fr->GetNbinsX()+1;

  // normalisation factor is the same for eta and pt r/f histograms: use eta
  //
  g_el_rr_tot = ( hist_el_eta_r_T->Integral(1,hist_el_eta_r_T->GetNbinsX()+1) ) / ( hist_el_eta_r_L->Integral(1,hist_el_eta_r_L->GetNbinsX()+1) );
  g_el_fr_tot = ( hist_el_eta_f_T->Integral(1,hist_el_eta_f_T->GetNbinsX()+1) ) / ( hist_el_eta_f_L->Integral(1,hist_el_eta_f_L->GetNbinsX()+1) );
  g_mu_rr_tot = ( hist_mu_eta_r_T->Integral(1,hist_mu_eta_r_T->GetNbinsX()+1) ) / ( hist_mu_eta_r_L->Integral(1,hist_mu_eta_r_L->GetNbinsX()+1) );
  g_mu_fr_tot = ( hist_mu_eta_f_T->Integral(1,hist_mu_eta_f_T->GetNbinsX()+1) ) / ( hist_mu_eta_f_L->Integral(1,hist_mu_eta_f_L->GetNbinsX()+1) );

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

  // Now converting rates to the efficiencies for the MM/FF
  //
  if ( g_verbose ) { Info("calc_weights()", "Rates = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  weights.at(0) = scaleRateToEfficiency(weights.at(0));
  weights.at(1) = scaleRateToEfficiency(weights.at(1));
  weights.at(2) = scaleRateToEfficiency(weights.at(2));

  if ( g_verbose ) { Info("calc_weights()", "MM/FF efficiency = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  return weights;
}

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
  MMW_out->at(0) = mm_weight;

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
    MMW_out->at(1) = ( calc_final_event_weight( region, f1.at(0), f2.at(0), r1up, r2up ) / mm_weight );
    // fdn syst
    //
    MMW_out->at(4) = ( calc_final_event_weight( region, f1dn, f2dn, r1.at(0), r2.at(0) ) / mm_weight );

    if ( (r1dn > f1.at(0)) && (r2dn > f2.at(0)) ) {
      // rdn syst
      //
      MMW_out->at(2) = ( calc_final_event_weight(region, f1.at(0), f2.at(0), r1dn, r2dn) / mm_weight );
    } else {
      if ( g_debug ) {
         Warning("recomputeMMW()", "Warning! Systematic MMWeight_rdn cannot be calculated because : \n r1dn = %f , r2dn = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", r1dn, r2dn,  f1.at(0), f2.at(0), lep_pt.at(0)/1e3, lep_eta.at(0), lep_pt.at(1)/1e3, lep_eta.at(1));
      }
    }

    if ( (r1.at(0) > f1up) && (r2.at(0) > f2up) ) {
      // fup syst
      //
      MMW_out->at(3) = ( calc_final_event_weight(region, f1up, f2up, r1.at(0),  r2.at(0)) / mm_weight );
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

void modifyttree_MM(std::string filename, std::string outfilename, std::string addWeight, std::string RR_dir, std::string doClosure,
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

  if ( addWeight.compare("NO") == 0 ) {

    // TO BE MODIFIED ACCORDINGLY TO YOUR NEEDS (name and type of the variables)
    //
    std::string in_eventNumber_name("eventNumber");
    std::string in_nlep_name("nlep");
    std::string in_isSS01_name("isSS01");
    std::string in_isTT_name("is_T_T");
    std::string in_isTL_name("is_T_AntiT");
    std::string in_isLT_name("is_AntiT_T");
    std::string in_isLL_name("is_AntiT_AntiT");
    std::string in_MMWeight_name("MMWeight");
    std::string in_lep_pt_name("lep_pt");
    std::string in_lep_eta_name("lep_caloCluster_eta");
    std::string in_lep_flavour_name("lep_flavour");
    std::string in_lep_isT_name("lep_isTightSelected");

    Long64_t               eventNumber_in; eventNumber_in = -1;
    Int_t                  nlep_in;        nlep_in = -1;
    Int_t                  isSS01_in;      isSS01_in = -1;
    Int_t 	           isTT_in;        isTT_in = -1;
    Int_t 	           isTL_in;        isTL_in = -1;
    Int_t 	           isLT_in;        isLT_in = -1;
    Int_t 	           isLL_in;        isLL_in = -1;
    std::vector<double>*   MMWeight_in;    MMWeight_in = 0;
    std::vector<float>*    lep_pt_in;      lep_pt_in = 0;
    std::vector<float>*    lep_eta_in;     lep_eta_in = 0;
    std::vector<int>*      lep_flavour_in; lep_flavour_in = 0;
    std::vector<int>*      lep_isT_in;     lep_isT_in = 0;

    // List of in branches
    //
    TBranch      *b_eventNumber = 0;      //!
    TBranch	 *b_nlep    = 0;          //!
    TBranch      *b_isSS01 = 0;           //!
    TBranch	 *b_is_T_T = 0;    	  //!
    TBranch	 *b_is_T_AntiT = 0;    	  //!
    TBranch	 *b_is_AntiT_T = 0;    	  //!
    TBranch	 *b_is_AntiT_AntiT = 0;   //!
    TBranch      *b_MMWeight = 0;         //!
    TBranch	 *b_lep_pt = 0;           //!
    TBranch	 *b_lep_eta = 0;          //!
    TBranch	 *b_lep_flavour = 0;      //!
    TBranch      *b_lep_isTightSelected = 0;   //!

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
    intree->SetBranchAddress(in_eventNumber_name.c_str(), &eventNumber_in, &b_eventNumber);
    intree->SetBranchAddress(in_nlep_name.c_str(), &nlep_in, &b_nlep);
    intree->SetBranchAddress(in_isSS01_name.c_str(), &isSS01_in, &b_isSS01);
    intree->SetBranchAddress(in_isTT_name.c_str(), &isTT_in, &b_is_T_T);
    intree->SetBranchAddress(in_isTL_name.c_str(), &isTL_in, &b_is_T_AntiT);
    intree->SetBranchAddress(in_isLT_name.c_str(), &isLT_in, &b_is_AntiT_T);
    intree->SetBranchAddress(in_isLL_name.c_str(), &isLL_in, &b_is_AntiT_AntiT);
    intree->SetBranchAddress(in_MMWeight_name.c_str(), &MMWeight_in, &b_MMWeight);
    intree->SetBranchAddress(in_lep_pt_name.c_str(), &lep_pt_in, &b_lep_pt);
    intree->SetBranchAddress(in_lep_eta_name.c_str(), &lep_eta_in, &b_lep_eta);
    intree->SetBranchAddress(in_lep_flavour_name.c_str(), &lep_flavour_in, &b_lep_flavour);
    intree->SetBranchAddress(in_lep_isT_name.c_str(), &lep_isT_in, &b_lep_isTightSelected);

    //Info("modifytree_MM()", "--> lep_pt before SetBranchAddress() %p\n", intree->GetBranch(in_lep_pt_name.c_str())->GetAddress());
    //Info("modifytree_MM()", "--> lep_pt after SetBranchAddress() %p\n", intree->GetBranch(in_lep_pt_name.c_str())->GetAddress());

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
      // Do it only for 2lepSS events
      //
      if ( nlep_in == 2 && isSS01_in == 1 ) {

        if ( g_debug ) {
          int idx_in(0);
          for ( const auto& itr : *MMWeight_in ) {
      	  Info("modifytree_MM()","\t\t IN MMWeight[%i] = %f", idx_in, itr );
      	  ++idx_in;
          }
        }

        int TT =  ( lep_isT_in->at(0) == 1 && lep_isT_in->at(1) == 1 );
        int TL =  ( lep_isT_in->at(0) == 1 && lep_isT_in->at(1) == 0 );
        int LT =  ( lep_isT_in->at(0) == 0 && lep_isT_in->at(1) == 1 );
        int LL =  ( lep_isT_in->at(0) == 0 && lep_isT_in->at(1) == 0 );

        recomputeMMW( MMWeight_out, TT, TL, LT, LL, *lep_pt_in, *lep_eta_in, *lep_flavour_in );

        if ( g_debug ) {
          int idx_out(0);
          for ( const auto& itr : *MMWeight_out ) {
      	  Info("modifytree_MM()","\t\t OUT MMWeight[%i] = %f", idx_out, itr );
      	  ++idx_out;
          }
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
    std::string in_lep_eta_0_name("lep_EtaBE2_0");
    std::string in_lep_eta_1_name("lep_EtaBE2_1");
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
    intree->SetBranchAddress(in_lep_eta_0_name.c_str(), &lep_eta_0_in, &b_lep_EtaBE2_0);
    intree->SetBranchAddress(in_lep_eta_1_name.c_str(), &lep_eta_1_in, &b_lep_EtaBE2_1);
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
          Info("modifytree_MM()","\t\t Default MMWeight[%i] = %f", idx, itr );
          ++idx;
        }
      }

      // and now, recompute the MM weights!
      // Do it only for 2lepSS events
      //
      if ( nlep_in == 2 && isSS01_in == 1 ) {

        int TT =  ( isTT_in );
        int TL =  ( isTL_in );
        int LT =  ( isLT_in );
        int LL =  ( isLL_in );

        // Create these vectors to avoid changing the signature of the function
        //
        std::vector<float> lep_pt_vec      = { lep_pt_0_in, lep_pt_1_in };
        std::vector<float> lep_eta_vec     = { lep_eta_0_in, lep_eta_1_in };
        std::vector<int> lep_flavour_vec   = { abs(static_cast<int>(lep_flavour_0_in)), abs(static_cast<int>(lep_flavour_1_in)) };

        recomputeMMW( &MMWeight_out, TT, TL, LT, LL, lep_pt_vec, lep_eta_vec, lep_flavour_vec );

        if ( g_debug ) {
          int idx_out(0);
          for ( const auto& itr : MMWeight_out ) {
	      Info("modifytree_MM()","\t\t OUT MMWeight[%i] = %f", idx_out, itr );
	      ++idx_out;
          }
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

