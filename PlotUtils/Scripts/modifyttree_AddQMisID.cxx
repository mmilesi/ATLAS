#include <vector>
#include <math.h>
#include <iostream>
#include <sstream>

#include "TError.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TTree.h"

/* ***************
/
/ Global variables
/
*************** */

bool g_debug(false);
bool g_verbose(false);

std::map< std::string, TH2D* > g_QMisID_hist_map;  

/* ***********************************
/
/ Function to get histograms from file
/
*********************************** */

template<typename T>
T* get_object( TFile& file, const std::string& name ) {
  T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
  if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
  return obj;
}

/* **********************************************
/
/ Function to get rates
/
********************************************** */

void getRates(std::string& input_path, std::string& filename_T, std::string& filename_AntiT)
{
  
  if ( input_path.back() != '/' ) { input_path += "/"; }

  std::string path_AntiT = input_path + filename_AntiT;
  std::string path_T     = input_path + filename_T;

  TFile *file_AntiT = TFile::Open(path_AntiT.c_str());
  TFile *file_T     = TFile::Open(path_T.c_str());
  
  if ( !file_AntiT->IsOpen() ) {
    SysError("getRates()", "Failed to open ROOT file from path: %s . Aborting", path_AntiT.c_str() );
    exit(-1);
  }
  if ( !file_T->IsOpen() ) {
    SysError("getRates()", "Failed to open ROOT file from path: %s . Aborting", path_T.c_str() );
    exit(-1);
  }  
  
  Info("readQMisIDRates()", "Successfully opened ROOT files with QMisID rates from path:\n AntiT --> %s \n T --> %s", path_AntiT.c_str(), path_T.c_str() );

  TH2D *hist_QMisID_AntiT = get_object<TH2D>( *file_AntiT, "Rates" );
  TH2D *hist_QMisID_T     = get_object<TH2D>( *file_T, "Rates" );

  // fill a map for later usage
  //
  g_QMisID_hist_map["AntiT"] = hist_QMisID_AntiT;
  g_QMisID_hist_map["T"]     = hist_QMisID_T;
 
}


/* ****************************************************************************************
/
/ This function checks in which 2D bin *this* electron falls, and read the rate accordingly
/
**************************************************************************************** */

void readRatesAndError( TH2D* rate_map, TH1D* proj_X, TH1D* proj_Y,
                        const float& x, const float& y,
			float& r, float& r_up, float& r_dn )
{

  float this_low_edge(-999.0),this_up_edge(-999.0);

  int xbin_nr(-1), ybin_nr(-1);

  // Loop over the projections, and keep track of the bin number where (x,y) is found
  //
  for ( int xbin = 0; xbin < proj_X->GetNbinsX()+1; ++xbin  ) {

    this_low_edge = proj_X->GetXaxis()->GetBinLowEdge(xbin);
    this_up_edge  = proj_X->GetXaxis()->GetBinLowEdge(xbin+1);

    if ( fabs(x) >= this_low_edge && fabs(x) < this_up_edge ) {
    
      if ( g_debug ) { Info("readRatesAndError()","\t\t x = %.2f found in %i-th bin", x, xbin ); }
      
      xbin_nr = proj_X->GetBin(xbin);
     
      break;
    }

  }
  for ( int ybin = 0; ybin < proj_Y->GetNbinsX()+1; ++ ybin ) {

    this_low_edge = proj_Y->GetXaxis()->GetBinLowEdge(ybin);
    this_up_edge  = proj_Y->GetXaxis()->GetBinLowEdge(ybin+1);

    if ( y >= this_low_edge && y < this_up_edge ) {
      
      if ( g_debug ) { Info("readRatesAndError()","\t\t y = %.2f found in %i-th bin", y, ybin ); }
      
      ybin_nr = proj_Y->GetBin(ybin);

      break;
    }

  }

  if ( g_debug ) { Info("readRatesAndError()","\t\t coordinates of efficiency bin = (%i,%i)", xbin_nr, ybin_nr ); }

  // Now get the NOMINAL rate via global bin number (x,y)

  r = rate_map->GetBinContent( rate_map->GetBin( xbin_nr, ybin_nr ) );

  // Get the UP and DOWN variations
  //
  // QUESTION: Why the hell ROOT has GetBinErrorUp and GetBinErrorLow for TH2 ??
  // They seem to give always the same result...
  //
  r_up = r + rate_map->GetBinErrorUp( rate_map->GetBin( xbin_nr, ybin_nr ) );
  r_dn = r - rate_map->GetBinErrorUp( rate_map->GetBin( xbin_nr, ybin_nr ) );
  r_dn = ( r_dn > 0.0 ) ? r_dn : 0.0;

}

/* ***********************************************************
/
/ Function to calculate QMisID weights and their unceratinties
/
*********************************************************** */

void QMisIDWeightCalculator (std::vector<float>* weights, 
                             const float& elA_eta, const float& elA_pt, const bool& elA_isT,
                             const float& elB_eta, const float& elB_pt, const bool& elB_isT )
{

  bool elA = ( fabs(elA_eta) < 2.5 && elA_pt >= 0.0 );
  bool elB = ( fabs(elB_eta) < 2.5 && elB_pt >= 0.0 );

  // If there are no electrons, return 
  //
  if ( !elA && !elB ) { return; }

  float rA(0.0), rA_up(0.0), rA_dn(0.0), rB(0.0), rB_up(0.0), rB_dn(0.0);

  // Get the 2D histogram from the map
  //
  TH2D* twoD_rates_T     = ( g_QMisID_hist_map.find("T")->second );
  TH2D* twoD_rates_AntiT = ( g_QMisID_hist_map.find("AntiT")->second );
  
  // Make X and Y projections of the twoD histogram with the rates
  //
  TH1D* proj_eta_T     = twoD_rates_T->ProjectionX("proj_eta_T");
  TH1D* proj_pt_T      = twoD_rates_T->ProjectionY("proj_pt_T");
  TH1D* proj_eta_AntiT = twoD_rates_AntiT->ProjectionX("proj_eta_AntiT");
  TH1D* proj_pt_AntiT  = twoD_rates_AntiT->ProjectionY("proj_pt_AntiT");

  // Look at elA first...
  //
  if ( elA_isT ) {
    readRatesAndError(twoD_rates_T, proj_eta_T, proj_pt_T, elA_eta, elA_pt, rA, rA_up, rA_dn);
  } else {
    readRatesAndError(twoD_rates_AntiT, proj_eta_AntiT, proj_pt_AntiT, elA_eta, elA_pt, rA, rA_up, rA_dn);
  }

  // .. and now at elB (if any...otherwise rB weights will be zero by default)
  //
  if ( elB ) {
    if ( elB_isT ) {
      readRatesAndError(twoD_rates_T, proj_eta_T, proj_pt_T, elB_eta, elB_pt, rB, rB_up, rB_dn);
    } else {
      readRatesAndError(twoD_rates_AntiT, proj_eta_AntiT, proj_pt_AntiT, elB_eta, elB_pt, rB, rB_up, rB_dn);
    }
  }

  if ( g_debug ) { 
    Info("QMisIDWeightCalculator()","\t rA = %f ( up = %f, dn = %f )", rA, rA_up, rA_dn ); 
    Info("QMisIDWeightCalculator()","\t rB = %f ( up = %f, dn = %f )", rB, rB_up, rB_dn ); 
  }

  // Finally, store the event weight + variations
  //
  weights->at(0) = ( rA + rB - 2.0 * rA * rB ) / ( 1.0 - rA - rB + 2.0 * rA * rB ) ;
  weights->at(1) = ( rA_up + rB_up - 2.0 * rA_up * rB_up ) / ( 1.0 - rA_up - rB_up + 2.0 * rA_up * rB_up );
  weights->at(2) = ( rA_dn + rB_dn - 2.0 * rA_dn * rB_dn ) / ( 1.0 - rA_dn - rB_dn + 2.0 * rA_dn * rB_dn );

}

/* ******************
/
/ The 'main' function
/
****************** */

void modifyttree_AddQMisID(std::string filename = "input.root", std::string  NENTRIES = "ALL", std::string treename= "physics", std::string newfilename = "output.root",std::string useGroupNTup = "")
{
  
  // This script loads a tree, clones it, removes a branch and substitutes it with another.
  // The branch can also have the same name and in this way you can change for example the type of the variable or the content.

  Info("modifytree()","Starting off...");

  //Get old file, old tree and set top branch address
  //
  TFile *oldfile = new TFile(filename.c_str());
  TTree *oldtree = dynamic_cast<TTree*>( oldfile->Get(treename.c_str()) );

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
  std::string old_eventNumber_name = ( useGroupNTup.empty() ) ? "eventNumber" : "EventNumber";
  std::string old_nel_name         = ( useGroupNTup.empty() ) ? "nel" : "nelectrons";
  std::string old_el_pt_name(""); 
  std::string old_el_eta_name("");
  std::string old_el_isT_name("");
  std::string old_lep0_ID_name("");
  std::string old_lep1_ID_name("");
  std::string old_lep0_pt_name(""); 
  std::string old_lep1_pt_name(""); 
  std::string old_lep0_eta_name("");
  std::string old_lep1_eta_name("");
  std::string old_lep0_isT_name("");
  std::string old_lep1_isT_name("");
  
  if ( useGroupNTup.empty() ) {
    std::string old_el_pt_name   = "el_pt";
    std::string old_el_eta_name  = "el_caloCluster_eta";
    std::string old_el_isT_name  = "el_isTightSelected";
  } else {
  
    old_lep0_ID_name  = "lep_ID_0"; 
    old_lep1_ID_name  = "lep_ID_1";   
    old_lep0_pt_name  = "lep_Pt_0"; 
    old_lep1_pt_name  = "lep_Pt_1"; 
    old_lep0_eta_name = "lep_EtaBE2_0";
    old_lep1_eta_name = "lep_EtaBE2_1";
    old_lep0_isT_name = "lep_isTightSelected_0";
    old_lep1_isT_name = "lep_isTightSelected_1";
  }

  ULong64_t             eventNumber_old; eventNumber_old = -1;
  Int_t                 nel_old;        nel_old = -1;
  std::vector<float>*   el_pt_old;      el_pt_old = 0;
  std::vector<float>*   el_eta_old;     el_eta_old = 0;
  std::vector<int>*     el_isT_old;     el_isT_old = 0;
  Float_t lep0_ID_old;  lep0_ID_old = 0.0; 
  Float_t lep1_ID_old; 	lep1_ID_old = 0.0;   
  Float_t lep0_pt_old;  lep0_pt_old = -1.0; 
  Float_t lep1_pt_old; 	lep1_pt_old = -1.0; 
  Float_t lep0_eta_old;	lep0_eta_old = -999.0;
  Float_t lep1_eta_old;	lep1_eta_old = -999.0;
  Char_t  lep0_isT_old;	lep0_isT_old = 0;
  Char_t  lep1_isT_old;	lep1_isT_old = 0;

  // List of old branches
  //
  TBranch        *b_eventNumber = 0;  //!
  TBranch	 *b_nel_old    = 0;   //!
  TBranch	 *b_el_pt_old = 0;    //!
  TBranch	 *b_el_eta_old = 0;   //!
  TBranch        *b_el_isT_old = 0;   //!
  TBranch        *b_lep0_ID_old = 0;    //!
  TBranch        *b_lep1_ID_old = 0;    //!
  TBranch        *b_lep0_pt_old = 0;    //!
  TBranch        *b_lep1_pt_old = 0;    //!
  TBranch        *b_lep0_eta_old = 0;   //!
  TBranch        *b_lep1_eta_old = 0;   //!
  TBranch        *b_lep0_isT_old = 0;   //!
  TBranch        *b_lep1_isT_old = 0;   //!

  // Before cloning input TTree, tell ROOT to process all the old branches,
  //
  oldtree->SetBranchStatus("*",1);

  // Create a new file + a clone of old tree in new file
  //
  TFile *newfile = new TFile(newfilename.c_str(),"RECREATE");
  TTree *newtree = oldtree->CloneTree(0); //clone only the structure

  // Get branches from input TTree
  //
  oldtree->SetBranchAddress(old_eventNumber_name.c_str(), &eventNumber_old, &b_eventNumber);
  oldtree->SetBranchAddress(old_nel_name.c_str(), &nel_old, &b_nel_old);
  
  if ( useGroupNTup.empty() ) {
    oldtree->SetBranchAddress(old_el_pt_name.c_str(),  &el_pt_old,  &b_el_pt_old);
    oldtree->SetBranchAddress(old_el_eta_name.c_str(), &el_eta_old, &b_el_eta_old);
    oldtree->SetBranchAddress(old_el_isT_name.c_str(), &el_isT_old, &b_el_isT_old);
  } else {
    oldtree->SetBranchAddress(old_lep0_ID_name.c_str(), &lep0_ID_old , &b_lep0_ID_old);
    oldtree->SetBranchAddress(old_lep1_ID_name.c_str(), &lep1_ID_old , &b_lep1_ID_old);
    oldtree->SetBranchAddress(old_lep0_pt_name.c_str(), &lep0_pt_old , &b_lep0_pt_old);
    oldtree->SetBranchAddress(old_lep1_pt_name.c_str(), &lep1_pt_old , &b_lep1_pt_old);
    oldtree->SetBranchAddress(old_lep0_eta_name.c_str(), &lep0_eta_old, &b_lep0_eta_old);
    oldtree->SetBranchAddress(old_lep1_eta_name.c_str(), &lep1_eta_old, &b_lep1_eta_old);
    oldtree->SetBranchAddress(old_lep0_isT_name.c_str(), &lep0_isT_old, &b_lep0_isT_old);
    oldtree->SetBranchAddress(old_lep1_isT_name.c_str(), &lep1_isT_old, &b_lep1_isT_old);  
  }
  
  // Set the "new" branches in the output TTree
  //
  std::vector<float>  QMisIDWeight_new; 
  newtree->Branch("QMisIDWeight", &QMisIDWeight_new);

  // read QMisID rates from ROOT histograms
  //
  std::string path("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/");
  
  //std::string filename_AntiT("QMisIDRates_Data_Loose.root");
  //std::string filename_T("QMisIDRates_Data_nominal_v4.root");
  std::string filename_AntiT("QMisIDRates_Data_antiTantiT_v030.root");
  std::string filename_T("QMisIDRates_Data_Nominal2_v030.root");

  Info("modifytree()","Reading QMisID rates from ROOT file(s)..");
  getRates(path,filename_T,filename_AntiT);

  // Loop over entries in TTree
  //
  Info("modifytree()","Begin loop on input tree entries...\n");
  int count_bad(0);
  Long64_t i = 0;
  for ( ; i < nentries; i++ ) {

    // Print out every N events to see where we are
    //
    if ( i > 0 && ( static_cast<int>(i) % 20000 == 0 ) ) { Info("modifytree()","\t Processed %lld entries",i); }

    oldtree->GetEntry(i);

    if ( g_debug ) { Info("modifytree()","\t Processing entry: %lld - eventNumber: %lli \n",i, eventNumber_old); }

    QMisIDWeight_new.assign(3,1.0);

    if ( g_debug ) {
      int idx(0);
      for ( const auto& itr : QMisIDWeight_new ) {
        Info("modifytree()","\t\t Default QMisIDWeight[%i] = %f", idx, itr );
        ++idx;
      }
    }

    float elA_eta(-999.0),elB_eta(-999.0);
    float elA_pt(-1.0),elB_pt(-1.0);
    bool  elA_isT(false),elB_isT(false);

    if ( nel_old > 0 ) {
      
      if ( useGroupNTup.empty() ) {
        elA_eta = el_eta_old->at(0);
        elA_pt  = el_pt_old->at(0) / 1e3;
        elA_isT = el_isT_old->at(0);
      } else {
	if ( abs(lep0_ID_old) == 11 ) {
	  elA_eta = lep0_eta_old;
          elA_pt  = lep0_pt_old/ 1e3;
          elA_isT = lep0_isT_old;
	} else if ( abs(lep0_ID_old) == 13 && abs(lep1_ID_old) == 11 ) {
	  elA_eta = lep1_eta_old;
          elA_pt  = lep1_pt_old/ 1e3;
          elA_isT = lep1_isT_old;
	}
      }
      
      if ( g_debug ) { Info("modifytree()","\t elA - pT = %.2f, eta = %.2f, isT = %i ", elA_pt, elA_eta, elA_isT ); }
    }
    if ( nel_old > 1 ) {
      
      if ( useGroupNTup.empty() ) {
        elB_eta = el_eta_old->at(1);
        elB_pt  = el_pt_old->at(1) / 1e3;
        elB_isT = el_isT_old->at(1);
      } else {
        elB_eta = lep1_eta_old;
        elB_pt  = lep1_pt_old/ 1e3;
        elB_isT = lep1_isT_old;
      }
      
      if ( g_debug ) { Info("modifytree()","\t elB - pT = %.2f, eta = %.2f, isT = %i ", elB_pt, elB_eta, elB_isT ); }
    } 
    
    QMisIDWeightCalculator( &QMisIDWeight_new, elA_eta, elA_pt, elA_isT, elB_eta, elB_pt, elB_isT );

    if ( g_debug ) {
      int idx_new(0);
      for ( const auto& itr : QMisIDWeight_new ) {
        Info("modifytree()","\t\t NEW QMisIDWeight[%i] = %f", idx_new, itr );
        ++idx_new;
      }
    }

    if ( g_debug ) { Info("modifytree()","\n\n"); }

    newtree->Fill();
    
    QMisIDWeight_new.clear();

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
