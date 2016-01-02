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


/* ******************
/ 
/ The 'main' function
/
****************** */

void debug_electronSF(std::string filename = "input.root", std::string  NENTRIES = "ALL", std::string treename= "physics")
{
  // This script loads a tree, clones it, removes a branch and substitutes it with another. 
  // The branch can also have the same name and in this way you can change for example the type of the variable or the content.

  Info("debug_electronSF()","Starting off...");

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
  std::string old_nel_name("nel");
  std::string old_el_pt_name("el_pt");
  std::string old_el_eta_name("el_eta"); 
  std::string old_LHTight_name("el_LHTight"); 
  std::string old_el_PIDEff_SF_LHTight_name("el_PIDEff_SF_LHTight");   
  std::string old_el_RecoEff_SF_name("el_RecoEff_SF");   
  
  Int_t                  eventNumber_old; eventNumber_old = -1;
  Int_t                  nel_old;        nel_old = -1;
  std::vector<float>*    el_pt_old;      el_pt_old = 0;
  std::vector<float>*    el_eta_old;     el_eta_old = 0;  
  std::vector<int>*      el_LHTight_old;     el_LHTight_old = 0;
  std::vector<vector<double> >* el_RecoEff_SF_old;  el_RecoEff_SF_old = 0;
  std::vector<vector<double> >* el_PIDEff_SF_LHTight_old;  el_PIDEff_SF_LHTight_old = 0;


  // List of old branches
  //
  TBranch        *b_eventNumber_old = 0;    //!
  TBranch	 *b_nel_old    = 0;         //!
  TBranch	 *b_el_pt_old = 0;          //!
  TBranch	 *b_el_eta_old = 0;         //!   
  TBranch        *b_el_LHTight_old = 0;     //!
  TBranch        *b_el_RecoEff_SF_old = 0;          //!
  TBranch	 *b_el_PIDEff_SF_LHTight_old = 0;   //!
  
  // Before cloning input TTree, tell ROOT to process all the old branches, 
  // except for the one you want to change
  //
  oldtree->SetBranchStatus("*",1);

  oldtree->SetBranchAddress(old_eventNumber_name.c_str(), &eventNumber_old, &b_eventNumber_old);
  oldtree->SetBranchAddress(old_nel_name.c_str(), &nel_old, &b_nel_old);
  oldtree->SetBranchAddress(old_el_pt_name.c_str(), &el_pt_old, &b_el_pt_old);
  oldtree->SetBranchAddress(old_el_eta_name.c_str(), &el_eta_old, &b_el_eta_old);
  oldtree->SetBranchAddress(old_LHTight_name.c_str(), &el_LHTight_old, &b_el_LHTight_old);
  oldtree->SetBranchAddress(old_el_RecoEff_SF_name.c_str(), &el_RecoEff_SF_old, &b_el_RecoEff_SF_old);
  oldtree->SetBranchAddress(old_el_PIDEff_SF_LHTight_name.c_str(), &el_PIDEff_SF_LHTight_old, &b_el_PIDEff_SF_LHTight_old);

  // Loop over entries in TTree
  //
  Info("debug_electronSF()","Begin loop on input tree entries...\n");
  Long64_t i = 0;
  for ( ; i < nentries; i++ ) {
  
    // Print out every N events to see where we are
    //
    if ( i > 0 && ( static_cast<int>(i) % 20000 == 0 ) ) { Info("debug_electronSF()","\t Processed %lld entries",i); }
  
    oldtree->GetEntry(i);
    
    if ( nel_old != 2 ) { continue; }
    
    if ( !(eventNumber_old == 664 || eventNumber_old == 1646 || eventNumber_old == 1650) ) { continue; }
    
    Info("debug_electronSF()","\t Processing entry: %lld - eventNumber: %i \n",i, eventNumber_old);
    
    if ( el_LHTight_old->at(0) == 1 &&  el_LHTight_old->at(1) == 1) {
    
      Info("debug_electronSF()","\t\t electron[0] (LHTight) - pT = %f \n", el_pt_old->at(0)/1e3 );
      
      for ( unsigned int sys = 0; sys < (el_PIDEff_SF_LHTight_old->at(0)).size(); sys++ ) {
         Info("debug_electronSF()","\t\t\t SF_LHTight[%i] = %f \n",sys, (el_PIDEff_SF_LHTight_old->at(0)).at(sys) );
      }
      
      Info("debug_electronSF()","\t\t electron[1] (LHTight) - pT = %f \n", el_pt_old->at(1)/1e3 );
      
      for ( unsigned int sys = 0; sys < (el_PIDEff_SF_LHTight_old->at(1)).size(); sys++ ) {
         Info("debug_electronSF()","\t\t\t SF_LHTight[%i] = %f \n",sys, (el_PIDEff_SF_LHTight_old->at(1)).at(sys) );
      }      
      
    }
  
  }
  
  Info("debug_electronSF()","End of loop!\n ---> total number of processed events: %lld ", i);

  delete oldfile;
  
}
