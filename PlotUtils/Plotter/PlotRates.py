 #!/usr/bin/python

import os
import sys
from array import array

from ROOT import gROOT, gDirectory, gStyle, TH1D, TH2D, TF1, TFile, TCanvas, TColor, TLegend, TLatex, TStyle, TPalletteAxis, TPaveText

gROOT.Reset()
gROOT.LoadMacro("Plotter/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True) 

def plotRateEff(homedir, extension="png",rate_or_eff="Efficiency"):

  if not ( rate_or_eff == "Efficiency" or rate_or_eff == "Rate" ) :    
      sys.exit("Error! Pass either \'Efficiency\' or \'Rate\'")

  filename = "Rates.root";  
   home_directory(HOME_DIRECTORY);
  //string home_directory = "~/PhD/ttH_MultiLeptons/RUN2/PlotUtils/common_ntuple_melbourne/OutputPlots_ChFlipRates/";
  //string home_directory = "~/PhD/ttH_MultiLeptons/RUN2/PlotUtils/common_ntuple_melbourne/OutputPlots_HFRates/"; 
  if ( home_directory.back() != '/' ) { home_directory += "/"; } 
  
  vector<string> lepton_flavours;
  lepton_flavours.push_back("El");
  lepton_flavours.push_back("Mu");

  vector<string> variables;
  variables.push_back("Eta");
  variables.push_back("Pt");

  vector<string> rates;
  rates.push_back("Real");
  rates.push_back("Fake");
  
  std::string data_type("observed");
  
  //****************************** 
  
  string extension(EXTENSION);

  string path = home_directory + filename;
 
  TFile *f = new TFile(path.c_str());
  if ( !f->IsOpen() ) {
     cout << "Error, file " << filename << " could not be opened" << endl;		    
     exit(-1);
  }      

  cout << "File contains the following histograms: " << endl;
  gDirectory->ls() ; 
       
  // loop over variables
  //
  for ( unsigned int iVar = 0; iVar < variables.size() ; ++iVar ) {
         
   
    cout << "Variable : " << variables.at(iVar) << endl;
    cout << "-------------------------------" << endl;    
 
    // loop over flavours
    //
    for ( unsigned int iFlav = 0; iFlav < lepton_flavours.size() ; ++iFlav ) {

      cout << "\t\t Lepton flavour : " << lepton_flavours.at(iFlav) << endl;
      cout << "-------------------------------" << endl;    

      TCanvas *canvas = new TCanvas();
      canvas = canvas; // get rid of the warning "unused variable 'c' "
   
      canvas->SetFrameFillColor(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);
  
      TLegend *legend = new TLegend(0.7,0.4,0.9,0.6); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      //string header = rate_or_eff + " : ";
      //legend->SetHeader(header.c_str());
      legend->AddEntry((TObject*)0, "", ""); // add an empty line
      legend->SetBorderSize(0);  // no border
      legend->SetFillColor(0);   // Legend background should be white
      legend->SetTextSize(0.05); // Increase entry font size! 
      legend->SetTextFont(42);   // Helvetica 

      TLatex* leg_ATLAS  = new TLatex();
      TLatex* leg_lumi   = new TLatex();
      leg_ATLAS->SetTextSize(0.05);
      leg_ATLAS->SetNDC();
      leg_lumi->SetTextSize(0.05);
      leg_lumi->SetNDC();

      // loop over rate types
      //
      for( unsigned int iRate = 0; iRate < rates.size() ; ++iRate ) {

         cout << "\t " << rate_or_eff  << " : " << rates.at(iRate) << endl;
         cout << "-------------------------------" << endl;    

         TH1D *h = new TH1D(); 
	 
	 string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + rates.at(iRate) + "_" + rate_or_eff + "_" + data_type;
	 
         gDirectory->GetObject(histname.c_str(), h);
 
         if ( !h ) {
           cout << "Error, could not get histogram " << histname << endl;
           exit(-1);
         }
 
         h->SetStats(kFALSE); // delete the stats box on the top right corner
         h->SetLineWidth(2);
         h->SetMarkerStyle(kFullCircle);
         h->SetMarkerSize(1.0);
	 
	 // For efficiency hist
	 //
	 if ( rate_or_eff == "Efficiency" ) {
	   h->GetYaxis()->SetRangeUser(0.0,1.0);
	 }
	 
	 if ( variables.at(iVar) == "Eta" )	{ h->GetXaxis()->SetTitle("Probe |#eta|");   }
	 else if ( variables.at(iVar) == "Pt" ) { h->GetXaxis()->SetTitle("Probe pT [GeV]"); }
	 
	 string flavour("");
	 if ( lepton_flavours.at(iFlav).find("El") != string::npos )      { flavour = "Electrons:"; }
	 else if ( lepton_flavours.at(iFlav).find("Mu") != string::npos ) { flavour = "Muons:"; }
         legend->SetHeader(flavour.c_str());
	 
	 string y_title = rate_or_eff;
         h->GetYaxis()->SetTitle(y_title.c_str()); 
	       
	 h->GetXaxis()->SetTitleOffset(1.0); 
	 h->GetYaxis()->SetTitleOffset(1.0); 
	 
         if ( iRate == 0 ) {
           h->SetLineColor(kRed);
           h->SetMarkerColor(kRed);
	   h->Draw("E0"); // E0 options draws error bars
         } else {
           h->SetLineColor(kBlue);
           h->SetMarkerColor(kBlue);
	   h->Draw("E0,SAME"); // E0 options draws error bars	   	 
	 }
        
	 string legend_entry = rates.at(iRate);
         legend->AddEntry(h, legend_entry.c_str(), "F");   
       
       } // loop over rate types

       legend->Draw();
       leg_ATLAS->DrawLatex(0.5,0.75,"#bf{#it{ATLAS}} Work In Progress");
       leg_lumi->DrawLatex(0.5,0.65,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}");

       string outputname = lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + rate_or_eff + "_" + data_type + "."; // the final period is important
       outputname += extension;

       canvas->SaveAs( outputname.c_str() );	     
       
       delete legend;
       delete leg_ATLAS;
       delete leg_lumi;

    } // loop over flavours

  } // loop over variables

}

//********************************************************

void PlotMaker2(const char *EXTENSION){
  
  gROOT->LoadMacro("AtlasStyle.C");
  gROOT->LoadMacro("AtlasUtils.C");  
  gROOT->SetStyle("ATLAS");	
  
  string filename = "Rates.root";  
  string home_directory = "~/PhD/ttH_MultiLeptons/RUN2/PlotUtils/common_ntuple_melbourne/OutputPlots_ttbarMMClosure/";
  
  vector<string> lepton_flavours;
  lepton_flavours.push_back("El");
  lepton_flavours.push_back("Mu");

  vector<string> variables;
  variables.push_back("Eta");
  variables.push_back("Pt");

  vector<string> rates;
  rates.push_back("Real");
  rates.push_back("Fake");
  
  //****************************** 
  
  string extension(EXTENSION);

  string path = home_directory + filename;
 
  TFile *f = new TFile(path.c_str());
  if ( !f->IsOpen() ) {
     cout << "Error, file " << filename << " could not be opened" << endl;		    
     exit(-1);
  }      

  cout << "File contains the following histograms: " << endl;
  gDirectory->ls() ; 
       
  // loop over variables
  //
  for ( unsigned int iVar = 0; iVar < variables.size() ; ++iVar ) {
         
   
    cout << "Variable : " << variables.at(iVar) << endl;
    cout << "-------------------------------" << endl;    
 
    // loop over rate types
    //
    for( unsigned int iRate = 0; iRate < rates.size() ; ++iRate ) {

      cout << "\t Rate : " << rates.at(iRate) << endl;
      cout << "-------------------------------" << endl;    


      for ( unsigned int iFlav = 0; iFlav < lepton_flavours.size() ; ++iFlav ) {
      
         cout << "\t\t Lepton flavour : " << lepton_flavours.at(iFlav) << endl;
         cout << "-------------------------------" << endl;    

     	 TCanvas *canvas = new TCanvas();
     	 canvas = canvas; // get rid of the warning "unused variable 'c' "
   
     	 canvas->SetFrameFillColor(0);
     	 canvas->SetFrameFillStyle(0);
     	 canvas->SetFrameBorderMode(0);
  
     	 TLegend *legend = new TLegend(0.2,0.7,0.4,0.85);
     	 legend->SetHeader("Lepton flavour");
     	 legend->SetBorderSize(0);  // no border
     	 legend->SetFillColor(0);   // Legend background should be white
     	 legend->SetTextSize(0.04); // Increase entry font size! 
     	 legend->SetTextFont(42);   // Helvetica 

         TH1D *h = new TH1D(); 
	 
	 string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + rates.at(iRate) + "_R_expected";
	 
         gDirectory->GetObject(histname.c_str(), h);
 
         if ( !h ) {
           cout << "Error, could not get histogram " << histname << endl;
           exit(-1);
         }
 
         h->SetStats(kFALSE); // delete the stats box on the top right corner
         h->SetLineWidth(1);
         h->SetMarkerStyle(kFullCircle);

         h->SetMarkerSize(0.7);
	 
	 if ( variables.at(iVar) == "Eta" )	{ h->GetXaxis()->SetTitle("Probe |#eta|"); }
	 else if ( variables.at(iVar) == "Pt" ) { h->GetXaxis()->SetTitle("Probe pT");  }
	 if ( rates.at(iRate) == "Real" )	{ h->GetYaxis()->SetTitle("Real rate");    }
	 else if ( rates.at(iRate) == "Fake" )  { h->GetYaxis()->SetTitle("Fake rate");    }
	       
	 h->GetXaxis()->SetTitleOffset(1.0); 
	 h->GetYaxis()->SetTitleOffset(1.0); 
	 
         if ( iFlav == 0 ) {
           h->SetLineColor(kRed);
           h->SetMarkerColor(kRed);
		 
         } else {
           h->SetLineColor(kBlue);
           h->SetMarkerColor(kBlue);	 
	 }
 
         h->Draw("E0"); // E0 options draws error bars
        
	 string legend_entry = lepton_flavours.at(iFlav);
         legend->AddEntry(h, legend_entry.c_str(), "F");   

         legend->Draw();
         
         string outputname = lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_" + rates.at(iRate) + "_R_expected."; // the period is inmportant
         outputname += extension;
         
         canvas->SaveAs( outputname.c_str() );        
       
       } // loop over flavours
       
    } // loop over rate types

  } // loop over variables
  
}
