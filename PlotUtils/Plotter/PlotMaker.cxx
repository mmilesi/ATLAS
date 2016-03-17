#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include "TPaveText.h"
#include "AtlasStyle.C"

using namespace std;

// for fancy 2-dim histograms!
void set_plot_style()
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}

void PlotRateEff(const char* EXTENSION, const char* HOME_DIRECTORY, const char* Rate_OR_EFF, const char* FLAV_COMP = "", const char* DATA_TYPE = "Data" ){

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  string Rate_or_eff(Rate_OR_EFF);
  if ( !( Rate_or_eff == "Efficiency" || Rate_or_eff == "Rate" ) ) {
     cout << "Error! Pass either 'Efficiency' or 'Rate' "<< endl;
     exit(-1);
  }

  string flav_comp(FLAV_COMP);
  if ( !( flav_comp == "" || flav_comp == "MuMu" || flav_comp == "ElEl" || flav_comp == "OF" ) ) {
     cout << "Error! Flavour composition not supported' "<< endl;
     exit(-1);
  }

  string data_type(DATA_TYPE);
  if ( !( data_type == "Data" || data_type == "MC" ) ) {
     cout << "Error! Data type not supported' "<< endl;
     exit(-1);
  }

  string filename(flav_comp + "Rates.root");
  //string filename(flav_comp + "AvgRates.root");
  string home_directory(HOME_DIRECTORY);
  if ( home_directory.back() != '/' ) { home_directory += "/"; }

  vector<string> lepton_flavours;
  lepton_flavours.push_back("El");
  lepton_flavours.push_back("Mu");

  vector<string> variables;
  variables.push_back("Eta");
  variables.push_back("Pt");
  variables.push_back("NJets");

  vector<string> Rates;
  Rates.push_back("Real");
  Rates.push_back("Fake");

  if ( data_type == "Data" )     data_type = "observed";
  else if ( data_type == "MC" )  data_type = "expected";

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

      if ( flav_comp == "ElEl" && lepton_flavours.at(iFlav) == "Mu" ) continue;
      if ( flav_comp == "MuMu" && lepton_flavours.at(iFlav) == "El" ) continue;

      TCanvas *canvas = new TCanvas();
      canvas = canvas; // get rid of the warning "unused variable 'c' "

      canvas->SetFrameFillColor(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend *legend = new TLegend(0.7,0.5,0.9,0.7); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      //string header = Rate_or_eff + " : ";
      //legend->SetHeader(header.c_str());
      legend->AddEntry((TObject*)0, "", ""); // add an empty line
      legend->SetBorderSize(0);  // no border
      legend->SetFillColor(0);   // Legend background should be white
      legend->SetTextSize(0.04); // Increase entry font size!
      legend->SetTextFont(42);   // Helvetica

      TLatex* leg_ATLAS  = new TLatex();
      TLatex* leg_lumi   = new TLatex();
      leg_ATLAS->SetTextSize(0.04);
      leg_ATLAS->SetNDC();
      leg_lumi->SetTextSize(0.04);
      leg_lumi->SetNDC();

      // loop over Rate types
      //
      for( unsigned int iRate = 0; iRate < Rates.size() ; ++iRate ) {

         cout << "\t " << Rate_or_eff  << " : " << Rates.at(iRate) << endl;
         cout << "-------------------------------" << endl;

         TH1D *h = new TH1D();

	 string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + Rate_or_eff + "_" + data_type;

         gDirectory->GetObject(histname.c_str(), h);

         if ( !h ) {
           cout << "Error, could not get histogram " << histname << endl;
           exit(-1);
         }

         h->SetStats(kFALSE); // delete the stats box on the top right corner
         h->SetLineWidth(2);
         h->SetMarkerStyle(kFullCircle);
         h->SetMarkerSize(1.0);

	 // For Efficiency hist
	 //
	 if ( Rate_or_eff == "Efficiency" ) {
	   h->GetYaxis()->SetRangeUser(0.0,1.0);
	 }

	 if ( variables.at(iVar) == "Eta" )	{ h->GetXaxis()->SetTitle("Probe |#eta|");   }
	 else if ( variables.at(iVar) == "Pt" ) { h->GetXaxis()->SetTitle("Probe pT [GeV]"); }
	 else if ( variables.at(iVar) == "NJets" ) { h->GetXaxis()->SetTitle("Jet multiplicity"); }

	 string flavour("");
	 if ( lepton_flavours.at(iFlav).find("El") != string::npos )      { flavour = "Electrons"; }
	 else if ( lepton_flavours.at(iFlav).find("Mu") != string::npos ) { flavour = "Muons"; }
         legend->SetHeader(flavour.c_str());

	 string y_title = Rate_or_eff;
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

	 string legend_entry = Rates.at(iRate);
         legend->AddEntry(h, legend_entry.c_str(), "F");

      } // loop over Rate types

      legend->Draw();
      leg_ATLAS->DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}");

      string prepend = ( flav_comp.empty() ) ? "" : (flav_comp + "_");
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + Rate_or_eff + "_" + data_type + "."; // the final period is important
      outputname += extension;

      canvas->SaveAs( outputname.c_str() );

      delete legend;
      delete leg_ATLAS;
      delete leg_lumi;

    } // loop over flavours

  } // loop over variables

}

//********************************************************

void PlotRateEff_DiffSamples( vector< pair<string,string> >& SAMPLE_LIST, 
			      const string& DATA_TYPE = "Data", 
			      const string& FLAV_COMP = "Inclusive",
			      const string& Rate_OR_EFF = "Efficiency",
			      const string& EXTENSION = "png" ) 
{

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  if ( !( Rate_OR_EFF == "Efficiency" || Rate_OR_EFF == "Rate" ) ) {
     cout << "Error! Pass either 'Efficiency' or 'Rate' "<< endl;
     exit(-1);
  }
  if ( !( FLAV_COMP == "Inclusive" || FLAV_COMP == "MuMu" || FLAV_COMP == "ElEl" || FLAV_COMP == "OF" ) ) {
     cout << "Error! Flavour composition not supported' "<< endl;
     exit(-1);
  }
  if ( !( DATA_TYPE == "Data" || DATA_TYPE == "MC" ) ) {
     cout << "Error! Data type not supported' "<< endl;
     exit(-1);
  }

  vector< pair< TFile*,string> > input_files;
  
  for ( auto& samp : SAMPLE_LIST ) {
    
    if ( samp.first.back() != '/' ) { samp.first += "/"; }
    
    string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : FLAV_COMP;
    string path = samp.first + prepend + "Rates.root"; // "AvgRates.root"

    TFile *f = TFile::Open(path.c_str());
    if ( !f->IsOpen() ) {
       cout << "Error, file " << path << " could not be opened" << endl;
       exit(-1);
    }
    pair<TFile*,string> this_pair = make_pair(f, samp.second);
    input_files.push_back(this_pair);
     
  }

  vector<string> lepton_flavours;
  lepton_flavours.push_back("El");
  lepton_flavours.push_back("Mu");

  vector<string> variables;
  variables.push_back("Eta");
  variables.push_back("Pt");
  //variables.push_back("NJets");

  vector<string> Rates;
  Rates.push_back("Real");
  Rates.push_back("Fake");

  string data_type("");
  if ( DATA_TYPE == "Data" )     { data_type = "observed"; }
  else if ( DATA_TYPE == "MC" )  { data_type = "expected"; }
   
  // loop over variables
  //
  for ( unsigned int iVar(0); iVar < variables.size() ; ++iVar ) {

    cout << "Variable : " << variables.at(iVar) << endl;
    cout << "-------------------------------" << endl;

    // loop over flavours
    //
    for ( unsigned int iFlav(0); iFlav < lepton_flavours.size() ; ++iFlav ) {

      cout << "\tLepton flavour : " << lepton_flavours.at(iFlav) << endl;
      cout << "-------------------------------" << endl;

      if ( FLAV_COMP == "ElEl" && lepton_flavours.at(iFlav) == "Mu" ) { continue; }
      if ( FLAV_COMP == "MuMu" && lepton_flavours.at(iFlav) == "El" ) { continue; }

      string flavour("");
      if ( lepton_flavours.at(iFlav).find("El") != string::npos )      { flavour = "Electrons"; }
      else if ( lepton_flavours.at(iFlav).find("Mu") != string::npos ) { flavour = "Muons"; }

      TCanvas *canvas = new TCanvas();
      canvas = canvas; // get rid of the warning "unused variable 'c' "

      canvas->SetFrameFillColor(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend *legend = new TLegend(0.68,0.6,0.925,0.8); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      legend->AddEntry((TObject*)0, "", ""); // add an empty line
      legend->SetBorderSize(0);  // no border
      legend->SetFillColor(0);   // Legend background should be white
      legend->SetTextSize(0.035); // Increase entry font size!
      legend->SetTextFont(42);   // Helvetica
      
      TLatex* leg_ATLAS  = new TLatex();
      TLatex* leg_lumi   = new TLatex();
      leg_ATLAS->SetTextSize(0.04);
      leg_ATLAS->SetNDC();
      leg_lumi->SetTextSize(0.04);
      leg_lumi->SetNDC();

      // loop over files
      //
      for ( unsigned int iFile(0); iFile < input_files.size(); ++iFile  ) {
        
	cout << "\t\tFile : " << (input_files.at(iFile).first)->GetName()  << endl;
        cout << "-------------------------------" << endl;
	
        // loop over Rate types
        //
        for( unsigned int iRate(0); iRate < Rates.size() ; ++iRate ) {

	  cout << "\t\t\t" << Rates.at(iRate) << " Rate/Efficiency " << endl;
          cout << "-------------------------------" << endl;

          TH1D *h(nullptr);

	  string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + Rate_OR_EFF  + "_" + data_type;

          (input_files.at(iFile).first)->GetObject(histname.c_str(), h);

          if ( !h ) {
            cout << "Error, could not get histogram " << histname << endl;
            exit(-1);
          }	
	  
          h->SetStats(kFALSE); // delete the stats box on the top right corner
          h->SetLineWidth(2);
          h->SetMarkerSize(1.0);

	  // For Efficiency hist
	  //
	  if ( Rate_OR_EFF == "Efficiency" ) {
	    h->GetYaxis()->SetRangeUser(0.0,1.0);
	  }
          
	  string title("");
	  if ( variables.at(iVar) == "Eta" )	    { title = "Probe |#eta| - " + flavour; }
	  else if ( variables.at(iVar) == "Pt" )    { title = "Probe pT [GeV] - " + flavour; }
	  else if ( variables.at(iVar) == "NJets" ) { title = "Jet multiplicity - " + flavour; }
          h->GetXaxis()->SetTitle(title.c_str());

          h->GetYaxis()->SetTitle(Rate_OR_EFF.c_str());
	  h->GetXaxis()->SetTitleOffset(1.0);
	  h->GetYaxis()->SetTitleOffset(1.0);

          switch (iFile) {
	    case 0:
	      h->SetLineStyle(1);
              h->SetMarkerStyle(kFullCircle);
	     break;
	    case 1:
	      h->SetLineStyle(3);
              h->SetMarkerStyle(kCircle);
	      break;
	    case 2:
	      h->SetLineStyle(6);
              h->SetMarkerStyle(kOpenTriangleUp);
	      break;
	    default:
	      h->SetLineStyle(1);
              h->SetMarkerStyle(kDot);
	      break;
	  }

          switch (iRate) {
	    case 0:
              h->SetLineColor(kRed);
              h->SetMarkerColor(kRed);
	     break;
	    default:
              h->SetLineColor(kBlue);
              h->SetMarkerColor(kBlue);
	      break;
	  }

          if ( iRate == 0 && iFile == 0 ) { h->Draw("E0"); } // E0 options draws error bars 
	  else                            { h->Draw("E0,SAME");}

	  string legend_entry = Rates.at(iRate) + " - " + input_files.at(iFile).second;
          legend->AddEntry(h, legend_entry.c_str(), "P");
          legend->AddEntry((TObject*)0, "", ""); 
	
	} // close loop over Rates
       
      } // close loop over files
       
      legend->Draw();
       
      leg_ATLAS->DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}");

      string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : ( FLAV_COMP + "_" );
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + Rate_OR_EFF + "_" + DATA_TYPE + "." + EXTENSION; 

      canvas->SaveAs( outputname.c_str() );

      delete legend;
      delete leg_ATLAS;
      delete leg_lumi;

    } // close loop on flavours
     
  } // close loop on variables

}

// ********************************************************

void PlotRateEff_DataVSMC( vector< pair<string,string> >& SAMPLE_LIST, 
			   const string& FLAV_COMP = "Inclusive",
			   const string& Rate_OR_EFF = "Efficiency",
			   const string& EXTENSION = "png" ) 
{

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  if ( !( Rate_OR_EFF == "Efficiency" || Rate_OR_EFF == "Rate" ) ) {
     cout << "Error! Pass either 'Efficiency' or 'Rate' "<< endl;
     exit(-1);
  }
  if ( !( FLAV_COMP == "Inclusive" || FLAV_COMP == "MuMu" || FLAV_COMP == "ElEl" || FLAV_COMP == "OF" ) ) {
     cout << "Error! Flavour composition not supported' "<< endl;
     exit(-1);
  }
 
  vector< pair< TFile*,string> > input_files;
  
  for ( auto& samp : SAMPLE_LIST ) {
  
    if ( samp.first.back() != '/' ) { samp.first += "/"; }
    
    string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : FLAV_COMP;
    string path = samp.first + prepend + "Rates.root"; // "AvgRates.root"    

    TFile *f = TFile::Open(path.c_str());
    if ( !f->IsOpen() ) {
       cout << "Error, file " << path << " could not be opened" << endl;
       exit(-1);
    }
    pair<TFile*,string> this_pair = make_pair(f, samp.second);
    input_files.push_back(this_pair);
     
  }

  vector<string> lepton_flavours;
  lepton_flavours.push_back("El");
  lepton_flavours.push_back("Mu");

  vector<string> variables;
  variables.push_back("Eta");
  variables.push_back("Pt");
  //variables.push_back("NJets");

  vector<string> Rates;
  Rates.push_back("Real");
  Rates.push_back("Fake");
   
  // loop over variables
  //
  for ( unsigned int iVar(0); iVar < variables.size() ; ++iVar ) {

    cout << "Variable : " << variables.at(iVar) << endl;
    cout << "-------------------------------" << endl;

    // loop over flavours
    //
    for ( unsigned int iFlav(0); iFlav < lepton_flavours.size() ; ++iFlav ) {

      cout << "\tLepton flavour : " << lepton_flavours.at(iFlav) << endl;
      cout << "-------------------------------" << endl;

      if ( FLAV_COMP == "ElEl" && lepton_flavours.at(iFlav) == "Mu" ) { continue; }
      if ( FLAV_COMP == "MuMu" && lepton_flavours.at(iFlav) == "El" ) { continue; }

      string flavour("");
      if ( lepton_flavours.at(iFlav).find("El") != string::npos )      { flavour = "Electrons"; }
      else if ( lepton_flavours.at(iFlav).find("Mu") != string::npos ) { flavour = "Muons"; }

      TCanvas *canvas = new TCanvas();
      canvas = canvas; // get rid of the warning "unused variable 'c' "

      canvas->SetFrameFillColor(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend *legend = new TLegend(0.6,0.6,0.89,0.8); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      legend->AddEntry((TObject*)0, "", ""); // add an empty line
      legend->SetBorderSize(0);  // no border
      legend->SetFillColor(0);   // Legend background should be white
      legend->SetTextSize(0.035); // Increase entry font size!
      legend->SetTextFont(42);   // Helvetica
 
      TLatex* leg_ATLAS  = new TLatex();
      TLatex* leg_lumi   = new TLatex();
      leg_ATLAS->SetTextSize(0.04);
      leg_ATLAS->SetNDC();
      leg_lumi->SetTextSize(0.04);
      leg_lumi->SetNDC();

      // loop over files
      //
      for ( unsigned int iFile(0); iFile < input_files.size(); ++iFile  ) {
        
	string filename((input_files.at(iFile).first)->GetName());
	cout << "\t\tFile : " << filename  << endl;
        cout << "-------------------------------" << endl;

       // loop over Rate types
        //
        for( unsigned int iRate(0); iRate < Rates.size() ; ++iRate ) {

	  cout << "\t\t\t" << Rates.at(iRate) << " Rate/Efficiency " << endl;
          cout << "-------------------------------" << endl;

          TH1D *h(nullptr);

          string data_type = ( input_files.at(iFile).second.find("Data") != string::npos ) ? "observed" : "expected";
	  
	  string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + Rate_OR_EFF  + "_" + data_type;

          (input_files.at(iFile).first)->GetObject(histname.c_str(), h);

          if ( !h ) {
            cout << "Error, could not get histogram " << histname << endl;
            exit(-1);
          }	
	  
          h->SetStats(kFALSE); // delete the stats box on the top right corner
          h->SetLineWidth(2);
          h->SetMarkerSize(1.0);

	  // For Efficiency hist
	  //
	  if ( Rate_OR_EFF == "Efficiency" ) {
	    h->GetYaxis()->SetRangeUser(0.0,1.0);
	  }

	  string title("");
	  if ( variables.at(iVar) == "Eta" )	    { title = "Probe |#eta| - " + flavour; }
	  else if ( variables.at(iVar) == "Pt" )    { title = "Probe pT [GeV] - " + flavour; }
	  else if ( variables.at(iVar) == "NJets" ) { title = "Jet multiplicity - " + flavour; }
          h->GetXaxis()->SetTitle(title.c_str());

          h->GetYaxis()->SetTitle(Rate_OR_EFF.c_str());
	  h->GetXaxis()->SetTitleOffset(1.0);
	  h->GetYaxis()->SetTitleOffset(1.0);

          switch (iFile) {
	    case 0:
	      h->SetLineStyle(1);
              h->SetMarkerStyle(kFullCircle);
	     break;
	    case 1:
	      h->SetLineStyle(3);
              h->SetMarkerStyle(kCircle);
	      break;
	    case 2:
	      h->SetLineStyle(6);
              h->SetMarkerStyle(kOpenTriangleUp);
	      break;
	    default:
	      h->SetLineStyle(1);
              h->SetMarkerStyle(kDot);
	      break;
	  }

          switch (iRate) {
	    case 0:
              h->SetLineColor(kBlack);
              h->SetMarkerColor(kBlack);
	     break;
	    default:
              h->SetLineColor(kOrange+7);
              h->SetMarkerColor(kOrange+7);
	      break;
	  }

          if ( iRate == 0 && iFile == 0 ) { h->Draw("E0"); } // E0 options draws error bars 
	  else                            { h->Draw("E0,SAME");}

	  string legend_entry = Rates.at(iRate) + " - " + input_files.at(iFile).second;
          legend->AddEntry(h, legend_entry.c_str(), "P");
          legend->AddEntry((TObject*)0, "", ""); 
	
	} // close loop over Rates
       
      } // close loop over files
       
      legend->Draw();
      leg_ATLAS->DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}");

      string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : ( FLAV_COMP + "_" );
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + Rate_OR_EFF + "_DataVSMC." + EXTENSION; 

      canvas->SaveAs( outputname.c_str() );

      delete legend;
      delete leg_ATLAS;
      delete leg_lumi;

    } // close loop on flavours
     
  } // close loop on variables

}

void execute_DiffSamples() {
  
  vector<pair<string,string> > vec;
  vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV/","Baseline"));
  vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIso_MCQMisID_Mllgt40GeV/","No Isolation"));
  vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIP_MCQMisID_Mllgt40GeV/","Relaxed IP"));
  
  PlotRateEff_DiffSamples(vec);

}

void execute_DataVSMC() {
  
  vector<pair<string,string> > vec;
  
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV/","Baseline - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_Baseline_Mllgt40GeV/","Baseline - MC t#bar{t}")); 
   
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIso_MCQMisID_Mllgt40GeV/","No Isolation - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIso_Mllgt40GeV/","No Isolation - MC t#bar{t}")); 
 
  vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIP_MCQMisID_Mllgt40GeV/","Relaxed IP - Data"));
  vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIP_Mllgt40GeV/","Relaxed IP - MC t#bar{t}")); 
  
  PlotRateEff_DataVSMC(vec);

}

