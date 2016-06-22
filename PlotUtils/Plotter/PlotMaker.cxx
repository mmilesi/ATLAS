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

void PlotRateEff( pair<string,string>& SAMPLE,
		  const string& DATA_TYPE = "Data",
		  const string& FLAV_COMP = "Inclusive",
		  const string& RATE_OR_EFF = "Efficiency",
		  const string& EXTENSION = "png" )
{

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  if ( !( RATE_OR_EFF == "Efficiency" || RATE_OR_EFF == "Rate" ) ) {
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

  if ( SAMPLE.first.back() != '/' ) { SAMPLE.first += "/"; }

  string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : FLAV_COMP;
  string path = SAMPLE.first + prepend + "Rates.root"; // "AvgRates.root"

  TFile *f = TFile::Open(path.c_str());
  if ( !f->IsOpen() ) {
     cout << "Error, file " << path << " could not be opened" << endl;
     exit(-1);
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
      legend->SetBorderSize(0);   // no border
      //legend->SetFillColor(0);   // Legend background should be white
      legend->SetFillStyle(0);    // Legend transparent background
      legend->SetTextSize(0.035); // Increase entry font size!
      legend->SetTextFont(42);    // Helvetica

      TLatex* leg_ATLAS  = new TLatex();
      TLatex* leg_lumi   = new TLatex();
      leg_ATLAS->SetTextSize(0.04);
      leg_ATLAS->SetNDC();
      leg_lumi->SetTextSize(0.04);
      leg_lumi->SetNDC();

      // loop over Rate types
      //
      for( unsigned int iRate(0); iRate < Rates.size() ; ++iRate ) {

        cout << "\t\t\t" << Rates.at(iRate) << " Rate/Efficiency " << endl;
        cout << "-------------------------------" << endl;

        TH1D *h(nullptr);

        string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + RATE_OR_EFF  + "_" + data_type;

        f->GetObject(histname.c_str(), h);

        if ( !h ) {
          cout << "Error, could not get histogram " << histname << endl;
          exit(-1);
        }

        h->SetStats(kFALSE); // delete the stats box on the top right corner
        h->SetLineWidth(2);
        h->SetMarkerSize(1.0);

        // For Efficiency hist
        //
        if ( RATE_OR_EFF == "Efficiency" ) {
          h->GetYaxis()->SetRangeUser(0.0,1.0);
        }

        string title("");
        if ( variables.at(iVar) == "Eta" )	  { title = "Probe |#eta| - " + flavour; }
        else if ( variables.at(iVar) == "Pt" )    { title = "Probe pT [GeV] - " + flavour; }
        else if ( variables.at(iVar) == "NJets" ) { title = "Jet multiplicity - " + flavour; }
        h->GetXaxis()->SetTitle(title.c_str());

        h->GetYaxis()->SetTitle(RATE_OR_EFF.c_str());
        h->GetXaxis()->SetTitleOffset(1.0);
        h->GetYaxis()->SetTitleOffset(1.0);

        h->SetLineStyle(1);
        h->SetMarkerStyle(kFullCircle);

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

        if ( iRate == 0 ) { h->Draw("E0"); } // E0 options draws error bars
        else		  { h->Draw("E0,SAME");}

        string legend_entry = Rates.at(iRate) + " - " + SAMPLE.second;
        legend->AddEntry(h, legend_entry.c_str(), "P");
        legend->AddEntry((TObject*)0, "", "");

      } // close loop over Rates

      legend->Draw();

      leg_ATLAS->DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.8 fb^{-1}");

      string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : ( FLAV_COMP + "_" );
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + RATE_OR_EFF + "_" + DATA_TYPE + "." + EXTENSION;

      canvas->SaveAs( outputname.c_str() );

      delete legend;
      delete leg_ATLAS;
      delete leg_lumi;

    } // close loop on flavours

  } // close loop on variables

}

//********************************************************

void PlotRateEff_DiffSamples( vector< pair<string,string> >& SAMPLE_LIST,
			      const string& DATA_TYPE = "Data",
			      const string& FLAV_COMP = "Inclusive",
			      const string& RATE_OR_EFF = "Efficiency",
			      const string& EXTENSION = "png" )
{

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  if ( !( RATE_OR_EFF == "Efficiency" || RATE_OR_EFF == "Rate" ) ) {
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

    //string path = samp.first;

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
  //variables.push_back("Eta");
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

      //TLegend *legend = new TLegend(0.68,0.6,0.925,0.8); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      TLegend *legend = new TLegend(0.45,0.5,0.925,0.8); // (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
      legend->AddEntry((TObject*)0, "", ""); // add an empty line
      legend->SetBorderSize(0);   // no border
      //legend->SetFillColor(0);  // Legend background should be white
      legend->SetFillStyle(0);    // Legend transparent background
      legend->SetTextSize(0.035); // Increase entry font size!
      legend->SetTextFont(42);    // Helvetica

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

	  string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + RATE_OR_EFF  + "_" + data_type;

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
	  if ( RATE_OR_EFF == "Efficiency" ) {
	    h->GetYaxis()->SetRangeUser(0.0,1.0);
	    //h->GetYaxis()->SetRangeUser(0.4,1.0);
	  }

	  string title("");
	  if ( variables.at(iVar) == "Eta" )	    { title = "Probe |#eta| - " + flavour; }
	  else if ( variables.at(iVar) == "Pt" )    { title = "Probe pT [GeV] - " + flavour; }
	  else if ( variables.at(iVar) == "NJets" ) { title = "Jet multiplicity - " + flavour; }
          h->GetXaxis()->SetTitle(title.c_str());

          h->GetYaxis()->SetTitle(RATE_OR_EFF.c_str());
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
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.8 fb^{-1}");

      string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : ( FLAV_COMP + "_" );
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + RATE_OR_EFF + "_" + DATA_TYPE + "." + EXTENSION;

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
			   const string& RATE_OR_EFF = "Efficiency",
			   const string& EXTENSION = "png" )
{

  // Use ATLAS style for plotting
  //
  SetAtlasStyle();

  if ( !( RATE_OR_EFF == "Efficiency" || RATE_OR_EFF == "Rate" ) ) {
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
      legend->SetBorderSize(0);   // no border
      //legend->SetFillColor(0);   // Legend background should be white
      legend->SetFillStyle(0);    // Legend transparent background
      legend->SetTextSize(0.035); // Increase entry font size!
      legend->SetTextFont(42);    // Helvetica

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

	  string histname = lepton_flavours.at(iFlav) + "_Probe" + variables.at(iVar) + "_" + Rates.at(iRate) + "_" + RATE_OR_EFF  + "_" + data_type;

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
	  if ( RATE_OR_EFF == "Efficiency" ) {
	    h->GetYaxis()->SetRangeUser(0.0,1.0);
	    //h->GetYaxis()->SetRangeUser(0.5,1.0);
	  }

	  string title("");
	  if ( variables.at(iVar) == "Eta" )	    { title = "Probe |#eta| - " + flavour; }
	  else if ( variables.at(iVar) == "Pt" )    { title = "Probe pT [GeV] - " + flavour; }
	  else if ( variables.at(iVar) == "NJets" ) { title = "Jet multiplicity - " + flavour; }
          h->GetXaxis()->SetTitle(title.c_str());

          h->GetYaxis()->SetTitle(RATE_OR_EFF.c_str());
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
	      // TEMP:
	      //h->SetLineStyle(1);
              //h->SetMarkerStyle(kFullCircle);
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
      leg_lumi->DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.8 fb^{-1}");

      string prepend = ( FLAV_COMP == "Inclusive" ) ? "" : ( FLAV_COMP + "_" );
      string outputname = prepend + lepton_flavours.at(iFlav) + "Probe" + variables.at(iVar) + "_RealFake" + "_" + RATE_OR_EFF + "_DataVSMC." + EXTENSION;

      canvas->SaveAs( outputname.c_str() );

      delete legend;
      delete leg_ATLAS;
      delete leg_lumi;

    } // close loop on flavours

  } // close loop on variables

}

void execute() {

  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV/","Data");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_v029_Baseline_Mllgt40GeV/","Baseline - MC t#bar{t}");

  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_25ns_v7_DDQMisID_NoZminCut/","Data");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_25ns_v7_NoZminCut/","MC t#bar{t}");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_25ns_v7_DDQMisID_ZminCut/","Data");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_25ns_v7_ZminCut/","MC t#bar{t}");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_v030_DDQMisID_NoZminCut/","Data");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_v030_NoZminCut/","MC t#bar{t}");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_25ns_v7_MCQMisID_NoZminCut_DataOnly/","Data - No sub.");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_25ns_v7_NoTruth/","MC t#bar{t} - No sub.");

  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/","Data - bkg sub.");
  //pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/","MC t#bar{t}");

  //pair<string,string> my_pair = make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake/","Data - bkg sub.");
  pair<string,string> my_pair = make_pair("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake/","MC t#bar{t}");

  //PlotRateEff(my_pair,"Data","Inclusive","Efficiency","png");
  PlotRateEff(my_pair,"MC","Inclusive","Efficiency","png");

}

void execute_DiffSamples() {

  vector<pair<string,string> > vec;

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/","Baseline"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIso_MCQMisID_Mllgt40GeV_AllElEtaCut/","No Isolation"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIP_MCQMisID_Mllgt40GeV_AllElEtaCut/","Relaxed IP"));

  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_Baseline_Mllgt40GeV_AllElEtaCut/","t#bar{t} -Baseline"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIso_Mllgt40GeV_AllElEtaCut/","t#bar{t} - No Isolation"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIP_Mllgt40GeV_AllElEtaCut/","t#bar{t} - Relaxed IP"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_ElTagEtaCut/","|#eta_{e}| < 1.37 on tag only"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/","|#eta_{e}| < 1.37 on all"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/Rates.root","OF+SF"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/MuMuRates.root","SF"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/ElElRates.root","SF"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/","MC QMisID"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_DDQMisID_Mllgt40GeV_AllElEtaCut/","DD QMisID"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_DDQMisID_NoZminCut/","Data - Group NTup v7"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v030_DDQMisID_NoZminCut/","Data - Melb NTup v30"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_NoZminCut/","MC t#bar{t} - Group NTup v7"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v030_NoZminCut/","MC t#bar{t} - Melb NTup v30"));


  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePrompt/","Data - OF+SF - No sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID/","Data - OF+SF - Prompt & QMisID Sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePrompt/","Data - OF - No sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID/","Data - OF - Prompt & QMisID Sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePrompt/","Data - SF - No sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID/","Data - SF - Prompt & QMisID Sub."));

  //vec.push_back(make_pair("../PlotVault/PLOTS_v030_25ns_v7/OutputPlots_MMRates_v030_DDQMisID_NoZminCut/","Data - OLR on loose"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v031_DDQmisID/","Data - OLR on tight"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePrompt/","Data - OF+SF - No sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID_BabarSub/Rates_NoSub/","Data - Babar sel. - No sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID/","Data - OF+SF - Prompt & QMisID Sub."));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FakeCRProbePromptORQMisID_BabarSub/Rates_Sub/","Data - Babar sel. - Prompt & QMisID Sub."));

  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_RealCRAllPromptChFlipVeto_OFplusSF_TTBar/","t#bar{t} - OF + SF"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_RealCRAllPromptChFlipVeto_OFplusSF_Zjets/","Z/#gamma* + jets - OF + SF"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_RealCRAllPromptChFlipVeto_OFplusSF_AllMC/","All MC - OF + SF"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_RealCRAllPromptChFlipVeto_ZjetsCR_Zjets/","Z/#gamma* + jets - Z CR"));

  vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection/Rates_NoSub/","Data - No sub."));
  vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection/Rates_ProbeSub/","Data - (!Prompt & QMisID) Sub."));

  PlotRateEff_DiffSamples(vec,"Data","Inclusive","Efficiency","png");
  //PlotRateEff_DiffSamples(vec,"Data","OF","Efficiency","png");
  //PlotRateEff_DiffSamples(vec,"Data","MuMu","Efficiency","png");
  //PlotRateEff_DiffSamples(vec, "MC","Inclusive","Efficiency","png");

}

void execute_DataVSMC() {

  vector<pair<string,string> > vec;

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_Baseline_MCQMisID_Mllgt40GeV_AllElEtaCut/","Baseline - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_Baseline_Mllgt40GeV_AllElEtaCut/","Baseline - MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIso_MCQMisID_Mllgt40GeV_AllElEtaCut/","No Isolation - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIso_Mllgt40GeV_AllElEtaCut/","No Isolation - MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v029_NoLepIP_MCQMisID_Mllgt40GeV_AllElEtaCut/","Relaxed IP - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v029_NoLepIP_Mllgt40GeV_AllElEtaCut/","Relaxed IP - MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_DDQMisID_NoZminCut/","Moriond sel. - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_NoZminCut/","Moriond sel. - MC t#bar{t}"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_DDQMisID_ZminCut/","Moriond sel. - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_ZminCut/","Moriond sel. - MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_v030_DDQMisID_NoZminCut/","Moriond sel. - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v030_NoZminCut/","Moriond sel. - MC t#bar{t} - TRUTH"));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_v030_NoZminCut_NoTruth/","Moriond sel. - MC t#bar{t}"));

  //vec.push_back(make_pair("../PLOTS_v030_25ns_v7/OutputPlots_MMRates_v030_DDQMisID_NoZminCut_TopEnrichedRealCR/","Moriond sel. - Data"));
  //vec.push_back(make_pair("../OutputPlots_MMRates_v030_DDQMisID_NoZminCut_TopEnrichedRealCR_RealCRProbeNonPromptOrChFlip/","Moriond sel. - Data (sub. !prompt, QMisID)"));
  //vec.push_back(make_pair("../PLOTS_v030_25ns_v7/OutputPlots_MMClosureRates_v030_NoZminCut_TopEnrichedRealCR/","Moriond sel. - MC t#bar{t} - prompt probe"));
  //vec.push_back(make_pair("../PLOTS_v030_25ns_v7/OutputPlots_MMClosureRates_v030_NoZminCut_TopEnrichedRealCR_NoTruth/","Moriond sel. - MC t#bar{t}"));

  //vec.push_back(make_pair("../PlotVault/PLOTS_v030_25ns_v7/OutputPlots_MMRates_v030_DDQMisID_NoZminCut/","Moriond sel. - Data - OF+SF"));
  //vec.push_back(make_pair("../PlotVault/PLOTS_v030_25ns_v7/OutputPlots_MMRates_v030_DDQMisID_NoZminCut_TopEnrichedRealCR/","Moriond sel. - Data - OF"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/","Data - bkg sub."));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/","MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake/","Data - bkg sub."));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake/","MC t#bar{t}"));

  vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v7_FinalSelection_DDQmisID/","Data - w/ sub."));
  vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v7_FinalSelection_DDQmisID/","MC t#bar{t}"));

  //vec.push_back(make_pair("../OutputPlots_MMRates_25ns_v14/","Data - w/ sub."));
  //vec.push_back(make_pair("../OutputPlots_MMClosureRates_25ns_v14/","MC t#bar{t}"));

  PlotRateEff_DataVSMC(vec,"Inclusive","Efficiency","png");

}

