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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "TError.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "AtlasStyle.h"
#include "TLine.h"
#include <Riostream.h>

TH1::SetDefaultSumw2() 

bool g_debug(true);

std::vector<TH1D*> g_histograms_pt;  

std::vector<std::string> g_charge = {"OS","SS"};
std::vector<std::string> g_flavour = {"ElEl","MuMu"};
std::vector<std::string> g_obs_selection = {"TT","TL","LL"};
std::vector<std::string> g_true_selection = {"RR","RF","FF"};

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
  
    path = ch + "_";
  
    for ( const auto& fl : g_flavour ) {
     
      path = path + fl + "_";

      for ( const auto& sl : g_obs_selection ) {
        path = path + sl;
	path = path + "/" + path + "LepPt.root";
	
	if ( g_debug ) { Info("getHists()","Reading histogram from:\t %s", path.c_str() ); }
	
	TFile *file = TFile::Open(path.c_str());
    	if ( !file->IsOpen() ) {
    	  SysError("getHists()", "Failed to open ROOT file from path: %s . Aborting", path.c_str() );
    	  exit(-1);
    	}
  
        TH1D *hist = get_object<TH1D>( *file, "observed" );
	
	if ( nbins == -1 ) { nbins = hist->GetNBinsX()+2; }
	
	// Do prompt and g_charge flip subtraction in SS
	//
        if  ( ch == "SS" ) {
          TH1D *hist_to_sub = get_object<TH1D>( *file, "expected" );
	  hist->Add(hist_to_sub, -1.0);
          for ( int i(0);  i < hist->GetNbinsX()+2; ++i ) {
            if ( hist->GetBinContent(i) < 0 ) { hist->SetBinContent(i,0.0); }	  
	  }
	}
	
	std::string new_name = ch + "_" + fl + "_" sl;
	
	if ( g_debug ) { Info("getHists()","Saving histogram w/ name:\t %s", new_name.c_str() ); }
	
	hist->SetName(new_name.c_str());
	hist->SetDirectory(0);	
	
	g_histograms_pt.push_back(hist);
      }
    
    }
  
  }
 
}

// -------------------------------------------------
// Define expected events with the matrix equation  
// -------------------------------------------------

double getExpected( const std::string& obs_selection, const double& r,  const double& f, const double& RR, const double& RF, const double& FF ) {

  double exp(-1.0);
  
  if ( obs_selection == "TT" ) {
  
    exp = r * r * RR + 2.0 * r * f * RF + f * f * FF; 
  
  } else if ( obs_selection == "TL" ) {
  
    exp = 2.0 * r * ( 1 - r ) * RR + 2.0 * ( r * ( 1 - f) + ( 1 - r ) * f ) * RF + 2.0 * f * ( 1 - f ) * FF;  
  
  } else if  ( obs_selection == "LL" ) {
    
    exp = ( 1 -r ) * ( 1 -r ) * RR + 2.0 * ( 1 -r ) * ( 1 - f ) * RF + ( 1 - f ) * ( 1 - f ) * FF; 
  
  }
  
  return exp;
}


void getParamIndex( std::string input_string, const int& ibin,  int& idx, const std::string& real_sel = "" ) {
  
  if ( input_string.find("r_eff_") != std::string::npos || input_string.find("f_eff_") != std::string::npos( ) ) { 
    
    if ( input_string.find("ElEl") != std::string::npos ) {
      input_string += "_el_pT";
    } else if ( input_string.find("MuMu") != std::string::npos ) {
      input_string += "_mu_pT";
    }
    
  } else {
    // Replace the last two characters of the input string
    //
    input_string.replace( input_string.end()-2, input_string.end(), real_sel);
  }
  
  // Add bin index to name
  //
  input_string += ( "_" + std::to_string(ibin) );

  // Get corresponding index in list of parameters
  //
  auto it = g_param_names.find(input_string);
  idx = std::distance(g_param_names.begin(),it);
  
  if ( g_debug ) { Info("getParamIndex()","Index of parameter with name:\t %s\n --> %i", input_string.c_str(), idx ); }

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
   int idx_FF(-1);
   int idx_reff(-1);
   int idx_feff(-1);
   
   double exp(-1);
   
   for ( auto hist : g_histograms_pt ) {

     if ( g_debug ) { Info("myLikelihood()","Getting observed yield from histogram:\t %s", hist->GetName() ); }
   
     std::string this_hist_name = hist->GetName();
     
     auto pos = this_hist_name.end()-2;
     std::string obs_selection = this_hist_name.substr(pos);
     if ( g_debug ) { Info("myLikelihood()","Observed selection:\t %s", obs_selection.c_str() ); }
     
     for ( int ibin(0); ibin < hist->GetNbinsX()+2; ++ibin  ) {

       obs = hist->GetBinContent(ibin);
       if ( g_debug ) { Info("myLikelihood()","\t bin %i (pT = [%.2f,%.2f]) - yield = %f.", hist->GetBinLowEdge(ibin), hist->GetBinLowEdge(ibin+1), obs ); }

       // Need to find correspondent set of parameters
       //
       getParamIndex( this_hist_name, ibin, idx_RR, real_sel = "RR" );
       getParamIndex( this_hist_name, ibin, idx_RF, real_sel = "RF" );
       getParamIndex( this_hist_name, ibin, idx_FF, real_sel = "FR" );

       getParamIndex( this_hist_name, ibin, idx_reff);
       getParamIndex( this_hist_name, ibin, idx_feff);

       exp = getExpected(obs_selection, par[idx_reff], par[idx_feff], par[idx_RR], par[idx_RF], par[idx_FF] );

       likelihood += ( obs * log( exp ) - ( exp ) );
     }
      
   }
  
   result = likelihood;
}


void setParameters( TMinuit *myFitter, const int& nbins ) {

   int idx_r_el_pt(0);
   int idx_r_mu_pt(1);
   int idx_f_el_pt(2);
   int idx_f_mu_pt(3);
   
   double start(0.0001);
   double step(1e-6);
   double up(1e4);
   double dn(0.0)
   int ierflg(0);
   
   
   for ( unsigned int ibin = 0; ibin < nbins; ++ibin ) {
   
     std::string r_eff_el_pt = "r_eff_el_pT_" + std::to_string(ibin);
     std::string r_eff_mu_pt = "r_eff_mu_pT_" + std::to_string(ibin);
     std::string f_eff_el_pt = "f_eff_el_pT_" + std::to_string(ibin);
     std::string f_eff_mu_pt = "f_eff_mu_pT_" + std::to_string(ibin);
     
     g_param_names.push_back(r_eff_el_pt);
     g_param_names.push_back(r_eff_mu_pt);
     g_param_names.push_back(f_eff_el_pt);
     g_param_names.push_back(f_eff_mu_pt);
   
   }
   
   for ( unsigned int icharge = 0; icharge < g_charge.size(); ++icharge ) {
     for ( unsigned int iflav = 0; iflav < g_flavour.size(); ++iflav ) {
       for ( unsigned int isel = 0; isel < g_true_selection.size(); ++isel ) {
	 for ( unsigned int ibin = 0; ibin < nbins; ++ibin ) {
	   std::string param_name = charge.at(icharge) + "_" + flavour.at(iflav) + "_" + true_selection.at(isel) + "_" + std::to_string(ibin);
	   g_param_names.push_back(param_name);
         }
       } // loop over selection
     } // loop over flavour composition
   } // loop over charge

   if ( g_debug ) { Info("setParameters()","Adding parameters to fit function..."); }
   int ipar(0);
   for ( auto& par : g_param_names ) {
     
     if ( par.find("r_eff_") != std::string::npos || par.find("f_eff_") != std::string::npos ) {
       up    = 0.99999;
     }
     
     if ( g_debug ) { std::cout << "\t " << par.c_str() << std::endl; }
     
     myFitter->mnparm( ipar, par.c_str(), start, step, dn, up, ierflg );
     ++ipar;
   }
   
}

void fitEff( std::string input_path = "" ) {

   int nbins(-1);
   
   // Read the input histograms, and set the number of bins
   //
   getHists( input_path, nbins ); 
   
   // Set the number of parameteres of the fit
   //
   int NFLAV   = g_flavour.size(); 
   int NCHARGE = g_charge.size(); 
   int NCOMP   = g_true_selection.size();   
   
   // Total number of parameters to be estimated in the fit
   //
   const int NPAR = nbins * NFLAV * NCHARGE * NCOMP;
  
   Double_t arglist[2*NPAR];
  
   Int_t ierflg = 0;

   // Create the TMinuit object
   //
   TMinuit *myFitter = new TMinuit(NPAR);

   // Set the fitting function
   //
   myFitter->SetFCN(myLikelihood);

   // Set the parameters of the fit 
   //
   setParameters(myFitter, nbins);

   arglist[0] = 1;
   myFitter->mnexcm("SET ERR",arglist,1,ierflg);
      
   arglist[0] = 500; // maximum number of iterations
   myFitter->mnexcm("MIGRAD",arglist,2,ierflg);
 
   // print statistic
   //
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
    
   myFitter->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   
   if ( icstat < 3 ) cout << "There is  a problem with Minuit! "<<endl;
 
   std::cout << std::endl;
   std::cout << " Minimum of function = " << amin <<std:: endl;
   std::cout << " Estimated vert. distance to min. = " << edm <<std:: endl;
   std::cout << " Number of variable parameters = " << nvpar << std::endl;
   std::cout << " Highest number of parameters defined by user = " << nparx << std::endl;
   std::cout << " Status of covariance matrix = " << icstat << std::endl;
   std::cout << std::endl;

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
   
   delete myFitter; myFitter = nullptr;
  
}






// ***********************************************************************************************


void AddName(TString name,TLatex* &pL,float sep,float alt);
void GetLatex(bool MC,TLatex* &pLat1, TString sample,bool mod);
void ModifyLabels(TH1D* &pH,float offset,float size);
void SetHistoStyle(TH1D* &pH,int Marker,int color,int line);

vector<string> ReturnLegendPt();

void SetGeneralCuts(int minJets,int maxJets,float minMET,float maxMET,int minVtx,int maxVtx,bool bkgSideBand,string leptonID);

extern void chi2Eta(int &npar, double *gin, double &f, double *par, Int_t iflag);

extern void chi2(int &npar, double *gin, double &f, double *par, Int_t iflag);

extern void chi2EtaPt(int &npar, double *gin, double &f, double *par, Int_t iflag);

 std::vector <TH2D*> HistosforLikelihood(2); 

 int NBinsEta;
 int NBinsPt;
 int NBinsPtTotal;
 int CrackBin;
 float * EtaBin;
 float * PtBin;

// Defining Z peak in the  SS and OS distri.
const float meanOS=90.36, meanSS=88.88;  
const float sigmaOS=2.974, sigmaSS=3.097;
float nsigma=4;

int GminJets;
int GmaxJets;
float GminMET;
float GmaxMET;
int GminVtx;
int GmaxVtx;
bool GSideBand;
string GleptonID;

double Vector4DSS[7][7][3][3];
double Vector4D[7][7][3][3];

 float lowOS,llowOS,upOS,uupOS;
 float lowSS,llowSS,upSS,uupSS; 

class MisID {

public:

  MisID();
  ~MisID();
  void SetBinning();
  void SetZwindow(float window);
  void ShowDefaultValues();
  void ComputeRates(bool isMC,bool PlotEta,bool PlotPt,bool MorePlots,bool OtherMethods,string syst);
  void RatesTruthMatching2D(TH2D *TrueLeptons, TH2D *BadCharge,TH2D* &pRate,TH1D** &ArrayTruthMatchEta,TH1D** &ArrayTruthMatchmedEta,TH1D** &ArrayAlphaEta,int bin,bool ttbar);
  double ComputeErrorCociente(double Num,double Den,double errorN,double errorD);
  void RatesElecPosi2D(TH2D *TrueLeptons,TH2D *BadCharge,TH2D* &pRate);
  void ComparElecPosi2D(bool isMC);
  void FillEvent(float L1EtaCl,float L2EtaCl,float L1Pt,float L2Pt,float mll, int L1Charge,int L2Charge,int L1isChFlip,int L2isChFlip );
  // void GetLatex(bool MC,TLatex* &pLat1, TString sample);
  
  float m_sigma;
  int m_statusMatrix;
 

  TH1D *pRateTagProbeEta;
  TH1D *pRateDirectExtractionEta;
  TH1D *pRateLikelihoodEta;
  TH1D *pRateTagProbePt;
  TH1D *pRateTruthMatchEta;
  TH1D *pRateTruthMatchPt;
  TH1D *pRateDirectExtractionPt;
  TH1D *pRateLikelihoodPt;
  TH2D *pRateLikelihoodEtaPt;
  TH1D **ArrayLikelihoodEta;
  TH1D **ArrayTruthMatchZEta;
  TH1D **ArrayTruthMatchZEtaMCTruth;	
  TH1D **ArrayTruthMatchZmedEta;
  TH1D **ArrayAlphaZEta;
  TH1D **ArrayTruthMatchZmedEtaMCTruth;
  TH1D **ArrayAlphaZEtaMCTruth;	
  TH1D **ArrayAlphaPt;
  TH1D **ArrayElecPosi;
  void CreateHistos();
	

	
private:

  void ResetValues(double Mll,int ChPl1l2,float WeightLum,float L1Eta,float L2Eta,int L1isChFlip,int L2isChFlip);

 void TagProbe(TH1D* &pTagProbeSS,TH1D* &pTagProbeOS, float L1,float L2); 
 void DirectExtraction(TH1D* &pDirectExSS,TH1D* &pDirectExOS,float L1,float L2);
 void DirectExtraction2D(TH2D* &pDE2DSS,TH2D* &pDE2DOS,float L1,float L2);
 void MatrixLikelihood(TH2D* &pLikelihoodSS,TH2D* &pLikelihoodOS,float L1,float L2);
 void TruthMatchingZ(TH1D* &TrueLeptons,TH1D* &BadCharge, float L1, float L2, float L1Charge, float L2Charge,string var,TH1D* &HistoBadDR, TH1D* &HistoAllDR); 
 void TruthMatchingZsingle(TH1D* &TrueElectrons,TH1D* &BadElectrons, float L1, float L2, float L1Charge, float L2Charge,TH1D* &TruePositrons,TH1D* &BadPositrons);
 void TruthMatchingZsingle2D(TH2D* &TrueElectrons,TH2D* &BadElectrons, float L1Eta, float L2Eta,float L1Pt,float L2Pt, float L1Charge, float L2Charge,TH2D* &TruePositrons,TH2D* &BadPositrons);
 void TruthMatchingZ2D(TH2D* &TrueLeptons,TH2D* &BadCharge, float L1Eta, float L2Eta,float L1Pt,float L2Pt, float L1Charge, float L2Charge);
 void TruthMatchingZ2DOld(TH2D* &TrueLeptons,TH2D* &BadCharge, float L1Eta, float L2Eta,float L1Pt,float L2Pt, float L1Charge, float L2Charge);

 void FillVectorsLikelihoodEtaPt(float L1Eta,float L2Eta,float L1Pt,float L2Pt);


 void RatesTagProbe(TH1D *pTagProbeSS,TH1D *pTagProbeOS, TH1D* &pRate,string dependence);
 void RatesDirectExtraction(TH1D *pDirectExSS,TH1D *pDirectExOS, TH1D* &pRate,string dependence);
 void RatesDE2D(TH2D *pSS,TH2D *pOS,TH2D* &pRates2D,string dependence);
 void RatesLikelihood(TH2D *pLikelihoodSS,TH2D *pLikelihoodOS,TH1D* &pRate,string dependence); 
  void RatesTruthMatching(TH1D *TrueLeptons, TH1D *BadCharge,TH1D* &pRate,TString sample);
  void RatesLikelihoodEtaPt(TH2D* &pRate,bool isMC);

 void PlotRates(TH1D *TP,TH1D *DE, TH1D *Likelihood,TH2D *DE2D,TH1D *pTM,bool MC,bool OtherMethods); 
 void PlotTruthMatching(TH1D *pRateTruthMatch,TH1D *pRateElecTruthMatch,TH1D *pRatePosiTruthMatch,bool isTTbar);
 void PlotRates2D(bool isMC,bool morePlots);

 void ShowMorePlots(bool PlotEta,bool PlotPt,bool OtherMethods,bool isMC);
 void ShowMoreOtherMethods(TH1D *pTagSS,TH1D *pTagOS,TH1D *pDESS,TH1D *pDEOS,TH2D *pDE2DSS,TH2D *pDE2DOS);
 void ShowMoreDefaultMethods(bool isMC,TH2D *pLikelihoodSS,TH2D *pLikelihoodOS,TH1D *pTrueLeptons,TH1D *pBadCharge,
                             TH1D *pTrueElectrons,TH1D *pBadElectrons,TH1D *pTruePositrons,TH1D *pBadPositrons);
  void ShowEventsMatrix();

 TH1D *pTagProbeSSEta;
 TH1D *pTagProbeOSEta;
 TH1D *pDirectExSSEta;
 TH1D *pDirectExOSEta;
 TH2D *pDE2DSSEta;
 TH2D *pDE2DOSEta;
 TH2D *pLikelihoodSSEta;
 TH2D *pLikelihoodOSEta;
 
 TH2D *pRateDE2DEta;  
 TH1D *pTrueLeptonsEta;
 TH1D *pBadChargeEta;
 TH1D *hBadDR;
 TH1D *hAllDR;
 TH1D *pTrueElectronsEta;
 TH1D *pBadElectronsEta;
 TH1D *pTruePositronsEta;
 TH1D *pBadPositronsEta;
 TH1D *pRateElecTruthMatchEta;
 TH1D *pRatePosiTruthMatchEta;

 TH1D *pTagProbeSSPt;
 TH1D *pTagProbeOSPt;
 TH1D *pDirectExSSPt;
 TH1D *pDirectExOSPt;
 TH2D *pDE2DSSPt;
 TH2D *pDE2DOSPt;
 TH2D *pLikelihoodSSPt;
 TH2D *pLikelihoodOSPt;

 TH2D *pRateDE2DPt;  
 TH1D *pTrueLeptonsPt;
 TH1D *pBadChargePt;
 TH1D *pTrueElectronsPt;
 TH1D *pBadElectronsPt;
 TH1D *pTruePositronsPt;
 TH1D *pBadPositronsPt;
 TH1D *pRateElecTruthMatchPt;
 TH1D *pRatePosiTruthMatchPt;

 TH2D *pTEventsTM2D;
 TH2D *pBEventsTM2D;
 TH2D *pTEventsTM2DMCTruth;
 TH2D *pBEventsTM2DMCTruth;	
 TH2D *pRatesTM2D;
 TH2D *pRatesTM2DMCTruth;	
 TH2D *pTElectrons;
 TH2D *pBElectrons;
 TH2D *pTPositrons;
 TH2D *pBPositrons;
 TH2D *pRatesElec2D;
 TH2D *pRatesPosi2D;

	
	
 TH1D *pEtaSSAll;
 TH1D *pEtaOSAll;
 TH1D *pEtaSSsig;
 TH1D *pEtaOSsig;

 TH1D *pPtSSAll;
 TH1D *pPtOSAll;
 TH1D *pPtSSsig;
 TH1D *pPtOSsig;

  vector<float> vec_L1EtaCl;
  vector<float> vec_L2EtaCl;
  vector<float> vec_mll;
  vector<float> vec_L1Pt;
  vector<float> vec_L2Pt;
  vector<int> vec_L1Charge;
  vector<int> vec_L2Charge;
  vector<int> vec_L1isChFlip;
  vector<int> vec_L2isChFlip;

  double m_Mll;
  int m_ChPl1l2; 
  float m_weight;
  float m_L1Eta,m_L2Eta;
  int m_L1isChFlip, m_L2isChFlip;
  int m_L1isTruthMatchedOther,m_L2isTruthMatchedOther;
  int m_L1isTruthMatchedElNoProdVtx,m_L2isTruthMatchedElNoProdVtx;

};


void MisIdCalc_run2()
{

 TH1::AddDirectory(kFALSE);
 TH2::AddDirectory(kFALSE);
 // TH1::SetDefaultSumw2(kTRUE);

  MisID ComputeMisID;

   ComputeMisID.SetBinning();
   ComputeMisID.SetZwindow(nsigma);

   ComputeMisID.ShowDefaultValues();

  bool isMC=true;
  bool PlotEta=true;
  bool PlotPt=true;
  bool MorePlots=false;
  bool OtherMethods=false;  
     
  bool RemoveFakes=false; // Remove fakes from Z peak;
  int bin=NBinsPt; 
  int minJets=0; // Minimum number jets required in the event. If RemoveFakes==true then jets=1;
  int maxJets=1000;
  float minMET=0;
  float maxMET=1000;
  int minVtx=1;
  int maxVtx=1000;
  bool bkgSideBand=true;
  string leptonID="tight";
  bool ttH=true;  
  string syst="default";
  // ##################  ttH ##################################  
  TString dataTTH="../../Merged_Melb15_HTopMultilep_test23_DxAOD/Data/physics_Egamma.root";
  TString FileMCttH="input_file.root";

  TString TTbarTTH="../../Merged_Melb15_HTopMultilep_test23_DxAOD/tops/PowhegPythia_P2011C_ttbar.root";

   TString File;
   TString inputTTbar;
   TString TTbarSyst;
   TString FileMC;
   TString OutputFile;
   TString OutputFile_MC;

   if (isMC==true) File=FileMCttH;
   else File=dataTTH;
   
   FileMC=FileMCttH;
   inputTTbar=TTbarTTH;
   TTbarSyst=TTbarTTH;  
   if (isMC==true) OutputFile="Rates-MC-TLH-";
   else OutputFile="Rates-Data-TLH-";
	 //}
   TFile *pInputFile = new TFile(File); 
	
   SetGeneralCuts(minJets,maxJets,minMET,maxMET,minVtx,maxVtx,bkgSideBand,leptonID);

   int nbEvts,test;
   int evt=0;
  float L2Charge,L1Charge;
  float L1Pt,L2Pt,L1Eta,L2Eta;
  float mll01;
  int ChPl1l2=0;
  int nel,nmuon;
  std::vector<float> *el_pt=0;
  std::vector<float> *el_eta=0;
  std::vector<float> *el_calo_eta=0;
  std::vector<float> *el_trkcharge=0;
  std::vector<int> *el_LHTight=0;
  std::vector<int> *el_isChFlip=0;
  std::vector<int> *el_isBrem=0;
  int NblJets;
  unsigned  int NbbJets;
  double Weight=-99.9;
  int ML2Tight;
  bool L1Tight=false;
  bool L2Tight=false;

  TTree *mTree;
  mTree=(TTree*)pInputFile->Get("physics");    
  if (isMC){               
    mTree->SetBranchAddress("el_isChFlip",&el_isChFlip);
    mTree->SetBranchAddress("el_isBrem",&el_isBrem);
  }
  mTree->SetBranchAddress("nel",&nel);
  mTree->SetBranchAddress("nmuon",&nmuon);
  mTree->SetBranchAddress("el_LHTight",&el_LHTight);
  mTree->SetBranchAddress("mll01",&mll01);
  mTree->SetBranchAddress("el_pt",&el_pt);
  mTree->SetBranchAddress("el_eta",&el_eta);
  mTree->SetBranchAddress("el_calo_eta",&el_calo_eta);
  mTree->SetBranchAddress("el_trkcharge",&el_trkcharge);
  mTree->SetBranchAddress("njets",&NblJets);
  mTree->SetBranchAddress("nMediumBjets",&NbbJets);
  
  nbEvts=mTree->GetEntries();

  string data="data";
  if (isMC) data="MC";

  cout<< "------->> Running on '" << data  << "' for '"<< GleptonID << "' leptons"  << endl; 
  std::cout<< "------->> Opening file '" << File <<"' with '" << nbEvts << "' events" <<std::endl;
  
  int TypeID=-99;
  int lepton=-99;
  if (leptonID=="tight") lepton=1;
  if (leptonID=="antitight") lepton=-1;
  if (leptonID=="loose") lepton=0;	
  ComputeMisID.CreateHistos();

  while(evt < nbEvts){
    mTree->GetEntry(evt);
    if (nel==2 && nmuon==0){
      
      Weight=1.;
      if ((el_LHTight[0][0])==1 && (el_LHTight[0][1])==1 && (NblJets+NbbJets>=minJets)){	
	test+=1;
	
	ChPl1l2=el_trkcharge[0][0]*el_trkcharge[0][1];
	if (el_LHTight[0][0]==1) {L1Tight=true;}
	if (el_LHTight[0][1]==1) {L2Tight=true;}
	if (el_LHTight[0][0]==1 && el_LHTight[0][1]==1){ML2Tight=1;}
	if ((ChPl1l2<0 && mll01/1000>meanOS-3*nsigma*sigmaOS && mll01/1000<meanOS+3*nsigma*sigmaOS) || (ChPl1l2>0 && mll01/1000>meanSS-3*nsigma*sigmaSS && mll01/1000<meanSS+3*nsigma*sigmaSS)){ 
	  if (L1Tight==true && L2Tight==true &&  ML2Tight==1) TypeID=1;
	  if (L1Tight==false && L2Tight==false) TypeID=-1;
	  if (lepton==0) TypeID=0;      		 
	  if (lepton==TypeID && (NblJets+NbbJets>=minJets) && (NblJets+NbbJets<=maxJets)){ 
	    ComputeMisID.FillEvent(el_calo_eta[0][0],el_calo_eta[0][1],el_pt[0][0],el_pt[0][1],mll01,el_trkcharge[0][0],el_trkcharge[0][1],el_isChFlip[0][0],el_isChFlip[0][1]);
	  }
	}
      }
      }   
    evt++;
  }

  ComputeMisID.ComputeRates(isMC,PlotEta,PlotPt,MorePlots,OtherMethods,syst);   
  
  return;
  
}

void SetGeneralCuts(int minJets,int maxJets,float minMET,float maxMET,int minVtx,int maxVtx,bool bkgSideBand,string leptonID){


  GminJets=minJets;
  GmaxJets=maxJets;
  GminMET=minMET;
  GmaxMET=maxMET;
  GminVtx=minVtx;
  GmaxVtx=maxVtx;  
  GSideBand=bkgSideBand;
  GleptonID=leptonID;

  return;

}

MisID::MisID()
{

}

MisID::~MisID()
{

}

void MisID::FillEvent(float L1EtaCl,float L2EtaCl,float L1Pt,float L2Pt,float mll, int L1Charge,int L2Charge,int L1isChFlip,int L2isChFlip ){
  vec_L1EtaCl.push_back(L1EtaCl);
  vec_L2EtaCl.push_back(L2EtaCl);
  vec_L1Pt.push_back(L1Pt/1000);
  vec_L2Pt.push_back(L2Pt/1000);
  vec_mll.push_back(mll/1000);
  vec_L1Charge.push_back(L1Charge);
  vec_L2Charge.push_back(L2Charge);
  vec_L1isChFlip.push_back(L1isChFlip);
  vec_L2isChFlip.push_back(L2isChFlip);
  return;

}


void MisID::SetBinning()
{

  // This function gives the general values of the binning in eta and pt 
 
 NBinsPtTotal=4;  // Number of total bins in Pt
 NBinsPt=3;    // Number of bins in Pt used for likelihood
 NBinsEta=7;   // Number of bins en Eta 7.

 CrackBin=4;   //Crack between between 1.37>|eta|<1.52; 

 EtaBin= new float [NBinsEta+1];

 EtaBin[0]=0.;   EtaBin[1]=0.6; EtaBin[2]=1.1;  EtaBin[3]=1.37;
 EtaBin[4]=1.52; EtaBin[5]=1.7; EtaBin[6]=2.3;  EtaBin[7]=2.6;  	

 PtBin= new float [NBinsPtTotal+1];

	PtBin[0]=0;  PtBin[1]=60;   PtBin[2]=90; PtBin[3]=130; PtBin[4]=1000; 
	
 
 return;

} // end SetBinning

void MisID::SetZwindow(float window)
{

 // Defining regions to extract the background.

 m_sigma=window;

 lowOS=meanOS-m_sigma*sigmaOS;
 llowOS=lowOS-2*m_sigma*sigmaOS;
 upOS=meanOS+m_sigma*sigmaOS;
 uupOS=upOS+2*m_sigma*sigmaOS;

 lowSS=meanSS-m_sigma*sigmaSS;
 llowSS=lowSS-2*m_sigma*sigmaSS;
 upSS=meanSS+m_sigma*sigmaSS;
 uupSS=upSS+2*m_sigma*sigmaSS;


  return;

} // end SetZwindow

void MisID::ShowDefaultValues()
{

  std::cout << "---------------------------------- Eta binning ----------------------------------" << std::endl;
  std::cout<< "There are '"<< NBinsEta << "' ETA bins:" << std::endl;

  for (int i=0; i < NBinsEta; i++)  std::cout<< "EtaBin " << i+1  << " = ["<< EtaBin[i] <<"," << EtaBin[i+1] << "], " << std::endl;
  
  std::cout << "---------------------------------- Pt binning ----------------------------------" << std::endl;
  std::cout<< "There are '"<< NBinsPt << "' Pt bins:" << std::endl;

  for (int i=0; i < NBinsPt; i++)  std::cout<< "PtBin " << i+1  << " = ["<< PtBin[i] <<"," << PtBin[i+1] << "], " << std::endl;

  std::cout << "---------------------------------- Background extracted for a Z window ---------------------------------- " << std::endl;
  
  std:: cout << "SS distribution: A = [" << llowSS << "," << lowSS <<"], B = [" << lowSS <<"," << upSS << "], C = [" << upSS << "," << uupSS << "]" << std::endl;
  std:: cout << "OS distribution: A = [" << llowOS << "," << lowOS <<"], B = [" << lowOS <<"," << upOS << "], C = [" << upOS << "," << uupOS << "]" << std::endl;
  

  return;
}

void MisID::CreateHistos()
{
   

    pTagProbeSSEta = new TH1D("TagProbeSSEta",";|#eta|; Number of events",NBinsEta,EtaBin);  
    pTagProbeOSEta = new TH1D("TagProbeOSEta",";|#eta|; Number of events",NBinsEta,EtaBin);

    pTagProbeSSEta->Sumw2();
    pTagProbeOSEta->Sumw2();

    pDirectExSSEta = new TH1D("DirectExtractionSSEta",";|#eta|; Number of events",NBinsEta,EtaBin);
    pDirectExOSEta = new TH1D("DirectExtractionOSEta",";|#eta|; Number of events",NBinsEta,EtaBin);

    pDirectExSSEta->Sumw2();
    pDirectExOSEta->Sumw2();

    pDE2DSSEta = new TH2D("DE2DSSEta",";|#eta|_{leadig lepton}; |#eta|_{sub-leading lepton};Number of SS events",NBinsEta,EtaBin,NBinsEta,EtaBin);
    pDE2DOSEta = new TH2D("DE2DOSEta",";|#eta|_{leadig lepton}; |#eta|_{sub-leading lepton};Number of OS events",NBinsEta,EtaBin,NBinsEta,EtaBin);

    pDE2DSSEta->Sumw2();
    pDE2DOSEta->Sumw2();

    pLikelihoodSSEta=new TH2D("LikelihoodSS",";|#eta|; |#eta|;Number of SS events",NBinsEta,EtaBin,NBinsEta,EtaBin);
    pLikelihoodOSEta=new TH2D("LikelihoodOS",";|#eta|; |#eta|;Number of OS events",NBinsEta,EtaBin,NBinsEta,EtaBin);

    pLikelihoodSSEta->Sumw2();
    pLikelihoodOSEta->Sumw2();

    pRateTagProbeEta= new TH1D("TagProbeEta",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
    pRateDirectExtractionEta= new TH1D("DirectExtraction",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
    pRateLikelihoodEta=new TH1D("Likelihood",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
    pRateDE2DEta = new TH2D("RateDE2D",";|#eta|_{leading lepton};|#eta|_{sub-leading lepton};#epsilon_{mis-id}",NBinsEta,EtaBin,NBinsEta,EtaBin);   

    hBadDR = new TH1D("NbBadElectronsDR",";dR(e,closest jet);Number of electrons",10,0,10);
    hAllDR = new TH1D("NbAllElectronsDR",";dR(e,closest jet);Number of electrons",10,0,10);

    pTrueLeptonsEta=new TH1D("TrueLeptonsEta",";|#eta|;Number of leptons",NBinsEta,EtaBin);
    pBadChargeEta=new TH1D("BadChargeEta",";|#eta|;Number of leptons",NBinsEta,EtaBin);

    pTrueLeptonsEta->Sumw2();
    pBadChargeEta->Sumw2();


    pTrueElectronsEta=new TH1D("TrueElectronsEta",";|#eta|;Number of electrons",NBinsEta,EtaBin);
    pBadElectronsEta=new TH1D("BadElectronsEta",";|#eta|;Number of electrons",NBinsEta,EtaBin);  
    pTruePositronsEta=new TH1D("TruePositronsEta",";|#eta|;Number of positrons",NBinsEta,EtaBin);
    pBadPositronsEta=new TH1D("BadPositronsEta",";|#eta|;Number of positrons",NBinsEta,EtaBin);  

    pTrueElectronsEta->Sumw2();
    pBadElectronsEta->Sumw2();  
    pTruePositronsEta->Sumw2();
    pBadPositronsEta->Sumw2();  


    pRateTruthMatchEta=new TH1D("RateEta",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
    pRateElecTruthMatchEta=new TH1D("RateEtaElec",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
    pRatePosiTruthMatchEta=new TH1D("RateEtaPosi",";|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);

  
     pTagProbeSSPt = new TH1D("TagProbeSSPt",";p_{T} [GeV]; Number of events",NBinsPtTotal,PtBin);  
     pTagProbeOSPt = new TH1D("TagProbeOSPt",";p_{T} [GeV]; Number of events",NBinsPtTotal,PtBin);

     pTagProbeSSPt->Sumw2();
     pTagProbeOSPt->Sumw2();

     pDirectExSSPt = new TH1D("DirectExtractionSSPt",";p_{T} [GeV]; Number of events",NBinsPtTotal,PtBin);
     pDirectExOSPt = new TH1D("DirectExtractionOSPt",";p_{T} [GeV]; Number of events",NBinsPtTotal,PtBin);
 
     pDirectExSSPt->Sumw2();
     pDirectExOSPt->Sumw2();


     pDE2DSSPt = new TH2D("DE2DSSPt",";Leading lepton p_{T} [GeV]; Sub-leading lepton p_{T} [GeV];Number of SS events",NBinsPtTotal,PtBin,NBinsPtTotal,PtBin);
     pDE2DOSPt = new TH2D("DE2DOSPt",";Leading lepton p_{T} [Gev]; Sub-leading lepton p_{T} [GeV];Number of OS events",NBinsPtTotal,PtBin,NBinsPtTotal,PtBin);

     pDE2DSSPt->Sumw2();
     pDE2DOSPt->Sumw2();

     pLikelihoodSSPt=new TH2D("LikelihoodSSPt",";Sub-leading lepton p_{T} [GeV]; Leading lepton p_{T} [GeV];Number of SS events",NBinsPtTotal,PtBin,NBinsPtTotal,PtBin);
     pLikelihoodOSPt=new TH2D("LikelihoodOSPt",";Sub-leading lepton p_{T} [GeV]; Leading lepton p_{T} [GeV];Number of OS events",NBinsPtTotal,PtBin,NBinsPtTotal,PtBin);
 
     pLikelihoodSSPt->Sumw2();
     pLikelihoodOSPt->Sumw2();

     pRateTagProbePt= new TH1D("TagProbePt",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);
     pRateDirectExtractionPt= new TH1D("DirectExtractionPt",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);
     pRateLikelihoodPt=new TH1D("LikelihoodPt",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);
     pRateDE2DPt = new TH2D("RateDE2DPt",";Leading p_{T} [GeV];Sub-leading p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin,NBinsPtTotal,PtBin);   
 
     pTrueLeptonsPt=new TH1D("TrueLeptonsPt",";p_{T} [GeV];Number of leptons",NBinsPtTotal,PtBin);
     pBadChargePt=new TH1D("BadChargePt",";p_{T} [GeV];Number of leptons",NBinsPtTotal,PtBin);

     pTrueLeptonsPt->Sumw2();
     pBadChargePt->Sumw2();

     pTrueElectronsPt=new TH1D("TrueElectronsPt",";p_{T} [GeV];Number of electrons",NBinsPtTotal,PtBin);
     pBadElectronsPt=new TH1D("BadElectronsPt",";p_{T} [GeV];Number of electrons",NBinsPtTotal,PtBin);  
     pTruePositronsPt=new TH1D("TruePositronsPt",";p_{T} [GeV];Number of electrons",NBinsPtTotal,PtBin);
     pBadPositronsPt=new TH1D("BadPositronsPt",";p_{T} [GeV];Number of electrons",NBinsPtTotal,PtBin);  

     pTrueElectronsPt->Sumw2();
     pBadElectronsPt->Sumw2();  
     pTruePositronsPt->Sumw2();
     pBadPositronsPt->Sumw2();  

     pRateTruthMatchPt=new TH1D("TruthMatchPt",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);   
     pRateElecTruthMatchPt=new TH1D("RatePtElec",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);
     pRatePosiTruthMatchPt=new TH1D("RatePtPosi",";p_{T} [GeV];#epsilon_{mis-id}",NBinsPtTotal,PtBin);  

     pRateLikelihoodEtaPt=new TH2D("LikelihoodEtaPt","Likelihood",NBinsEta,EtaBin,NBinsPtTotal,PtBin);

     ArrayLikelihoodEta= new TH1D*[NBinsPt];
 
     pTEventsTM2D=new TH2D("tTM2D","ElectronsT TM2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
     pBEventsTM2D=new TH2D("bTM2D","ElectronsB TM2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
	pTEventsTM2DMCTruth=new TH2D("tTM2D","ElectronsT TM2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
	pBEventsTM2DMCTruth=new TH2D("bTM2D","ElectronsB TM2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);

	
	
     pRatesTM2D=new TH2D("TM2D","Truth-matching 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
	 pRatesTM2DMCTruth=new TH2D("TM2DMCTruth","Truth-matching 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);

     pTElectrons=new TH2D("Elec2D","Electrons 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
     pBElectrons=new TH2D("ElecB2D","ElectronsB 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
     pTPositrons=new TH2D("Posi2D","Positrons 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
     pBPositrons=new TH2D("PosiB2D","PositronsB 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
	 pRatesElec2D=new TH2D("RatesElec2D","RatesElec 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);
	 pRatesPosi2D=new TH2D("RatesPosi2D","RatesPosi 2D",NBinsEta,EtaBin,NBinsPtTotal,PtBin);

     pTEventsTM2D->Sumw2();
     pBEventsTM2D->Sumw2();
	 pTEventsTM2DMCTruth->Sumw2();
	 pBEventsTM2DMCTruth->Sumw2();
	
	pTElectrons->Sumw2();
	pBElectrons->Sumw2();
    pTPositrons->Sumw2();
	pBPositrons->Sumw2();
	

     ArrayTruthMatchZEta= new TH1D*[NBinsPtTotal];
	 ArrayTruthMatchZEtaMCTruth= new TH1D*[NBinsPtTotal];
     ArrayTruthMatchZmedEta= new TH1D*[NBinsPtTotal];
     ArrayAlphaZEta= new TH1D*[NBinsPtTotal];

	 ArrayTruthMatchZmedEtaMCTruth= new TH1D*[NBinsPtTotal];
	 ArrayAlphaZEtaMCTruth= new TH1D*[NBinsPtTotal];

	 ArrayAlphaPt= new TH1D*[NBinsEta];
     ArrayElecPosi= new TH1D*[NBinsPtTotal];
 
    for (int j=0; j<NBinsPtTotal; j++){
       ArrayTruthMatchZEta[j] = new TH1D(Form("TMZEta%d",j),"Histogram;|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
	   ArrayTruthMatchZEtaMCTruth[j] = new TH1D(Form("MCTMZEta%d",j),"Histogram;|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
       ArrayTruthMatchZmedEta[j] = new TH1D(Form("TMZmedEta%d",j),"Histogram;|#eta|;#epsilon(|#eta|,p_{T})/#epsilon(|#eta|)",NBinsEta,EtaBin);
       ArrayAlphaZEta[j]= new TH1D(Form("AlphaZEta%d",j),"Histogram;|#eta|;#alpha_{Z}(|#eta|,p_{T})",NBinsEta,EtaBin);
	   ArrayTruthMatchZmedEtaMCTruth[j] = new TH1D(Form("TMZmedEtaMC%d",j),"Histogram;|#eta|;#epsilon(|#eta|,p_{T})/#epsilon(|#eta|)",NBinsEta,EtaBin);
	   ArrayAlphaZEtaMCTruth[j]= new TH1D(Form("AlphaZEtaMC%d",j),"Histogram;|#eta|;#alpha_{Z}(|#eta|,p_{T})",NBinsEta,EtaBin);
	
       ArrayElecPosi[j]= new TH1D(Form("ElecPosi%d",j),"Histogram;|#eta|;#epsilon_{mis-id,e^{-}}/#epsilon_{mis-id,e^{+}}",NBinsEta,EtaBin);
     }

    for (int i=0; i<NBinsEta; i++) ArrayAlphaPt[i]= new TH1D(Form("AlphattbarEta%d",i),"Histogram;|#eta|;#alpha_{t#bar{t}}(|#eta|,p_{T})",NBinsPtTotal,PtBin);


   

     pEtaSSAll = new TH1D("pSSEta",";Leading lepton #eta;Number of same-sign events",NBinsEta,EtaBin);
     pEtaOSAll = new TH1D("pOSEta",";Leading lepton #eta;Number of opposite-sign events",NBinsEta,EtaBin);

     pEtaSSsig = new TH1D("pSSEtasig",";Leading lepton #eta;Number of same-sign events",NBinsEta,EtaBin);
     pEtaOSsig = new TH1D("pOSEtasig",";Leading lepton #eta;Number of opposite-sing events",NBinsEta,EtaBin);

     pPtSSAll = new TH1D("pSSPt",";Leading lepton p_{T};Number of same-sign events",NBinsPt,PtBin);
     pPtOSAll = new TH1D("pOSPt",";Leading lepton p_{T};Number of opposite-sign events",NBinsPt,PtBin); 
 
     pPtSSsig = new TH1D("pSSPtsig",";Leading lepton p_{T};Number of same-sign events",NBinsPt,PtBin);
     pPtOSsig = new TH1D("pOSPtsig",";Leading lepton p_{T};Number of opposite-sign events",NBinsPt,PtBin); 


    return;

}

void MisID::ResetValues(double Mll,int ChPl1l2,float WeightLum, float L1Eta,float L2Eta,int L1isChFlip,int L2isChFlip)
{

  m_Mll=Mll;
  m_ChPl1l2=ChPl1l2;
  m_weight=WeightLum;
  m_L1Eta=L1Eta;
  m_L2Eta=L2Eta;
  m_L1isChFlip=L1isChFlip;
  m_L2isChFlip=L2isChFlip;

  return;

} // end ResetValues

void MisID::ComputeRates(bool isMC,bool PlotEta,bool PlotPt,bool MorePlots,bool OtherMethods,string syst)
{
   
  for (int i=0;i<7; i++){
    for (int j=0;j<7; j++){
      for(int k=0;k<3; k++){
	for(int m=0;m<3; m++){
	  Vector4DSS[i][j][k][m]=0;
	  Vector4D[i][j][k][m]=0;
	}
      }
    }
  }
  double Weight;
  int count=0;
  int ChPl1l2=0;
  int L1Charge;
  int L2Charge;
  double L1EtaCl;
  double L2EtaCl;
  double L1Pt;
  double L2Pt;
  while(count < vec_L2Charge.size()){
    // if (count%10000==0){
    //   cout<<"Value of l1eta vector "<<vec_L1EtaCl[count]<<endl;
    //   cout<<"Value of l2eta vector "<<vec_L2EtaCl[count]<<endl;
    //   cout<<"Value of l1pt vector "<<vec_L1Pt[count]<<endl;
    //   cout<<"Value of l2pt vector "<<vec_L2Pt[count]<<endl;
    //   cout<<"Value of l1charge vector "<<vec_L1Charge[count]<<endl;
    //   cout<<"Value of l2charge vector "<<vec_L2Charge[count]<<endl;
    //   cout<<"Value of l1ischflip vector "<<vec_L1isChFlip[count]<<endl;
    //   cout<<"Value of l2ischflip vector "<<vec_L2isChFlip[count]<<endl;
    //   cout<<"Value of mll vector "<<vec_mll[count]<<endl;
    // }
    ChPl1l2=vec_L1Charge[count]*vec_L2Charge[count];
    L1EtaCl=vec_L1EtaCl[count];
    L2EtaCl=vec_L2EtaCl[count];
    L1Pt=vec_L1Pt[count];
    L2Pt=vec_L2Pt[count];
    L1Charge=vec_L1Charge[count];
    L2Charge=vec_L2Charge[count];
    Weight=1.;
    if (isMC){
      ResetValues(vec_mll[count],ChPl1l2,Weight,fabs(L1EtaCl),fabs(L2EtaCl),vec_L1isChFlip[count],vec_L2isChFlip[count]);
    }
    else{
      ResetValues(vec_mll[count],ChPl1l2,Weight,fabs(L1EtaCl),fabs(L2EtaCl),0,0);
    }
    TagProbe(pTagProbeSSEta,pTagProbeOSEta,fabs(L1EtaCl),fabs(L2EtaCl));         
    DirectExtraction(pDirectExSSEta,pDirectExOSEta,fabs(L1EtaCl),fabs(L2EtaCl));
    DirectExtraction2D(pDE2DSSEta,pDE2DOSEta,fabs(L1EtaCl),fabs(L2EtaCl));
    MatrixLikelihood(pLikelihoodSSEta,pLikelihoodOSEta,fabs(L1EtaCl),fabs(L2EtaCl)); 
    TagProbe(pTagProbeSSPt,pTagProbeOSPt,L1Pt,L2Pt);         
    DirectExtraction(pDirectExSSPt,pDirectExOSPt,L1Pt,L2Pt);  
    DirectExtraction2D(pDE2DSSPt,pDE2DOSPt,L1Pt,L2Pt);    
    MatrixLikelihood(pLikelihoodSSPt,pLikelihoodOSPt,L1Pt,L2Pt); 
    
    if (L1Pt<PtBin[NBinsPt] && L2Pt<PtBin[NBinsPt]) FillVectorsLikelihoodEtaPt(fabs(L1EtaCl),fabs(L2EtaCl),L1Pt,L2Pt);	  
    if (isMC){
      TruthMatchingZ(pTrueLeptonsEta,pBadChargeEta, fabs(L1EtaCl), fabs(L2EtaCl), L1Charge, L2Charge,"eta",hBadDR,hAllDR);
      TruthMatchingZsingle(pTrueElectronsEta,pBadElectronsEta,fabs(L1EtaCl),fabs(L2EtaCl),L1Charge,L2Charge,pTruePositronsEta,pBadPositronsEta);
      TruthMatchingZ(pTrueLeptonsPt,pBadChargePt, L1Pt,L2Pt, L1Charge, L2Charge,"pt",hBadDR,hAllDR);
      TruthMatchingZsingle(pTrueElectronsPt,pBadElectronsPt,L1Pt,L2Pt,L1Charge,L2Charge,pTruePositronsPt,pBadPositronsPt);
      TruthMatchingZ2D(pTEventsTM2D,pBEventsTM2D,fabs(L1EtaCl),fabs(L2EtaCl),L1Pt,L2Pt,L1Charge,L2Charge);
      TruthMatchingZ2DOld(pTEventsTM2DMCTruth,pBEventsTM2DMCTruth,fabs(L1EtaCl),fabs(L2EtaCl),L1Pt,L2Pt,L1Charge,L2Charge);		 
      TruthMatchingZsingle2D(pTElectrons,pBElectrons,fabs(L1EtaCl), fabs(L2EtaCl),L1Pt,L2Pt,L1Charge,L2Charge,pTPositrons,pBPositrons);
      }
    count=count+1;
  }

  RatesTagProbe(pTagProbeSSEta,pTagProbeOSEta,pRateTagProbeEta,"eta"); 
  RatesDirectExtraction(pDirectExSSEta,pDirectExOSEta,pRateDirectExtractionEta,"eta"); 
  RatesDE2D(pDE2DSSEta,pDE2DOSEta,pRateDE2DEta,"eta"); 
  RatesLikelihood(pLikelihoodSSEta,pLikelihoodOSEta,pRateLikelihoodEta,"eta");
  ShowEventsMatrix();
  RatesTagProbe(pTagProbeSSPt,pTagProbeOSPt,pRateTagProbePt,"pt");
  RatesDirectExtraction(pDirectExSSPt,pDirectExOSPt,pRateDirectExtractionPt,"pt");
  RatesDE2D(pDE2DSSPt,pDE2DOSPt,pRateDE2DPt,"pt");
  RatesLikelihood(pLikelihoodSSPt,pLikelihoodOSPt,pRateLikelihoodPt,"pt");
  RatesLikelihoodEtaPt(pRateLikelihoodEtaPt,isMC);
  
  if (isMC){
    RatesTruthMatching(pTrueLeptonsEta,pBadChargeEta,pRateTruthMatchEta,"Leptons");
    RatesTruthMatching(pTrueElectronsEta,pBadElectronsEta,pRateElecTruthMatchEta,"Electrons");
    RatesTruthMatching(pTruePositronsEta,pBadPositronsEta,pRatePosiTruthMatchEta,"Positrons");
    RatesTruthMatching(pTrueLeptonsPt,pBadChargePt,pRateTruthMatchPt,"Leptons");
    RatesTruthMatching(pTrueElectronsPt,pBadElectronsPt,pRateElecTruthMatchPt,"Electrons");
    RatesTruthMatching(pTruePositronsPt,pBadPositronsPt,pRatePosiTruthMatchPt,"Positrons");
    RatesTruthMatching2D(pTEventsTM2D,pBEventsTM2D,pRatesTM2D,ArrayTruthMatchZEta,ArrayTruthMatchZmedEta,ArrayAlphaZEta,1,false);
    RatesTruthMatching2D(pTEventsTM2DMCTruth,pBEventsTM2DMCTruth,pRatesTM2DMCTruth,ArrayTruthMatchZEtaMCTruth,ArrayTruthMatchZmedEtaMCTruth,ArrayAlphaZEtaMCTruth,1,false);  
    RatesElecPosi2D(pTElectrons,pBElectrons,pRatesElec2D);
    RatesElecPosi2D(pTPositrons,pBPositrons,pRatesPosi2D);
    
    if (PlotEta) PlotTruthMatching(pRateTruthMatchEta,pRateElecTruthMatchEta,pRatePosiTruthMatchEta,false);
    if (PlotPt)  PlotTruthMatching(pRateTruthMatchPt,pRateElecTruthMatchPt,pRatePosiTruthMatchPt,false);
    ComparElecPosi2D(isMC);
  }
  if (PlotEta) PlotRates(pRateTagProbeEta,pRateDirectExtractionEta,pRateLikelihoodEta,pRateDE2DEta,pRateTruthMatchEta,isMC,OtherMethods);
  if (PlotPt)  PlotRates(pRateTagProbePt,pRateDirectExtractionPt,pRateLikelihoodPt,pRateDE2DPt,pRateTruthMatchPt,isMC,OtherMethods);
  PlotRates2D(isMC,MorePlots);
  if (MorePlots) ShowMorePlots(PlotEta,PlotPt,OtherMethods,isMC);
  return;

  } // end ComputeRates

void MisID::TagProbe(TH1D* &pTagProbeSS,TH1D* &pTagProbeOS, float L1,float L2)
{
 
 float tag=pTagProbeSS->GetBinLowEdge(2);

if ((L1<tag) || (L2<tag))
{
  if (m_ChPl1l2>0){
    
         if ((L1<tag) && (L2<tag))
	   {
	     if (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS))) {pTagProbeSS->Fill(L1,-0.5*m_weight);
	     }
	     if ((m_Mll>lowSS) && (m_Mll<upSS)) pTagProbeSS->Fill(L1,m_weight);     
	 } // end both tag

	 if ((L1<tag) && (L2>tag))
	 {   
	   if (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS))) pTagProbeSS->Fill(L2,-0.5*m_weight);
	   if ((m_Mll>lowSS) && (m_Mll<upSS))pTagProbeSS->Fill(L2,m_weight);
	 } // end L1 tag, L2 proble;

	 if ((L1>tag) && (L2<tag))
         {   
	   if (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS))) pTagProbeSS->Fill(L1,-0.5*m_weight);
	   if ((m_Mll>lowSS) && (m_Mll<upSS)) pTagProbeSS->Fill(L1,m_weight);
	 } // end L1 probe, L2 Tag

   } // end same sign.

  if (m_ChPl1l2<0){

         if ((L1<tag) && (L2<tag))
	 {   
	     if (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS))) pTagProbeOS->Fill(L1,-0.5*m_weight);
	     if ((m_Mll>lowOS) && (m_Mll<upOS)) pTagProbeOS->Fill(L1,m_weight);     
	 } // end both tag

	 if ((L1<tag) && (L2>tag))
	 {   
	     if (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS))) pTagProbeOS->Fill(L2,-0.5*m_weight);
	     if ((m_Mll>lowOS) && (m_Mll<upOS)) pTagProbeOS->Fill(L2,m_weight);
	 } // end L1 tag, L2 proble;

	 if ((L1>tag) && (L2<tag))
         {   
	     if (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS))) pTagProbeOS->Fill(L1,-0.5*m_weight);
	     if ((m_Mll>lowOS) && (m_Mll<upOS)) pTagProbeOS->Fill(L1,m_weight);
	 } // end L2 tag, L1 probe

  } // end opposite sign.

 } // end if tag

 return;


} // end function tag&probeGeneral 

void MisID::DirectExtraction(TH1D* &pDirectExSS,TH1D* &pDirectExOS,float L1,float L2)

{

int nbins=pDirectExSS->GetNbinsX();


int bin1=-1,bin2=-1;


for(int i=1; i<=nbins; i++){
    if((L1>pDirectExSS->GetBinLowEdge(i)) && (L1< pDirectExSS->GetBinLowEdge(i+1))) bin1 = i;
    if((L2>pDirectExSS->GetBinLowEdge(i)) && (L2< pDirectExSS->GetBinLowEdge(i+1))) bin2 = i;
 }



if(bin1==bin2){

  if (m_ChPl1l2>0){
    if (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS))) {pDirectExSS->Fill(L1,-0.5*m_weight);
    }
    if ((m_Mll>lowSS) && (m_Mll<upSS)) {pDirectExSS->Fill(L1,m_weight);
    }
  }
  
  if (m_ChPl1l2<0){
    if (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS))) {pDirectExOS->Fill(L1,-0.5*m_weight);
    }
    if ((m_Mll>lowOS) && (m_Mll<upOS)) {pDirectExOS->Fill(L1,m_weight);
    }
  }

}


 return;

} // end function Direct Extraction

void MisID::DirectExtraction2D(TH2D* &pDE2DSS,TH2D* &pDE2DOS,float L1,float L2)
{

  if (m_ChPl1l2>0){
    if (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS))) {pDE2DSS->Fill(L1,L2,-0.5*m_weight);
    }
    if ((m_Mll>lowSS) && (m_Mll<upSS)) {pDE2DSS->Fill(L1,L2,m_weight);
    }
  }
  
  if (m_ChPl1l2<0){
    if (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS))) {pDE2DOS->Fill(L1,L2,-0.5*m_weight);
    }
    if ((m_Mll>lowOS) && (m_Mll<upOS)) {pDE2DOS->Fill(L1,L2,m_weight);
    }
  }


  return;

}

void MisID::MatrixLikelihood(TH2D* &pLikelihoodSS,TH2D* &pLikelihoodOS,float L1,float L2)
{

 float low=-1,up=-1;

 low = L1;
 up = L2;


 if (low>up){
     low=L2;
     up=L1;
 }

  if (m_ChPl1l2>0){
    if (GSideBand==true && (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS)))) {pLikelihoodSS->Fill(low,up,-0.5*m_weight);
    }
    if ((m_Mll>lowSS) && (m_Mll<upSS)) {pLikelihoodSS->Fill(low,up,m_weight);
    }
  }
  
  if (m_ChPl1l2<0){
    if (GSideBand==true && (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS)))){ pLikelihoodOS->Fill(low,up,-0.5*m_weight);
  }
    if ((m_Mll>lowOS) && (m_Mll<upOS)){ pLikelihoodOS->Fill(low,up,m_weight);
    }
  }


  return;

}

void MisID:: FillVectorsLikelihoodEtaPt(float L1Eta,float L2Eta,float L1Pt,float L2Pt)
{

int Eta1=-1,Eta2=-1,Pt1=-1,Pt2=-1;

 for (int i=0;i<NBinsEta; i++)
 {
   if ((L1Eta>EtaBin[i]) && (L1Eta<=EtaBin[i+1]))  Eta1=i;
   if ((L2Eta>EtaBin[i]) && (L2Eta<=EtaBin[i+1]))  Eta2=i;
 }

for (int j=0;j<NBinsPt; j++)
 {
   if ((L1Pt>PtBin[j]) && (L1Pt<=PtBin[j+1]))  Pt1=j;
   if ((L2Pt>PtBin[j]) && (L2Pt<=PtBin[j+1]))  Pt2=j;
 }



int lowE=Eta1;
int upE=Eta2;
int lowPt=Pt1;
int upPt=Pt2;

if (lowE>upE){
  lowE=Eta2;
  upE=Eta1;
  lowPt=Pt2;
  upPt=Pt1;
}

 if ((Eta1==-1) || (Eta2==-1) || (Pt1==-1) || (Pt2==-1)) std::cout << "Something is wrong with this event: Check eta and pT" << std::endl;

 if (m_ChPl1l2>0){
   
    if ( GSideBand==true && (((m_Mll>llowSS) && (m_Mll< lowSS)) || ((m_Mll>upSS) && (m_Mll<uupSS)))){
         Vector4DSS[lowE][upE][lowPt][upPt]=Vector4DSS[lowE][upE][lowPt][upPt]-0.5*m_weight;
         pEtaSSsig->Fill(L1Eta,-0.5*m_weight);
	 pPtSSsig->Fill(L1Pt,-0.5*m_weight);
    } 
     if ((m_Mll>lowSS) && (m_Mll<upSS)){
         Vector4DSS[lowE][upE][lowPt][upPt]=Vector4DSS[lowE][upE][lowPt][upPt]+m_weight;
         pEtaSSsig->Fill(L1Eta,m_weight);
	 pPtSSsig->Fill(L1Pt,m_weight);
	 pEtaSSAll->Fill(L1Eta,m_weight);
         pPtSSAll->Fill(L1Pt,m_weight);
    }
 }

if (m_ChPl1l2<0){
 
   if ( GSideBand==true && (((m_Mll>llowOS) && (m_Mll< lowOS)) || ((m_Mll>upOS) && (m_Mll<uupOS)))){
       Vector4D[lowE][upE][lowPt][upPt]=Vector4D[lowE][upE][lowPt][upPt]-0.5*m_weight;
      
       pEtaOSsig->Fill(L1Eta,-0.5*m_weight);
       pPtOSsig->Fill(L1Pt,-0.5*m_weight);
    }
    if ((m_Mll>lowOS) && (m_Mll<upOS)){
      Vector4D[lowE][upE][lowPt][upPt]=Vector4D[lowE][upE][lowPt][upPt]+m_weight;
      pEtaOSsig->Fill(L1Eta,m_weight);
      pPtOSsig->Fill(L1Pt,m_weight);
      pEtaOSAll->Fill(L1Eta,m_weight);
      pPtOSAll->Fill(L1Pt,m_weight);
    }
 }

  return;
}
void MisID::TruthMatchingZ(TH1D* &TrueLeptons,TH1D* &BadCharge, float L1, float L2, float L1Charge, float L2Charge,string var,TH1D* &HistoBadDR, TH1D* &HistoAllDR)
{

	double peso=1000000;
	
	if ((m_ChPl1l2>0 && m_Mll>lowSS && m_Mll<upSS) || (m_ChPl1l2<0 && m_Mll>lowOS && m_Mll<upOS)) peso=m_weight;
	else peso=-0.5*m_weight;	
	
	if (m_L1isChFlip==1 ){ 
	  BadCharge->Fill(L1,peso);
	}        
	TrueLeptons->Fill(L1,peso);

	if (m_L2isChFlip==1 ){ 
	  BadCharge->Fill(L2,peso);
	}        
	TrueLeptons->Fill(L2,peso);

   

  return;


 } // end method TruthMatching


void MisID::TruthMatchingZsingle(TH1D* &TrueElectrons,TH1D* &BadElectrons, float L1, float L2, float L1Charge, float L2Charge,TH1D* &TruePositrons,TH1D* &BadPositrons)
{

	
  double peso=m_weight;
	
  if ((m_ChPl1l2>0 && m_Mll>lowSS && m_Mll<upSS) || (m_ChPl1l2<0 && m_Mll>lowOS && m_Mll<upOS)) peso=m_weight;
  else peso=-0.5*m_weight;
  
  if (L1Charge==-1){
    if (m_L1isChFlip==1){ 
      BadPositrons->Fill(L1,peso);
      TruePositrons->Fill(L1,peso); 
    }
    else {
      TrueElectrons->Fill(L1,peso);    
    }
  }
    
  if (L1Charge==1){
    if (m_L1isChFlip==1)
      {
	BadElectrons->Fill(L1,peso);
	TrueElectrons->Fill(L1,peso); 
      }
    else{
      TruePositrons->Fill(L1,peso);    
    }
  }
  
  if (L2Charge==-1){
    if (m_L2isChFlip==1 )
      {
	BadPositrons->Fill(L2,peso);
	TruePositrons->Fill(L2,peso);
      }
    else {
      TrueElectrons->Fill(L2,peso);    
    }
  }
    
  if (L2Charge==1){
    if (m_L2isChFlip==1)
      {
	BadElectrons->Fill(L2,peso);
	TrueElectrons->Fill(L2,peso);
      }
    else {
      TruePositrons->Fill(L2,peso);    
    }
  }
  
  return;
}

void MisID::TruthMatchingZsingle2D(TH2D* &TrueElectrons,TH2D* &BadElectrons, float L1Eta, float L2Eta,float L1Pt,float L2Pt, 
							float L1Charge, float L2Charge,TH2D* &TruePositrons,TH2D* &BadPositrons)
{
	
  double peso=m_weight;
  
  if ((m_ChPl1l2>0 && m_Mll>lowSS && m_Mll<upSS) || (m_ChPl1l2<0 && m_Mll>lowOS && m_Mll<upOS)) peso=m_weight;
  else peso=-0.5*m_weight;
  if (L1Charge==-1){				
    if (m_L1isChFlip==1 )
      {
	BadPositrons->Fill(L1Eta,L1Pt,peso);
	TruePositrons->Fill(L1Eta,L1Pt,peso); 
      }
    else {					
      TrueElectrons->Fill(L1Eta,L1Pt,peso);    
    }
  }
    
  if (L1Charge==1){				
    if (m_L1isChFlip==1)
      {
	BadElectrons->Fill(L1Eta,L1Pt,peso);
	TrueElectrons->Fill(L1Eta,L1Pt,peso); 
      }
    else{					
      TruePositrons->Fill(L1Eta,L1Pt,peso);    
    }
  }
  

  
  if (L2Charge==-1){				
    if (m_L2isChFlip==1 )
      {
	BadPositrons->Fill(L2Eta,L2Pt,peso);
	TruePositrons->Fill(L2Eta,L2Pt,peso);					
      }
    else {
      TrueElectrons->Fill(L2Eta,L2Pt,peso);    
    }
  }
  
  if (L2Charge==1){				
    if (m_L2isChFlip==1)
      {
	BadElectrons->Fill(L2Eta,L2Pt,peso);
	TrueElectrons->Fill(L2Eta,L2Pt,peso);			
      }
    else {
      TruePositrons->Fill(L2Eta,L2Pt,peso);    
    }
  }			
  

	return;
}

void MisID::TruthMatchingZ2DOld(TH2D* &TrueLeptons,TH2D* &BadCharge, float L1Eta, float L2Eta,float L1Pt,float L2Pt, float L1Charge, float L2Charge)
{

  double peso=m_weight;
	
  if ((m_ChPl1l2>0 && m_Mll>lowSS && m_Mll<upSS) || (m_ChPl1l2<0 && m_Mll>lowOS && m_Mll<upOS)) peso=m_weight;
  else peso=-0.5*m_weight;


    if (m_L1isChFlip==1) BadCharge->Fill(L1Eta,L1Pt,peso);
    TrueLeptons->Fill(L1Eta,L1Pt,peso);	
	
    if (m_L2isChFlip==1)BadCharge->Fill(L2Eta,L2Pt,peso);		
    TrueLeptons->Fill(L2Eta,L2Pt,peso);

  return;

 }
	
void MisID::TruthMatchingZ2D(TH2D* &TrueLeptons,TH2D* &BadCharge, float L1Eta, float L2Eta,float L1Pt,float L2Pt, float L1Charge, float L2Charge)
{

  double peso=m_weight;
	
  if ((m_ChPl1l2>0 && m_Mll>lowSS && m_Mll<upSS) || (m_ChPl1l2<0 && m_Mll>lowOS && m_Mll<upOS)) peso=m_weight;
  else peso=-0.5*m_weight;

  if (m_L1isChFlip==1 )BadCharge->Fill(L1Eta,L1Pt,peso);	      
  TrueLeptons->Fill(L1Eta,L1Pt,peso);
  
  if (m_L2isChFlip==1 )BadCharge->Fill(L2Eta,L2Pt,peso);
  TrueLeptons->Fill(L2Eta,L2Pt,peso);
  	
  return;

 } // end method TruthMatching2D

void MisID::RatesTagProbe(TH1D *pTagProbeSS,TH1D *pTagProbeOS, TH1D* &pRate,string dependence)
{


  int Crack=-1;
  
  if (dependence=="eta") Crack=CrackBin;

  int nbins=pTagProbeSS->GetNbinsX();
 
  double deltai=-1,Rate0=-1,Rate=-1;

  double ss=-1,os=-1;
  
  double error=-1;

  for (int i=1;i<=nbins;i++){
    
    ss = pTagProbeSS->GetBinContent(i);
    os = pTagProbeOS->GetBinContent(i);
    deltai=ss*1.0/(ss+os);

    if (i==1){
      Rate=deltai/2;
      Rate0=Rate;
    }
    if (i>1){
      Rate=deltai-Rate0;
    }
    
    error=sqrt((Rate*(1-Rate))/(ss+os));
    if (i!=Crack){
    pRate->SetBinContent(i,Rate);
    pRate->SetBinError(i,error);
    }
  }
 

  return;
}

void MisID::RatesDirectExtraction(TH1D *pDirectExSS,TH1D *pDirectExOS, TH1D* &pRate,string dependence)
{

int nbins=pDirectExSS->GetNbinsX();
  
int Crack=-1;
   
   if (dependence=="eta") Crack=CrackBin;
 
  double Rate=-1;

  double error=-1;

  double ss=-1,os=-1;


  for (int i=1;i<=nbins;i++){
    
    ss = pDirectExSS->GetBinContent(i);
    os = pDirectExOS->GetBinContent(i);

    Rate=ss*0.5/(ss+os);

    error=sqrt((Rate*(1-Rate))/(ss+os));
   
    if (i!=Crack){
    pRate->SetBinContent(i,Rate);
    pRate->SetBinError(i,error);
    }
  }
 
  
  return;

}


void MisID::RatesDE2D(TH2D *pSS,TH2D *pOS,TH2D* &pRates2D,string dependence){


  int nbins=pSS->GetNbinsX();
  double ss=0,os=0;
  double rate=0,error=0;

  int Crack=-1;
 
  if (dependence=="eta") Crack=CrackBin;
 
  for(int i=1; i<=nbins; i++){
    for(int j=1; j<=nbins; j++){
                
      ss=pSS->GetBinContent(i,j);
      os=pOS->GetBinContent(i,j);   
  
       rate = ss*0.5/(os+ss); 
       error=sqrt((rate*(1-rate))/(ss+os));

         
      
       if ((i!=Crack) && (j!=Crack) && (rate>=0) && (error>=0)){
	 pRates2D->SetBinContent(i,j,rate);
          pRates2D->SetBinError(i,j,error);
       }
       else {
         pRates2D->SetBinContent(i,j,0);
         pRates2D->SetBinError(i,j,0);
       }

    }
  }


  return;

}

void MisID::RatesLikelihood(TH2D *pLikelihoodSS,TH2D *pLikelihoodOS,TH1D* &pRate,string dependence)
{

int max=pLikelihoodSS->GetNbinsX(); 

if (dependence=="eta") max=max-1;

 std::cout <<"-------------------  Starting minimization in likelihood for '" << dependence << "' --------------------------------------------------" << std::endl;


   const int NPAR=max;
   Double_t vstart[NPAR];
   Double_t step[NPAR];
   TMinuit *pMinuit = new TMinuit(NPAR);

  Double_t arglist[2*NPAR];
  Int_t ierflg = 0;

  int Crack=NBinsPtTotal+1;

 if (dependence=="eta"){
   Crack=CrackBin;
   pMinuit->SetFCN(chi2Eta);
 }
 else {
   pMinuit->SetFCN(chi2);
 }  
    HistosforLikelihood[0]=pLikelihoodSS;
    HistosforLikelihood[1]=pLikelihoodOS;
    
   arglist[0] = 1;
   pMinuit->mnexcm("SET ERR",arglist,1,ierflg);
 
 
   for (int i=0;i<NPAR;++i){
      vstart[i]=0.0001;
      step[i]=1e-6;
      const TString name=TString::Format("epsilon%d",i+1);
     pMinuit->mnparm(i,name,vstart[i],step[i],0.00001,0.9999,ierflg);
   }


 
    // determination of parametres
    arglist[0] = 500;
    arglist[1] = 1.;
    pMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
 
    // print statistic
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    pMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
 
    if(icstat < 3) cout<<"There is  a problem with Minuit! "<<endl;
 
   std::cout << std::endl;
   std::cout << " Minimum fonction = " << amin <<std:: endl;
   std::cout << " Estimated vert. distance to min. = " << edm <<std:: endl;
   std::cout << " Number of variable parameters = " << nvpar << std::endl;
   std::cout << " Highest number of parameters defined by user = " << nparx << std::endl;
   std::cout << " Status of covariance matrix = " << icstat << std::endl;
   std::cout << std::endl;
 
  double pvalue, perror, plow,phigh;
   TString pname;
 
   int lu=0;
 
   if(amin < 99999) {
     for (int i=0;i<max;i++)
       {
         lu=i;
         if ((i>=Crack-1)) lu=i+1;
 
         pMinuit->mnpout(i,pname,pvalue,perror,plow,phigh,ierflg);
         pRate->SetBinContent(lu+1,pvalue);
         pRate->SetBinError(lu+1,perror);
 
       }
   }

  HistosforLikelihood.clear();

  delete pMinuit;

  return;

}


void chi2Eta(int &npar, double *gin, double &f, double *par, Int_t iflag)
{

  int max=HistosforLikelihood[0]->GetNbinsX();

   int l=0,m=0;

   float ss=0,All=0;
 
   double fonction = 0;

   for(int i=0; i<max-1; i++){
     for(int j=i; j<max-1; j++){
  
        l=i; m=j;
        
        if (i>=CrackBin-1) l=i+1;
  
        if (j>=CrackBin-1) m=j+1;
       
         ss=HistosforLikelihood[0]->GetBinContent(l+1,m+1);
         All=ss+HistosforLikelihood[1]->GetBinContent(l+1,m+1);
	 if (All > 0) fonction += log(All*(par[i]+par[j]))*ss - All*(par[i]+par[j]);
     }
   }

  
   f = -1*fonction;

   return;

}

void chi2(int &npar, double *gin, double &f, double *par, Int_t iflag)
{

  int max=HistosforLikelihood[0]->GetNbinsX();

 
   float ss=0,All=0;
 
   double fonction = 0;
   for(int i=0; i<max; i++){
     for(int j=i; j<max; j++){
  
         ss=HistosforLikelihood[0]->GetBinContent(i+1,j+1);
         All=ss+HistosforLikelihood[1]->GetBinContent(i+1,j+1);
	 if (All > 0) fonction += log(All*(par[i]+par[j]))*ss - All*(par[i]+par[j]);
     }
   }
  
   f = -1*fonction;

   return;

}


void MisID:: RatesTruthMatching(TH1D *TrueLeptons, TH1D *BadCharge,TH1D* &pRate,TString sample)
{


  double GlobalNum=0.;
  double GlobalDen=0.;

 int nb=pRate->GetNbinsX();

   double rate=-1.;
  
   double errorNum=-1,errorDen=-1,error=-1;

   for(int j=1; j<=nb; j++){
   
    double tr=TrueLeptons->GetBinContent(j);
    double bd=BadCharge->GetBinContent(j);
    
    GlobalNum=GlobalNum+bd;
    GlobalDen=GlobalDen+tr;
     errorNum=BadCharge->GetBinError(j);
     errorDen=TrueLeptons->GetBinError(j);
     rate=bd/tr;
     
     error=TEfficiency::ClopperPearson(static_cast<int>(tr+0.5),static_cast<int>(bd+0.5),0.68,true)-rate;
     if (rate>=0){
       pRate->SetBinContent(j,rate); 
       pRate->SetBinError(j,error);
     }
    }
    
   
   std::cout<< "-------------------------------- Global Mis-id for '" << sample << "' ------------------------" << std::endl;
   if (nb==NBinsEta)  std::cout<< "eta: Global Mis-id = " << GlobalNum/GlobalDen << std::endl;  
   else std::cout<< "pT: Global Mis-id = " << GlobalNum/GlobalDen << std::endl; 
  //******************************************************************

      return;

}

void MisID::RatesTruthMatching2D(TH2D *TrueLeptons, TH2D *BadCharge,TH2D* &pRate,TH1D** &ArrayTruthMatchEta,TH1D** &ArrayTruthMatchmedEta,TH1D** &ArrayAlphaEta,int bin,bool ttbar)
{

  std::cout<< "--------->> Computing mis-id rates using truth-matching in 2D " << std::endl;

    int nX=pRate->GetNbinsX();
    int nY=pRate->GetNbinsY();
  

      for(int j=1; j<=nY; j++){
         for (int i=1; i<=nX; i++){
          
	   

            double tr=TrueLeptons->GetBinContent(i,j);
            double bd=BadCharge->GetBinContent(i,j);
  
            double rate=bd/tr;

            double error=TEfficiency::ClopperPearson(static_cast<int>(tr+0.5),static_cast<int>(bd+0.5),0.68,true)-rate;

	    double numMed=0;  double denMed=0;

	    for (int k=1;k<=nY;k++){
                numMed=numMed+BadCharge->GetBinContent(i,k);
                denMed=denMed+TrueLeptons->GetBinContent(i,k);
            }

            double Med=numMed/denMed;
 
            double errorMed=TEfficiency::ClopperPearson(static_cast<int>(denMed+0.5),static_cast<int>(numMed+0.5),0.68,true)-Med;

	    double errorM=ComputeErrorCociente(rate,Med,error,errorMed);

            double rate0 = BadCharge->GetBinContent(i,bin)/TrueLeptons->GetBinContent(i,bin);
            double error0= TEfficiency::ClopperPearson(static_cast<int>(TrueLeptons->GetBinContent(i,bin)+0.5),static_cast<int>(BadCharge->GetBinContent(i,bin)+0.5),0.68,true)-rate0;

            double errorAlpha=ComputeErrorCociente(rate,rate0,error,error0);

           if (rate>0){
			   pRate->SetBinContent(i,j,rate);
			   pRate->SetBinError(i,j,error); 
              ArrayTruthMatchEta[j-1]->SetBinContent(i,rate);
              ArrayTruthMatchEta[j-1]->SetBinError(i,error);
	      ArrayTruthMatchmedEta[j-1]->SetBinContent(i,rate/Med);
	      ArrayTruthMatchmedEta[j-1]->SetBinError(i,errorM);
	      ArrayAlphaEta[j-1]->SetBinContent(i,rate/rate0);
	      ArrayAlphaEta[j-1]->SetBinError(i,errorAlpha);
	      
	     if (ttbar){
	        ArrayAlphaPt[i-1]->SetBinContent(j,rate/rate0);
                ArrayAlphaPt[i-1]->SetBinError(j,errorAlpha);
	     }

           }


	   

  
         } // end for eta bin
      } // end for pt bin
	
	
	for (int k=1; k<=NBinsPtTotal; k++){
		for (int j=1; j<=NBinsEta; j++){
			
			double alpha=0;
			double errorAl=0;
			alpha=ArrayAlphaEta[k-1]->GetBinContent(j);
			if (alpha==0 && j!=CrackBin){
				if (j==1){
					alpha=ArrayAlphaEta[k-1]->GetBinContent(2);
					errorAl=ArrayAlphaEta[k-1]->GetBinError(2);
				}	
				else if (j==NBinsEta){
					alpha=ArrayAlphaEta[k-1]->GetBinContent(NBinsEta-1);
					errorAl=ArrayAlphaEta[k-1]->GetBinError(NBinsEta-1);
				}	
				else {
					alpha=0.5*(ArrayAlphaEta[k-1]->GetBinContent(j-1)+ArrayAlphaEta[k-1]->GetBinContent(j+1));
				    errorAl= 0.5*(ArrayAlphaEta[k-1]->GetBinError(j-1)+ArrayAlphaEta[k-1]->GetBinError(j+1));
				}	
				
				ArrayAlphaEta[k-1]->SetBinContent(j,alpha);
				ArrayAlphaEta[k-1]->SetBinError(j,errorAl);
				   
			}	
			
		}
	}	



  return;

 }

void MisID::RatesElecPosi2D(TH2D *TrueLeptons,TH2D *BadCharge,TH2D* &pRate)
{

   int nX=pRate->GetNbinsX();
    int nY=pRate->GetNbinsY();
  
	
      for(int j=1; j<=nY; j++){
         for (int i=1; i<=nX; i++){
          
            double tr=TrueLeptons->GetBinContent(i,j);
            double bd=BadCharge->GetBinContent(i,j);
  
            double rate=bd/tr;

            double errorUp=TEfficiency::ClopperPearson(static_cast<int>(tr+0.5),static_cast<int>(bd+0.5),0.68,true)-rate;
			double errorDown=rate-TEfficiency::ClopperPearson(static_cast<int>(tr+0.5),static_cast<int>(bd+0.5),0.68,false);
			 
			 double error=errorUp;
			 if (errorUp<errorDown) error=errorDown;
           if (rate>=0){
              pRate->SetBinContent(i,j,rate);
              pRate->SetBinError(i,j,error);
           }

         } // end for eta bin
      } // end for pt bin


      TCanvas *pC = new TCanvas;
      pC->Divide(1,2);
      pC->cd(1);
      BadCharge->Draw("LEGO");
      pC->cd(2);
      TrueLeptons->Draw("LEGO");
      


return;
}
void MisID::RatesLikelihoodEtaPt(TH2D* &pRate,bool isMC)
{

m_statusMatrix=-99;
	
std::cout <<"----------------------- Starting minimization in likelihood 2D  -------------------------------------------------" << std::endl;

int max=NBinsEta-1; // 

   const int NPAR=max*NBinsPt;
   Double_t vstart[NPAR];
   Double_t step[NPAR];
   TMinuit *pMinuit = new TMinuit(NPAR);

  Double_t arglist[2*NPAR];
  Int_t ierflg = 1;

  int Crack=-99;
  Crack=CrackBin;

   pMinuit->SetFCN(chi2EtaPt);
 

   arglist[0] = 1;
   pMinuit->mnexcm("SET ERR",arglist,1,ierflg);


   vstart[0]=0.0005; // vstart[1]=0.001;
  
   for (int j=1; j<NPAR;++j){
     vstart[j]=vstart[j-1]+0.005;
     vstart[1]=0.001;
   }


   for (int i=0;i<NPAR;++i){
     
      step[i]=1e-6;
      const TString name=TString::Format("Ptepsilon%d",i+1);
      pMinuit->mnparm(i,name,vstart[i],step[i],0.000005,0.9,ierflg);
   }



    // determination of parametres
    arglist[0] = 500;
    arglist[1] = 1.;
    pMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
   
     // print statistic
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    pMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    if(icstat < 3) cout<<"There is  a problem with Minuit! "<<endl;

   std::cout << std::endl;
   std::cout << " Minimum fonction = " << amin <<std:: endl;
   std::cout << " Estimated vert. distance to min. = " << edm <<std:: endl;
   std::cout << " Number of variable parameters = " << nvpar << std::endl;
   std::cout << " Highest number of parameters defined by user = " << nparx << std::endl;
   std::cout << " Status of covariance matrix = " << icstat << std::endl;
   std::cout << std::endl;
	
	m_statusMatrix=icstat;	

  double pvalue, perror, plow,phigh;
   TString pname;

   int lu=0;



 for (int j=0; j<NBinsPt; j++) ArrayLikelihoodEta[j]= new TH1D(Form("LikelihoodEta%d",j),"Histogram;|#eta|;#epsilon_{mis-id}",NBinsEta,EtaBin);
  

  if(amin < 99999) {

    for (int i=0; i<max; i++){
      for (int j=0; j<NBinsPt; j++){
        
           lu=i;
            if ((i>=Crack-1)) lu=i+1;
            pMinuit->mnpout(NBinsPt*i+j,pname,pvalue,perror,plow,phigh,ierflg);

            ArrayLikelihoodEta[j]->SetBinContent(lu+1,pvalue);
	    ArrayLikelihoodEta[j]->SetBinError(lu+1,perror);

 
            pRate->SetBinContent(lu+1,j+1,pvalue);
            pRate->SetBinError(lu+1,j+1,perror);

	    cout<<"For bin "<<j+1<<" in eta and bin "<<lu+1<<" in pt we have as rates "<<pvalue<<endl;
      }
    }
  }


  TCanvas *pR = new TCanvas;
  pR->cd();
  pR->cd()->SetLogy();
  pR->cd()->SetLogz();
  pRate->GetYaxis()->SetTitle("p_{T} [GeV]");
  pRate->GetXaxis()->SetTitle("|#eta|");

  pRate->Draw("LEGO2");

   

  TLatex *pL= new TLatex;
  GetLatex(isMC,pL,"Z+jets",false); 

  delete pMinuit;

  return;

}

void chi2EtaPt(int &npar, double *gin, double &f, double *par, Int_t iflag)
{

  int maxEta=NBinsEta;
  int maxPt=NBinsPt;

   int l=0,m=0;

   float ss=0,All=0;

   int par1=-1,par2=-1;

   double fonction = 0;
   for(int i=0; i<maxEta-1; i++){
     for(int j=i; j<maxEta-1; j++){
        for (int pp=0;pp<maxPt; pp++){    
          for (int tt=0; tt<maxPt; tt++) {

                l=i; m=j;

        	if (i>=CrackBin-1)  l=i+1;
           
        	if (j>=CrackBin-1)  m=j+1;

          	ss=Vector4DSS[l][m][pp][tt];
          	All=ss+Vector4D[l][m][pp][tt];

                par1=maxPt*i+pp;
                par2=maxPt*j+tt;

 		if (ss!=0){                  
       		  fonction += log(All*(par[par1]+par[par2]))*ss - All*(par[par1]+par[par2]);
 		 }
 		else{
 		  fonction += - All*(par[par1]+par[par2]);
 		 }
     }
   }
  }
 }

   f=-1*fonction;


   return;

}


void MisID::PlotRates(TH1D *TP,TH1D *DE, TH1D *Likelihood,TH2D *DE2D,TH1D *pTM,bool MC,bool OtherMethods)
{

  Likelihood->GetYaxis()->SetRangeUser(0.00007, 1.11);
  Likelihood->GetXaxis()->SetLimits(0,2.55);
  TP->SetMarkerStyle(23);
  TP->SetLineColor(8);
  TP->SetMarkerColor(8);

  DE->SetMarkerStyle(22);
  DE->SetLineColor(2);
  DE->SetMarkerColor(2);

  pTM->SetMarkerColor(1);
  pTM->SetMarkerStyle(24);
  pTM->SetLineColor(1);

  Likelihood->SetMarkerStyle(22);
  Likelihood->SetMarkerColor(4);
  Likelihood->SetLineColor(4);

 TCanvas *pC = new TCanvas;
 
 TPad *pad1;

 if (!MC) {

   pad1 = new TPad("pad1","pad1",0,0,1,1);
 }
 else {

   Likelihood->GetYaxis()->SetTitleSize(0.08);
   Likelihood->GetYaxis()->SetTitleOffset(0.8);

   pad1 = new TPad("pad1","pad1",0,0.35,1,1);
   pad1->SetBottomMargin(0);
 }


   pad1->Draw();
   pad1->cd();
   pad1->cd()->SetLogy();
   Likelihood->Draw("P");

  TLegend *leg1=new TLegend(0.68,0.90,0.93,0.78);
  leg1->SetFillColor(0);
  leg1->AddEntry(Likelihood,"Likelihood","lp"); 


  if (OtherMethods){
      TP->Draw("p,same");
      DE->Draw("p,same");
      leg1->AddEntry(TP,"Tag&Probe","lp");
      leg1->AddEntry(DE,"Direct Extraction","lp");
   }

   TLatex *pL = new TLatex;
   TString sample="";

   if (MC){      
      pTM->Draw("p,same");
      leg1->AddEntry(pTM,"Truth-matching");  
      sample="Z+jets";
    }
 
   bool mod=false;
   if (MC==true) mod=true; 

   GetLatex(MC,pL,sample,mod); 

 leg1->Draw("same");       
   
 pC->cd();

 if (MC){
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->Draw();
   pad2->cd();
    
   TH1D* pRatio = (TH1D*)Likelihood->Clone();
   
   pRatio->GetYaxis()->SetTitle("Likelihood/Truth-matching");
   pRatio->Divide(pTM);
   pRatio->GetYaxis()->SetRangeUser(0.5,1.45);
   pRatio->GetXaxis()->SetTitleSize(0.11);
   pRatio->GetXaxis()->SetTitleOffset(0.8);
   pRatio->GetXaxis()->SetLabelSize(0.07);
   pRatio->GetYaxis()->SetLabelSize(0.07);
   pRatio->Draw("ep");

	 TLine *pli= new TLine;
	 pli->SetLineWidth(2);
	 pli->DrawLine(0,1,pRatio->GetBinLowEdge(pRatio->GetNbinsX()+1),1);
	 
   pC->cd();

   TCanvas *pD = new TCanvas;
   pTM->GetYaxis()->SetRangeUser(0.00007, 1.11);
   pD->cd();
   pD->cd()->SetLogy();
   pTM->Draw();
   

   TLatex *pll = new TLatex;
   GetLatex(MC,pll,sample,false); 


 }


 if (OtherMethods){
   TCanvas *pp = new TCanvas;
   pp->cd();
   pp->cd()->SetLogz();
   DE2D->GetZaxis()->SetRangeUser(0.0001,1.01);
   DE2D->Draw("LEGO2");
   TLatex *pLat4=new TLatex();
   AddName("Direct extraction 2D",pLat4,0.21,0.75);
 }


  return;

} 
void MisID::ComparElecPosi2D(bool isMC)
{

  vector <string> full1(NBinsPtTotal); 

 
  full1=ReturnLegendPt();

  TH2D* pRatio = (TH2D*)pRatesElec2D->Clone();
  pRatio->Divide(pRatesPosi2D);
  
	TCanvas *pR = new TCanvas;
	pR->Divide(1,2);
	pR->cd(1);
	pRatesElec2D->Draw("LEGO");
	pR->cd(2);
	pRatesPosi2D->Draw("LEGO");
	
	
	TCanvas *pC = new TCanvas;
	pC->cd();
	
	TLegend *leg1=new TLegend(0.71,0.19,0.91,0.34);
    leg1->SetFillColor(0);		
	
	ArrayElecPosi[0]->GetYaxis()->SetRangeUser(0,4);
	
  for (int i=1;i<=pRatio->GetNbinsY();i++){
   for (int j=1;j<=pRatio->GetNbinsX();j++){
     ArrayElecPosi[i-1]->SetBinContent(j,pRatio->GetBinContent(j,i));
	 ArrayElecPosi[i-1]->SetBinError(j,pRatio->GetBinError(j,i));   
   }

	  SetHistoStyle(ArrayElecPosi[i-1],19+i,i,i);
	  leg1->AddEntry(ArrayElecPosi[i-1],full1[i-1].c_str(),"lp");
	  if (i==1) ArrayElecPosi[i-1]->Draw();
	  else  ArrayElecPosi[i-1]->Draw("same");
  
  }
	
	TLine *pli = new TLine;
	pli->SetLineWidth(2);
	pli->DrawLine(0,1,EtaBin[NBinsEta],1);

	leg1->Draw("same");

	TLatex *pL= new TLatex;
	GetLatex(isMC,pL,"Z+jets",false);
	
	
	
return;
}

void MisID::PlotRates2D(bool isMC,bool morePlots)
{

     vector <string> full1(NBinsPtTotal);

     full1=ReturnLegendPt();
 
    ArrayLikelihoodEta[0]->GetYaxis()->SetRangeUser(0.0001,1.01);
    ArrayTruthMatchZmedEta[0]->GetYaxis()->SetRangeUser(0,5);
    ArrayAlphaZEta[0]->GetYaxis()->SetRangeUser(0,5);

    TLegend *leg1=new TLegend(0.71,0.19,0.91,0.34);
    leg1->SetHeader("Likelihood");    
    leg1->SetFillColor(0);

    TLegend *leg2=new TLegend(0.71,0.34,0.91,0.49);
    leg2->SetHeader("Truth-matching");    
    leg2->SetFillColor(0);


    TCanvas *pp = new TCanvas;

    TPad *pad1;

    if (isMC==false) pad1 = new TPad("pad1","pad1",0,0,1,1);

    else {
      ArrayLikelihoodEta[0]->GetYaxis()->SetTitleSize(0.08);
      ArrayLikelihoodEta[0]->GetYaxis()->SetTitleOffset(0.8);

      pad1 = new TPad("pad1","pad1",0,0.35,1,1);
      pad1->SetBottomMargin(0);
    }

  
    pad1->Draw();
    pad1->cd();
    pad1->cd()->SetLogy();


  for (int j=0;j<NBinsPt; j++){
      std::stringstream bin1, bin2;
      SetHistoStyle(ArrayLikelihoodEta[j],20+j,j+1,1);
      SetHistoStyle(ArrayTruthMatchZEta[j],24+j,j+1,2);     

      bin1 << PtBin[j];
      bin2 << PtBin[j+1];      
      string full="p_{T} #in ["+bin1.str()+","+bin2.str()+"] GeV";

      if (j==0) ArrayLikelihoodEta[j]->DrawCopy();
      else ArrayLikelihoodEta[j]->DrawCopy("same");

      if (isMC) ArrayTruthMatchZEta[j]->Draw("same");
      leg1->AddEntry(ArrayLikelihoodEta[j],full.c_str(),"lp");  
      leg2->AddEntry(ArrayTruthMatchZEta[j],full.c_str(),"lp");  

    }

  TLatex *pLt = new TLatex;

  float mod=false;
  if (isMC==true) mod=true;
 
  GetLatex(isMC,pLt,"Z+jets",mod); 

  if (isMC==true) leg2->Draw("same");
    
  leg1->Draw("same");
   

  pp->cd();

  if (isMC){
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->Draw();
   pad2->cd();

   TH1D **pRatio = new TH1D*[NBinsPt];

   for (int i=0; i<NBinsPt; i++)
   {
     pRatio[i]=(TH1D*)ArrayLikelihoodEta[i]->Clone();
     pRatio[i]->Divide(ArrayTruthMatchZEta[i]);

     if (i==0){
         pRatio[i]->GetYaxis()->SetTitle("Likelihood/Truth-matching");
         pRatio[i]->GetYaxis()->SetRangeUser(0,3);
         pRatio[i]->GetYaxis()->SetTitleSize(0.08);
	 pRatio[i]->GetXaxis()->SetTitleSize(0.08);
         pRatio[i]->GetXaxis()->SetTitleOffset(0.8);
         pRatio[i]->GetXaxis()->SetLabelSize(0.07);
         pRatio[i]->GetYaxis()->SetLabelSize(0.07);
         pRatio[i]->Draw("ep");
     }
     else pRatio[i]->Draw("ep,same");

   }


   pp->cd();


   TCanvas *pp1 = new TCanvas;
   pp1->cd();

   TLegend *legd=new TLegend(0.63,0.90,0.93,0.75);
   legd->SetFillColor(0);

    for (int i=0; i<NBinsPt; i++)
   {
    
     if (i==0){
         pRatio[i]->GetYaxis()->SetTitle("Likelihood/Truth-matching");
         pRatio[i]->GetYaxis()->SetTitleSize(0.05);
	 pRatio[i]->GetXaxis()->SetTitleSize(0.05);
         pRatio[i]->GetXaxis()->SetTitleOffset(1.4);
         pRatio[i]->GetYaxis()->SetTitleOffset(1.4);
         pRatio[i]->GetXaxis()->SetLabelSize(0.05);
         pRatio[i]->GetYaxis()->SetLabelSize(0.05);
         pRatio[i]->Draw("ep");
     }
     else pRatio[i]->Draw("ep,same");
      
     legd->AddEntry(pRatio[i],full1[i].c_str(),"lp");

     legd->Draw("same");

     TLatex *pLt1 = new TLatex;
 
     GetLatex(isMC,pLt1,"Z+jets",false); 

   }



  }


 if (morePlots && isMC){

     TLegend *leg3=new TLegend(0.63,0.90,0.93,0.7);
     leg3->SetFillColor(0);
     for (int j=0;j<NBinsPt; j++){
      std::stringstream bin1, bin2;
      SetHistoStyle(ArrayTruthMatchZmedEta[j],20+j,j+1,1);
      SetHistoStyle(ArrayAlphaZEta[j],20+j,j+1,1);     

      bin1 << PtBin[j];
      bin2 << PtBin[j+1];      
      string full="p_{T} #in ["+bin1.str()+","+bin2.str()+"] GeV";

      leg3->AddEntry(ArrayTruthMatchZmedEta[j],full.c_str(),"lp");  
     }



      TCanvas *pC = new TCanvas;
      pC->cd();
      
      for (int j=0;j<NBinsPt; j++){
         if (j==0) ArrayTruthMatchZmedEta[j]->Draw();
         else ArrayTruthMatchZmedEta[j]->Draw("same");
      }
     
      leg3->Draw("same");

       TLatex *pt1 = new TLatex;
       GetLatex(true,pt1,"Z+jets",false);


      TCanvas *pCC = new TCanvas;
      pCC->cd();
      
      for (int j=0;j<NBinsPt; j++){
         if (j==0) ArrayAlphaZEta[j]->Draw();
         else ArrayAlphaZEta[j]->Draw("same");
      }
      leg3->Draw("same");
  
       TLatex *pt2 = new TLatex;
       GetLatex(true,pt2,"Z+jets",false);

  }


  return;

}

void MisID::ShowMoreOtherMethods(TH1D *pTagSS,TH1D *pTagOS,TH1D *pDESS,TH1D *pDEOS,TH2D *pDE2DSS,TH2D *pDE2DOS)
{


     TCanvas *pC1 = new TCanvas;

     TLegend *leg1=new TLegend(0.66,0.92,0.93,0.72);
     leg1->SetFillColor(0);
     leg1->AddEntry(pTagOS,"Opposite-sign","lp");
     leg1->AddEntry(pTagSS,"Same-sign","lp");

     pC1->Divide(2,2);
     pC1->cd(1);
     pC1->cd(1)->SetLogy();
     pTagOS->GetYaxis()->SetRangeUser(1,100000000);
     pTagOS->Draw();
     pTagSS->SetLineColor(2);
     pTagSS->Draw("same");
    
     TLatex *pLat2=new TLatex();
     AddName("Tag&Probe",pLat2,0.23,0.85);
     leg1->Draw("same");

     pC1->cd(2);
     pC1->cd(2)->SetLogy();
     pDEOS->GetYaxis()->SetRangeUser(1,100000000);
     pDEOS->Draw();
     pDESS->SetLineColor(2);
     pDESS->Draw("same");
     TLatex *pLat3=new TLatex();
     AddName("Direct extraction",pLat3,0.23,0.85);
     leg1->Draw("same");

     pC1->cd(3);
     pDE2DSS->Draw("LEGO2");

     TLatex *pLat4=new TLatex();
     AddName("Direct extraction 2D SS",pLat4,0.21,0.75);

     pC1->cd(4);
     pDE2DOS->Draw("LEGO2");
      TLatex *pLat5 =new TLatex();
      AddName("Direct extraction 2D OS",pLat5,0.21,0.75);

      return;

}

void MisID::ShowMoreDefaultMethods(bool isMC,TH2D *pLikelihoodSS,TH2D *pLikelihoodOS,TH1D *pTrueLeptons,TH1D *pBadCharge,
                                   TH1D *pTrueElectrons,TH1D *pBadElectrons,TH1D *pTruePositrons,TH1D *pBadPositrons){

   TCanvas *pC3 = new TCanvas;
   if (!isMC){
      pC3->Divide(1,2);
   }
   else pC3->Divide(2,3);

   pC3->cd(1);
  
   if (pLikelihoodSS->GetNbinsX()==NBinsPt) pC3->cd(1)->SetLogz();
   pLikelihoodSS->Draw("LEGO2");  
   
   TLatex *pL1=new TLatex();
   AddName("Likelihood",pL1,0.23,0.85);  

   pC3->cd(2);
   if (pLikelihoodOS->GetNbinsX()==NBinsPt) pC3->cd(2)->SetLogz();
   pLikelihoodOS->Draw("LEGO2"); 
   TLatex *pL2=new TLatex();
   AddName("Likelihood",pL2,0.23,0.85);  

   if (isMC){

     TLegend *leg2=new TLegend(0.66,0.92,0.93,0.72);

     leg2->SetFillColor(0);
     leg2->AddEntry(pTrueLeptons,"True leptons","lp");
     leg2->AddEntry(pBadCharge,"Mis-id leptons","lp");

     pC3->cd(5);
     pC3->cd(5)->SetLogy();
     pTrueLeptons->GetYaxis()->SetRangeUser(1,100000000);
     pTrueLeptons->Draw();
     pBadCharge->SetLineColor(2);
     pBadCharge->Draw("SAME");
     leg2->Draw("same");
     TLatex *pL3=new TLatex();
     AddName("Truth-matching: leptons",pL3,0.23,0.85); 

     TLegend *leg3=new TLegend(0.66,0.92,0.93,0.72);
     leg3->SetFillColor(0);
     leg3->AddEntry(pTrueElectrons,"True electrons","lp");
     leg3->AddEntry(pBadElectrons,"Mis-id electrons","lp");

     pC3->cd(3);
     pC3->cd(3)->SetLogy();
     pTrueElectrons->GetYaxis()->SetRangeUser(1,100000000);
     pTrueElectrons->Draw();
     pBadElectrons->SetLineColor(2);
     pBadElectrons->Draw("SAME");
     leg3->Draw("same");
     TLatex *pL4=new TLatex();
     AddName("Truth-matching: electrons",pL4,0.23,0.85); 

     TLegend *leg4=new TLegend(0.66,0.92,0.93,0.72);
     leg4->SetFillColor(0);
     leg4->AddEntry(pTruePositrons,"True positrons","lp");
     leg4->AddEntry(pBadPositrons,"Mis-id positrons","lp");

     pC3->cd(4);
     pC3->cd(4)->SetLogy();
     pTruePositrons->GetYaxis()->SetRangeUser(1,100000000);
     pTruePositrons->Draw();
     pBadPositrons->SetLineColor(2);
     pBadPositrons->Draw("SAME");
     leg4->Draw("same");
     TLatex *pL5=new TLatex();
     AddName("Truth-matching: positrons",pL5,0.23,0.85); 

    TLegend *leg5=new TLegend(0.66,0.92,0.93,0.72);
     leg5->SetFillColor(0);
     leg5->AddEntry(pTruePositrons,"True positrons","lp");
     leg5->AddEntry(pBadPositrons,"Mis-id positrons","lp");
     leg5->AddEntry(pTrueElectrons,"True electrons","lp");
     leg5->AddEntry(pBadElectrons,"Mis-id electrons","lp");

     pC3->cd(6);
     pC3->cd(6)->SetLogy();
     pTruePositrons->SetLineStyle(2);
     TH1D* pTruePositronsClone = (TH1D*)pTruePositrons->Clone();
     pTruePositronsClone->GetYaxis()->SetTitle("Number of leptons");
     pTruePositronsClone->Draw();
     pBadPositrons->SetLineStyle(2);
     pBadPositrons->SetLineColor(2);
     pTrueElectrons->Draw("same");
     pBadElectrons->Draw("same");
     pBadPositrons->Draw("SAME");
     leg5->Draw("same");
     TLatex *pL6=new TLatex();
     AddName("Truth-matching: positrons+electrons",pL6,0.23,0.85); 



   }


  return;
}

void MisID::ShowMorePlots(bool PlotEta,bool PlotPt,bool OtherMethods,bool isMC){


 if (isMC==true){

   TCanvas *pC = new TCanvas;
   pC->cd();

  TH1F* pRatiodR = (TH1F*)hBadDR->Clone();
     pRatiodR->GetYaxis()->SetTitle("#epsilon_{mis-id}");
     pRatiodR->Divide(hAllDR);
     pRatiodR->Draw();
 
   TLatex *pL = new TLatex;
   GetLatex(true,pL,"Z+jets",false);
 }


   if (OtherMethods){

     if (PlotEta) ShowMoreOtherMethods(pTagProbeSSEta,pTagProbeOSEta,pDirectExSSEta,pDirectExOSEta,pDE2DSSEta,pDE2DOSEta);
     if (PlotPt)  ShowMoreOtherMethods(pTagProbeSSPt,pTagProbeOSPt,pDirectExSSPt,pDirectExOSPt,pDE2DSSPt,pDE2DOSPt);
   }

   if (PlotEta) ShowMoreDefaultMethods(isMC,pLikelihoodSSEta,pLikelihoodOSEta,pTrueLeptonsEta,pBadChargeEta,pTrueElectronsEta,pBadElectronsEta,pTruePositronsEta,pBadPositronsEta);
   if (PlotPt) ShowMoreDefaultMethods(isMC,pLikelihoodSSPt,pLikelihoodOSPt,pTrueLeptonsPt,pBadChargePt,pTrueElectronsPt,pBadElectronsPt,pTruePositronsPt,pBadPositronsPt);


   if (isMC) {
    
      TCanvas *pp = new TCanvas;
    pp->cd();
    pp->cd()->SetLogy();
    ArrayTruthMatchZEta[0]->GetYaxis()->SetRangeUser(0.0001,1.01);
   
    TLegend *leg2=new TLegend(0.63,0.38,0.9,0.18);
    leg2->SetFillColor(0);
    leg2->SetHeader("Truth-matching");    


     for (int j=0;j<NBinsPtTotal; j++){
      std::stringstream bin1, bin2;
      
      SetHistoStyle(ArrayTruthMatchZEta[j],24+j,j+1,2);     

      bin1 << PtBin[j];
      bin2 << PtBin[j+1];      
      string full="p_{T} #in ["+bin1.str()+","+bin2.str()+"] GeV";

      if (j==0) ArrayTruthMatchZEta[j]->Draw();
      else ArrayTruthMatchZEta[j]->Draw("same");
     
      leg2->AddEntry(ArrayTruthMatchZEta[j],full.c_str(),"lp");  
     }

    leg2->Draw("same");
  
     TLatex *pL1 = new TLatex;
     GetLatex(true,pL1,"Z+jets",false);
          
   }
  

   TCanvas *pep = new TCanvas;

   pep->Divide(2,2);
   pep->cd(1);

  
   TH1D* pEtaSSbkg = (TH1D*)pEtaSSAll->Clone();

    pEtaSSbkg->Add(pEtaSSsig,-1);

    pEtaSSAll->SetFillColor(4);
    pEtaSSbkg->SetFillColor(2);
    pEtaSSAll->SetLineColor(4);
    pEtaSSbkg->SetLineColor(2);
      

    pEtaSSAll->Draw();
    pEtaSSbkg->Draw("same");

   TLegend *leg1=new TLegend(0.67,0.7,0.92,0.90);
   leg1->SetFillColor(0); 
   
   leg1->AddEntry(pEtaSSAll,"Signal","lp");
   leg1->AddEntry(pEtaSSbkg,"Background","lp");

   pEtaSSAll->Draw("sameaxis");
   leg1->Draw("same");  
 
    pep->cd(3);
    TH1D* pEtaOSbkg = (TH1D*)pEtaOSAll->Clone();

    pEtaOSbkg->Add(pEtaOSsig,-1);

    pEtaOSAll->SetFillColor(4);
    pEtaOSbkg->SetFillColor(2);
    pEtaOSAll->SetLineColor(4);
    pEtaOSbkg->SetLineColor(2);
   
    pEtaOSAll->Draw();
    pEtaOSbkg->Draw("same");
    pEtaOSAll->Draw("sameaxis");

    leg1->Draw("same");

   pep->cd(2);


   TH1D* pPtSSbkg = (TH1D*)pPtSSAll->Clone();

    pPtSSbkg->Add(pPtSSsig,-1);

    pPtSSAll->SetFillColor(4);
    pPtSSbkg->SetFillColor(2);
    pPtSSAll->SetLineColor(4);
    pPtSSbkg->SetLineColor(2);
   
    pPtSSAll->Draw();
    pPtSSbkg->Draw("same");
    pPtSSAll->Draw("sameaxis");

   leg1->Draw("same");
 
    pep->cd(4);
    TH1D* pPtOSbkg = (TH1D*)pPtOSAll->Clone();

    pPtOSbkg->Add(pPtOSsig,-1);

    pPtOSAll->SetFillColor(4);
    pPtOSbkg->SetFillColor(2);
    pPtOSAll->SetLineColor(4);
    pPtOSbkg->SetLineColor(2);
   
    pPtOSAll->Draw();
    pPtOSbkg->Draw("same");
    pPtOSAll->Draw("sameaxis");
   leg1->Draw("same");



  return;

}

void AddName(TString name,TLatex* &pL,float sep,float alt)
{
  pL->SetTextSize(0.03);

  pL->SetNDC();
  pL->DrawLatex(sep,alt,name);


  return;

}
void MisID::PlotTruthMatching(TH1D *pRateTruthMatch,TH1D *pRateElecTruthMatch,TH1D *pRatePosiTruthMatch,bool isTTbar){

  TCanvas *pC = new TCanvas;
  pC->cd();
  pC->cd()->SetLogy();
  pRateElecTruthMatch->GetYaxis()->SetRangeUser(0.00009, 1.11);

  SetHistoStyle(pRateElecTruthMatch,22,3,2);
  SetHistoStyle(pRatePosiTruthMatch,23,2,2);
  SetHistoStyle(pRateTruthMatch,24,1,1);

  pRateElecTruthMatch->Draw("P");
  pRatePosiTruthMatch->Draw("P,same");
  pRateTruthMatch->Draw("P,same");

  TLegend *leg1=new TLegend(0.7,0.2,0.95,0.40);
  leg1->SetFillColor(0); 

  leg1->SetHeader("Truth-matching"); 
  leg1->AddEntry(pRateElecTruthMatch,"Electrons","lp"); 
  leg1->AddEntry(pRatePosiTruthMatch,"Positrons","lp"); 
  leg1->AddEntry(pRateTruthMatch,"Leptons","lp"); 
  leg1->Draw("same");

  TString sample="Z+jets";
   if (isTTbar==true) sample="t#bar{t}";
  
   TLatex *pL = new TLatex;
   GetLatex(true,pL,sample,false);


   TCanvas *pC1= new TCanvas;
   pC1->cd();

   TH1D* pRatioEP = (TH1D*)pRateElecTruthMatch->Clone();
     pRatioEP->GetYaxis()->SetTitle("#epsilon_{mis-id,e^{-}}/#epsilon_{mis-id,e^{+}}");
     pRatioEP->Divide(pRatePosiTruthMatch);
     pRatioEP->GetYaxis()->SetRangeUser(0,2);
     pRatioEP->Draw(); 

   TLatex *pL1 = new TLatex;
   GetLatex(true,pL1,sample,false);

   TLine *pP = new TLine;
   pP->SetLineWidth(2);
   pP->DrawLine(pRateElecTruthMatch->GetBinLowEdge(1),1,pRateElecTruthMatch->GetBinLowEdge(pRateElecTruthMatch->GetNbinsX()+1),1);
   


  return;
}

void SetHistoStyle(TH1D* &pH,int Marker,int color,int line)
{
  if (color==4) color=8;
  if (color==3) color=4;
  if (color==5) color=kViolet;

  pH->SetMarkerColor(color);
  pH->SetLineColor(color);
  pH->SetMarkerStyle(Marker);
  pH->SetLineStyle(line);

  return;
}

double MisID:: ComputeErrorCociente(double Num,double Den,double errorN,double errorD){

  double error=TMath::Sqrt(pow(errorN/Den,2)+pow(errorD*Num/pow(Den,2),2));

  return error;
}

void MisID:: ShowEventsMatrix()
{

 int max=HistosforLikelihood[0]->GetNbinsX();

   int l=0,m=0;

  
 
   for(int i=0; i<max-1; i++){
     for(int j=i; j<max-1; j++){
  
        double ss=0,All=0;

           l=i; m=j;
        if (i>=CrackBin-1){
            l=i+1;
            m=j+1;
        }

        if (j>=CrackBin-1){
            m=j+1;
       }
 
         ss=HistosforLikelihood[0]->GetBinContent(l+1,m+1);
         All=ss+HistosforLikelihood[1]->GetBinContent(l+1,m+1);
	 std::cout<< "(i,j)=(" << l+1 << "," << m+1 << ")--> ss=" << ss << ", All=" << All << std::endl;       

     }
   }

  int maxEta=NBinsEta;
  int ll=0,mm=0;

  

   for(int i=0; i<maxEta-1; i++){
     for(int j=i; j<maxEta-1; j++){
       
                double sss=0,Alll=0;

                ll=i; mm=j;
        	if (i>=CrackBin-1){
            	   ll=i+1;
            	   mm=j+1;
        	}

        	if (j>=CrackBin-1){
            	   mm=j+1;
       		}

         	sss=Vector4DSS[ll][mm][0][0]+Vector4DSS[ll][mm][0][1]+Vector4DSS[ll][mm][1][0]+Vector4DSS[ll][mm][1][1];
         	Alll=sss+Vector4D[ll][mm][0][0]+Vector4D[ll][mm][0][1]+Vector4D[ll][mm][1][0]+Vector4D[ll][mm][1][1];

		 std::cout<< "(i,j)=(" << ll+1 << "," << mm+1 << ")--> sss=" << sss << ", Alll=" << Alll << std::endl;  
		 std::cout<< "cero=" << Vector4DSS[ll][mm][1][0] << ", zero=" << Vector4D[ll][mm][1][0] << std::endl;
  }
 }
   return;


}


void GetLatex(bool MC,TLatex* &pLat1, TString sample, bool mod)
{
  float f=1.0;  
  float m=1.0;


  if (mod==true){
     f=1.4;
     m=1.3;
  }
  

    pLat1->SetNDC();
    pLat1->SetTextColor(1);
    pLat1->SetTextSize(f*0.05);
    pLat1->SetTextFont(72);

    float sep=0.178;
    float alt=0.85;
    float alts=m*0.07;
   pLat1->DrawLatex(sep,alt,"ATLAS");

      pLat1->SetTextFont(42);
      pLat1->SetTextSize(f*0.04);
     
       TString label="#font[42]{Internal}";
       float deltaX=0.115*696*gPad->GetWh()/(472*gPad->GetWw());
       pLat1->DrawLatex(sep+deltaX,alt,label);
       if (!MC) pLat1->DrawLatex(sep,alt-alts, "#intLdt =5.4  fb^{-1}");
       
       else  pLat1->DrawLatex(sep,alt-alts, "Simulation "+sample);      
    

       pLat1->DrawLatex(sep,alt-alts-alts, "#sqrt{s} = 8 TeV");

       return;

}

void ModifyLabels(TH1D* &pH,float offset,float size)
{

  pH->GetYaxis()->SetTitleOffset(offset);
  pH->GetYaxis()->SetTitleSize(size);

  pH->GetXaxis()->SetTitleOffset(offset);
  pH->GetXaxis()->SetTitleSize(size);


  return;
}


vector<string> ReturnLegendPt()
{

 
  vector <string> Legend(NBinsPtTotal);

    for (int j=0;j<NBinsPtTotal; j++){
      std::stringstream bin1, bin2;
       
      bin1 << PtBin[j];
      bin2 << PtBin[j+1];      
      string full="p_{T} #in ["+bin1.str()+","+bin2.str()+"] GeV";

      Legend[j]=full;
  }

  return Legend;
}
