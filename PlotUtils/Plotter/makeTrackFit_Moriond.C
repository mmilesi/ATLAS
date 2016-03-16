//==================================================================
//
//
//
//==================================================================
#include <TROOT.h>
#include <THStack.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TCanvas.h>
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TMatrixD.h"
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include "TEventList.h"
#include "TLorentzVector.h"
#include <vector>
#include "TVector2.h"
#include "TVector3.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TStyle.h>
#include <stdlib.h>
#include <fstream>
#include "TH1F.h"
#include "TEventList.h"
#include "TObjArray.h"
#include "TLine.h"
#include <iomanip>
#include "TSystem.h"
#include "TMinuit.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TRandom.h"
#include "TMath.h"
#include "utils.h"


int NDATA;
double DATA[10000];

double sumTAU; double dTAU[64];
double sumQCD; double dQCD[64];
double sumHIG; double dHIG[64];
double sumELW; double dELW[64];
double sumELW_SS; double dELW_SS[64];
double sumDATA; double trackDATA[64];
double sumQCDstat; double statQCD[64];
double sumTAUstat; double statTAU[64];
double sumQCD2; double dQCD2[64];
double sumTAU_pileuplow;  double dTAU_pileuplow[64];
double sumTAU_pileuphigh; double dTAU_pileuphigh[64];

double shapeUncertainty = 0.05;
double signalAcceptanceUncertainty = 0.10;

const int iRandomSeed = 65539;   // Seed for flat random number.
int category;

//=========================
// Log Likelihood Function
//=========================

float stepFunction(float x){
    if(x>0.) return 1.;
    else return 0.;
    
}

void trackFit(Int_t &npar, Double_t *gin, Double_t &f,Double_t *par, Int_t iflag) {
    
    double xqcd = par[0];
    double xelw = par[1]*sumELW/sumDATA;
    double xtau = 1.0-xqcd-xelw;
    
    
    f = 1.0e30;
    if (xtau<0.0) { f = 1.0e30; return; }
    if (xqcd<0.0) { f = 1.0e30; return; }
    if (xtau>1.0) { f = 1.0e30; return; }
    if (xqcd>1.0) { f = 1.0e30; return; }
    
    double xshape = 0.0;
    
    for (int i = 0; i < 64; i++) {
        if (trackDATA[i]<1) { continue; }
        if (dTAU[i]+dQCD[i]+dELW[i]==0.0) { continue; }
        double wtau = sumDATA*xtau*dTAU[i];
        double wqcd = sumDATA*xqcd*dQCD[i];
        double welw = sumDATA*xelw*dELW[i];
        
		// QCD shape syst:
        wqcd += sumDATA*xqcd*par[2]*(dQCD2[i]-dQCD[i]);
		// Ztt shape syst:
        wtau += sumDATA*xtau*(  par[3]*(dTAU_pileuphigh[i]-dTAU[i])*stepFunction(par[3]) +   par[3]*(dTAU[i]-dTAU_pileuplow[i])*stepFunction(-1.*par[3]) );
        //wtau += sumDATA*xtau*par[4]*(dTAU_pileuphigh[i]-dTAU[i]);
        
        if (i==9) {
            if (statQCD[i]>0) {
                double xqcdexp = par[5];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==10) {
            if (statQCD[i]>0) {
                double xqcdexp = par[6];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }  
        }
        if (i==11) {
            if (statQCD[i]>0) {
                double xqcdexp = par[7];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==17) {
            if (statQCD[i]>0) {
                double xqcdexp = par[8];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==18) {
            if (statQCD[i]>0) {
                double xqcdexp = par[9];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==19) {
            if (statQCD[i]>0) {
                double xqcdexp = par[10];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==25) {
            if (statQCD[i]>0) {
                double xqcdexp = par[11];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==26) {
            if (statQCD[i]>0) {
                double xqcdexp = par[12];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==27) {
            if (statQCD[i]>0) {
                double xqcdexp = par[13];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }
        }
        if (i==12) {
            if (statQCD[i]>0) {
                double xqcdexp = par[14];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==20) {
            if (statQCD[i]>0) {
                double xqcdexp = par[15];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==28) {
            if (statQCD[i]>0) {
                double xqcdexp = par[16];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==33) {
            if (statQCD[i]>0) {
                double xqcdexp = par[17];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==34) {
            if (statQCD[i]>0) {
                double xqcdexp = par[18];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==35) {
            if (statQCD[i]>0) {
                double xqcdexp = par[19];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        if (i==36) {
            if (statQCD[i]>0) {
                double xqcdexp = par[20];
                double xqcdpois = 1.0*statQCD[i]*TMath::Log(1.0*statQCD[i]/xqcdexp); 
                xshape += (xqcdpois+xqcdexp-1.0*statQCD[i]);
                wqcd   *= xqcdexp/statQCD[i];
            }            
        }
        
        double xexpall = wtau+wqcd+welw;
        double xexpois = 1.0*trackDATA[i]*TMath::Log(1.0*trackDATA[i]/xexpall);
        xshape += (xexpois+xexpall-1.0*trackDATA[i]);
    }
    
    double normELWK   = TMath::Gaus(par[1],1.0,0.2);
    
    f = 2.0*(xshape-TMath::Log(normELWK));
    
}

int makeTrackFit(
                 
                 TH2F* hist2D_kTtrack_data,
                 TH2F* hist2D_kTtrack_h130,
                 TH2F* hist2D_kTtrack_ztt,
                 TH2F* hist2D_kTtrack_ztt_pileuplow,
                 TH2F* hist2D_kTtrack_ztt_pileuphigh,
                 TH2F* hist2D_kTtrack_qcd,
                 TH2F* hist2D_kTtrack_qcdSyst,
                 TH2F* hist2D_kTtrack_qcdstat,
                 TH2F* hist2D_kTtrack_elwk,
                 TH2F* hist2D_kTtrack_elwk_SS,
                 TH2F* hist2D_kTtrack_ztt_SS,
                 float& QCD_norm,
                 float& QCD_err,
                 float& QCD_errL,
                 float& QCD_errR,
                 float& Ztt_norm,
                 float& Ztt_err, 
                 float& EW_norm,
                 float& EW_err,
                 float& H_norm,
                 float& QCD_norm13,
                 float& QCD_err13, 
                 float& Ztt_norm13,
                 float& Ztt_err13, 
                 float& EW_norm13,
                 float& H_norm13,
                 int region,
                 bool mode=false
                 
                 ){
    
    int imode = 0; // 0:read data  1:make psedo-data
    if(mode) imode = 1;
    category = region;
    TRandom rand(iRandomSeed);     // Initialization
    
    
    std::cout<<hist2D_kTtrack_data->GetName()<<"\t"<<hist2D_kTtrack_data->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_h130->GetName()<<"\t"<<hist2D_kTtrack_h130->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_ztt->GetName()<<"\t"<<hist2D_kTtrack_ztt->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_ztt_pileuplow->GetName()<<"\t"<<hist2D_kTtrack_ztt_pileuplow->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_ztt_pileuphigh->GetName()<<"\t"<<hist2D_kTtrack_ztt_pileuphigh->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_qcd->GetName()<<"\t"<<hist2D_kTtrack_qcd->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_qcdSyst->GetName()<<"\t"<<hist2D_kTtrack_qcdSyst->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_qcdstat->GetName()<<"\t"<<hist2D_kTtrack_qcdstat->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_elwk->GetName()<<"\t"<<hist2D_kTtrack_elwk->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_elwk_SS->GetName()<<"\t"<<hist2D_kTtrack_elwk_SS->Integral()<<std::endl;
    std::cout<<hist2D_kTtrack_ztt_SS->GetName()<<"\t"<<hist2D_kTtrack_ztt_SS->Integral()<<std::endl;

    
    TH1F hist1D_track_TAU("hist1D_track_TAU","",64,0,64);
    TH1F hist1D_track_QCD("hist1D_track_QCD","",64,0,64);
    TH1F hist1D_track_ELW("hist1D_track_ELW","",64,0,64);
    
    //=====================================
    // Make PDF for track multiplicity fit
    //=====================================
    sumDATA = 0.0;
    sumTAU = 0.0;
    sumQCD = 0.0;
    sumHIG = 0.0;
    sumELW = 0.0;
    sumELW_SS = 0.0;
    sumQCD2 = 0.0;
    sumTAU_pileuplow = 0.0;
    sumTAU_pileuphigh = 0.0;
    sumQCDstat = 0.0;
    sumTAUstat = 0.0;

    
    /*for(int i=0;i<64;i++){
     dTAU[i]=0.0;
     dQCD[i]=0.0;
     dHIG[i]=0.0;
     dELW[i]=0.0;
     dELW_SS[i]=0.0;
     trackDATA[i]=0.0;
     statQCD[i]=0.0;
     dQCD2[i]=0.0;
     dTAU_pileuplow[i]=0.0;
     dTAU_pileuphigh[i]=0.0;
     }*/
    
    
    //hist2D_kTtrack_qcd->Scale(hist2D_kTtrack_qcdstat->Integral()/hist2D_kTtrack_qcd->Integral());
    //hist2D_kTtrack_qcdSyst->Scale(hist2D_kTtrack_qcdstat->Integral()/hist2D_kTtrack_qcdSyst->Integral());

    
    int NbinX    = hist2D_kTtrack_ztt->GetNbinsX();
    int NbinY    = hist2D_kTtrack_ztt->GetNbinsY();
    for (int i=0; i<NbinX; i++) {
        for (int j=0; j<NbinY; j++) {
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)>0)  { sumTAU += hist2D_kTtrack_ztt->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)<=0)  { sumTAU += 0.000001; }

            //if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)>0)  { sumTAU += 5.*hist2D_kTtrack_h130->GetBinContent(i+1,j+1); } // SIGNAL INJECTION
            
            
            if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)>0)  { sumQCD += hist2D_kTtrack_qcd->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)<=0)  { sumQCD += 0.000001; }

            // if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)>0)  { sumQCD += 0.2*hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1); }
            
            if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)>0) { sumQCD2 += hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)<=0)  { sumQCD2 += 0.000001; }
            
            // if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)>0)  { sumQCD += 0.2*hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1); }
            
            if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)>0) { sumHIG += hist2D_kTtrack_h130->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)<=0)  { sumHIG += 0.000001; }

            if (hist2D_kTtrack_elwk->GetBinContent(i+1,j+1)>0) { sumELW += hist2D_kTtrack_elwk->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_elwk->GetBinContent(i+1,j+1)<=0)  { sumELW += 0.000001; }

            if (hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)>0) { sumELW_SS += hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)<=0)  { sumELW_SS += 0.000001; }


            
            if (hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1)>0)  { sumTAU_pileuplow += hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1)<=0)  { sumTAU_pileuplow += 0.000001; }

            if (hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1)>0) { sumTAU_pileuphigh += hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1)<=0)  { sumTAU_pileuphigh += 0.000001; }

            
            sumQCDstat += hist2D_kTtrack_qcdstat->GetBinContent(i+1,j+1);

            
        }
    }
    for (int i=0; i<NbinX; i++) {
        for (int j=0; j<NbinY; j++) {
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)>0)  { dTAU[NbinY*i+j] = hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)/sumTAU; }
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)<=0)  { dTAU[NbinY*i+j] = 0.000001/sumTAU; }
            
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)>0)  { statTAU[NbinY*i+j] = hist2D_kTtrack_ztt->GetBinError(i+1,j+1)/hist2D_kTtrack_ztt->GetBinContent(i+1,j+1); }
            if (hist2D_kTtrack_ztt->GetBinContent(i+1,j+1)<=0)  { statTAU[NbinY*i+j] = 0.0; }
            
            //if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)>0)  { dTAU[NbinY*i+j] += 5.*hist2D_kTtrack_h130->GetBinContent(i+1,j+1)/sumTAU; }  // SIGNAL INJECTION
            
            if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)>0)  { dQCD[NbinY*i+j] = hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)/sumQCD; }
            if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)<=0)  { dQCD[NbinY*i+j] = 0.000001/sumQCD; }

            //if (hist2D_kTtrack_qcd->GetBinContent(i+1,j+1)>0)  { dQCD[NbinY*i+j] += 0.2*hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)/sumQCD; }
            
            if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)>0) { dQCD2[NbinY*i+j] = hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)/sumQCD2; }
            if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)<=0)  { dQCD2[NbinY*i+j] = 0.000001/sumQCD2; }
            
            //if (hist2D_kTtrack_qcdSyst->GetBinContent(i+1,j+1)>0)  { dQCD2[NbinY*i+j] += 0.2*hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)/sumQCD2; }

            
            if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)>0) { dHIG[NbinY*i+j] = hist2D_kTtrack_h130->GetBinContent(i+1,j+1)/sumHIG; }
            if (hist2D_kTtrack_h130->GetBinContent(i+1,j+1)<=0)  { dHIG[NbinY*i+j] = 0.000001/sumHIG; }

            if (hist2D_kTtrack_elwk->GetBinContent(i+1,j+1)>0) { dELW[NbinY*i+j] = hist2D_kTtrack_elwk->GetBinContent(i+1,j+1)/sumELW; }
            if (hist2D_kTtrack_elwk->GetBinContent(i+1,j+1)<=0)  { dELW[NbinY*i+j] = 0.000001/sumELW; }

            if (hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)>0) { dELW_SS[NbinY*i+j] = hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)/sumELW_SS; }
            if (hist2D_kTtrack_elwk_SS->GetBinContent(i+1,j+1)<=0)  { dELW_SS[NbinY*i+j] = 0.000001/sumELW_SS; }
            
            if (hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1)>0)  { dTAU_pileuplow[NbinY*i+j] = hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1)/sumTAU_pileuplow; }
            if (hist2D_kTtrack_ztt_pileuplow->GetBinContent(i+1,j+1)<=0)  { dTAU_pileuplow[NbinY*i+j] = 0.000001/sumTAU_pileuplow; }

            if (hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1)>0) { dTAU_pileuphigh[NbinY*i+j] = hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1)/sumTAU_pileuphigh; }
            if (hist2D_kTtrack_ztt_pileuphigh->GetBinContent(i+1,j+1)<=0)  { dTAU_pileuphigh[NbinY*i+j] = 0.000001/sumTAU_pileuphigh; }

            
            statQCD[NbinY*i+j] = hist2D_kTtrack_qcdstat->GetBinContent(i+1,j+1);
            
            hist1D_track_TAU.SetBinContent(NbinY*i+j+1,dTAU[NbinY*i+j]);
            hist1D_track_QCD.SetBinContent(NbinY*i+j+1,dQCD[NbinY*i+j]);
            hist1D_track_ELW.SetBinContent(NbinY*i+j+1,dELW[NbinY*i+j]);
        }
    }
    
    
    
    
    //====================================
    // Read real DATA or make pseudo-data
    //====================================
    sumDATA = 0;
    for (int i=0; i<NbinX; i++) {
        for (int j=0; j<NbinY; j++) {
            trackDATA[NbinY*i+j] = 0.0;
        }
    }
    
    // Read data
    if (imode==false) {
        for (int i=0; i<NbinX; i++) {
            for (int j=0; j<NbinY; j++) {
                trackDATA[NbinY*i+j] = hist2D_kTtrack_data->GetBinContent(i+1,j+1);
                //trackDATA[NbinY*i+j] += 5.*hist2D_kTtrack_h130->GetBinContent(i+1,j+1); // SIGNAL INJECTION!!!
                sumDATA += trackDATA[NbinY*i+j];
            }
        }
    }
    
    // Or make pseudo-data
    if (imode==true) {
        sumTAU*=1;
        sumELW*=1;
        sumQCD*=1;
        
        
        for (int i=0; i<TMath::Nint(sumTAU); i++)  { sumDATA++; trackDATA[TMath::Nint(hist1D_track_TAU.GetRandom()-0.5)]++; }
        for (int i=0; i<TMath::Nint(sumQCD); i++)  { sumDATA++; trackDATA[TMath::Nint(hist1D_track_QCD.GetRandom()-0.5)]++; }
        for (int i=0; i<TMath::Nint(sumELW); i++)  { sumDATA++; trackDATA[TMath::Nint(hist1D_track_ELW.GetRandom()-0.5)]++; }
    }
    
    
    //==============
    // Do track fit
    //==============
    double fitQCD,errQCD;
    double xelwk,xelwker;
    double qcdsh,qcdsher;
    double pllw,pllwer;
    double plhg,plhger;
    double cov11,cov11er;
    double cov12,cov12er;
    double cov13,cov13er;
    double cov14,cov14er;
    double cov21,cov21er;
    double cov22,cov22er;
    double cov23,cov23er;
    double cov24,cov24er;
    double cov31,cov31er;
    double cov32,cov32er;
    double cov33,cov33er;
    double cov34,cov34er;
    double cov41,cov41er;
    double cov42,cov42er;
    double cov43,cov43er;
    double cov44,cov44er;

        
    double effTauMediumErrR,effTauMediumErrL;
    
    double xLog1,xLog2,edm1,edm2,errdef;
    int nvpar,nparx,icstat;
    double eparab;
    double globcc;
    
    int ierflg = 0;
    double arglist[21];
    double vstart[21];
    double step[21];
    
    vstart[0] = sumQCDstat/(sumTAU+sumQCDstat+sumELW);
    vstart[1] = 1.0;
    vstart[2] = 0.0;
    vstart[3] = 0.0;
    vstart[4] = 0.0;
    vstart[5] = statQCD[9];
    vstart[6] = statQCD[10];
    vstart[7] = statQCD[11];
    vstart[8] = statQCD[17];
    vstart[9] = statQCD[18];
    vstart[10] = statQCD[19];
    vstart[11] = statQCD[25];
    vstart[12] = statQCD[26];
    vstart[13] = statQCD[27];
    vstart[14] = statQCD[12];
    vstart[15] = statQCD[20];
    vstart[16] = statQCD[28];
    vstart[17] = statQCD[33];
    vstart[18] = statQCD[34];
    vstart[19] = statQCD[35];
    vstart[20] = statQCD[36];
    
    
    for (int i=0; i<21; i++) { step[i] = 0.001; }
    
    TMinuit trkMinuit1(21);
    trkMinuit1.SetFCN(trackFit);
    trkMinuit1.mnparm(0,"Fraction of QCD         ",vstart[0], step[0],  0.,0.,ierflg);
    trkMinuit1.mnparm(1,"Electroweak             ",vstart[1], step[1],  0.,0.,ierflg);
    trkMinuit1.mnparm(2,"qcd syst                ",vstart[2], step[2],  0.,0.,ierflg);
    trkMinuit1.mnparm(3,"tau low pu              ",vstart[3], step[3],  0.,0.,ierflg);
    trkMinuit1.mnparm(4,"tau high pu             ",vstart[4], step[4],  0.,0.,ierflg);    
    trkMinuit1.mnparm(5,"stat uncertainty (1,1) ",vstart[5],step[5], 0.,0.,ierflg);
    trkMinuit1.mnparm(6,"stat uncertainty (1,2) ",vstart[6],step[6], 0.,0.,ierflg);
    trkMinuit1.mnparm(7,"stat uncertainty (1,3) ",vstart[7],step[7], 0.,0.,ierflg);
    trkMinuit1.mnparm(8,"stat uncertainty (2,1) ",vstart[8],step[8], 0.,0.,ierflg);
    trkMinuit1.mnparm(9,"stat uncertainty (2,2) ",vstart[9],step[9], 0.,0.,ierflg);
    trkMinuit1.mnparm(10,"stat uncertainty (2,3) ",vstart[10],step[10], 0.,0.,ierflg);
    trkMinuit1.mnparm(11,"stat uncertainty (3,1) ",vstart[11],step[11], 0.,0.,ierflg);
    trkMinuit1.mnparm(12,"stat uncertainty (3,2) ",vstart[12],step[12], 0.,0.,ierflg);
    trkMinuit1.mnparm(13,"stat uncertainty (3,3) ",vstart[13],step[13], 0.,0.,ierflg);
    trkMinuit1.mnparm(14,"stat uncertainty (1,4) ",vstart[14],step[14], 0.,0.,ierflg);
    trkMinuit1.mnparm(15,"stat uncertainty (2,4) ",vstart[15],step[15], 0.,0.,ierflg);
    trkMinuit1.mnparm(16,"stat uncertainty (3,4) ",vstart[16],step[16], 0.,0.,ierflg);
    trkMinuit1.mnparm(17,"stat uncertainty (4,1) ",vstart[17],step[17], 0.,0.,ierflg);
    trkMinuit1.mnparm(18,"stat uncertainty (4,2) ",vstart[18],step[18], 0.,0.,ierflg);
    trkMinuit1.mnparm(19,"stat uncertainty (4,3) ",vstart[19],step[19], 0.,0.,ierflg);
    trkMinuit1.mnparm(20,"stat uncertainty (4,4) ",vstart[20],step[20], 0.,0.,ierflg);
    
    //arglist[0] =  4; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
      arglist[0] =  5; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  6; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  7; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  8; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  9; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  10; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  11; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  12; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  13; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  14; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  15; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  16; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  17; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  18; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  19; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  20; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);
//    arglist[0] =  21; trkMinuit1.mnexcm("FIX", arglist ,1,ierflg);


    
    
    arglist[0] = 0.0; trkMinuit1.mnexcm("MINI", arglist ,0,ierflg);
    arglist[0] = 1000; trkMinuit1.mnexcm("MIGRAD", arglist ,0,ierflg);
    
    if (!trkMinuit1.fCstatu.Contains("CONVERGED")) {
        std::cout<<"No convergence at fitting"<<std::endl;
        std::cout<<"Minuit return string:"<<trkMinuit1.fCstatu.Data();
        return -1;
    }
    
    arglist[0] = 0.0; arglist[1] = 1.0; trkMinuit1.mnexcm("MINOS", arglist ,2,ierflg);
    
    if (!trkMinuit1.fCstatu.Contains("SUCCESSFUL")) {
        std::cout<<"No convergence at fitting"<<std::endl;
        std::cout<<"Minuit return string:"<<trkMinuit1.fCstatu.Data();
        return -1;
    }
    
    trkMinuit1.mnerrs(0,effTauMediumErrR,effTauMediumErrL,errQCD,globcc);
    std::cout<<"QCD err: L "<<effTauMediumErrL<<" R "<<effTauMediumErrR<<" center "<<errQCD<<std::endl;
    
    trkMinuit1.GetParameter(0,fitQCD,errQCD);
    trkMinuit1.GetParameter(1,xelwk,xelwker);
    trkMinuit1.GetParameter(2,qcdsh,qcdsher);
    trkMinuit1.GetParameter(3,pllw,pllwer);
    trkMinuit1.GetParameter(4,plhg,plhger);
    trkMinuit1.GetParameter(5,cov11,cov11er);
    trkMinuit1.GetParameter(6,cov12,cov12er);
    trkMinuit1.GetParameter(7,cov13,cov13er);
    trkMinuit1.GetParameter(8,cov21,cov21er);
    trkMinuit1.GetParameter(9,cov22,cov22er);
    trkMinuit1.GetParameter(10,cov23,cov23er);
    trkMinuit1.GetParameter(11,cov31,cov31er);
    trkMinuit1.GetParameter(12,cov32,cov32er);
    trkMinuit1.GetParameter(13,cov33,cov33er);
    trkMinuit1.GetParameter(14,cov14,cov14er);
    trkMinuit1.GetParameter(15,cov24,cov24er);
    trkMinuit1.GetParameter(16,cov34,cov34er);
    trkMinuit1.GetParameter(17,cov41,cov41er);
    trkMinuit1.GetParameter(18,cov42,cov42er);
    trkMinuit1.GetParameter(19,cov43,cov43er);
    trkMinuit1.GetParameter(20,cov44,cov44er);
    
    
    trkMinuit1.mnstat(xLog1,edm1,errdef,nvpar,nparx,icstat);
    trkMinuit1.mnexcm("STOP",arglist,0,ierflg);
    trkMinuit1.Delete();
    
    QCD_errL = effTauMediumErrL;
    QCD_errR = effTauMediumErrR;
    
    double sumModel; double dModel[64]; double oldQCD[64]; double oldTAU[64];
    
    sumModel = sumDATA;
    for(int i=0;i<64;i++){
        oldQCD[i] = dQCD[i];
        oldTAU[i] = dTAU[i];
    }
    
    if (statQCD[9]>0)  { dQCD[9]  *= cov11/statQCD[9];  }
    if (statQCD[10]>0) { dQCD[10] *= cov12/statQCD[10]; }
    if (statQCD[11]>0) { dQCD[11] *= cov13/statQCD[11]; }
    if (statQCD[17]>0) { dQCD[17] *= cov21/statQCD[17]; }
    if (statQCD[18]>0) { dQCD[18] *= cov22/statQCD[18]; }
    if (statQCD[19]>0) { dQCD[19] *= cov23/statQCD[19]; }
    if (statQCD[25]>0) { dQCD[25] *= cov31/statQCD[25]; }
    if (statQCD[26]>0) { dQCD[26] *= cov32/statQCD[26]; }
    if (statQCD[27]>0) { dQCD[27] *= cov33/statQCD[27]; }
    if (statQCD[12]>0) { dQCD[12] *= cov14/statQCD[12]; }
    if (statQCD[20]>0) { dQCD[20] *= cov24/statQCD[20]; }
    if (statQCD[28]>0) { dQCD[28] *= cov34/statQCD[28]; }
    if (statQCD[33]>0) { dQCD[33] *= cov41/statQCD[33]; }
    if (statQCD[34]>0) { dQCD[34] *= cov42/statQCD[34]; }
    if (statQCD[35]>0) { dQCD[35] *= cov43/statQCD[35]; }
    if (statQCD[36]>0) { dQCD[36] *= cov44/statQCD[36]; }

    
    for (int i=0; i<64; i++) {
        dQCD[i] += qcdsh*(dQCD2[i]-dQCD[i]);
		dTAU[i] += pllw*(dTAU_pileuphigh[i]-dTAU[i])*stepFunction(pllw) +   pllw*(dTAU[i]-dTAU_pileuplow[i])*stepFunction(-1.*pllw);
        //dTAU[i] += pllw*(dTAU_pileuplow[i]-dTAU[i]);
        //dTAU[i] += plhg*(dTAU_pileuphigh[i]-dTAU[i]);
    }
    
    for(int i=0;i<64;i++){
        
        dModel[i] = dQCD[i]*sumDATA*fitQCD+dELW[i]*xelwk*sumELW+dTAU[i]*(sumDATA-sumDATA*fitQCD-xelwk*sumELW);
    }
    
    
    if(!mode){
        // Draw best fit:
        
        TH1F histo_data("histo_data","histo_data",16,0.5,16.5);
        TH1F histo_model("histo_model","histo_model",16,0.5,16.5);
        TH1F histo_QCD("histo_QCD","histo_QCD",16,0.5,16.5);
        TH1F histo_QCDsyst("histo_QCDsyst","histo_QCDsyst",16,0.5,16.5);
        TH1F histo_TAU("histo_TAU","histo_TAU",16,0.5,16.5);
        TH1F histo_EW("histo_EW","histo_EW",16,0.5,16.5);
        TH1F histo_bestQCD("histo_bestQCD","histo_bestQCD",16,0.5,16.5);
        TH1F histo_bestTAU("histo_bestTAU","histo_bestTAU",16,0.5,16.5);
        
//        TH1F histo_data("histo_data","histo_data",25,0.5,25.5);
//        TH1F histo_model("histo_model","histo_model",25,0.5,25.5);
//        TH1F histo_QCD("histo_QCD","histo_QCD",25,0.5,25.5);
//        TH1F histo_QCDsyst("histo_QCDsyst","histo_QCDsyst",25,0.5,25.5);
//        TH1F histo_TAU("histo_TAU","histo_TAU",25,0.5,25.5);
//        TH1F histo_EW("histo_EW","histo_EW",25,0.5,25.5);
//        TH1F histo_bestQCD("histo_bestQCD","histo_bestQCD",25,0.5,25.5);
//        TH1F histo_bestTAU("histo_bestTAU","histo_bestTAU",25,0.5,25.5);
        
        
        for(int i=0;i<64;i++){
            if(i!=9 && i!=10 && i!=11 && i!=12 && i!=17 && i!=18 && i!=19 && i!=20 && i!=25 && i!=26 && i!=27 && i!=28 && i!=33 && i!=34 && i!=35 && i!=36 ) continue;
            //if(i<9 || i==14 || i==15 || i==16 || i==22 || i==23 || i==24 || i==30 || i==31 || i==32 || i==38 || i==39 || i==40 || i>46) continue;

            int bin=0;
            if(i==9)  bin =1;
            if(i==10) bin =2;
            if(i==11) bin =3;
            if(i==12) bin =4;
            if(i==17) bin =5;
            if(i==18) bin =6;
            if(i==19) bin =7;
            if(i==20) bin =8;
            if(i==25) bin =9;
            if(i==26) bin =10;
            if(i==27) bin =11;
            if(i==28) bin =12;
            if(i==33) bin =13;
            if(i==34) bin =14;
            if(i==35) bin =15;
            if(i==36) bin =16;
            
//            if(i==9)  bin =1;
//            if(i==10) bin =2;
//            if(i==11) bin =3;
//            if(i==12) bin =4;
//            if(i==13) bin =5;
//            if(i==17) bin =6;
//            if(i==18) bin =7;
//            if(i==19) bin =8;
//            if(i==20) bin =9;
//            if(i==21) bin =10;
//            if(i==25) bin =11;
//            if(i==26) bin =12;
//            if(i==27) bin =13;
//            if(i==28) bin =14;
//            if(i==29) bin =15;
//            if(i==33) bin =16;
//            if(i==34) bin =17;
//            if(i==35) bin =18;
//            if(i==36) bin =19;
//            if(i==37) bin =20;
//            if(i==41) bin =21;
//            if(i==42) bin =22;
//            if(i==43) bin =23;
//            if(i==44) bin =24;
//            if(i==45) bin =25;
            
            histo_data.SetBinContent(bin, trackDATA[i]);
            histo_model.SetBinContent(bin, dModel[i]);
            histo_QCD.SetBinContent(bin, oldQCD[i]*sumDATA*fitQCD);
            histo_QCDsyst.SetBinContent(bin, dQCD2[i]*sumDATA*fitQCD);
            histo_bestQCD.SetBinContent(bin, dQCD[i]*sumDATA*fitQCD);
            histo_TAU.SetBinContent(bin, oldTAU[i]*(sumDATA-sumDATA*fitQCD-xelwk*sumELW));
            histo_bestTAU.SetBinContent(bin, dTAU[i]*(sumDATA-sumDATA*fitQCD-xelwk*sumELW));
            histo_EW.SetBinContent(bin, dELW[i]*xelwk*sumELW);
        }
        
        gStyle->SetTextSize(0.05);
        //gStyle->SetTitleYOffset(1.4);
        
        histo_model.SetLineColor(kGray+2);
        histo_model.SetLineWidth(2);
        histo_model.SetLineStyle(2);
        
        histo_QCD.SetLineColor(kRed);
        histo_QCD.SetLineWidth(2);
        histo_QCD.SetLineStyle(1);
        
        histo_QCDsyst.SetLineColor(kBlue);
        histo_QCDsyst.SetLineWidth(2);
        histo_QCDsyst.SetLineStyle(1);
        
        histo_bestQCD.SetLineColor(kRed);
        histo_bestQCD.SetLineWidth(1);
        histo_bestQCD.SetLineStyle(2);
        
        histo_TAU.SetLineColor(kBlue);
        histo_TAU.SetLineWidth(2);
        histo_TAU.SetLineStyle(1);
        
        histo_bestTAU.SetLineColor(kBlue);
        histo_bestTAU.SetLineWidth(1);
        histo_bestTAU.SetLineStyle(2);
        
        histo_EW.SetLineColor(kGreen);
        histo_EW.SetLineWidth(2);
        histo_EW.SetLineStyle(1);
        
        float max =0.;
        if(histo_data.GetMaximum()>max) max = histo_data.GetMaximum();
        if(histo_model.GetMaximum()>max) max = histo_model.GetMaximum();
        if(histo_QCD.GetMaximum()>max) max = histo_QCD.GetMaximum();
        if(histo_bestQCD.GetMaximum()>max) max = histo_bestQCD.GetMaximum();
        if(histo_TAU.GetMaximum()>max) max = histo_TAU.GetMaximum();
        if(histo_bestTAU.GetMaximum()>max) max = histo_bestTAU.GetMaximum();
        if(histo_EW.GetMaximum()>max) max = histo_EW.GetMaximum();
        
        histo_model.SetMaximum(max*1.1);
        histo_model.GetYaxis()->SetTitle("Events");
        histo_model.GetXaxis()->SetTitle("Number of tracks (lead,sublead)");
        histo_model.GetYaxis()->SetTitleSize(0.06);
        histo_model.GetYaxis()->SetTitleOffset(1.2);
        histo_model.GetYaxis()->SetLabelSize(0.06);
        histo_model.GetXaxis()->SetTitleSize(0.06);
        histo_model.GetXaxis()->SetLabelSize(0.06);
        
        histo_model.GetXaxis()->SetBinLabel(1, "(1,1)");
        histo_model.GetXaxis()->SetBinLabel(2, "(1,2)");
        histo_model.GetXaxis()->SetBinLabel(3, "(1,3)");
        histo_model.GetXaxis()->SetBinLabel(4, "(1,4)");
        histo_model.GetXaxis()->SetBinLabel(5, "(2,1)");
        histo_model.GetXaxis()->SetBinLabel(6, "(2,2)");
        histo_model.GetXaxis()->SetBinLabel(7, "(2,3)");
        histo_model.GetXaxis()->SetBinLabel(8, "(2,4)");
        histo_model.GetXaxis()->SetBinLabel(9, "(3,1)");
        histo_model.GetXaxis()->SetBinLabel(10, "(3,2)");
        histo_model.GetXaxis()->SetBinLabel(11, "(3,3)");
        histo_model.GetXaxis()->SetBinLabel(12, "(3,4)");
        histo_model.GetXaxis()->SetBinLabel(13, "(4,1)");
        histo_model.GetXaxis()->SetBinLabel(14, "(4,2)");
        histo_model.GetXaxis()->SetBinLabel(15, "(4,3)");
        histo_model.GetXaxis()->SetBinLabel(16, "(4,4)");
        
//        histo_model.GetXaxis()->LabelsOption("v");
//        histo_model.GetXaxis()->SetBinLabel(1, "(<2#pi^{0},<2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(2, "(<2#pi^{0},2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(3, "(<2#pi^{0},2)");
//        histo_model.GetXaxis()->SetBinLabel(4, "(<2#pi^{0},3)");
//        histo_model.GetXaxis()->SetBinLabel(5, "(<2#pi^{0},4)");
//        histo_model.GetXaxis()->SetBinLabel(6, "(2#pi^{0},<2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(7, "(2#pi^{0},2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(8, "(2#pi^{0},2)");
//        histo_model.GetXaxis()->SetBinLabel(9, "(2#pi^{0},3)");
//        histo_model.GetXaxis()->SetBinLabel(10, "(2#pi^{0},4)");
//        histo_model.GetXaxis()->SetBinLabel(11, "(2,<2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(12, "(2,2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(13, "(2,2)");
//        histo_model.GetXaxis()->SetBinLabel(14, "(2,3)");
//        histo_model.GetXaxis()->SetBinLabel(15, "(2,4)");
//        histo_model.GetXaxis()->SetBinLabel(16, "(3,<2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(17, "(3,2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(18, "(3,2)");
//        histo_model.GetXaxis()->SetBinLabel(19, "(3,3)");
//        histo_model.GetXaxis()->SetBinLabel(20, "(3,4)");
//        histo_model.GetXaxis()->SetBinLabel(21, "(4,<2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(22, "(4,2#pi^{0})");
//        histo_model.GetXaxis()->SetBinLabel(23, "(4,2)");
//        histo_model.GetXaxis()->SetBinLabel(24, "(4,3)");
//        histo_model.GetXaxis()->SetBinLabel(25, "(4,4)");
        
        histo_data.GetYaxis()->SetTitle("Events");
        histo_data.GetXaxis()->SetTitle("Number of tracks (lead,sublead)");
        histo_data.GetYaxis()->SetTitleSize(0.06);
        histo_data.GetYaxis()->SetTitleOffset(1.2);
        histo_data.GetYaxis()->SetLabelSize(0.06);
        histo_data.GetXaxis()->SetTitleSize(0.06);
        histo_data.GetXaxis()->SetLabelSize(0.06);
        
        histo_data.GetXaxis()->SetBinLabel(1, "(1,1)");
        histo_data.GetXaxis()->SetBinLabel(2, "(1,2)");
        histo_data.GetXaxis()->SetBinLabel(3, "(1,3)");
        histo_data.GetXaxis()->SetBinLabel(4, "(1,4)");
        histo_data.GetXaxis()->SetBinLabel(5, "(2,1)");
        histo_data.GetXaxis()->SetBinLabel(6, "(2,2)");
        histo_data.GetXaxis()->SetBinLabel(7, "(2,3)");
        histo_data.GetXaxis()->SetBinLabel(8, "(2,4)");
        histo_data.GetXaxis()->SetBinLabel(9, "(3,1)");
        histo_data.GetXaxis()->SetBinLabel(10, "(3,2)");
        histo_data.GetXaxis()->SetBinLabel(11, "(3,3)");
        histo_data.GetXaxis()->SetBinLabel(12, "(3,4)");
        histo_data.GetXaxis()->SetBinLabel(13, "(4,1)");
        histo_data.GetXaxis()->SetBinLabel(14, "(4,2)");
        histo_data.GetXaxis()->SetBinLabel(15, "(4,3)");
        histo_data.GetXaxis()->SetBinLabel(16, "(4,4)");
        
//        histo_data.GetXaxis()->LabelsOption("v");
//        histo_data.GetXaxis()->SetBinLabel(1, "(<2#pi^{0},<2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(2, "(<2#pi^{0},2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(3, "(<2#pi^{0},2)");
//        histo_data.GetXaxis()->SetBinLabel(4, "(<2#pi^{0},3)");
//        histo_data.GetXaxis()->SetBinLabel(5, "(<2#pi^{0},4)");
//        histo_data.GetXaxis()->SetBinLabel(6, "(2#pi^{0},<2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(7, "(2#pi^{0},2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(8, "(2#pi^{0},2)");
//        histo_data.GetXaxis()->SetBinLabel(9, "(2#pi^{0},3)");
//        histo_data.GetXaxis()->SetBinLabel(10, "(2#pi^{0},4)");
//        histo_data.GetXaxis()->SetBinLabel(11, "(2,<2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(12, "(2,2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(13, "(2,2)");
//        histo_data.GetXaxis()->SetBinLabel(14, "(2,3)");
//        histo_data.GetXaxis()->SetBinLabel(15, "(2,4)");
//        histo_data.GetXaxis()->SetBinLabel(16, "(3,<2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(17, "(3,2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(18, "(3,2)");
//        histo_data.GetXaxis()->SetBinLabel(19, "(3,3)");
//        histo_data.GetXaxis()->SetBinLabel(20, "(3,4)");
//        histo_data.GetXaxis()->SetBinLabel(21, "(4,<2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(22, "(4,2#pi^{0})");
//        histo_data.GetXaxis()->SetBinLabel(23, "(4,2)");
//        histo_data.GetXaxis()->SetBinLabel(24, "(4,3)");
//        histo_data.GetXaxis()->SetBinLabel(25, "(4,4)");
        
        
        TCanvas canvas("can","can",800,600);
        
        canvas.Divide(1,2);
        
        double t_size = 0.06;
        double t_off = 0.8;
        canvas.cd(1);
        gPad->SetPad( .005, .3525,.995,.995);
        gPad->SetBottomMargin(0.001);
        gPad->SetRightMargin(0.10); 
        double lenght_1 = gPad->XtoPixel(gPad->GetX2());
        double height_1 = gPad->YtoPixel(gPad->GetY1());        
        float chi2 = 0.;
        for(int i=0;i<histo_data.GetNbinsX();i++){
            float totBKg = histo_model.GetBinContent(i+1);
            float totData = histo_data.GetBinContent(i+1);
            if(totBKg || totData){ 
                if(TMath::Nint(totBKg)) chi2+=(totData-TMath::Nint(totBKg))*(totData-TMath::Nint(totBKg))/(TMath::Nint(totBKg)+totBKg*totBKg*errQCD*errQCD);
                if(!TMath::Nint(totBKg)) chi2+=(totData-TMath::Nint(totBKg))*(totData-TMath::Nint(totBKg))/(1.);
                
            }
            //std::cout<<i<<" "<<totData<<" "<<totBKg<<" "<<chi2<<std::endl;
        }
        
        canvas.cd(2);
        gPad->SetPad( .005, .005,.995,.3525);
        gPad->SetTopMargin(0.001);
        gPad->SetBottomMargin(0.30); 
        double lenght_2 = gPad->XtoPixel(gPad->GetX2());
        double height_2 = gPad->YtoPixel(gPad->GetY1());
        
        //TPad pad("pad","pad",0.,0.25,1.,1.);
        //pad.Draw();
        canvas.cd(1);
        
        histo_model.Draw("HIST");
        histo_bestQCD.Draw("SAME HIST");
        histo_bestTAU.Draw("SAME HIST");
        histo_EW.Draw("SAME HIST");
        histo_TAU.Draw("SAME HIST");
        histo_QCD.Draw("SAME HIST");
        histo_data.Draw("SAME P");
        
        TLegend* leg = new TLegend(0.55,0.5,0.85,0.9);
        leg->SetBorderSize(0.);
        leg->SetFillColor(0);
        
        leg->AddEntry(&histo_data,"#bf{data}","P");
        leg->AddEntry(&histo_model,"#bf{Fitted Model}","L");
        leg->AddEntry(&histo_QCD,"#bf{Input SS template}","L");
        leg->AddEntry(&histo_bestQCD,"#bf{Fitted SS template}","L");
        leg->AddEntry(&histo_TAU,"#bf{Input di-tau template}","L");
        leg->AddEntry(&histo_bestTAU,"#bf{Fitted di-tau template}","L");
        leg->AddEntry(&histo_EW,"#bf{EW}","L");
        
        std::stringstream label1,label2,label3,label7;
        float inputZ = sumTAU;
        float fittedZ = sumDATA-sumDATA*fitQCD-xelwk*sumELW;
        float diff = 0.;
        if(fittedZ) diff = (fittedZ-sumTAU)/fittedZ;
        label1 << "f_{Multi-jet}: "<<std::setprecision(2)<<fitQCD<<" #pm "<<std::setprecision(2)<<errQCD;
        label2 << "#chi^{2}/ndof: "<<std::setprecision(2)<<chi2<<"/"<<histo_data.GetNbinsX()-1;
        label7 << "Z (fit-MC)/fit: "<<std::setprecision(2)<<diff;
        if(category==0) label3<<"Z-Normalization";
        if(category==1) label3<<"Preselection";
        if(category==2) label3<<"VBF";
        if(category==3) label3<<"Boosted";
        if(category==4) label3<<"VBF CR w/ MW";
        if(category==5) label3<<"VBF CR w/o MW";
        if(category==6) label3<<"Boosted CR w/ MW";
        if(category==7) label3<<"Boosted CR w/o MW";
        if(category==8) label3<<"Z-Normalization";
        if(category==9) label3<<"Z-Normalization";
        if(category==10) label3<<"Z-Normalization";
        //label2 << "Estimated - Fitted EW: ("<<EW_fit<<" #pm "<<EW_fit_err<<") %";
        //label3 << "H contamination: "<<Higgs_contamination<<" % of fitted Z";
		
		std::stringstream lumi_string;
		lumi_string << "#int Ldt = 20.3 fb^{-1}, #sqrt{s}=8 TeV";

        
        leg->Draw();
        TLatex l1; //l.SetTextAlign(12);
        l1.SetTextSize(0.05); 
        l1.SetNDC();
        l1.SetTextColor(1);
        l1.DrawLatex(0.23,0.85,lumi_string.str().c_str());
		l1.DrawLatex(0.23,0.77,label3.str().c_str());
		l1.DrawLatex(0.23,0.70,label2.str().c_str());
        l1.DrawLatex(0.23,0.63,label1.str().c_str());
        l1.DrawLatex(0.23,0.56,label7.str().c_str());

        
        //ATLAS_LABEL(0.27,0.85,1);myText(0.35,0.85,1,"Internal");
        
        
        canvas.cd(2);
        //gPad->SetGridy();
        
        std::string histo_name(histo_data.GetName());
        histo_name += "_ratio";
        TH1F* tmp_tot_histo = (TH1F*)histo_data.Clone(histo_name.c_str());
        tmp_tot_histo->Sumw2();
        tmp_tot_histo->Divide(&histo_model);
        int nbins = tmp_tot_histo->GetNbinsX();
        for(int i=0;i<nbins+2;i++){
            float totBKg = histo_model.GetBinContent(i);
            float errorB = sqrt(TMath::Nint(totBKg)+totBKg*totBKg*errQCD*errQCD);
            float totData = histo_data.GetBinContent(i);
            float errorA = sqrt(TMath::Nint(totData));
            float error = sqrt(errorA*errorA/(totBKg*totBKg)+ totData*totData*errorB*errorB/(totBKg*totBKg*totBKg*totBKg));
            tmp_tot_histo->SetBinError(i,error);
        }
        
        tmp_tot_histo->SetMarkerColor(1);
        tmp_tot_histo->SetMarkerStyle(20);
        tmp_tot_histo->SetLineColor(1);
        tmp_tot_histo->SetLineStyle(1);
        tmp_tot_histo->SetMinimum(0.01);
        tmp_tot_histo->SetMaximum(1.99);

        
        tmp_tot_histo->GetYaxis()->CenterTitle(true);
        tmp_tot_histo->GetYaxis()->SetNdivisions(605,kTRUE);
        tmp_tot_histo->GetYaxis()->SetLabelSize(t_size*height_1/height_2);
        tmp_tot_histo->GetYaxis()->SetTitleSize(t_size*height_1/height_2);
        tmp_tot_histo->GetYaxis()->SetTitleOffset(t_off*height_2/height_1);
        
        tmp_tot_histo->GetXaxis()->SetLabelSize(t_size*height_1/height_2);
        tmp_tot_histo->GetXaxis()->SetTitleSize(t_size*height_1/height_2);
        tmp_tot_histo->GetXaxis()->SetLabelOffset(t_off*0.05*height_2/height_1);
        
        TLine line(tmp_tot_histo->GetXaxis()->GetXmin(),1.,tmp_tot_histo->GetXaxis()->GetXmax(),1.);
        line.SetLineColor(2);
        line.SetLineStyle(2);
        line.SetLineWidth(2);
        tmp_tot_histo->Draw("ep");
        
        line.Draw();
        tmp_tot_histo->GetYaxis()->SetTitle("Data/Fit");
        //tmp_tot_histo->GetYaxis()->SetTitleSize(0.08);
        tmp_tot_histo->GetXaxis()->SetTitleOffset(1.4);
        
        canvas.cd();
        
        std::string canvas_name("trackFit_plots/Fit_plots/bestFit_");
        if(category==0) canvas_name = canvas_name +"0.eps";
        if(category==1) canvas_name = canvas_name +"1.eps";
        if(category==2) canvas_name = canvas_name +"2.eps";
        if(category==3) canvas_name = canvas_name +"3.eps";
        if(category==4) canvas_name = canvas_name +"4.eps";
        if(category==5) canvas_name = canvas_name +"5.eps";
        if(category==6) canvas_name = canvas_name +"6.eps";
        if(category==7) canvas_name = canvas_name +"7.eps";
        if(category==8) canvas_name = canvas_name +"8.eps";
        if(category==9) canvas_name = canvas_name +"9.eps";
        if(category==10) canvas_name = canvas_name +"10.eps";
        
        canvas.Print(canvas_name.c_str());
        
        
        max =0.;
        if(histo_QCD.GetMaximum()>max) max = histo_QCD.GetMaximum();
        if(histo_bestQCD.GetMaximum()>max) max = histo_bestQCD.GetMaximum();
        if(histo_QCDsyst.GetMaximum()>max) max = histo_QCDsyst.GetMaximum();
        
        histo_QCD.SetMaximum(max*1.1);
        histo_QCD.GetYaxis()->SetTitle("Events");
        histo_QCD.GetXaxis()->SetTitle("Number of tracks (lead,sublead)");
        histo_QCD.GetYaxis()->SetTitleSize(0.06);
        histo_QCD.GetYaxis()->SetTitleOffset(1.2);
        histo_QCD.GetYaxis()->SetLabelSize(0.06);
        histo_QCD.GetXaxis()->SetTitleSize(0.06);
        histo_QCD.GetXaxis()->SetLabelSize(0.06);
        
        histo_QCD.GetXaxis()->SetBinLabel(1, "(1,1)");
        histo_QCD.GetXaxis()->SetBinLabel(2, "(1,2)");
        histo_QCD.GetXaxis()->SetBinLabel(3, "(1,3)");
        histo_QCD.GetXaxis()->SetBinLabel(4, "(1,4)");
        histo_QCD.GetXaxis()->SetBinLabel(5, "(2,1)");
        histo_QCD.GetXaxis()->SetBinLabel(6, "(2,2)");
        histo_QCD.GetXaxis()->SetBinLabel(7, "(2,3)");
        histo_QCD.GetXaxis()->SetBinLabel(8, "(2,4)");
        histo_QCD.GetXaxis()->SetBinLabel(9, "(3,1)");
        histo_QCD.GetXaxis()->SetBinLabel(10, "(3,2)");
        histo_QCD.GetXaxis()->SetBinLabel(11, "(3,3)");
        histo_QCD.GetXaxis()->SetBinLabel(12, "(3,4)");
        histo_QCD.GetXaxis()->SetBinLabel(13, "(4,1)");
        histo_QCD.GetXaxis()->SetBinLabel(14, "(4,2)");
        histo_QCD.GetXaxis()->SetBinLabel(15, "(4,3)");
        histo_QCD.GetXaxis()->SetBinLabel(16, "(4,4)");
		
		for(unsigned int bin = 0; bin<histo_QCD.GetNbinsX()+2; bin++){
			float content = histo_QCDsyst.GetBinContent(bin);
			content *= 1./histo_QCD.GetBinContent(bin);
			std::cout<< histo_QCD.GetXaxis()->GetBinLabel(bin)<<" QCD shape syst: "<<content<<std::endl;
		}
        
//        histo_QCD.GetXaxis()->LabelsOption("v");
//        histo_QCD.GetXaxis()->SetBinLabel(1, "(<2#pi^{0},<2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(2, "(<2#pi^{0},2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(3, "(<2#pi^{0},2)");
//        histo_QCD.GetXaxis()->SetBinLabel(4, "(<2#pi^{0},3)");
//        histo_QCD.GetXaxis()->SetBinLabel(5, "(<2#pi^{0},4)");
//        histo_QCD.GetXaxis()->SetBinLabel(6, "(2#pi^{0},<2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(7, "(2#pi^{0},2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(8, "(2#pi^{0},2)");
//        histo_QCD.GetXaxis()->SetBinLabel(9, "(2#pi^{0},3)");
//        histo_QCD.GetXaxis()->SetBinLabel(10, "(2#pi^{0},4)");
//        histo_QCD.GetXaxis()->SetBinLabel(11, "(2,<2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(12, "(2,2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(13, "(2,2)");
//        histo_QCD.GetXaxis()->SetBinLabel(14, "(2,3)");
//        histo_QCD.GetXaxis()->SetBinLabel(15, "(2,4)");
//        histo_QCD.GetXaxis()->SetBinLabel(16, "(3,<2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(17, "(3,2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(18, "(3,2)");
//        histo_QCD.GetXaxis()->SetBinLabel(19, "(3,3)");
//        histo_QCD.GetXaxis()->SetBinLabel(20, "(3,4)");
//        histo_QCD.GetXaxis()->SetBinLabel(21, "(4,<2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(22, "(4,2#pi^{0})");
//        histo_QCD.GetXaxis()->SetBinLabel(23, "(4,2)");
//        histo_QCD.GetXaxis()->SetBinLabel(24, "(4,3)");
//        histo_QCD.GetXaxis()->SetBinLabel(25, "(4,4)");
        
        TCanvas canvas1("can1","can1",800,600);
        
        histo_QCD.Draw("HIST");
        histo_bestQCD.Draw("SAME HIST");
        histo_QCDsyst.Draw("SAME HIST");
        
        TLegend* leg1 = new TLegend(0.5,0.7,0.85,0.9);
        leg1->SetBorderSize(0.);
        leg1->SetFillColor(0);
    
        leg1->AddEntry(&histo_QCD,"Input SS template","L");
        leg1->AddEntry(&histo_bestQCD,"Fitted SS template","L");
        leg1->AddEntry(&histo_QCDsyst,"Input OS template","L");
        
        std::stringstream label4,label5,label6;

        if(category==0) label4<<"Z-Normalization";
        if(category==1) label4<<"Preselection";
        if(category==2) label4<<"VBF";
        if(category==3) label4<<"Boosted";
        if(category==4) label4<<"VBF CR w/ MW";
        if(category==5) label4<<"VBF CR w/o MW";
        if(category==6) label4<<"Boosted CR w/ MW";
        if(category==7) label4<<"Boosted CR w/o MW";
        if(category==8) label4<<"Z-Normalization";
        if(category==9) label4<<"Z-Normalization";
        if(category==10) label4<<"Z-Normalization";
        //label2 << "Estimated - Fitted EW: ("<<EW_fit<<" #pm "<<EW_fit_err<<") %";
        //label3 << "H contamination: "<<Higgs_contamination<<" % of fitted Z";
        
        leg1->Draw();

        l1.DrawLatex(0.35,0.51,label4.str().c_str());
        
        float totQCDss  = hist2D_kTtrack_qcd->Integral();
        float totEWss   = hist2D_kTtrack_elwk_SS->Integral();
        float totZttss  = hist2D_kTtrack_ztt_SS->Integral();
        
        float EWcontSS  = 0.;
        float ZttcontSS = 0.;
        float HcontOS   = 0.;
        if(totQCDss) EWcontSS = totEWss/(totEWss+totZttss+totQCDss);
        if(totQCDss) ZttcontSS = totZttss/(totZttss+totEWss+totQCDss);

        label5 << "EW cont in SS: "<<EWcontSS;
        label6 << "Ztt cont in SS: "<<ZttcontSS;
        l1.DrawLatex(0.35,0.44,label5.str().c_str());
        l1.DrawLatex(0.35,0.37,label6.str().c_str());

            
        std::string canvas_name1("trackFit_plots/Fit_plots/QCDsyst_");
        if(category==0) canvas_name1 = canvas_name1 +"0.eps";
        if(category==1) canvas_name1 = canvas_name1 +"1.eps";
        if(category==2) canvas_name1 = canvas_name1 +"2.eps";
        if(category==3) canvas_name1 = canvas_name1 +"3.eps";
        if(category==4) canvas_name1 = canvas_name1 +"4.eps";
        if(category==5) canvas_name1 = canvas_name1 +"5.eps";
        if(category==6) canvas_name1 = canvas_name1 +"6.eps";
        if(category==7) canvas_name1 = canvas_name1 +"7.eps";
        if(category==8) canvas_name1 = canvas_name1 +"8.eps";
        if(category==9) canvas_name1 = canvas_name1 +"9.eps";
        if(category==10) canvas_name1 = canvas_name1 +"10.eps";
        
        canvas1.Print(canvas_name1.c_str());

        
    }
    
    
//    double fQCD13 = dQCD[9];
//    double fTAU13 = dTAU[9];
//    double fHIG13 = dHIG[9];
//    double fELW13 = dELW[9];
    
    
    double fQCD13 = dQCD[9]+dQCD[11]+dQCD[25]+dQCD[27];
    double fTAU13 = dTAU[9]+dTAU[11]+dTAU[25]+dTAU[27];
    double fHIG13 = dHIG[9]+dHIG[11]+dHIG[25]+dHIG[27];
    double fELW13 = dELW[9]+dELW[11]+dELW[25]+dELW[27];

//    double fQCD13 = dQCD[9]+dQCD[10]+dQCD[12]+dQCD[17]+dQCD[18]+dQCD[20]+dQCD[33]+dQCD[34]+dQCD[36];
//    double fTAU13 = dTAU[9]+dTAU[10]+dTAU[12]+dTAU[17]+dTAU[18]+dTAU[20]+dTAU[33]+dTAU[34]+dTAU[36];
//    double fHIG13 = dHIG[9]+dHIG[10]+dHIG[12]+dHIG[17]+dHIG[18]+dHIG[20]+dHIG[33]+dHIG[34]+dHIG[36];
//    double fELW13 = dELW[9]+dELW[10]+dELW[12]+dELW[17]+dELW[18]+dELW[20]+dELW[33]+dELW[34]+dELW[36];
    
    // linearity correction (nominal lumi, no syst):
//    if(category==2) {
//        fitQCD = (fitQCD - 0.0831)/0.917;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.917*0.917)+0.0067*0.0067/(0.917*0.917)+(fitQCD-0.0831)*(fitQCD-0.0831)*0.013*0.013/(0.917*0.917*0.917*0.917));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.917*0.917)+0.0067*0.0067/(0.917*0.917)+(fitQCD-0.0831)*(fitQCD-0.0831)*0.013*0.013/(0.917*0.917*0.917*0.917));
//        errQCD = sqrt(errQCD*errQCD/(0.917*0.917)+0.0067*0.0067/(0.917*0.917)+(fitQCD-0.0831)*(fitQCD-0.0831)*0.013*0.013/(0.917*0.917*0.917*0.917));
//    }
//    if(category==3) {
//        fitQCD = (fitQCD - 0.0191)/0.978;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.978*0.978)+0.00083*0.00083/(0.978*0.978)+(fitQCD-0.0191)*(fitQCD-0.0191)*0.0021*0.0021/(0.978*0.978*0.978*0.978));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.978*0.978)+0.00083*0.00083/(0.978*0.978)+(fitQCD-0.0191)*(fitQCD-0.0191)*0.0021*0.0021/(0.978*0.978*0.978*0.978));
//        errQCD = sqrt(errQCD*errQCD/(0.978*0.978)+0.00083*0.00083/(0.978*0.978)+(fitQCD-0.0191)*(fitQCD-0.0191)*0.0021*0.0021/(0.978*0.978*0.978*0.978));
//    }
//    if(category==4) {
//        fitQCD = (fitQCD - 0.0248)/0.980;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.980*0.980)+0.0023*0.0023/(0.980*0.980)+(fitQCD-0.0248)*(fitQCD-0.0248)*0.0033*0.0033/(0.980*0.980*0.980*0.980));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.980*0.980)+0.0023*0.0023/(0.980*0.980)+(fitQCD-0.0248)*(fitQCD-0.0248)*0.0033*0.0033/(0.980*0.980*0.980*0.980));
//        errQCD = sqrt(errQCD*errQCD/(0.980*0.980)+0.0023*0.0023/(0.980*0.980)+(fitQCD-0.0248)*(fitQCD-0.0248)*0.0033*0.0033/(0.980*0.980*0.980*0.980));
//    }
//    if(category==5) {
//        fitQCD = (fitQCD - 0.0203)/0.986;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.986*0.986)+0.0019*0.0019/(0.986*0.986)+(fitQCD-0.0203)*(fitQCD-0.0203)*0.0027*0.0027/(0.986*0.986*0.986*0.986));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.986*0.986)+0.0019*0.0019/(0.986*0.986)+(fitQCD-0.0203)*(fitQCD-0.0203)*0.0027*0.0027/(0.986*0.986*0.986*0.986));
//        errQCD = sqrt(errQCD*errQCD/(0.986*0.986)+0.0019*0.0019/(0.986*0.986)+(fitQCD-0.0203)*(fitQCD-0.0203)*0.0027*0.0027/(0.986*0.986*0.986*0.986));
//    }
//////    
//    if(category==6) {
//        fitQCD = (fitQCD - 0.0132)/0.991;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.991*0.991)+0.0013*0.0013/(0.991*0.991)+(fitQCD-0.0132)*(fitQCD-0.0132)*0.0019*0.0019/(0.991*0.991*0.991*0.991));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.991*0.991)+0.0013*0.0013/(0.991*0.991)+(fitQCD-0.0132)*(fitQCD-0.0132)*0.0019*0.0019/(0.991*0.991*0.991*0.991));
//        errQCD = sqrt(errQCD*errQCD/(0.991*0.991)+0.0013*0.0013/(0.991*0.991)+(fitQCD-0.0132)*(fitQCD-0.0132)*0.0019*0.0019/(0.991*0.991*0.991*0.991));
//    }
//////    
//    if(category==7) {
//        fitQCD = (fitQCD - 0.00564)/0.995;
//        QCD_errL = -sqrt(QCD_errL*QCD_errL/(0.995*0.995)+0.0011*0.0011/(0.995*0.995)+(fitQCD-0.00576)*(fitQCD-0.00576)*0.0015*0.0015/(0.995*0.995*0.995*0.995));
//        QCD_errR = sqrt(QCD_errR*QCD_errR/(0.995*0.995)+0.0011*0.0011/(0.995*0.995)+(fitQCD-0.00576)*(fitQCD-0.00576)*0.0015*0.0015/(0.995*0.995*0.995*0.995));
//        errQCD = sqrt(errQCD*errQCD/(0.995*0.995)+0.0011*0.0011/(0.995*0.995)+(fitQCD-0.00576)*(fitQCD-0.00576)*0.0015*0.0015/(0.995*0.995*0.995*0.995));
//    }
    
    
    //  double obsQCD = sumDATA*fitQCD*fQCD13;
    double obsQCD = sumDATA*fitQCD;
    double obsQCDerr = sumDATA*errQCD;
    double obsH   = sumHIG;
    double obsEWK = xelwk*sumELW;
    double errEWK = xelwker*sumELW;
    double obsZ   = sumDATA-obsQCD-obsEWK;
    //double errZ   = sqrt(sumFitErr*sumFitErr*(1-fitQCD)*(1-fitQCD)+sumFit*sumFit*errQCD*errQCD+errEWK*errEWK);
    double errZ   = sumDATA*errQCD;
    
    printf("\n");
    printf("Inputs QCD: %9.4f \n",sumQCD);
    printf("Inputs ELW: %9.4f \n",sumELW);
    printf("Inputs Ztt: %9.4f \n",sumTAU);
    
    
    printf("Fitted QCD fraction:  %9.4f +- %9.4f\n",fitQCD,errQCD);
    printf("---------------------------------\n");
    printf("Fitted QCD:  %9.4f +- %9.4f\n",obsQCD,obsQCDerr);
    printf("Fitted Z:    %9.4f +- %9.4f\n",obsZ,errZ);
    printf("Fitted ELWK: %9.4f +- %9.4f\n",obsEWK,errEWK);
    printf("Tot data:    %9.4f\n",sumDATA);
    printf("Obs data:    %9.4f\n",sumDATA);
    printf("fQCD13          :  %9.4f\n",fQCD13);
    printf("fTAU13          :  %9.4f\n",fTAU13);
    
    
    QCD_norm = obsQCD;
    QCD_err = obsQCDerr;
    QCD_errL = obsQCD*QCD_errL/fitQCD;
    QCD_errR = obsQCD*QCD_errR/fitQCD;
    QCD_norm13 = obsQCD*fQCD13;
    QCD_err13 = QCD_err*fQCD13;
    
    Ztt_norm = obsZ;
    Ztt_err = errZ;
    Ztt_norm13 = obsZ*fTAU13;
    Ztt_err13 = Ztt_err*fTAU13;
    EW_norm = obsEWK;
    EW_err = errEWK;
    H_norm = obsH;
    EW_norm13 = obsEWK*fELW13;
    H_norm13 = obsH*fHIG13;
    
    
     return 0;

}

