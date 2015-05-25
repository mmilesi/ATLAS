#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TApplication.h"

using namespace std;

void drawFromTree();

int main () {

  //TApplication *app =new TApplication("app",0,0);

  drawFromTree();
  
  //app->Run();
  
  return 0;
}

void drawFromTree()
{
   TFile *f = new TFile("~/xAOD_vs_DxAOD/DxAODdata/user.mmilesi.HTopMultilep.test12_DxAOD.00205113.physics_Egamma_tree.root.22165149/user.mmilesi.5142513._000001.tree.root");
   TTree* tree = (TTree*)f->Get("physics");
  
   TCanvas *c = new TCanvas();
   c = c;
   TH1F *h_event_number = new TH1F("event_number", "event_number", 45000000, 0, 45000000);

   int evtsel_event_number;
   tree->SetBranchAddress("evtsel_event_number",&evtsel_event_number);

   Int_t nevent = tree->GetEntries();   
   for (Int_t i=0;i<nevent;i++) {
      tree->GetEntry(i);
      h_event_number->Fill(evtsel_event_number);
   }
   h_event_number->Draw();
   c->SaveAs("event_number.png");
}
