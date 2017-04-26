#include "TLorentzVector.h"
#include "TMath.h"

Float_t mTLepMET( Float_t lep_pt, Float_t lep_eta, Float_t lep_phi, Float_t MET_et, Float_t MET_phi)
{
    TLorentzVector lep, MET;

    lep.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,0.00);
    MET.SetPtEtaPhiM(MET_et,0.0,MET_phi,0.0);

    Float_t mT = 2.0 * lep.Pt() * MET.Pt() * ( 1.0 - TMath::Cos(lep.DeltaPhi(MET)) );
    return TMath::Sqrt(mT);
}
