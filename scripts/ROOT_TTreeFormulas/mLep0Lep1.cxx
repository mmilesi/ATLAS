#include "TLorentzVector.h"

Float_t mLep0Lep1( Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2 )
{
    TLorentzVector lep1, lep2, sum;

    lep1.SetPtEtaPhiM(pt1,eta1,phi1,0.0);
    lep2.SetPtEtaPhiM(pt2,eta2,phi2,0.0);

    sum = lep1 + lep2;
    return sum.M();
}
