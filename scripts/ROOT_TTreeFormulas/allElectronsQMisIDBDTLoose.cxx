#include "TMath.h"

Int_t allElectronsQMisIDBDTLoose( Int_t dilep_type, Int_t lep0_flavour, Int_t lep1_flavour, Float_t lep0_chargeIDBDTLoose, Float_t lep1_chargeIDBDTLoose, Float_t scoreBDT )
{
    if ( dilep_type == 1 ) { return 1; }
    else if ( dilep_type == 2 ) {
	if ( TMath::Abs(lep0_flavour) == 11 && lep0_chargeIDBDTLoose > scoreBDT ) { return 1; }
	if ( TMath::Abs(lep1_flavour) == 11 && lep1_chargeIDBDTLoose > scoreBDT ) { return 1; }
    } else if ( dilep_type == 3 ) {
	if ( lep0_chargeIDBDTLoose > scoreBDT && lep1_chargeIDBDTLoose > scoreBDT ) { return 1; }
    }

    return 0;
}
