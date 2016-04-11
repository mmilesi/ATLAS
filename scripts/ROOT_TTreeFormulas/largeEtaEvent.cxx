#include "TMath.h"

Int_t largeEtaEvent( Int_t nelectrons, Int_t lep0_flavour, Int_t lep1_flavour, Float_t lep0_etaBE2, Float_t lep1_etaBE2 )
{
  if ( nelectrons == 1 ) {
    if ( TMath::Abs(lep0_flavour) == 11 && TMath::Abs( lep0_etaBE2 ) >= 1.37 ) { return 1; }
    if ( TMath::Abs(lep1_flavour) == 11 && TMath::Abs( lep1_etaBE2 ) >= 1.37 ) { return 1; }
  } else if ( nelectrons == 2 ) {
    if ( TMath::Abs( lep0_etaBE2 ) >= 1.37 || TMath::Abs( lep1_etaBE2 ) >= 1.37 ) { return 1; }
  }
  return 0;
}
