#include "TMath.h"

Int_t passBabar( Int_t dilep_type, Int_t lep_probe_flavour, Int_t lep0_flavour )
{
  // OF event should be removed if:
  // -) probe is electron, and the leading lepton in the event
  // -) probe is muon, and the leading lepton in the event
  if ( dilep_type == 2 ) {
    if ( TMath::Abs( lep_probe_flavour ) == 11 && TMath::Abs( lep0_flavour ) == 11 ) { return 0; }
    if ( TMath::Abs( lep_probe_flavour ) == 13 && TMath::Abs( lep0_flavour ) == 13 ) { return 0; }
  }

  // Event passes by default if not OF
  return 1;
}
