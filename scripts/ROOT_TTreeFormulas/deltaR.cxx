#include "TVector2.h"
#include "TMath.h"

Float_t deltaR( Float_t lepID1, Float_t eta1, Float_t etaBE21, Float_t phi1, Float_t lepID2, Float_t eta2, Float_t etaBE22, Float_t phi2 )
{
    Float_t my_eta1 = ( TMath::Abs(lepID1) == 11 ) ? etaBE21 : eta1;
    Float_t my_eta2 = ( TMath::Abs(lepID2) == 11 ) ? etaBE22 : eta2;
    
    Float_t deta = my_eta1 - my_eta2;
    Float_t dphi = static_cast<Float_t>( TVector2::Phi_mpi_pi( phi1 - phi2 ) );
    return TMath::Sqrt( deta*deta + dphi*dphi );

}
