#include "TVector2.h"
#include "TMath.h"

Float_t deltaR( Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2 )
{

    Float_t deta = eta1 - eta2;
    Float_t dphi = static_cast<Float_t>( TVector2::Phi_mpi_pi( phi1 - phi2 ) );
    return TMath::Sqrt( deta*deta + dphi*dphi );

}
