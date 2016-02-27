#include "TVector2.h"
#include "TMath.h"

float deltaR( const float& eta1, const float& phi1, const float& eta2, const float& phi2 )
{

    float deta = eta1 - eta2;
    float dphi = static_cast<float>( TVector2::Phi_mpi_pi( phi1 - phi2 ) );
    return TMath::Sqrt( deta*deta + dphi*dphi );

}
