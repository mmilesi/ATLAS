#include "TVector2.h"

Float_t deltaPhi( Float_t phi1, Float_t phi2 )
{
    return static_cast<Float_t>( TVector2::Phi_mpi_pi( phi1 - phi2 ) );
}
