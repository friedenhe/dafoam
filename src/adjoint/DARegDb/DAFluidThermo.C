/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAFluidThermo::DAFluidThermo(
    const fvMesh& mesh,
    fluidThermo& fluidThermo)
    : regIOobject(
        IOobject(
            "DAFluidThermo", // always use DAFluidThermo for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      fluidThermo_(fluidThermo)
{
}

DAFluidThermo::~DAFluidThermo()
{
}

// this is a virtual function for regIOobject
bool DAFluidThermo::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
