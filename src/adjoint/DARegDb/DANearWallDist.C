/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DANearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DANearWallDist::DANearWallDist(
    const fvMesh& mesh,
    nearWallDist& nearWallDist)
    : regIOobject(
        IOobject(
            "DANearWallDist", // always use DANearWallDist for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      nearWallDist_(nearWallDist)
{
}

DANearWallDist::~DANearWallDist()
{
}

// this is a virtual function for regIOobject
bool DANearWallDist::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
