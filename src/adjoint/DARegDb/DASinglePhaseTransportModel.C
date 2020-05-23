/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASinglePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DASinglePhaseTransportModel::DASinglePhaseTransportModel(
    const fvMesh& mesh,
    singlePhaseTransportModel& singlePhaseTransportModel)
    : regIOobject(
        IOobject(
            "DASinglePhaseTransportModel", // always use DASinglePhaseTransportModel for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      singlePhaseTransportModel_(singlePhaseTransportModel)
{
}

DASinglePhaseTransportModel::~DASinglePhaseTransportModel()
{
}

// this is a virtual function for regIOobject
bool DASinglePhaseTransportModel::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
