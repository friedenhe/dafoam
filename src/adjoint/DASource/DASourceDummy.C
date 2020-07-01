/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASourceDummy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASourceDummy, 0);
addToRunTimeSelectionTable(DASource, DASourceDummy, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASourceDummy::DASourceDummy(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DASource(modelType, mesh, daOption, daModel, daIndex)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DASourceDummy::calcSource(volVectorField& source)
{
    /*
    Description:
        compute the actuator disk source term
    */

    // do nothing for the dummy source
}

} // End namespace Foam

// ************************************************************************* //
