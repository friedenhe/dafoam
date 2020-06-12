/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAJacConDummy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAJacConDummy, 0);
addToRunTimeSelectionTable(DAJacCon, DAJacConDummy, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAJacConDummy::DAJacConDummy(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAJacCon(modelType, mesh, daOption, daModel, daIndex)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAJacConDummy::setupJacCon(const dictionary& options)
{
}

void DAJacConDummy::initializeJacCon(const dictionary& options)
{
}

} // End namespace Foam

// ************************************************************************* //
