/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAResidualRhoSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualRhoSimpleFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualRhoSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualRhoSimpleFoam::DAResidualRhoSimpleFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualRhoSimpleFoam::calcResiduals(const dictionary& options)
{
}

void DAResidualRhoSimpleFoam::updateIntermediateVariables()
{
}

/// update the boundary condition for all the states in the selected solver
void DAResidualRhoSimpleFoam::correctBoundaryConditions()
{
}

} // End namespace Foam

// ************************************************************************* //
