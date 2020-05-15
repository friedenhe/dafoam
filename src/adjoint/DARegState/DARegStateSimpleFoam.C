/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DARegStateSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DARegStateSimpleFoam, 0);
addToRunTimeSelectionTable(DARegState, DARegStateSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DARegStateSimpleFoam::DARegStateSimpleFoam(
    const fvMesh& mesh)
    : DARegState(mesh)
{
    // Register state variables
    // NOTE: *State Lists are for each solver and *ModelState lists are for selected physical
    // models at runtime, such as turbulence model and radiation models.
    // For *State lists, register the primitive variables for this solver. For model variables,
    // register specific names to *ModelStates, refer to details in model's DA*.H header files,
    // For example, register "nut" for RANS turbulence model to volScalarModelStates,
    // refer to DATurbModel.H. Then in DATurbModel child class, we will call setModelState() to modify
    // "nut" based on selected turbulence model, for SA model, setModelState will just replace "nut"
    // with "nuTilda", for SST model, it will replace "nut" with "k" and append "omega" to
    // volScalarModelStates. In other words, *State lists will not change for this solver while
    // *ModelState lists will be modified based on the selected models at runtime.

    regStates_["volScalarField"].append("p");
    regStates_["volScalarField"].append("nut");
    regStates_["volVectorField"].append("U");
    regStates_["surfaceScalarField"].append("phi");

    Info << "regStates: " << regStates_ << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
