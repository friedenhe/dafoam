/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DARegStateRhoSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DARegStateRhoSimpleFoam, 0);
addToRunTimeSelectionTable(DARegState, DARegStateRhoSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DARegStateRhoSimpleFoam::DARegStateRhoSimpleFoam(
    const fvMesh& mesh)
    : DARegState(mesh)
{
    // Register state variables
    // NOTE:
    // For model variables, such as turbulence model, register specific names
    // For example, register "nut" to modelStates for RANS turbulence models,
    // Then, we will call correctModelStates(regStates_["modelStates"]) to modify
    // "nut" based on the selected turbulence model. For example, for SA model, 
    // correctModelStates will just replace "nut" with "nuTilda", for SST model,
    // it will replace "nut" with "k" and append "omega" to modelStates. 
    // In other words, the model variables will be modified based on the selected 
    // models at runtime.

    regStates_["volScalarStates"].append("p");
    regStates_["volScalarStates"].append("T");
    regStates_["modelStates"].append("nut");
    regStates_["volVectorStates"].append("U");
    regStates_["surfaceScalarStates"].append("phi");

    // correct the names for model states based on the selected physical model at runtime
    this->correctModelStates(regStates_["modelStates"]);

    Info << "regStates: " << regStates_ << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
