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
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DARegState(mesh, daOption, daModel)
{
    // Register the names of state variables
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
    regStates_["modelStates"].append("nut");
    regStates_["volVectorStates"].append("U");
    regStates_["surfaceScalarStates"].append("phi");

    // correct the names for model states based on the selected physical model at runtime
    daModel.correctModelStates(regStates_["modelStates"]);

    Info << "regStates: " << regStates_ << endl;

    /* 
    Adjoint state connectivity info, numbers denote the level of connectivity
    N/A means this state does not connect to the corrsponding residual 

                  U      p     nut    phi
     URes         2      1      1      0
     pRes         3      2      2      1
     phiRes       2      1      1      0

    ******************************** NOTE 1 **********************************
    One does not need to specify connectivity for each physical model, set the 
    connectivity for original variables instead. For example, for turbulence models,
    set nut. Then, how is nut connected to the other turbulence states will be 
    set in the DAModel class. This is done by calling correctAdjStateResidualModelCon. 
    For example, for SA model we just replace nut with nuTilda, for SST model, we need 
    to add extract connectivity since nut depends on grad(U), k, and omega. We need
    to do this for other pyhsical models such as radiation models.
    **************************************************************************

    ******************************** NOTE 2 **********************************
    Do not specify physical model connectivity here, because they will be added
    by calling addAdjModelResidualCon. For example, for the SA turbulence
    model, it will add the nuTildaRes to adjStateResidualConInfo_ and setup
    its connectivity automatically.
    **************************************************************************

    */

    adjStateResidualConInfo_.set(
        "URes",
        {
            {"U", "p", "nut", "phi"}, // lv0
            {"U", "p", "nut"}, // lv1
            {"U"} // lv2
        });

    adjStateResidualConInfo_.set(
        "pRes",
        {
            {"U", "p", "nut", "phi"}, // lv0
            {"U", "p", "nut", "phi"}, // lv1
            {"U", "p", "nut"}, // lv2
            {"U"} // lv3
        });

    adjStateResidualConInfo_.set(
        "phiRes",
        {
            {"U", "p", "nut", "phi"}, // lv0
            {"U", "p", "nut"}, // lv1
            {"U"}, // lv2
        });

    // need to correct connectivity for physical models for each residual
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["URes"]);
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["pRes"]);
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["phiRes"]);

    // add physical model residual connectivity
    daModel.addAdjModelResidualCon(adjStateResidualConInfo_);

    Info << adjStateResidualConInfo_ << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
