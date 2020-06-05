/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAJacConRhoSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAJacConRhoSimpleFoam, 0);
addToRunTimeSelectionTable(DAJacCon, DAJacConRhoSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAJacConRhoSimpleFoam::DAJacConRhoSimpleFoam(
    const fvMesh& mesh)
    : DAJacCon(mesh)
{
    /* 
    Adjoint state connectivity info, numbers denote the level of connectivity
    N/A means this state does not connect to the corrsponding residual 

                 U      T      p     nut    phi
    URes         2      2      1      1      0
    TRes         2      2      2      1      0
    pRes         3      2      2      2      1
    phiRes       2      2      1      1      0

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
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "T"} // lv2
        });

    adjStateResidualConInfo_.set(
        "TRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "p", "T"} // lv2
        });

    adjStateResidualConInfo_.set(
        "pRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut", "phi"}, // lv1
            {"U", "p", "T", "nut"}, // lv2
            {"U"} // lv3
        });

    adjStateResidualConInfo_.set(
        "phiRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "T"}, // lv2
        });

    const DAModel& daModel = mesh.thisDb().lookupObject<DAModel>("DAModel");

    // need to correct connectivity for physical models for each residual
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["URes"]);
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["TRes"]);
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["pRes"]);
    daModel.correctAdjStateResidualModelCon(adjStateResidualConInfo_["phiRes"]);

    // add physical model residual connectivity
    daModel.addAdjModelResidualCon(adjStateResidualConInfo_);

    //Info<<adjStateResidualConInfo_<<endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
