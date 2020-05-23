/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAModel::DAModel(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DAModel", // always use DAModel for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh)
{
    // check whether we have register any physical models
    hasTurbulenceModel_ = mesh.thisDb().foundObject<DATurbulenceModel>("DATurbulenceModel");
    hasRadiationModel_ = mesh.thisDb().foundObject<DARadiationModel>("DARadiationModel");
}

DAModel::~DAModel()
{
}

// this is a virtual function for regIOobject
bool DAModel::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

void DAModel::correctModelStates(wordList& modelStates)
{
    /*
    Update the name in modelStates based on the selected physical model at runtime

    Example:
    -------
    In DARegState, if the modelStates reads {"nut"}, then for the SA model,
    calling correctModelStates(wordList& modelStates) for the SA model will give 
    modelStates={"nuTilda"}, while calling correctModelStates(wordList& modelStates)
    for the SST model will give modelStates={"k","omega"}. We don't udpate the names
    for the radiation model becasue users are supposed to set modelStates={"G"}
    */

    if (hasTurbulenceModel_)
    {
        DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel"));
        daTurb.correctModelStates(modelStates);
    }

    // correct radiation regState
    if (hasTurbulenceModel_)
    {
        // correct nothing because we should have register G for modelStates
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
