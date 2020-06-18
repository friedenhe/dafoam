/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAModel::DAModel(
    const fvMesh& mesh,
    const DAOption& daOption)
    : mesh_(mesh),
      daOption_(daOption)
{
    // check whether we have register any physical models
    hasTurbulenceModel_ = mesh.thisDb().foundObject<DATurbulenceModel>("DATurbulenceModel");
    hasRadiationModel_ = mesh.thisDb().foundObject<DARadiationModel>("DARadiationModel");
}

DAModel::~DAModel()
{
}

void DAModel::correctModelStates(wordList& modelStates) const
{
    /*
    Update the name in modelStates based on the selected physical model at runtime

    Example:
    -------
    In DAStateInfo, if the modelStates reads:
    
    modelStates = {"nut"}
    
    then for the SA model, calling correctModelStates(wordList& modelStates) will give:

    modelStates={"nuTilda"}
    
    while calling correctModelStates(wordList& modelStates) for the SST model will give 
    
    modelStates={"k","omega"}
    
    We don't udpate the names for the radiation model becasue users are 
    supposed to set modelStates={"G"}
    */

    if (hasTurbulenceModel_)
    {
        const DATurbulenceModel& daTurb =
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
        daTurb.correctModelStates(modelStates);
    }

    // correct radiation regState
    if (hasRadiationModel_)
    {
        // correct nothing because we should have register G for modelStates
    }
}

void DAModel::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    /*
    Update the original variable connectivity for the adjoint state 
    residuals in stateCon. Basically, we modify/add state variables based on the
    original model variables defined in stateCon.

    Input:
    -----
    stateResCon: the connectivity levels for a state residual, defined in DAJacCon

    Example:
    -------

    If stateCon reads:
    stateCon=
    {
        {"U", "p", "nut"},
        {"p"}
    }

    For the SA turbulence model, calling this function for will get a new stateCon
    stateCon=
    {
        {"U", "p", "nuTilda"},
        {"p"}
    }

    For the SST turbulence model, calling this function will give
    stateCon=
    {
        {"U", "p", "k", "omega"},
        {"p", "U"}
    }
    ***NOTE***: we add a extra level of U connectivity because nut is 
    related to grad(U), k, and omega in SST!
    */

    // correct turbulence model states
    if (hasTurbulenceModel_)
    {
        const DATurbulenceModel& daTurb =
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
        daTurb.correctStateResidualModelCon(stateCon);
    }

    // correct radiation model states
    if (hasRadiationModel_)
    {
        // correct nothing because we should have register G for modelStates
    }
}

/// add the model residual connectivity to stateCon
void DAModel::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    /*
    Add the connectivity levels for all physical model residuals to allCon

    Input:
    -----
    allCon: the connectivity levels for all state residual, defined in DAJacCon

    Example:
    -------

    If stateCon reads:
    allCon=
    {
        "URes":
        {
           {"U", "p", "nut"},
           {"p"}
        }
    }

    For the SA turbulence model, calling this function for will get a new stateCon,
    something like this:
    allCon=
    {
        "URes":
        {
           {"U", "p", "nuTilda"},
           {"p"}
        },
        "nuTildaRes": 
        {
            {"U", "phi", "nuTilda"},
            {"U"}
        }
    }

    */

    // add turbulence model state residuals
    if (hasTurbulenceModel_)
    {
        const DATurbulenceModel& daTurb =
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
        daTurb.addModelResidualCon(allCon);
    }

    // add radiation model state residuals
    if (hasRadiationModel_)
    {
        const DARadiationModel& daRadiation =
            mesh_.thisDb().lookupObject<DARadiationModel>("DARadiationModel");
        daRadiation.addModelResidualCon(allCon);
    }
}

const DATurbulenceModel& DAModel::getDATurbulenceModel() const
{
    if (!hasTurbulenceModel_)
    {
        FatalErrorIn("getDATurbulenceModel") << "DATurbulenceModel not found in mesh.thisDb()! "
                                             << abort(FatalError);
    }
    return mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
}

/// calculate the residuals for model state variables
void DAModel::calcResiduals(const dictionary& options)
{
    if (hasTurbulenceModel_)
    {
        DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel"));
        daTurb.calcResiduals(options);
    }

    if (hasRadiationModel_)
    {
    }
}

void DAModel::correctBoundaryConditions()
{
    if (hasTurbulenceModel_)
    {
        DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel"));
        daTurb.correctBoundaryConditions();
    }

    if (hasRadiationModel_)
    {
    }
}

void DAModel::updateIntermediateVariables()
{
    if (hasTurbulenceModel_)
    {
        DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(
            mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel"));
        daTurb.updateIntermediateVariables();
    }

    if (hasRadiationModel_)
    {
    }
}

#ifdef CompressibleFlow
const fluidThermo& DAModel::getThermo() const
{
    if (!hasTurbulenceModel_)
    {
        FatalErrorIn("DATurbulence not found!") << abort(FatalError);
    }

    const DATurbulenceModel& daTurb = mesh_.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
    return daTurb.getThermo();
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
