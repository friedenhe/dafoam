/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAObjFuncStateErrorNorm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAObjFuncStateErrorNorm, 0);
addToRunTimeSelectionTable(DAObjFunc, DAObjFuncStateErrorNorm, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAObjFuncStateErrorNorm::DAObjFuncStateErrorNorm(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAResidual& daResidual,
    const word objFuncName,
    const word objFuncPart,
    const dictionary& objFuncDict)
    : DAObjFunc(
        mesh,
        daOption,
        daModel,
        daIndex,
        daResidual,
        objFuncName,
        objFuncPart,
        objFuncDict)
{

    // Assign type, this is common for all objectives
    objFuncDict_.readEntry<word>("type", objFuncType_);

    stateName_ = objFuncDict_.getWord("stateName");
    stateRefName_ = objFuncDict_.getWord("stateRefName");
    stateType_ = objFuncDict_.getWord("stateType");
    argVarName_ = objFuncDict_.getWord("argVarName");
    argVarRefVal_ = objFuncDict_.getScalar("argVarRefVal");
    argVarCoeff_ = objFuncDict_.getScalar("argVarCoeff");
    scale_ = objFuncDict_.getScalar("scale");

    // setup the connectivity, this is needed in Foam::DAJacCondFdW
    // this objFunc only depends on the state variable at the zero level cell
    objFuncConInfo_ = {
        {stateName_}}; // level 0
}

/// calculate the value of objective function
void DAObjFuncStateErrorNorm::calcObjFunc(
    const labelList& objFuncFaceSources,
    const labelList& objFuncCellSources,
    scalarList& objFuncFaceValues,
    scalarList& objFuncCellValues,
    scalar& objFuncValue)
{
    /*
    Description:
        Calculate the stateErrorNorm
        f = L2Norm( state-stateRef ) + argVarCoeff_ * L2Norm( argVar - argVarRef )

    Input:
        objFuncFaceSources: List of face source (index) for this objective
    
        objFuncCellSources: List of cell source (index) for this objective

    Output:
        objFuncFaceValues: the discrete value of objective for each face source (index). 
        This  will be used for computing df/dw in the adjoint.
    
        objFuncCellValues: the discrete value of objective on each cell source (index). 
        This will be used for computing df/dw in the adjoint.
    
        objFuncValue: the sum of objective, reduced across all processsors and scaled by "scale"
    */

    // initialize to zero
    forAll(objFuncFaceValues, idxI)
    {
        objFuncCellValues[idxI] = 0.0;
    }
    // initialize objFunValue
    objFuncValue = 0.0;

    const objectRegistry& db = mesh_.thisDb();
    if (stateType_ == "scalar")
    {
        const volScalarField& state = db.lookupObject<volScalarField>(stateName_);
        const volScalarField& stateRef = db.lookupObject<volScalarField>(stateRefName_);
        const volScalarField& argVar = db.lookupObject<volScalarField>(argVarName_);
        forAll(objFuncCellSources, idxI)
        {
            const label& cellI = objFuncCellSources[idxI];
            objFuncCellValues[idxI] =
                scale_ * (sqr(state[cellI] - stateRef[cellI]) + argVarCoeff_ * sqr(argVar[cellI] - argVarRefVal_));
            objFuncValue += objFuncCellValues[idxI];
        }
    }
    else if (stateType_ == "vector")
    {
        const volVectorField& state = db.lookupObject<volVectorField>(stateName_);
        const volVectorField& stateRef = db.lookupObject<volVectorField>(stateRefName_);
        const volScalarField& argVar = db.lookupObject<volScalarField>(argVarName_);
        forAll(objFuncCellSources, idxI)
        {
            const label& cellI = objFuncCellSources[idxI];
            objFuncCellValues[idxI] =
                scale_ * (sqr(mag(state[cellI] - stateRef[cellI])) + argVarCoeff_ * sqr(argVar[cellI] - argVarRefVal_));
            objFuncValue += objFuncCellValues[idxI];
        }
    }

    // need to reduce the sum of force across all processors
    reduce(objFuncValue, sumOp<scalar>());

    return;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
