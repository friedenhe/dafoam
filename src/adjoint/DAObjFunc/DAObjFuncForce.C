/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAObjFuncForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAObjFuncForce, 0);
addToRunTimeSelectionTable(DAObjFunc, DAObjFuncForce, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAObjFuncForce::DAObjFuncForce(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const word objFuncName,
    const word objFuncPart,
    const dictionary& objFuncDict)
    : DAObjFunc(
        mesh,
        daOption,
        daModel,
        objFuncName,
        objFuncPart,
        objFuncDict),
      daTurb_(daModel.getDATurbulenceModel()),
      daIndexPtr_(nullptr)
{

    objFuncDict_.readEntry<word>("type", objFuncType_);

    scalarList dir;
    objFuncDict_.readEntry<scalarList>("direction", dir);
    forceDir_[0] = dir[0];
    forceDir_[1] = dir[1];
    forceDir_[2] = dir[2];

    objFuncDict_.readEntry<scalar>("scale", scale_);

    // initialize daIndex
    daIndexPtr_.reset(new DAIndex(mesh, daOption, daModel));
    
}

/// calculate the value of objective function
void DAObjFuncForce::calcObjFunc(
    const labelList& objFuncFaceSources,
    const labelList& objFuncCellSources,
    scalarList& objFuncFaceValues,
    scalarList& objFuncCellValues,
    scalar& objFuncValue)
{
    /*
    Calculate the force which consist of pressure and viscous components.
    This code is modified based on:
    src/functionObjects/forcces/forces.C

    Input:
    -----
    objFuncFaceSources: List of face source (index) for this objective

    objFuncCellSources: List of cell source (index) for this objective

    Output:
    ------
    objFuncFaceValues: the discrete value of objective for each face source (index). 
    This  will be used for computing df/dw in the adjoint.

    objFuncCellValues: the discrete value of objective on each cell source (index). 
    This will be used for computing df/dw in the adjoint.

    objFuncValue: the sum of objective, reduced across all processsors and scaled by "scale"
    */

    // initialize faceValues to zero
    forAll(objFuncFaceValues, idxI)
    {
        objFuncFaceValues[idxI]=0.0;
    }
    // initialize objFunValue
    objFuncValue = 0.0;

    const objectRegistry& db = mesh_.thisDb();
    const volScalarField& p = db.lookupObject<volScalarField>("p");

    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = daTurb_.devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();

    // calculate discrete force for each objFuncFace
    forAll(objFuncFaceSources, idxI)
    {
        const label& objFuncFaceI = objFuncFaceSources[idxI];
        label bFaceI = objFuncFaceI - daIndexPtr_->nLocalInternalFaces;
        const label patchI = daIndexPtr_->bFacePatchI[bFaceI];
        const label faceI = daIndexPtr_->bFaceFaceI[bFaceI];

        // normal force
        vector fN(Sfb[patchI][faceI] * p.boundaryField()[patchI][faceI]);
        // tangential force
        vector fT(Sfb[patchI][faceI] & devRhoReffb[patchI][faceI]);
        // project the force to forceDir
        objFuncFaceValues[idxI] = scale_ * ((fN + fT) & forceDir_);

        objFuncValue += objFuncFaceValues[idxI];
    }

    // need to reduce the sum of force across all processors
    reduce(objFuncValue, sumOp<scalar>());

    return;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
