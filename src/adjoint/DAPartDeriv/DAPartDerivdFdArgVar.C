/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAPartDerivdFdArgVar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAPartDerivdFdArgVar, 0);
addToRunTimeSelectionTable(DAPartDeriv, DAPartDerivdFdArgVar, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPartDerivdFdArgVar::DAPartDerivdFdArgVar(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAJacCon& daJacCon,
    const DAResidual& daResidual)
    : DAPartDeriv(modelType,
                  mesh,
                  daOption,
                  daModel,
                  daIndex,
                  daJacCon,
                  daResidual)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAPartDerivdFdArgVar::initializePartDerivMat(
    const dictionary& options,
    Mat* jacMat)
{
    /*
    Description:
        Initialize jacMat
    
    Input:
        options. This is not used
    */

    // create dFdArgVar
    label nCells = mesh_.nCells();
    MatCreate(PETSC_COMM_WORLD, jacMat);
    MatSetSizes(
        *jacMat,
        PETSC_DECIDE,
        nCells,
        1,
        PETSC_DETERMINE);
    MatSetFromOptions(*jacMat);
    MatMPIAIJSetPreallocation(*jacMat, nCells, NULL, nCells, NULL);
    MatSeqAIJSetPreallocation(*jacMat, nCells, NULL);
    //MatSetOption(jacMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*jacMat);
    MatZeroEntries(*jacMat);
    Info << "Partial derivative matrix created. " << mesh_.time().elapsedClockTime() << " s" << endl;
}

void DAPartDerivdFdArgVar::calcPartDerivMat(
    const dictionary& options,
    const Vec xvVec,
    const Vec wVec,
    Mat jacMat)
{
    /*
    Description:
        Compute jacMat. We use the analytical method
    
    Input:
        options.objFuncSubDictPart: the objFunc subDict, obtained from DAOption

        options.objFuncName: the name of the objective

        options.objFuncPart: the part of the objective

        xvVec: the volume mesh coordinate vector

        wVec: the state variable vector
    
    Output:
        jacMat: the partial derivative matrix dFdArgVar to compute
    */

    word objFuncName, objFuncPart;
    dictionary objFuncSubDictPart = options.subDict("objFuncSubDictPart");
    scalar argVarCoeff = objFuncSubDictPart.getScalar("argVarCoeff");
    word argVarName = objFuncSubDictPart.getWord("argVarName");
    const objectRegistry& db = mesh_.thisDb();
    const volScalarField& argVar = db.lookupObject<volScalarField>(argVarName);

    // zero all the matrices
    MatZeroEntries(jacMat);

    forAll(mesh_.cells(), cellI)
    {
        label globalCellI = daIndex_.getGlobalCellIndex(cellI);
        PetscScalar partDeriv = 2.0 * argVarCoeff * (argVar[cellI] - 1.0);
        MatSetValue(jacMat, 0, globalCellI, partDeriv, INSERT_VALUES);
    }

    MatAssemblyBegin(jacMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacMat, MAT_FINAL_ASSEMBLY);
}

} // End namespace Foam

// ************************************************************************* //
