/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAPartDerivdRdArgVar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAPartDerivdRdArgVar, 0);
addToRunTimeSelectionTable(DAPartDeriv, DAPartDerivdRdArgVar, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPartDerivdRdArgVar::DAPartDerivdRdArgVar(
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

void DAPartDerivdRdArgVar::initializePartDerivMat(
    const dictionary& options,
    Mat* jacMat)
{
    /*
    Description:
        Initialize jacMat
    
    Input:
        options. this is not used
    */

    // now initialize the memory for the jacobian itself
    label localSize = daIndex_.nLocalAdjointStates;
    label nCells = mesh_.nCells();

    // create dRdArgVar
    MatCreate(PETSC_COMM_WORLD, jacMat);
    MatSetSizes(
        *jacMat,
        localSize,
        nCells,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*jacMat);
    MatMPIAIJSetPreallocation(*jacMat, 1, NULL, 1, NULL);
    MatSeqAIJSetPreallocation(*jacMat, 1, NULL);
    //MatSetOption(jacMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*jacMat);
    MatZeroEntries(*jacMat);
    Info << "Partial derivative matrix created. " << mesh_.time().elapsedClockTime() << " s" << endl;
}

void DAPartDerivdRdArgVar::calcPartDerivMat(
    const dictionary& options,
    const Vec xvVec,
    const Vec wVec,
    Mat jacMat)
{
    /*
    Description:
        Compute jacMat. We use analytical method
    
    Input:

        options.isPC: whether to compute the jacMat for preconditioner

        xvVec: the volume mesh coordinate vector

        wVec: the state variable vector
    
    Output:
        jacMat: the partial derivative matrix dRdArgVar to compute
    */

    // zero all the matrices
    MatZeroEntries(jacMat);

    scalarList prodTerm(mesh_.nCells());
    daModel_.getTurbProdTerm(prodTerm);

    // dRdBetaSA only has diagonal component for the turbulence residual
    forAll(mesh_.cells(), cellI)
    {
        label globalIndex = daIndex_.getGlobalAdjointStateIndex("nuTilda", cellI);
        scalar val = prodTerm[cellI];
        MatSetValue(jacMat, globalIndex, cellI, val, INSERT_VALUES);
    }

    MatAssemblyBegin(jacMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacMat, MAT_FINAL_ASSEMBLY);
}

} // End namespace Foam

// ************************************************************************* //
