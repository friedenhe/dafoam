/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAPartDerivdRdW.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAPartDerivdRdW, 0);
addToRunTimeSelectionTable(DAPartDeriv, DAPartDerivdRdW, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPartDerivdRdW::DAPartDerivdRdW(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAJacCon& daJacCon)
    : DAPartDeriv(modelType, mesh, daOption, daModel, daIndex, daJacCon)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAPartDerivdRdW::initializePartDeriv(
    Mat* jacMat,
    const dictionary& options)
{
    label transposed = 0;
    options.readEntry<label>("transposed", transposed);

    // now initialize the memory for the jacobian itself
    label localSize = daIndex_.nLocalAdjointStates;

    // create dRdWTPC
    MatCreate(PETSC_COMM_WORLD, jacMat);
    MatSetSizes(
        *jacMat,
        localSize,
        localSize,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*jacMat);
    if (daOption_.getOption<label>("useColoring"))
    {
        daJacCon_.preallocatedRdW(*jacMat, transposed);
        //MatSetOption(jacMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    else
    {
        MatMPIAIJSetPreallocation(*jacMat, 2000, NULL, 2000, NULL);
        MatSeqAIJSetPreallocation(*jacMat, 2000, NULL);
        MatSetOption(*jacMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
    MatSetUp(*jacMat); // create dRdWTPC
    Info << "Partial deriative matrix created. " << this->getRunTime() << " s" << endl;
}

void DAPartDerivdRdW::calcPartDeriv(const dictionary& options)
{
    
}

} // End namespace Foam

// ************************************************************************* //
