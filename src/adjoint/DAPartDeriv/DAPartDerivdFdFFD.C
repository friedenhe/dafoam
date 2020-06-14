/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAPartDerivdFdFFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAPartDerivdFdFFD, 0);
addToRunTimeSelectionTable(DAPartDeriv, DAPartDerivdFdFFD, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPartDerivdFdFFD::DAPartDerivdFdFFD(
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

void DAPartDerivdFdFFD::initializePartDerivMat(
    const dictionary& options,
    Mat* jacMat)
{
    label nDesignVars = options.getLabel("nDesignVars");

    // create dFdFFD
    MatCreate(PETSC_COMM_WORLD, jacMat);
    MatSetSizes(
        *jacMat,
        PETSC_DECIDE,
        PETSC_DECIDE,
        1,
        nDesignVars);
    MatSetFromOptions(*jacMat);
    MatMPIAIJSetPreallocation(*jacMat, nDesignVars, NULL, nDesignVars, NULL);
    MatSeqAIJSetPreallocation(*jacMat, nDesignVars, NULL);
    //MatSetOption(jacMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*jacMat);
    MatZeroEntries(*jacMat);
    Info << "Partial deriative matrix created. " << mesh_.time().elapsedClockTime() << " s" << endl;
}

void DAPartDerivdFdFFD::calcPartDerivMat(
    const dictionary& options,
    const Vec xvVec,
    const Vec wVec,
    Mat jacMat)
{
    // for dFdFFD, we do brute force finite-difference

    label nDesignVars = options.getLabel("nDesignVars");

    word objFuncName, objFuncPart;
    dictionary objFuncSubDictPart = options.subDict("objFuncSubDictPart");
    options.readEntry<word>("objFuncName", objFuncName);
    options.readEntry<word>("objFuncPart", objFuncPart);

    autoPtr<DAObjFunc> daObjFunc(
        DAObjFunc::New(
            mesh_,
            daOption_,
            daModel_,
            daIndex_,
            daResidual_,
            objFuncName,
            objFuncPart,
            objFuncSubDictPart)
            .ptr());

    // zero all the matrices
    MatZeroEntries(jacMat);

    dictionary mOptions;
    mOptions.set("updateState", 1);
    mOptions.set("updateMesh", 1);
    scalar fRef = daObjFunc->masterFunction(mOptions, xvVec, wVec);

    scalar delta = daOption_.getOption<scalar>("adjEpsDerivFFD");
    scalar rDelta = 1.0 / delta;

    Vec xvVecNew;
    VecDuplicate(xvVec, &xvVecNew);
    VecZeroEntries(xvVecNew);

    for (label i = 0; i < nDesignVars; i++)
    {

        // perturb FFD
        VecZeroEntries(xvVecNew);
        MatGetColumnVector(dXvdFFDMat_, xvVecNew, i);
        VecAXPY(xvVecNew, 1.0, xvVec);

        // compute object
        scalar fNew = daObjFunc->masterFunction(mOptions, xvVecNew, wVec);

        // no need to reset FFD here

        scalar partDeriv = (fNew - fRef) * rDelta;

        MatSetValue(jacMat, 0, i, partDeriv, INSERT_VALUES);
    }

    // reset the perturbation to the original pointsField in OpenFOAM
    daObjFunc->masterFunction(mOptions, xvVec, wVec);

    MatAssemblyBegin(jacMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacMat, MAT_FINAL_ASSEMBLY);
}

} // End namespace Foam

// ************************************************************************* //
