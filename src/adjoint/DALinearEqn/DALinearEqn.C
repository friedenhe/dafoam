/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DALinearEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DALinearEqn::DALinearEqn(
    const fvMesh& mesh,
    const DAOption& daOption)
    : mesh_(mesh),
      daOption_(daOption)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DALinearEqn::createMLRKSP(
    const dictionary options,
    const Mat jacMat,
    const Mat jacPCMat,
    KSP* genksp)
{

    PC MLRMasterPC, MLRGlobalPC;
    PC MLRsubpc;
    KSP MLRMasterPCKSP;
    KSP* MLRsubksp;
    // ASM Preconditioner variables
    PetscInt MLRoverlap; // width of subdomain overlap
    PetscInt MLRnlocal, MLRfirst; // number of local subblocks, first local subblock

    // Create linear solver context
    KSPCreate(PETSC_COMM_WORLD, genksp);

    // Set operators. Here the matrix that defines the linear
    // system also serves as the preconditioning matrix.
    KSPSetOperators(*genksp, jacMat, jacPCMat);

    // This code sets up the supplied kspObject in the following
    // specific fashion.
    //
    // The hierarchy of the setup is:
    //  kspObject --> Supplied KSP object
    //  |
    //  --> master_PC --> Preconditioner type set to KSP
    //      |
    //      --> master_PC_KSP --> KSP type set to Richardson with 'globalPreConIts'
    //          |
    //           --> globalPC --> PC type set to 'globalPCType'
    //               |            Usually Additive Schwartz and overlap is set
    //               |            with 'ASMOverlap'. Use 0 to get BlockJacobi
    //               |
    //               --> subKSP --> KSP type set to Richardon with 'LocalPreConIts'
    //                   |
    //                   --> subPC -->  PC type set to 'localPCType'.
    //                                  Usually ILU. 'localFillLevel' is
    //                                  set and 'localMatrixOrder' is used.
    //
    // Note that if globalPreConIts=1 then maser_PC_KSP is NOT created and master_PC=globalPC
    // and if localPreConIts=1 then subKSP is set to preOnly.

    // First, KSPSetFromOptions MUST be called
    KSPSetFromOptions(*genksp);

    // Set GMRES
    // Set the type of solver to GMRES
    KSPType kspObjectType = KSPGMRES;

    KSPSetType(*genksp, kspObjectType);
    // Set the gmres restart
    PetscInt restartGMRES = readLabel(options.lookup("GMRESRestart"));

    KSPGMRESSetRestart(*genksp, restartGMRES);
    // Set the GMRES refinement type
    KSPGMRESSetCGSRefinementType(*genksp, KSP_GMRES_CGS_REFINE_IFNEEDED);

    // Set the preconditioner side
    KSPSetPCSide(*genksp, PC_RIGHT);

    // Set global and local PC iters
    PetscInt globalPreConIts = readLabel(options.lookup("GlobalPCIters"));

    // Since there is an extraneous matMult required when using the
    // richardson precondtiter with only 1 iteration, only use it when we need
    // to do more than 1 iteration.
    if (globalPreConIts > 1)
    {
        // Extract preconditioning context for main KSP solver: (MLRMasterPC)
        KSPGetPC(*genksp, &MLRMasterPC);

        // Set the type of MLRMasterPC to ksp. This lets us do multiple
        // iterations of preconditioner application
        PCSetType(MLRMasterPC, PCKSP);

        // Get the ksp context from MLRMasterPC which is the actual preconditioner:
        PCKSPGetKSP(MLRMasterPC, &MLRMasterPCKSP);

        // MLRMasterPCKSP type will always be of type richardson. If the
        // number  of iterations is set to 1, this ksp object is transparent.
        KSPSetType(MLRMasterPCKSP, KSPRICHARDSON);

        // Important to set the norm-type to None for efficiency.
        KSPSetNormType(MLRMasterPCKSP, KSP_NORM_NONE);

        // Do one iteration of the outer ksp preconditioners. Note the
        // tolerances are unsued since we have set KSP_NORM_NONE
        KSPSetTolerances(MLRMasterPCKSP, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, globalPreConIts);

        // Get the 'preconditioner for MLRMasterPCKSP, called 'MLRGlobalPC'. This
        // preconditioner is potentially run multiple times.
        KSPGetPC(MLRMasterPCKSP, &MLRGlobalPC);
    }
    else
    {
        // Just pull out the pc-object if we are not using kspRichardson
        KSPGetPC(*genksp, &MLRGlobalPC);
    }

    // Set the type of 'MLRGlobalPC'. This will almost always be additive schwartz
    PCSetType(MLRGlobalPC, PCASM);

    // Set the overlap required
    MLRoverlap = readLabel(options.lookup("ASMOverlap"));
    PCASMSetOverlap(MLRGlobalPC, MLRoverlap);

    //label KSPCalcEigen = readLabel(options.lookup("KSPCalcEigen"));
    //if (KSPCalcEigen)
    //{
    //    KSPSetComputeEigenvalues(*genksp, PETSC_TRUE);
    //}

    //Setup the main ksp context before extracting the subdomains
    KSPSetUp(*genksp);

    // Extract the ksp objects for each subdomain
    PCASMGetSubKSP(MLRGlobalPC, &MLRnlocal, &MLRfirst, &MLRsubksp);

    //Loop over the local blocks, setting various KSP options
    //for each block.
    PetscInt localPreConIts = readLabel(options.lookup("LocalPCIters"));
    word matOrdering = word(options.lookup("JacMatReOrdering"));
    PetscInt localFillLevel = readLabel(options.lookup("PCFillLevel"));
    for (PetscInt i = 0; i < MLRnlocal; i++)
    {
        // Since there is an extraneous matMult required when using the
        // richardson precondtiter with only 1 iteration, only use it we need
        // to do more than 1 iteration.
        if (localPreConIts > 1)
        {
            // This 'subksp' object will ALSO be of type richardson so we can do
            // multiple iterations on the sub-domains
            KSPSetType(MLRsubksp[i], KSPRICHARDSON);

            // Set the number of iterations to do on local blocks. Tolerances are ignored.
            KSPSetTolerances(MLRsubksp[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, localPreConIts);

            // Again, norm_type is NONE since we don't want to check error
            KSPSetNormType(MLRsubksp[i], KSP_NORM_NONE);
        }
        else
        {
            KSPSetType(MLRsubksp[i], KSPPREONLY);
        }

        // Extract the preconditioner for subksp object.
        KSPGetPC(MLRsubksp[i], &MLRsubpc);

        // The subpc type will almost always be ILU
        PCType localPCType = PCILU;
        PCSetType(MLRsubpc, localPCType);

        // Set PC factor
        PCFactorSetPivotInBlocks(MLRsubpc, PETSC_TRUE);
        PCFactorSetShiftType(MLRsubpc, MAT_SHIFT_NONZERO);
        PCFactorSetShiftAmount(MLRsubpc, PETSC_DECIDE);

        // Setup the matrix ordering for the subpc object:
        // 'natural':'natural',
        // 'rcm':'rcm',
        // 'nested dissection':'nd' (default),
        // 'one way dissection':'1wd',
        // 'quotient minimum degree':'qmd',
        MatOrderingType localMatrixOrdering;
        if (matOrdering == "natural")
        {
            localMatrixOrdering = MATORDERINGNATURAL;
        }
        else if (matOrdering == "nd")
        {
            localMatrixOrdering = MATORDERINGND;
        }
        else if (matOrdering == "rcm")
        {
            localMatrixOrdering = MATORDERINGRCM;
        }
        else if (matOrdering == "1wd")
        {
            localMatrixOrdering = MATORDERING1WD;
        }
        else if (matOrdering == "qmd")
        {
            localMatrixOrdering = MATORDERINGQMD;
        }
        else
        {
            Info << "matOrdering not known. Using default: nested dissection" << endl;
            localMatrixOrdering = MATORDERINGND;
        }
        PCFactorSetMatOrderingType(MLRsubpc, localMatrixOrdering);

        // Set the ILU parameters
        PCFactorSetLevels(MLRsubpc, localFillLevel);
    }

    // Set the norm to unpreconditioned
    KSPSetNormType(*genksp, KSP_NORM_UNPRECONDITIONED);
    // Setup monitor if necessary:
    if (readLabel(options.lookup("printInfo")))
    {
        KSPMonitorSet(*genksp, myKSPMonitor, this, 0);
    }

    PetscInt maxIts = readLabel(options.lookup("GMRESMaxIters"));
    PetscScalar rtol, atol;
    rtol = readScalar(options.lookup("GMRESRelTol"));
    atol = readScalar(options.lookup("GMRESAbsTol"));
    KSPSetTolerances(*genksp, rtol, atol, PETSC_DEFAULT, maxIts);

    if (readLabel(options.lookup("printInfo")))
    {
        Info << "Solver Type: " << kspObjectType << endl;
        Info << "GMRES Restart: " << restartGMRES << endl;
        Info << "ASM Overlap: " << MLRoverlap << endl;
        Info << "Global PC Iters: " << globalPreConIts << endl;
        Info << "Local PC Iters: " << localPreConIts << endl;
        Info << "Mat ReOrdering: " << matOrdering << endl;
        Info << "ILU PC Fill Level: " << localFillLevel << endl;
        Info << "GMRES Max Iterations: " << maxIts << endl;
        Info << "GMRES Relative Tolerance: " << rtol << endl;
        Info << "GMRES Absolute Tolerance: " << atol << endl;
    }
}

label DALinearEqn::solveLinearEqn(
    const KSP ksp,
    const Vec rhsVec,
    Vec solVec)
{
    Info << "Solving Linear Euqation... "
         << this->getRunTime() << " s" << endl;

    //Solve adjoint
    VecZeroEntries(solVec);

    KSPSolve(ksp, rhsVec, solVec);

    //Print convergence information
    label its;
    KSPGetIterationNumber(ksp, &its);
    PetscScalar finalResNorm;
    KSPGetResidualNorm(ksp, &finalResNorm);
    PetscPrintf(
        PETSC_COMM_WORLD,
        "Main iteration %D KSP Residual norm %14.12e %d s \n",
        its,
        finalResNorm,
        this->getRunTime());
    PetscPrintf(PETSC_COMM_WORLD, "Total iterations %D\n", its);

    VecAssemblyBegin(solVec);
    VecAssemblyEnd(solVec);

    Info << "Solving Lineq Equation... Completed! "
         << this->getRunTime() << " s" << endl;
        
    
    return 0;
}

PetscErrorCode DALinearEqn::myKSPMonitor(
    KSP ksp,
    PetscInt n,
    PetscReal rnorm,
    void* ctx)
{

    /*
      Write the solution vector and residual norm to stdout.
      - PetscPrintf() handles output for multiprocessor jobs
      by printing from only one processor in the communicator.
      - The parallel viewer PETSC_VIEWER_STDOUT_WORLD handles
      data from multiple processors so that the output
      is not jumbled.
    */

    DALinearEqn* daLinearEqn = (DALinearEqn*)ctx;

    PetscInt printFrequency = 10; // residual print frequency
    if (n % printFrequency == 0)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Main iteration %D KSP Residual norm %14.12e %d s\n",
            n,
            rnorm,
            daLinearEqn->getRunTime());
    }
    return 0;
}

label DALinearEqn::getRunTime()
{
    return mesh_.time().elapsedClockTime();
}

} // End namespace Foam

// ************************************************************************* //
