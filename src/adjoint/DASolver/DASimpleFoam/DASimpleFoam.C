/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASimpleFoam, 0);
addToRunTimeSelectionTable(DASolver, DASimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASimpleFoam::DASimpleFoam(
    char* argsAll,
    PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      simplePtr_(nullptr),
      pPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      laminarTransportPtr_(nullptr),
      turbulencePtr_(nullptr)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DASimpleFoam::initSolver()
{
    /*
    Description:
        Initialize variables for DASolver
    */
    
    Info << "Initializing fields for DASimpleFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
#include "createSimpleControlPython.H"
#include "createFieldsSimple.H"
#include "createAdjointIncompressible.H"
    // initialize checkMesh
    daCheckMeshPtr_.reset(new DACheckMesh(runTime, mesh));
}

label DASimpleFoam::solvePrimal(
    const Vec xvVec,
    Vec wVec)
{
    /*
    Description:
        Call the primal solver to get converged state variables

    Input:
        xvVec: a vector that contains all volume mesh coordinates

    Output:
        wVec: state variable vector
    */

#include "createRefsSimple.H"

    turbulencePtr_->validate();

    Info << "\nStarting time loop\n"
         << endl;

    // deform the mesh based on the xvVec
    this->pointVec2OFMesh(xvVec);

    // check mesh quality
    label meshOK = this->checkMesh();

    if (!meshOK)
    {
        return 1;
    }

    label nSolverIters = 1;
    //while (simple.loop()) // using simple.loop() will have seg fault in parallel
    while (this->loop(runTime))
    {
        if (nSolverIters % 100 == 0 || nSolverIters == 1)
        {
            Info << "Time = " << runTime.timeName() << nl << endl;
        }

        p.storePrevIter();

        // --- Pressure-velocity SIMPLE corrector
        {
#include "UEqnSimple.H"
#include "pEqnSimple.H"
        }

        laminarTransport.correct();
        daTurbulenceModelPtr_->correct();

        if (nSolverIters % 100 == 0 || nSolverIters == 1)
        {
            this->printAllObjFuncs();

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        runTime.write();

        nSolverIters++;
    }

    this->calcPrimalResidualStatistics("print");

    // primal converged, assign the OpenFoam fields to the state vec wVec
    this->ofField2StateVec(wVec);

    // write the mesh to files
    mesh.write();

    Info << "End\n"
         << endl;

    return 0;
}

} // End namespace Foam

// ************************************************************************* //
