/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DARhoSimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DARhoSimpleFoam, 0);
addToRunTimeSelectionTable(DASolver, DARhoSimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DARhoSimpleFoam::DARhoSimpleFoam(
    char* argsAll,
    PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      simplePtr_(nullptr),
      pThermoPtr_(nullptr),
      pPtr_(nullptr),
      rhoPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      pressureControlPtr_(nullptr),
      turbulencePtr_(nullptr),
      initialMass_(dimensionedScalar("initialMass", dimensionSet(1, 0, 0, 0, 0, 0, 0), 0.0))
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DARhoSimpleFoam::initSolver()
{
    Info << "Initializing fields for DARhoSimpleFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
    argList& args = argsPtr_();
#include "createSimpleControlPython.H"
#include "createFieldsRhoSimple.H"
#include "createAdjointRhoSimple.H"
    // initialize checkMesh
    daCheckMeshPtr_.reset(new DACheckMesh(runTime, mesh));
}

label DARhoSimpleFoam::solvePrimal(
    const Vec xvVec,
    Vec wVec)
{
    /*
    Call the primal solver to get converged state variables

    Input:
    -----
    xvVec: a vector that contains all volume mesh coordinates

    Output:
    ------
    wVec: state variable vector
    */

#include "createRefsRhoSimple.H"

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
        rho.storePrevIter();

        // Pressure-velocity SIMPLE corrector
#include "UEqnRhoSimple.H"
#include "EEqnRhoSimple.H"
#include "pEqnRhoSimple.H"

        turbulencePtr_->correct();

        if (nSolverIters % 100 == 0 || nSolverIters == 1)
        {
            this->printAllObjFuncs();
            
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        nSolverIters++;
    }

    // primal converged, assign the OpenFoam fields to the state vec wVec
    this->ofField2StateVec(wVec);

    Info << "End\n"
         << endl;

    return 0;
}
void DARhoSimpleFoam::solveAdjoint()
{
}
void DARhoSimpleFoam::calcTotalDerivs()
{
}

} // End namespace Foam

// ************************************************************************* //
