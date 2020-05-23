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
      runTimePtr_(nullptr),
      meshPtr_(nullptr),
      simplePtr_(nullptr),
      pPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      laminarTransportPtr_(nullptr),
      turbulencePtr_(nullptr),
      daUtilPtr_(nullptr),
      daOptionPtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      daModelPtr_(nullptr),
      daRegStatePtr_(nullptr),
      daIndexPtr_(nullptr)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DASimpleFoam::initSolver()
{
    Info << "Initializing solvers" << endl;
#include "setArgs.H"
#include "setRootCasePython.H"
#include "createTimePython.H"
#include "createMeshPython.H"
#include "createSimpleControlPython.H"
#include "createFieldsSimple.H"
#include "createAdjointSimple.H"
}

void DASimpleFoam::solvePrimal()
{
#include "createRefsSimple.H"

    turbulencePtr_->validate();

    Info << "\nStarting time loop\n"
         << endl;

    //while (simple.loop()) // using simple.loop() will have seg fault in parallel
    while (this->loop(runTime))
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        p.storePrevIter();

        // --- Pressure-velocity SIMPLE corrector
        {
#include "UEqnSimple.H"
#include "pEqnSimple.H"
        }

        laminarTransport.correct();
        turbulencePtr_->correct();

        //runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;
}
void DASimpleFoam::solveAdjoint()
{
}
void DASimpleFoam::calcTotalDerivs()
{
}

/// basically, we call DAIndex::getGlobalXvIndex
label DASimpleFoam::getGlobalXvIndex(
    const label idxPoint,
    const label idxCoord) const
{
    return daIndexPtr_->getGlobalXvIndex(idxPoint, idxCoord);
}

} // End namespace Foam

// ************************************************************************* //
