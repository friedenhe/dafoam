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
      runTimePtr_(nullptr),
      meshPtr_(nullptr),
      simplePtr_(nullptr),
      pThermoPtr_(nullptr),
      pPtr_(nullptr),
      rhoPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      pressureControlPtr_(nullptr),
      turbulencePtr_(nullptr),
      initialMass_(dimensionedScalar("initialMass", dimensionSet(1, 0, 0, 0, 0, 0, 0), 0.0)),
      daUtilPtr_(nullptr),
      daOptionPtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      daModelPtr_(nullptr),
      daRegStatePtr_(nullptr),
      daIndexPtr_(nullptr)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DARhoSimpleFoam::initSolver()
{
    Info << "Initializing solvers" << endl;
#include "setArgs.H"
#include "setRootCasePython.H"
#include "createTimePython.H"
#include "createMeshPython.H"
#include "createSimpleControlPython.H"
#include "createFields.H"
#include "createAdjoint.H"
}

void DARhoSimpleFoam::solvePrimal()
{
#include "createRefs.H"

    turbulencePtr_->validate();

    Info << "\nStarting time loop\n"
         << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
#include "UEqn.H"
#include "EEqn.H"
#include "pEqn.H"

        turbulencePtr_->correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n"
         << endl;
}
void DARhoSimpleFoam::solveAdjoint()
{
}
void DARhoSimpleFoam::calcTotalDerivs()
{
}

/// basically, we call DAIndex::getGlobalXvIndex
label DARhoSimpleFoam::getGlobalXvIndex(
    const label idxPoint,
    const label idxCoord) const
{
    return daIndexPtr_->getGlobalXvIndex(idxPoint, idxCoord);
}

} // End namespace Foam

// ************************************************************************* //
