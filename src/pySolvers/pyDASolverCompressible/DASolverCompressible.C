/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/
#include "DASolverCompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DASolverCompressible::DASolverCompressible(
    char* argsAll,
    PyObject* pyOptions)
    : argsAll_(argsAll),
      pyOptions_(pyOptions),
      DASolverPtr_(nullptr)
{
    DASolverPtr_.reset(DASolver::New(argsAll, pyOptions));
}

DASolverCompressible::~DASolverCompressible()
{
}

void DASolverCompressible::initSolver()
{
    DASolverPtr_->initSolver();
    return;
}

void DASolverCompressible::solvePrimal()
{
    DASolverPtr_->solvePrimal();
    return;
}

void DASolverCompressible::solveAdjoint()
{
    DASolverPtr_->solveAdjoint();
    return;
}

void DASolverCompressible::calcTotalDerivs()
{
    DASolverPtr_->calcTotalDerivs();
    return;
}

label DASolverCompressible::getGlobalXvIndex(
    const label idxPoint,
    const label idxCoord)
{
    return DASolverPtr_->getGlobalXvIndex(idxPoint, idxCoord);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
