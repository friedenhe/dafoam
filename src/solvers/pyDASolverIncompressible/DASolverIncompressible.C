/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/
#include "DASolverIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DASolverIncompressible::DASolverIncompressible(
    char* argsAll,
    PyObject* pyOptions)
    : argsAll_(argsAll),
      pyOptions_(pyOptions),
      DASolverPtr_(nullptr)
{
    DASolverPtr_.reset(DASolver::New(argsAll, pyOptions));
}

DASolverIncompressible::~DASolverIncompressible()
{
}

void DASolverIncompressible::initSolver()
{
    DASolverPtr_->initSolver();
    return;
}

void DASolverIncompressible::solvePrimal()
{
    DASolverPtr_->solvePrimal();
    return;
}

void DASolverIncompressible::solveAdjoint()
{
    DASolverPtr_->solveAdjoint();
    return;
}

void DASolverIncompressible::calcTotalDerivs()
{
    DASolverPtr_->calcTotalDerivs();
    return;
}

label DASolverIncompressible::getGlobalXvIndex(
    const label idxPoint,
    const label idxCoord)
{
    return DASolverPtr_->getGlobalXvIndex(idxPoint, idxCoord);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
