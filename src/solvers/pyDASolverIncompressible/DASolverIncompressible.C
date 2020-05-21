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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
