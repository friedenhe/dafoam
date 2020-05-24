/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Python wrapper library for DASolverCompressible

\*---------------------------------------------------------------------------*/

#ifndef DASolverCompressible_H
#define DASolverCompressible_H

#include <petscksp.h>
#include "Python.h"
#include "DASolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASolverCompressible Declaration
\*---------------------------------------------------------------------------*/

class DASolverCompressible
{

private:
    /// Disallow default bitwise copy construct
    DASolverCompressible(const DASolverCompressible&);

    /// Disallow default bitwise assignment
    void operator=(const DASolverCompressible&);

    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

    /// DASolver
    autoPtr<DASolver> DASolverPtr_;

public:
    // Constructors

    /// Construct from components
    DASolverCompressible(
        char* argsAll,
        PyObject* pyOptions);

    /// Destructor
    virtual ~DASolverCompressible();

    /// initialize fields and variables
    void initSolver()
    {
        DASolverPtr_->initSolver();
    }

    /// solve the primal equations
    void solvePrimal()
    {
        DASolverPtr_->solvePrimal();
    }

    /// solve the adjoint linear equations
    void solveAdjoint()
    {
        DASolverPtr_->solveAdjoint();
    }

    /// compute the total derivatives
    void calcTotalDerivs()
    {
        DASolverPtr_->calcTotalDerivs();
    }

    /// basically, we call DAIndex::getGlobalXvIndex
    label getGlobalXvIndex(
        const label idxPoint,
        const label idxCoord) const
    {
        return DASolverPtr_->getGlobalXvIndex(idxPoint, idxCoord);
    }
    
    /// basically, we call DASolver::ofField2StateVec
    void ofField2StateVec(Vec stateVec) const
    {
        return DASolverPtr_->ofField2StateVec(stateVec);
    }

    /// basically, we call DASolver::stateVec2OFField
    void stateVec2OFField(const Vec stateVec) const
    {
        return DASolverPtr_->stateVec2OFField(stateVec);
    }

    /// return the number of local adjoint states
    label getNLocalAdjointStates() const
    {
        return DASolverPtr_->getNLocalAdjointStates();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
