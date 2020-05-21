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
    void initSolver();

    /// solve the primal equations
    void solvePrimal();

    /// solve the adjoint linear equations
    void solveAdjoint();

    /// compute the total derivatives
    void calcTotalDerivs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
