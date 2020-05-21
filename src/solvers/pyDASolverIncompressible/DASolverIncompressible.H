/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Python wrapper library for DASolverIncompressible

\*---------------------------------------------------------------------------*/

#ifndef DASolverIncompressible_H
#define DASolverIncompressible_H

#include <petscksp.h>
#include "Python.h"
#include "DASolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASolverIncompressible Declaration
\*---------------------------------------------------------------------------*/

class DASolverIncompressible
{

private:
    /// Disallow default bitwise copy construct
    DASolverIncompressible(const DASolverIncompressible&);

    /// Disallow default bitwise assignment
    void operator=(const DASolverIncompressible&);

    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

    /// DASolver
    autoPtr<DASolver> DASolverPtr_;

public:
    // Constructors

    /// Construct from components
    DASolverIncompressible(
        char* argsAll,
        PyObject* pyOptions);

    /// Destructor
    virtual ~DASolverIncompressible();

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
