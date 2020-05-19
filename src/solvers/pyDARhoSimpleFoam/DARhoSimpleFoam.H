/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Python wrapper library for DARhoSimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DARhoSimpleFoam_H
#define DARhoSimpleFoam_H

#include <petscksp.h>
#include "Python.h"
#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "DAUtility.H"
#include "DAOption.H"
#include "DARegState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DARhoSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DARhoSimpleFoam
{

private:
    /// Disallow default bitwise copy construct
    DARhoSimpleFoam(const DARhoSimpleFoam&);

    /// Disallow default bitwise assignment
    void operator=(const DARhoSimpleFoam&);

    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

public:
    // Constructors

    /// Construct from components
    DARhoSimpleFoam(
        char* argsAll,
        PyObject* pyOptions);

    /// Destructor
    virtual ~DARhoSimpleFoam();

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
