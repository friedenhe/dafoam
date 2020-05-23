/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Solver class that calls primal and adjoint solvers,
        and compute the total derivatives

\*---------------------------------------------------------------------------*/

#ifndef DASolver_H
#define DASolver_H

#include <petscksp.h>
#include "Python.h"
#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "DAUtility.H"
#include "DAOption.H"
#include "DARegState.H"
#include "DAModel.H"
#include "DATurbulenceModel.H"
#include "DAIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DASolver Declaration
\*---------------------------------------------------------------------------*/

class DASolver
{

private:
    /// Disallow default bitwise copy construct
    DASolver(const DASolver&);

    /// Disallow default bitwise assignment
    void operator=(const DASolver&);

protected:
    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

public:
    /// Runtime type information
    TypeName("DASolver");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DASolver,
        dictionary,
        (
            char* argsAll,
            PyObject* pyOptions),
        (argsAll, pyOptions));

    // Constructors

    //- Construct from components
    DASolver(
        char* argsAll,
        PyObject* pyOptions);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DASolver> New(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DASolver()
    {
    }

    // Member functions

    /// initialize fields and variables
    virtual void initSolver() = 0;

    /// solve the primal equations
    virtual void solvePrimal() = 0;

    /// solve the adjoint linear equations
    virtual void solveAdjoint() = 0;

    /// compute the total derivatives
    virtual void calcTotalDerivs() = 0;

    /// basically, we call DAIndex::getGlobalXvIndex
    virtual label getGlobalXvIndex(
        const label idxPoint,
        const label idxCoord) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
