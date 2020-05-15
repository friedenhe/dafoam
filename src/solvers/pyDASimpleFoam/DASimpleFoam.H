/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Python wrapper library for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DASimpleFoam_H
#define DASimpleFoam_H

#include <petscksp.h>
#include "Python.h"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "DAUtility.H"
#include "DAOption.H"
#include "DARegState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DASimpleFoam
{

private:
    /// Disallow default bitwise copy construct
    DASimpleFoam(const DASimpleFoam&);

    /// Disallow default bitwise assignment
    void operator=(const DASimpleFoam&);

    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

protected:
    /// runTime pointer
    autoPtr<Time> runTimePtr_;

    /// fvMesh pointer
    autoPtr<fvMesh> meshPtr_;

    /// simple pointer
    autoPtr<simpleControl> simplePtr_;

    /// pressure field pointer
    autoPtr<volScalarField> pPtr_;

    /// velocity field pointer
    autoPtr<volVectorField> UPtr_;

    /// surface flux field pointer
    autoPtr<surfaceScalarField> phiPtr_;

    /// laminar transport properties pointer
    autoPtr<singlePhaseTransportModel> laminarTransportPtr_;

    /// turbulence pointer
    autoPtr<incompressible::turbulenceModel> turbulencePtr_;

    /// continuity error
    scalar cumulativeContErr_ = 0.0;

    /// pressure referefence cell id
    label pRefCell_ = 0;

    /// pressure reference value
    scalar pRefValue_ = 0.0;

public:
    // Constructors

    /// Construct from components
    DASimpleFoam(
        char* argsAll,
        PyObject* pyOptions);

    /// Destructor
    virtual ~DASimpleFoam();

    /// initialize all the state and intermedia variables
    void init();

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
