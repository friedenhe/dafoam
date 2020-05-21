/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DARhoSimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DARhoSimpleFoam_H
#define DARhoSimpleFoam_H

#include "DASolver.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DARhoSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DARhoSimpleFoam
    : public DASolver
{

protected:
    /// runTime pointer
    autoPtr<Time> runTimePtr_;

    /// fvMesh pointer
    autoPtr<fvMesh> meshPtr_;

    /// simple pointer
    autoPtr<simpleControl> simplePtr_;

    /// thermal model
    autoPtr<fluidThermo> pThermoPtr_;

    /// pressure field pointer
    autoPtr<volScalarField> pPtr_;

    /// density field pointer
    autoPtr<volScalarField> rhoPtr_;

    /// velocity field pointer
    autoPtr<volVectorField> UPtr_;

    /// surface flux field pointer
    autoPtr<surfaceScalarField> phiPtr_;

    /// pressure control pointer
    autoPtr<pressureControl> pressureControlPtr_;

    /// turbulence pointer
    autoPtr<compressible::turbulenceModel> turbulencePtr_;

    /// initial mass
    dimensionedScalar initialMass_;

    /// continuity error
    scalar cumulativeContErr_ = 0.0;

    /// DAUtility pointer
    autoPtr<DAUtility> daUtilPtr_;

    /// DAOption pointer
    autoPtr<DAOption> daOptionPtr_;

public:
    TypeName("DARhoSimpleFoam");
    // Constructors

    //- Construct from components
    DARhoSimpleFoam(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DARhoSimpleFoam()
    {
    }

    /// initialize fields and variables
    virtual void initSolver();

    /// solve the primal equations
    virtual void solvePrimal();

    /// solve the adjoint linear equations
    virtual void solveAdjoint();

    /// compute the total derivatives
    virtual void calcTotalDerivs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
