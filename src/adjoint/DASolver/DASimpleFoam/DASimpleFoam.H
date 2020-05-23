/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DASimpleFoam_H
#define DASimpleFoam_H

#include "DASolver.H"
#include "addToRunTimeSelectionTable.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DASimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DASimpleFoam
    : public DASolver
{

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

    /// DAUtility pointer
    autoPtr<DAUtility> daUtilPtr_;

    /// DAOption pointer
    autoPtr<DAOption> daOptionPtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// DAModel pointer
    autoPtr<DAModel> daModelPtr_;

    /// DARegState pointer
    autoPtr<DARegState> daRegStatePtr_;

    /// DAIndex pointer
    autoPtr<DAIndex> daIndexPtr_;
    
public:
    TypeName("DASimpleFoam");
    // Constructors

    //- Construct from components
    DASimpleFoam(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DASimpleFoam()
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

    /// basically, we call DAIndex::getGlobalXvIndex
    virtual label getGlobalXvIndex(
        const label idxPoint,
        const label idxCoord) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
