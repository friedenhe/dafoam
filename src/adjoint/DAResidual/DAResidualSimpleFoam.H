/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DAResidualSimpleFoam_H
#define DAResidualSimpleFoam_H

#include "DAResidual.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleControl.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAResidualSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAResidualSimpleFoam
    : public DAResidual
{
protected:
    /// \name These are state variables, reference state variables, flow residuals, and partial derivatives
    //@{
    volVectorField& U_;
    volVectorField URef_;
    volVectorField URes_;
    volVectorField UResRef_;
    volVectorField UResPartDeriv_;

    volScalarField& p_;
    volScalarField pRef_;
    volScalarField pRes_;
    volScalarField pResRef_;
    volScalarField pResPartDeriv_;

    surfaceScalarField& phi_;
    surfaceScalarField phiRef_;
    surfaceScalarField phiRes_;
    surfaceScalarField phiResRef_;
    surfaceScalarField phiResPartDeriv_;
    //@}

    DATurbulenceModel& daTurb_;

    simpleControl simple_;

public:
    TypeName("DASimpleFoam");
    // Constructors

    //- Construct from components
    DAResidualSimpleFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAResidualSimpleFoam()
    {
    }

    // Members

    virtual void calcResiduals(const dictionary& options);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
