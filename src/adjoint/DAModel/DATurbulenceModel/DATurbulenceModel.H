/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Augmented turbulence model for the adjoint method, including residual 
    calculation functions, etc

    NOTE 1:
    Instead of inheriting from the OpenFOAM turbulence implementation, in 
    RASModel's children, we re-write all the correspndong functions for each 
    turbulence model. This is to avoid using template classes and template 
    functions for all the other classes in DAFoam. The downside is that we 
    need to update all the RASModel's children when upgrading to a new version 
    of OpenFOAM. Hopefully, the turbulence model part does not change too much 
    from version to version so the modification will be minimal. 

    NOTE 2:
    This function will lookup turbulence model object in fvMesh, so make sure
    you initialize a turbulence model BEFORE initialzing DATurbulenceModel

\*---------------------------------------------------------------------------*/

#ifndef DATurbulenceModel_H
#define DATurbulenceModel_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "fvc.H"
#include "fvm.H"
#include "surfaceFields.H"
#include "geometricOneField.H"
#include "wallDist.H"
#include "nearWallDist.H"
#include "zeroGradientFvPatchField.H"
#include "DAOption.H"
#include "DARegState.H"
#include "DARegDbNearWallDist.H"

#ifdef IncompressibleFlow
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "DARegDbTurbulenceModelIncompressible.H"
#include "DARegDbSinglePhaseTransportModel.H"
#endif

#ifdef CompressibleFlow
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "DARegDbTurbulenceModelCompressible.H"
#include "DARegDbFluidThermo.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DATurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class DATurbulenceModel
    : public regIOobject
{

private:
    //- Disallow default bitwise copy construct
    DATurbulenceModel(const DATurbulenceModel&);

    //- Disallow default bitwise assignment
    void operator=(const DATurbulenceModel&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// turbulence viscosity
    volScalarField& nut_;

    /// velocity
    volVectorField& U_;

    /// face flux
    surfaceScalarField& phi_;

    /// phase field
    volScalarField phase_;

    /// phase*phi*density field
    surfaceScalarField& phaseRhoPhi_;

#ifdef IncompressibleFlow
    /// a DARegDb object for looking up laminarTransport_
    const DARegDbSinglePhaseTransportModel& daRegDbTransport_;
    /// laminar viscosity for incompressible flow
    const singlePhaseTransportModel& laminarTransport_;
    /// a DARegDb object for looking up turbulence_
    const DARegDbTurbulenceModelIncompressible& daRegDbTurbIncomp_;
    /// turbulence model object
    const incompressible::turbulenceModel& turbulence_;
    /// density field
    volScalarField rho_;
#endif

#ifdef CompressibleFlow
    /// a DARegDb object for looking up thermo_
    const DARegDbFluidThermo& daRegDbThermo_;
    /// thermo model for compressible flow
    const fluidThermo& thermo_;
    /// a DARegDb object for looking up turbulence_
    const DARegDbTurbulenceModelCompressible& daRegDbTurbComp_;
    /// turbulence model object
    const compressible::turbulenceModel& turbulence_;
    /// density field
    volScalarField& rho_;
#endif

    /// turbulence model property dict
    IOdictionary turbDict_;

    /// turbulence model parameters dict
    dictionary coeffDict_;

    /// Lower limit of k
    dimensionedScalar kMin_;

    /// Lower limit of epsilon
    dimensionedScalar epsilonMin_;

    /// Lower limit for omega
    dimensionedScalar omegaMin_;

    /// Lower limit for nuTilda
    dimensionedScalar nuTildaMin_;

    /// Prandtl number
    scalar Pr_;

public:
    //- Runtime type information
    TypeName("DATurbulenceModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DATurbulenceModel,
        dictionary,
        (const fvMesh& mesh),
        (mesh));

    // Constructors

    //- Construct from components
    DATurbulenceModel(const fvMesh& mesh);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DATurbulenceModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~DATurbulenceModel()
    {
    }

    // Members

    /// update wall distance for d_. Note: y_ will be automatically updated in mesh_ object
    void correctWallDist();

    /// update nut based on other turbulence variables and update the BCs
    virtual void updateNut() = 0;

    /// update the turbulence state for DARegState::regStates_
    virtual void correctModelStates(wordList& modelStates) const = 0;

    /// update turbulence variable boundary values
    virtual void correctTurbBoundaryConditions() = 0;

    /// dev terms
    tmp<volSymmTensorField> devRhoReff() const;

    /// divDev terms
    tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U);

    /// divDev terms
    tmp<fvVectorMatrix> divDevReff(volVectorField& U);

    /// return effective viscosity
    tmp<volScalarField> nuEff() const;

    /// get the nut field
    tmp<volScalarField> getNut()
    {
        return nut_;
    }

    /// return effective thermal diffusivity
    tmp<volScalarField> alphaEff();

    /// get the nu field
    tmp<volScalarField> getNu() const;

    tmp<volScalarField> getAlpha() const;

    /// get the density field
    tmp<volScalarField> getRho()
    {
        return rho_;
    }

    /// get the phase field
    tmp<volScalarField> getPhase()
    {
        return phase_;
    }

    /// get mu
    tmp<Foam::volScalarField> getMu() const;

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
