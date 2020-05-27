/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Augmented turbulence model for the adjoint method, including residual 
    calculation functions, etc

    NOTE:
    Instead of inheriting from the OpenFOAM turbulence implementation, in 
    RASModel's children, we re-write all the correspndong functions for each 
    turbulence model. This is to avoid using template classes and template 
    functions for all the other classes in DAFoam. The downside is that we 
    need to update all the RASModel's children when upgrading to a new version 
    of OpenFOAM. Hopefully, the turbulence model part does not change too much 
    from version to version so the modification will be minimal. 

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
#include "turbulenceModel.H"
#include "nearWallDist.H"
#include "zeroGradientFvPatchField.H"
#include "DAOption.H"
#include "DARegState.H"

#ifdef IncompressibleFlow
#include "singlePhaseTransportModel.H"
#endif

#ifdef CompressibleFlow
#include "fluidThermo.H"
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
    /// laminar viscosity for incompressible flow
    singlePhaseTransportModel& laminarTransport_;
    /// density field
    volScalarField rho_;
#endif

#ifdef CompressibleFlow
    /// thermo model for compressible flow
    fluidThermo& thermo_;
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
    virtual void correctTurbBoundaryConditions()=0;

    /// dev terms
    tmp<volSymmTensorField> devRhoReff();
    
    /// divDev terms
    tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U);
    
    /// divDev terms
    tmp<fvVectorMatrix> divDevReff(volVectorField& U);
    
    /// return effective viscosity
    tmp<volScalarField> nuEff();

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
