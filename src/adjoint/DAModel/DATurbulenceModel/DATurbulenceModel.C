/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DATurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DATurbulenceModel, 0);
defineRunTimeSelectionTable(DATurbulenceModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DATurbulenceModel::DATurbulenceModel(
    const fvMesh& mesh,
    const DAOption& daOption)
    : regIOobject(
        IOobject(
            "DATurbulenceModel",
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      daOption_(daOption),
      nut_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("nut"))),
      U_(const_cast<volVectorField&>(
          mesh.thisDb().lookupObject<volVectorField>("U"))),
      phi_(const_cast<surfaceScalarField&>(
          mesh.thisDb().lookupObject<surfaceScalarField>("phi"))),
      phase_(
          IOobject(
              "phase",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false),
          mesh,
          dimensionedScalar("phase", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
          zeroGradientFvPatchScalarField::typeName),
      phaseRhoPhi_(const_cast<surfaceScalarField&>(
          mesh.thisDb().lookupObject<surfaceScalarField>("phi"))),
#ifdef IncompressibleFlow
      laminarTransportPtr_(nullptr),
      turbulencePtr_(nullptr),
      // for incompressible, we use uniform one field for rho
      rho_(
          IOobject(
              "rho",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false),
          mesh,
          dimensionedScalar("rho", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
          zeroGradientFvPatchScalarField::typeName),
#endif
#ifdef CompressibleFlow
      thermoPtr_(nullptr),
      turbulencePtr_(nullptr),
      // for compressible flow, we lookup rho in fvMesh
      rho_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("rho"))),
#endif
      turbDict_(
          IOobject(
              "turbulenceProperties",
              mesh.time().constant(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE,
              false)),
      coeffDict_(turbDict_.subDict("RAS")),
      kMin_(
          dimensioned<scalar>::lookupOrAddToDict(
              "kMin",
              coeffDict_,
              sqr(dimVelocity),
              SMALL)),
      epsilonMin_(
          dimensioned<scalar>::lookupOrAddToDict(
              "epsilonMin",
              coeffDict_,
              kMin_.dimensions() / dimTime,
              SMALL)),
      omegaMin_(
          dimensioned<scalar>::lookupOrAddToDict(
              "omegaMin",
              coeffDict_,
              dimless / dimTime,
              SMALL)),
      nuTildaMin_(
          dimensioned<scalar>::lookupOrAddToDict(
              "nuTildaMin",
              coeffDict_,
              nut_.dimensions(),
              SMALL))
{
    // Now we need to initialize other variables
#ifdef IncompressibleFlow

    // initialize the Prandtl number
    IOdictionary transportProperties(
        IOobject(
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false));
    Pr_ = readScalar(transportProperties.lookup("Pr"));

    // initialize laminarTransportPtr_
    // first check if DARegDbSinglePhaseTransportModel has been registered
    // to fvMesh, if not, print an error
    if (mesh.thisDb().foundObject<DARegDbSinglePhaseTransportModel>(
            "DARegDbSinglePhaseTransportModel"))
    {
        const DARegDbSinglePhaseTransportModel& regDb =
            mesh.thisDb().lookupObject<DARegDbSinglePhaseTransportModel>(
                "DARegDbSinglePhaseTransportModel");
        singlePhaseTransportModel& regDbObj =
            const_cast<singlePhaseTransportModel&>(regDb.getObject());
        laminarTransportPtr_.reset(&regDbObj);
    }
    else
    {
        FatalErrorIn(
            "DARegDbSinglePhaseTransportModel not found in mesh.thisDb()!"
            "Please use DARegDb to register it before initializing DATurbulenceModel!"
            "Check src/adjoint/DARegDb/DARegDbSinglePhaseTransportModel.H")
            << exit(FatalError);
    }

    // initialize turbulencePtr_
    // first check if DARegDbTurbulenceModelIncompressible has been registered
    // to fvMesh, if not, print an error
    if (mesh.thisDb().foundObject<DARegDbTurbulenceModelIncompressible>(
            "DARegDbTurbulenceModelIncompressible"))
    {
        const DARegDbTurbulenceModelIncompressible& regDb =
            mesh.thisDb().lookupObject<DARegDbTurbulenceModelIncompressible>(
                "DARegDbTurbulenceModelIncompressible");
        incompressible::turbulenceModel& regDbObj =
            const_cast<incompressible::turbulenceModel&>(regDb.getObject());
        turbulencePtr_.reset(&regDbObj);
    }
    else
    {
        FatalErrorIn(
            "DARegDbTurbulenceModelIncompressible not found in mesh.thisDb()!"
            "Please use DARegDb to register it before initializing DATurbulenceModel!"
            "Check src/adjoint/DARegDb/DARegDbTurbulenceModelIncompressible.H")
            << exit(FatalError);
    }

#endif

#ifdef CompressibleFlow

    // initialize the Prandtl number
    IOdictionary thermophysicalProperties(
        IOobject(
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false));
    Pr_ = readScalar(
        thermophysicalProperties.subDict("mixture").subDict("transport").lookup("Pr"));

    // initialize thermoPtr_
    // first check if DARegDbFluidThermo has been registered
    // to fvMesh, if not, print an error
    if (mesh.thisDb().foundObject<DARegDbFluidThermo>("DARegDbFluidThermo"))
    {
        const DARegDbFluidThermo& regDb =
            mesh.thisDb().lookupObject<DARegDbFluidThermo>("DARegDbFluidThermo");
        fluidThermo& regDbObj = const_cast<fluidThermo&>(regDb.getObject());
        thermoPtr_.reset(&regDbObj);
    }
    else
    {
        FatalErrorIn(
            "DARegDbSinglePhaseTransportModel not found in mesh.thisDb()!"
            "Please use DARegDb to register it before initializing DATurbulenceModel!"
            "Check src/adjoint/DARegDb/DARegDbSinglePhaseTransportModel.H")
            << exit(FatalError);
    }

    // initialize turbulencePtr_
    // first check if DARegDbTurbulenceModelCompressible has been registered
    // to fvMesh, if not, print an error
    if (mesh.thisDb().foundObject<DARegDbTurbulenceModelCompressible>(
            "DARegDbTurbulenceModelCompressible"))
    {
        const DARegDbTurbulenceModelCompressible& regDb =
            mesh.thisDb().lookupObject<DARegDbTurbulenceModelCompressible>(
                "DARegDbTurbulenceModelCompressible");
        compressible::turbulenceModel& regDbObj =
            const_cast<compressible::turbulenceModel&>(regDb.getObject());
        turbulencePtr_.reset(&regDbObj);
    }
    else
    {
        FatalErrorIn(
            "DARegDbTurbulenceModelCompressible not found in mesh.thisDb()!"
            "Please use DARegDb to register it before initializing DATurbulenceModel!"
            "Check src/adjoint/DARegDb/DARegDbTurbulenceModelCompressible.H")
            << exit(FatalError);
    }

#endif
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DATurbulenceModel> DATurbulenceModel::New(
    const fvMesh& mesh,
    const DAOption& daOption)
{
    // look up the solver name
    word modelType = daOption.getOption<word>("turbulenceModel");

    Info << "Selecting " << modelType << " for DATurbulenceModel" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DATurbulenceModel::New"
            "("
            "    const fvMesh&,"
            "    const DAOption&"
            ")")
            << "Unknown DATurbulenceModel type "
            << modelType << nl << nl
            << "Valid DATurbulenceModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<DATurbulenceModel>(
        cstrIter()(mesh, daOption));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// this is a virtual function for regIOobject
bool DATurbulenceModel::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

tmp<volScalarField> DATurbulenceModel::nuEff() const
{

    return tmp<volScalarField>(
        new volScalarField(
            "nuEff",
            this->getNu() + nut_));
}

tmp<volScalarField> DATurbulenceModel::alphaEff()
{

#ifdef IncompressibleFlow
    const volScalarField& alphat = mesh_.thisDb().lookupObject<volScalarField>("alphat");
    return tmp<volScalarField>(
        new volScalarField(
            "alphaEff",
            this->getAlpha() + alphat));
#endif

#ifdef CompressibleFlow
    const volScalarField& alphat = mesh_.thisDb().lookupObject<volScalarField>("alphat");
    return tmp<volScalarField>(
        new volScalarField(
            "alphaEff",
            thermoPtr_->alphaEff(alphat)));
#endif
}

tmp<volScalarField> DATurbulenceModel::getNu() const
{

#ifdef IncompressibleFlow
    return laminarTransportPtr_->nu();
#endif

#ifdef CompressibleFlow
    return thermoPtr_->mu() / rho_;
#endif
}

tmp<volScalarField> DATurbulenceModel::getAlpha() const
{
    return this->getNu() / Pr_;
}

tmp<Foam::volScalarField> DATurbulenceModel::getMu() const
{

#ifdef CompressibleFlow
    return thermoPtr_->mu();
#else
    FatalErrorIn("flowCondition not valid!") << abort(FatalError);
    return nut_;
#endif
}

tmp<volSymmTensorField> DATurbulenceModel::devRhoReff() const
{

    return tmp<volSymmTensorField>(
        new volSymmTensorField(
            IOobject(
                IOobject::groupName("devRhoReff", U_.group()),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            (-phase_ * rho_ * nuEff()) * dev(twoSymm(fvc::grad(U_)))));
}

tmp<fvVectorMatrix> DATurbulenceModel::divDevRhoReff(
    volVectorField& U)
{

#ifdef IncompressibleFlow
    word divScheme = "div((nuEff*dev2(T(grad(U)))))";
#endif
#ifdef CompressibleFlow
    word divScheme = "div(((rho*nuEff)*dev2(T(grad(U)))))";
#endif

    volScalarField& phase = phase_;
    volScalarField& rho = rho_;

    return (
        -fvm::laplacian(phase * rho * nuEff(), U)
        - fvc::div((phase * rho * nuEff()) * dev2(T(fvc::grad(U))), divScheme));
}

tmp<fvVectorMatrix> DATurbulenceModel::divDevReff(
    volVectorField& U)
{
    return divDevRhoReff(U);
}

void DATurbulenceModel::correctWallDist()
{
    //d_.correct();
    // need to correct turbulence boundary conditions
    // this is because when the near wall distance changes, the nut, omega, epsilon at wall
    // may change if you use wall functions
    this->correctTurbBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
