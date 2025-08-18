
Info << "Reading field p_rgh\n"
     << endl;
p_rghPtr_.reset(
    new volScalarField(
        IOobject(
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volScalarField& p_rgh = p_rghPtr_();

Info << "Reading field U\n"
     << endl;
UPtr_.reset(
    new volVectorField(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volVectorField& U = UPtr_();

#include "createPhiPython.H"
phi.write();

Info << "Reading transportProperties\n"
     << endl;

mixturePtr_.reset(new immiscibleIncompressibleTwoPhaseMixture(U, phi));
immiscibleIncompressibleTwoPhaseMixture& mixture = mixturePtr_();

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U)
rhoPtr_.reset(
    new volScalarField(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        alpha1* rho1 + alpha2 * rho2));
volScalarField& rho = rhoPtr_();
rho.oldTime();
rho.write();

// Mass flux
rhoPhiPtr_.reset(
    new surfaceScalarField(
        IOobject(
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        fvc::interpolate(rho) * phi));

// Construct incompressible turbulence model
turbulencePtr_.reset(
    incompressible::turbulenceModel::New(U, phi, mixture));

#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info << "Calculating field g.h\n"
     << endl;
dimensionedScalar ghRef(
    mag(g.value()) > SMALL
        ? g & (cmptMag(g.value()) / mag(g.value())) * hRef
        : dimensionedScalar("ghRef", g.dimensions() * dimLength, scalar(0)));
ghPtr_.reset(new volScalarField("gh", (g & mesh.C()) - ghRef));
ghfPtr_.reset(new surfaceScalarField("ghf", (g & mesh.Cf()) - ghRef));
volScalarField& gh = ghPtr_();

pPtr_.reset(
    new volScalarField(
        IOobject(
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        p_rgh + rho * gh));
volScalarField& p = pPtr_();

pRefCell_ = 0;
pRefValue_ = 0.0;
setRefCell(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell_,
    pRefValue_);

if (p_rgh.needReference())
{
    p += dimensionedScalar(
        "p",
        p.dimensions(),
        pRefValue_ - getRefCellValue(p, pRefCell_));
    p_rgh = p - rho * gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// initialize the advector pointer
Info << "Initializing isoAdvection" << endl;
advectorPtr_.set(
    new isoAdvectionDF(alpha1, phi, U));