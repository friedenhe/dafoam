/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
        Calculating lift per surface area. This will be used to plot the 
        spanwise lift distribution

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{

    argList::addOption(
        "patchNames",
        "'(inlet)'",
        "List of patch names to compute");

    argList::addOption(
        "liftDir",
        "y",
        "Lift direction, can be either y or z");

    argList::addOption(
        "time",
        "1000",
        "Tme instance to compute");

    argList::addOption(
        "aoa",
        "1.5",
        "angle of attack in degree");

#include "setRootCase.H"
#include "createTime.H"

    scalar time;
    if (args.optionFound("time"))
    {
        time = readScalar(args.optionLookup("time")());
    }
    else
    {
        Info << "time not set! Exit." << endl;
        return 1;
    }
    runTime.setTime(time, 0);

#include "createMesh.H"
#include "createFields.H"

    List<wordRe> patchNames;
    if (args.optionFound("patchNames"))
    {
        patchNames = wordReList(args.optionLookup("patchNames")());
    }
    else
    {
        Info << "patchNames not set! Exit." << endl;
        return 1;
    }

    scalar aoa;
    if (args.optionFound("aoa"))
    {
        aoa = readScalar(args.optionLookup("aoa")());
    }
    else
    {
        Info << "aoa not set! Exit." << endl;
        return 1;
    }
    scalar aoaRad = aoa * Foam::constant::mathematical::pi / 180.0;

    word liftDir1 = "y";
    if (args.optionFound("liftDir"))
    {
        liftDir1 = word(args.optionLookup("liftDir")());
    }
    else
    {
        Info << "liftDir not set! Using default y." << endl;
    }

    vector liftDir(vector::zero);
    liftDir[0] = -Foam::sin(aoaRad);

    if (liftDir1 == "y")
    {
        liftDir[1] = Foam::cos(aoaRad);
    }
    else if (liftDir1 == "z")
    {
        liftDir[2] = Foam::cos(aoaRad);
    }
    else
    {
        Info << "liftDir can be either y or z" << endl;
    }

    Info << "Computing liftPerS: patches = " << patchNames << " aoa = " << aoa
         << " lift direction = " << liftDir << " time = " << time << endl;

    volScalarField liftPerS(
        IOobject(
            "liftPerS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("liftPerS", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
        fixedValueFvPatchScalarField::typeName);

    // this code is pulled from:
    // src/functionObjects/forcces/forces.C
    // modified slightly
    vector forces(vector::zero);

    const surfaceVectorField::Boundary& Sfb = mesh.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh.magSf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = turbulence->devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();

    scalar totalLift = 0;
    forAll(patchNames, cI)
    {
        // get the patch id label
        label patchI = mesh.boundaryMesh().findPatchID(patchNames[cI]);
        // create a shorter handle for the boundary patch
        const fvPatch& patch = mesh.boundary()[patchI];
        // normal force
        vectorField fN(
            Sfb[patchI] * p.boundaryField()[patchI]);
        // tangential force
        vectorField fT(Sfb[patchI] & devRhoReffb[patchI]);
        // sum them up
        forAll(patch, faceI)
        {
            forces.x() = fN[faceI].x() + fT[faceI].x();
            forces.y() = fN[faceI].y() + fT[faceI].y();
            forces.z() = fN[faceI].z() + fT[faceI].z();
            scalar lift = forces & liftDir;
            totalLift += lift;
            liftPerS.boundaryFieldRef()[patchI][faceI] = lift / magSfb[patchI][faceI];
        }
    }
    liftPerS.write();

    Info << "totalLift: " << totalLift << endl;

    Info << "Computing liftPerS.... Completed!" << endl;

    return 0;
}

// ************************************************************************* //
