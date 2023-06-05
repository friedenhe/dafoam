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
#include "fvOptions.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volSymmTensorField> devRhoReff(
    fvMesh& mesh,
    volScalarField& rho,
    volScalarField& nuEff,
    volVectorField& U)
{
    return tmp<volSymmTensorField>(
        new volSymmTensorField(
            IOobject(
                IOobject::groupName("devRhoReff", U.group()),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            (-rho * nuEff) * dev(twoSymm(fvc::grad(U)))));
}

int main(int argc, char* argv[])
{

    // define the input arguments
    argList::addOption(
        "patchNames",
        "'(inlet)'",
        "List of patch names to compute");

    argList::addOption(
        "nSections",
        "50",
        "Number of spanwise sections");

    argList::addOption(
        "time",
        "1000",
        "Tme instance to compute");

    argList::addOption(
        "spanAxisIdx",
        "2",
        "Spanwise axis, can be 0, 1, or 2 for x, y, or z axes, respectively");

#include "setRootCase.H"
#include "createTime.H"

    // read the arguments and change the runTime
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
    //Info << "Setting time to " << time << endl;
    runTime.setTime(time, 0);

#include "createMesh.H"

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

    label nSections = 20;
    if (args.optionFound("nSections"))
    {
        nSections = readLabel(args.optionLookup("nSections")());
    }
    else
    {
        Info << "nSections not set! Using default 20" << endl;
    }

    label spanAxisIdx = 2;
    if (args.optionFound("spanAxisIdx"))
    {
        spanAxisIdx = readLabel(args.optionLookup("spanAxisIdx")());
    }
    else
    {
        Info << "spanAxisIdx not set! Using default 2" << endl;
    }

    Info << "Computing spanwise force " << endl
         << "patches = " << patchNames << " time = " << time << " nSections = " << nSections << " spanAxisIdx = " << spanAxisIdx << endl;

    // Now we can read the flow field for the given time
    Info << "Reading field p" << endl;
    volScalarField p(
        IOobject(
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);

    Info << "Reading field U" << endl;
    volVectorField U(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);

    Info << "Reading field nut" << endl;
    volScalarField nut(
        IOobject(
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);

    autoPtr<volScalarField> rhoPtr;
    autoPtr<volScalarField> nuPtr;

    if (p.dimensions() == dimensionSet(0, 2, -2, 0, 0, 0, 0))
    {
        // incompressible
        rhoPtr.reset(new volScalarField(
            IOobject(
                "rho",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("rho", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
            "zeroGradient"));

        IOdictionary transportProperties(
            IOobject(
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));

        scalar nuVal = transportProperties.getScalar("nu");

        nuPtr.reset(new volScalarField(
            IOobject(
                "nu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("nu", nut.dimensions(), nuVal),
            "zeroGradient"));
    }
    else if (p.dimensions() == dimensionSet(1, -1, -2, 0, 0, 0, 0))
    {
        // compressible
        rhoPtr.reset(new volScalarField(
            IOobject(
                "rho",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE),
            mesh));

        IOdictionary thermophysicalProperties(
            IOobject(
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));

        scalar muVal = thermophysicalProperties.subDict("mixture").subDict("transport").getScalar("mu");

        nuPtr.reset(new volScalarField(
            IOobject(
                "nu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("nu", nut.dimensions(), 1.0),
            "zeroGradient"));

        forAll(nuPtr(), cellI)
        {
            nuPtr()[cellI] = muVal / rhoPtr()[cellI];
        }
        nuPtr().correctBoundaryConditions();
    }
    else
    {
        FatalErrorIn("")
            << "p dimension not valid!"
            << abort(FatalError);
    }

    volScalarField nuEff = nut + nuPtr();

    // now we can find the min/max spanwise location
    scalar minS = 1e10, maxS = -1e10;
    const surfaceVectorField::Boundary& CfB = mesh.Cf().boundaryField();
    forAll(patchNames, cI)
    {
        label patchI = mesh.boundaryMesh().findPatchID(patchNames[cI]);
        forAll(mesh.boundaryMesh()[patchI], faceI)
        {
            if (CfB[patchI][faceI][spanAxisIdx] > maxS)
            {
                maxS = CfB[patchI][faceI][spanAxisIdx];
            }
            if (CfB[patchI][faceI][spanAxisIdx] < minS)
            {
                minS = CfB[patchI][faceI][spanAxisIdx];
            }
        }
    }
    reduce(maxS, maxOp<scalar>());
    reduce(minS, minOp<scalar>());
    scalar span = maxS - minS + 1e-8;

    scalarList spanPoints(nSections + 1);
    scalarList spanCenters(nSections);

    scalar deltaZ = span / nSections;
    forAll(spanPoints, idxI)
    {
        spanPoints[idxI] = minS - 0.5e-8 + deltaZ * idxI;
    }
    forAll(spanCenters, idxI)
    {
        spanCenters[idxI] = minS - 0.5e-8 + 0.5 * deltaZ + deltaZ * idxI;
    }

    // this code is pulled from:
    // src/functionObjects/forcces/forces.C
    // modified slightly

    List<vector> spanwiseForce(nSections, vector::zero);

    const surfaceVectorField::Boundary& Sfb = mesh.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff(mesh, rhoPtr(), nuEff, U);
    const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();

    vector totalForce = vector::zero;
    forAll(patchNames, cI)
    {
        // get the patch id label
        label patchI = mesh.boundaryMesh().findPatchID(patchNames[cI]);

        // create a shorter handle for the boundary patch
        const fvPatch& patch = mesh.boundary()[patchI];
        // normal force
        vectorField fN(Sfb[patchI] * p.boundaryField()[patchI]);
        // tangential force
        vectorField fT(Sfb[patchI] & devRhoReffb[patchI]);
        // sum them up
        forAll(patch, faceI)
        {

            const vector& meshCf = mesh.Cf().boundaryField()[patchI][faceI];
            // const scalar& meshS = mesh.magSf().boundaryField()[patchI][faceI];

            vector force = fN[faceI] + fT[faceI];

            totalForce += force;

            scalar spanwiseDist = meshCf[spanAxisIdx];

            label found = 0;
            forAll(spanwiseForce, i)
            {
                if (spanwiseDist > spanPoints[i] && spanwiseDist < spanPoints[i + 1])
                {
                    spanwiseForce[i] += force;
                    found += 1;
                }
            }
            if (found != 1)
            {
                Info << "found " << found << " spanwiseDist " << spanwiseDist << endl;
            }
        }
    }

    reduce(totalForce, sumOp<vector>());

    Info << "totalForce: " << totalForce << endl;

    forAll(spanwiseForce, idxI)
    {
        reduce(spanwiseForce[idxI], sumOp<vector>());
    }

    forAll(spanwiseForce, idxI)
    {
        for (label i = 0; i < 3; i++)
        {
            spanwiseForce[idxI][i] /= totalForce[i];
        }
    }

    OFstream forceOut("spanwiseForce.txt");

    forAll(spanwiseForce, idxI)
    {
        forceOut << spanCenters[idxI] << " " << spanwiseForce[idxI][0] << " " << spanwiseForce[idxI][1] << " " << spanwiseForce[idxI][2] << endl;
    }

    Info << "Computing spanwise force Done!" << endl;

    return 0;
}

// ************************************************************************* //
