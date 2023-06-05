/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

\*---------------------------------------------------------------------------*/

#include "DAPostSpanwiseLift.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAPostSpanwiseLift, 0);
addToRunTimeSelectionTable(DAPost, DAPostSpanwiseLift, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPostSpanwiseLift::DAPostSpanwiseLift(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAPost(modelType, mesh, daOption, daModel)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAPostSpanwiseLift::run()
{
    dictionary postDict = daOption_.getAllOptions().subDict("postProcessing");

    wordList patchNames;

    postDict.readEntry<wordList>("patchNames", patchNames);

    scalar aoa = postDict.getScalar("aoa");
    scalar aoaRad = aoa * Foam::constant::mathematical::pi / 180.0;

    label dirMode = postDict.getLabel("dirMode");

    label flowAxisIndex;
    label normalAxisIndex;
    label spanAxisIndex;

    if (dirMode == 1)
    {
        flowAxisIndex = 0;
        normalAxisIndex = 1;
        spanAxisIndex = 2;
    }
    else if (dirMode == 2)
    {
        flowAxisIndex = 0;
        normalAxisIndex = 2;
        spanAxisIndex = 1;
    }

    label nSections = postDict.getLabel("nSections");

    scalar span = postDict.getScalar("span");

    vector liftDir(vector::zero);
    liftDir[flowAxisIndex] = -Foam::sin(aoaRad);
    liftDir[normalAxisIndex] = Foam::cos(aoaRad);

    scalarList spanCoords(nSections + 1);

    scalar deltaZ = span / nSections;
    forAll(spanCoords, idxI)
    {
        spanCoords[idxI] = deltaZ * idxI;
    }

    DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(daModel_.getDATurbulenceModel());

    const volScalarField& p = mesh_.thisDb().lookupObject<volScalarField>("p");

    vector forces(vector::zero);

    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = daTurb.devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb = tdevRhoReff().boundaryField();

    scalarList spanwiseLift(nSections, 0.0);
    scalarList spanwiseArea(nSections, 0.0);

    forAll(patchNames, cI)
    {
        // get the patch id label
        label patchI = mesh_.boundaryMesh().findPatchID(patchNames[cI]);
        // create a shorter handle for the boundary patch
        const fvPatch& patch = mesh_.boundary()[patchI];
        // normal force
        vectorField fN(Sfb[patchI] * p.boundaryField()[patchI]);
        // tangential force
        vectorField fT(Sfb[patchI] & devRhoReffb[patchI]);
        // sum them up
        forAll(patch, faceI)
        {
            forces = fN[faceI] + fT[faceI];
            scalar lift = forces & liftDir;

            vector meshCf = mesh_.Cf().boundaryField()[patchI][faceI];
            scalar meshS = mesh_.magSf().boundaryField()[patchI][faceI];

            scalar spanwiseDist = meshCf[spanAxisIndex];

            for (label i = 0; i < spanCoords.size() - 1; i++)
            {
                if (spanwiseDist > spanCoords[i] && spanwiseDist < spanCoords[i + 1])
                {
                    spanwiseLift[i] += spanwiseDist;
                    spanwiseArea[i] += meshS;
                }
            }
        }
    }

    forAll(spanwiseArea, idxI)
    {
        reduce(spanwiseArea[idxI], sumOp<scalar>());
    }

    forAll(spanwiseLift, idxI)
    {
        spanwiseLift[idxI] /= spanwiseArea[idxI];
    }
}

} // End namespace Foam

// ************************************************************************* //
