/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        A modified version of CMULES from
        src/finiteVolume/fvMatrices/solvers/MULES

License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MULESDF.H"

namespace Foam
{
defineTypeNameAndDebug(MULESDF, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::MULESDF::MULESDF(
    fvMesh& mesh)
    : regIOobject(
        IOobject(
            "MULESDF",
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            ))
{
}

void Foam::MULESDF::correct(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsiCorr,
    const scalar psiMax,
    const scalar psiMin)
{
    correct(
        geometricOneField(),
        psi,
        phi,
        phiPsiCorr,
        zeroField(),
        zeroField(),
        psiMax,
        psiMin);
}

void Foam::MULESDF::explicitSolve(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin)
{
    addProfiling(solve, "MULESDF::explicitSolve");

    explicitSolve(
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(),
        zeroField(),
        psiMax,
        psiMin);
}

void Foam::MULESDF::limitSum(UPtrList<scalarField>& phiPsiCorrs)
{
    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumPos = 0;
        scalar sumNeg = 0;

        for (int phasei = 0; phasei < phiPsiCorrs.size(); phasei++)
        {
            if (phiPsiCorrs[phasei][facei] > 0)
            {
                sumPos += phiPsiCorrs[phasei][facei];
            }
            else
            {
                sumNeg += phiPsiCorrs[phasei][facei];
            }
        }

        scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > VSMALL)
        {
            scalar lambda = -sumNeg / sumPos;

            for (int phasei = 0; phasei < phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] > 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -VSMALL)
        {
            scalar lambda = -sumPos / sumNeg;

            for (int phasei = 0; phasei < phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}

// ************************************************************************* //
