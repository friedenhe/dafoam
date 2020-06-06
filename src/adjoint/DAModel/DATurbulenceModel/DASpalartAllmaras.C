/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASpalartAllmaras.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DASpalartAllmaras, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DASpalartAllmaras, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASpalartAllmaras::DASpalartAllmaras(
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(mesh, daOption)
{
}

void DASpalartAllmaras::correctModelStates(wordList& modelStates) const
{
    // replace nut with nuTilda
    forAll(modelStates, idxI)
    {
        word stateName = modelStates[idxI];
        if (stateName == "nut")
        {
            modelStates[idxI] = "nuTilda";
        }
    }
}

/// update nut based on other turbulence variables and update the BCs
void DASpalartAllmaras::updateNut()
{
}

/// update turbulence variable boundary values
void DASpalartAllmaras::correctTurbBoundaryConditions()
{
}

void DASpalartAllmaras::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    // update the original variable connectivity for the adjoint state residuals in stateCon
    // For SA model just replace nut with nuTilda
    forAll(stateCon, idxI)
    {
        forAll(stateCon[idxI], idxJ)
        {
            word conStateName = stateCon[idxI][idxJ];
            if (conStateName == "nut")
            {
                stateCon[idxI][idxJ] = "nuTilda";
            }
        }
    }
}

void DASpalartAllmaras::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    // add the SA model residual connectivity to stateCon

    word pName;

    if (mesh_.thisDb().foundObject<volScalarField>("p"))
    {
        pName = "p";
    }
    else if (mesh_.thisDb().foundObject<volScalarField>("p_rgh"))
    {
        pName = "p_rgh";
    }
    else
    {
        FatalErrorIn(
            "Neither p nor p_rgh was found in mesh.thisDb()!"
            "addModelResidualCon failed to setup turbulence residuals!")
            << exit(FatalError);
    }

#ifdef IncompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "nuTilda", "phi"}, // lv0
            {"U", "nuTilda"}, // lv1
            {"nuTilda"} // lv2
        });
#endif

#ifdef CompressibleFlow
    allCon.set(
        "nuTildaRes",
        {
            {"U", "T", pName, "nuTilda", "phi"}, // lv0
            {"U", "T", pName, "nuTilda"}, // lv1
            {"T", pName, "nuTilda"} // lv2
        });
#endif
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
