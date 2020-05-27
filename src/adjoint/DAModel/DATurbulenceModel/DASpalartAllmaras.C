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

DASpalartAllmaras::DASpalartAllmaras(const fvMesh& mesh)
    : DATurbulenceModel(mesh)
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
