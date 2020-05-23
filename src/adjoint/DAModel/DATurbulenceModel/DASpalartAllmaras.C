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

void DASpalartAllmaras::correctNut()
{
    Info << "In SA correctNut" << endl;
}

void DASpalartAllmaras::correctModelStates(wordList& modelStates)
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
