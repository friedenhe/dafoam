/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DARadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DARadiationModel, 0);
defineRunTimeSelectionTable(DARadiationModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DARadiationModel::DARadiationModel(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DARadiationModel",
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh)
{
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DARadiationModel> DARadiationModel::New(const fvMesh& mesh)
{
    // look up the solver name 
    const DAOption& daOption = mesh.thisDb().lookupObject<DAOption>("DAOption");
    word solverName = daOption.getOption<word>("radiationModel");

    Info << "Selecting " << solverName << " for DARadiationModel" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DARadiationModel::New"
            "("
            "    const fvMesh&"
            ")")
            << "Unknown DARadiationModel type "
            << solverName << nl << nl
            << "Valid DARadiationModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<DARadiationModel>(
        cstrIter()(mesh));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// this is a virtual function for regIOobject
bool DARadiationModel::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
