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

DATurbulenceModel::DATurbulenceModel(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DATurbulenceModel",
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

autoPtr<DATurbulenceModel> DATurbulenceModel::New(const fvMesh& mesh)
{
    // look up the solver name 
    const DAOption& daOption = mesh.thisDb().lookupObject<DAOption>("DAOption");
    word solverName = daOption.getOption<word>("turbulenceModel");

    Info << "Selecting " << solverName << " for DATurbulenceModel" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DATurbulenceModel::New"
            "("
            "    const fvMesh&"
            ")")
            << "Unknown DATurbulenceModel type "
            << solverName << nl << nl
            << "Valid DATurbulenceModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<DATurbulenceModel>(
        cstrIter()(mesh));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// this is a virtual function for regIOobject
bool DATurbulenceModel::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
