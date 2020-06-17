/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAResidual.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DAResidual, 0);
defineRunTimeSelectionTable(DAResidual, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidual::DAResidual(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel),
      daIndex_(daIndex),
      daField_(mesh, daOption, daModel, daIndex)
{
    /*
    Description:
        Construct from Foam::fvMesh
    Input:
        modelType: the type of model
        mesh: a fvMesh object
        daOption: DAOption object
        daModel: DAModel object
    */
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DAResidual> DAResidual::New(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
{
    // standard setup for runtime selectable classes

    Info << "Selecting " << modelType << " for DAResidual" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DAResidual::New"
            "("
            "    const word,"
            "    const fvMesh&,"
            "    const DAOption&,"
            "    const DAModel&,"
            "    const DAIndex&"
            ")")
            << "Unknown DAResidual type "
            << modelType << nl << nl
            << "Valid DAResidual types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DAResidual>(
        cstrIter()(modelType, mesh, daOption, daModel, daIndex));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAResidual::masterFunction(
    const dictionary& options,
    const Vec xvVec,
    const Vec wVec,
    Vec resVec)
{
    // the master function that compute the residual vector given the state and point vectors

    VecZeroEntries(resVec);

    DAModel& daModel = const_cast<DAModel&>(daModel_);

    label updateState = options.getLabel("updateState");

    label updateMesh = options.getLabel("updateMesh");

    label setResVec = options.getLabel("setResVec");

    if (updateMesh)
    {
        daField_.pointVec2OFMesh(xvVec);
    }

    if (updateState)
    {
        daField_.stateVec2OFField(wVec);

        // now update intermediate states and boundry conditions
        this->correctBoundaryConditions();
        this->updateIntermediateVariables();
        daModel.correctBoundaryConditions();
        daModel.updateIntermediateVariables();
    }

    this->calcResiduals(options);
    daModel.calcResiduals(options);

    if (setResVec)
    {
        // asssign the openfoam residual field to resVec
        daField_.ofResField2ResVec(resVec);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
