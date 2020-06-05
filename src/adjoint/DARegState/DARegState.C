/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DARegState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DARegState, 0);
defineRunTimeSelectionTable(DARegState, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DARegState::DARegState(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel)
{
    /*
    Description:
        Construct from Foam::fvMesh
    Input:
        mesh: a fvMesh object
        daOption: DAOption object
        daModel: DAModel object
    */

    // initialize regStates
    regStates_.set("volScalarStates", {});
    regStates_.set("volVectorStates", {});
    regStates_.set("surfaceScalarStates", {});
    regStates_.set("modelStates", {});

    //Info<<regStates<<endl;
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DARegState> DARegState::New(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
{
    // standard setup for runtime selectable classes

    // look up the solver name defined in system/DADict
    word modelType = daOption.getOption<word>("modelType");

    Info << "Selecting " << modelType << " for DARegState" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DARegState::New"
            "("
            "    const fvMesh&,"
            "    const DAOption&,"
            "    const DAModel&"
            ")")
            << "Unknown DARegState type "
            << modelType << nl << nl
            << "Valid DARegState types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DARegState>(
        cstrIter()(mesh, daOption, daModel));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
