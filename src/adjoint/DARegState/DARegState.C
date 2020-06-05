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

DARegState::DARegState(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DARegState",
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh)
{
    /*
    Description:
        Construct from Foam::fvMesh
    Input:
        mesh: a fvMesh object
    */

    // initialize regStates
    regStates_.set("volScalarStates", {});
    regStates_.set("volVectorStates", {});
    regStates_.set("surfaceScalarStates", {});
    regStates_.set("modelStates", {});

    //Info<<regStates<<endl;
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DARegState> DARegState::New(const fvMesh& mesh)
{
    // standard setup for runtime selectable classes

    // look up the solver name defined in system/DADict
    const DAOption& daOption = mesh.thisDb().lookupObject<DAOption>("DAOption");
    word solverName = daOption.getOption<word>("solverName");

    Info << "Selecting " << solverName << " for DARegState" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DARegState::New"
            "("
            "    const fvMesh&"
            ")")
            << "Unknown DARegState type "
            << solverName << nl << nl
            << "Valid DARegState types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DARegState>(
        cstrIter()(mesh));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// this is a virtual function for regIOobject
bool DARegState::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

void DARegState::correctModelStates(wordList& modelStates)
{
    const DAModel& daModel = mesh_.thisDb().lookupObject<DAModel>("DAModel");
    daModel.correctModelStates(modelStates);
    return;
}

/// return the name of pressure field, it can be either p or p_rgh
word DARegState::getPName() const
{
    word pName = "p";
    forAll(regStates_.toc(), idxI)
    {
        word key  = regStates_.toc()[idxI];
        forAll(regStates_[key], idxJ)
        {
            word stateName = regStates_[key][idxJ];
            if(stateName == "p_rgh")
            {
                pName = "p_rgh";
            }
        }
    }
    return pName;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
