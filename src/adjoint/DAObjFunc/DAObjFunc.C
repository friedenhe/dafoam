/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAObjFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DAObjFunc, 0);
defineRunTimeSelectionTable(DAObjFunc, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAObjFunc::DAObjFunc(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const word objFuncName,
    const word objFuncPart,
    const dictionary& objFuncDict)
    : mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel),
      daIndex_(daIndex),
      objFuncName_(objFuncName),
      objFuncPart_(objFuncPart),
      objFuncDict_(objFuncDict)
{
    /*
    Description:
        Construct from Foam::fvMesh
    Input:
        mesh: a fvMesh object

        objFuncName: the name of the objective function prescribed by users
        NOTE: this is different from the "type" attribute in the objFuncDict_

        objFuncDict: a dictionary that contains information for this objective
    */

    this->calcObjFuncSources(objFuncFaceSources_, objFuncCellSources_);

    // initialize
    objFuncFaceValues_.setSize(objFuncFaceSources_.size());
    forAll(objFuncFaceValues_, idxI)
    {
        objFuncFaceValues_[idxI] = 0.0;
    }

    objFuncCellValues_.setSize(objFuncCellSources_.size());
    forAll(objFuncCellValues_, idxI)
    {
        objFuncCellValues_[idxI] = 0.0;
    }
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DAObjFunc> DAObjFunc::New(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const word objFuncName,
    const word objFuncPart,
    const dictionary& objFuncDict)
{
    // standard setup for runtime selectable classes

    // look up the solver name
    word modelType;
    objFuncDict.readEntry<word>("type", modelType);

    Info << "Selecting type: " << modelType << " for DAObjFunc. Name: " << objFuncName
         << " part: " << objFuncPart << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DAObjFunc::New"
            "("
            "    const fvMesh&,"
            "    const DAOption&,"
            "    const DAModel&,"
            "    const DAIndex&,"
            "    const word,"
            "    const word,"
            "    const dictionary&"
            ")")
            << "Unknown DAObjFunc type "
            << modelType << nl << nl
            << "Valid DAObjFunc types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DAObjFunc>(
        cstrIter()(mesh, daOption, daModel, daIndex, objFuncName, objFuncPart, objFuncDict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAObjFunc::calcObjFuncSources(
    labelList& faceSources,
    labelList& cellSources)
{
    /*
    Compute the face and cell sources for the objective function.
    A typical objFunc dictionary reads:

    {
        "type": "force",
        "source": "patchToFace",
        "patches": ["walls", "wallsbump"],
        "scale": 0.5,
        "addToAdjoint": False
    }

    */

    word objSource;
    objFuncDict_.readEntry("source", objSource);
    if (objSource == "patchToFace")
    {
        // create a topoSet
        autoPtr<topoSet> currentSet(
            topoSet::New(
                "faceSet",
                mesh_,
                "set0",
                IOobject::NO_READ));
        // create the source
        autoPtr<topoSetSource> sourceSet(
            topoSetSource::New(objSource, mesh_, objFuncDict_));

        // add the sourceSet to topoSet
        sourceSet().applyToSet(topoSetSource::NEW, currentSet());
        // get the face index from currentSet, we need to use
        // this special for loop
        for (const label i : currentSet())
        {
            faceSources.append(i);
        }
    }
    else if (objSource == "boxToCell")
    {
        // create a topoSet
        autoPtr<topoSet> currentSet(
            topoSet::New(
                "cellSet",
                mesh_,
                "set0",
                IOobject::NO_READ));
        // we need to change the min and max because they need to
        // be of type point; however, we can't parse point type
        // in pyDict, we need to change them here.
        dictionary objFuncTmp = objFuncDict_;
        scalarList boxMin;
        scalarList boxMax;
        objFuncDict_.readEntry("min", boxMin);
        objFuncDict_.readEntry("max", boxMax);

        point boxMin1;
        point boxMax1;
        boxMin1[0] = boxMin[0];
        boxMin1[1] = boxMin[1];
        boxMin1[2] = boxMin[2];
        boxMax1[0] = boxMax[0];
        boxMax1[1] = boxMax[1];
        boxMax1[2] = boxMax[2];

        objFuncTmp.set("min", boxMin1);
        objFuncTmp.set("max", boxMax1);

        // create the source
        autoPtr<topoSetSource> sourceSet(
            topoSetSource::New(objSource, mesh_, objFuncTmp));

        // add the sourceSet to topoSet
        sourceSet().applyToSet(topoSetSource::NEW, currentSet());
        // get the face index from currentSet, we need to use
        // this special for loop
        for (const label i : currentSet())
        {
            cellSources.append(i);
        }
    }
    else
    {
        FatalErrorIn("calcObjFuncSources") << "source: " << objSource << " not supported!"
                                           << "Options are: patchToFace, boxToCell!"
                                           << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
