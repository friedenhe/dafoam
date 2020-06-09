/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DASolver, 0);
defineRunTimeSelectionTable(DASolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASolver::DASolver(
    char* argsAll,
    PyObject* pyOptions)
    : argsAll_(argsAll),
      pyOptions_(pyOptions),
      argsPtr_(nullptr),
      runTimePtr_(nullptr),
      meshPtr_(nullptr),
      daOptionPtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      daModelPtr_(nullptr),
      daIndexPtr_(nullptr),
      daFieldPtr_(nullptr),
      daCheckMeshPtr_(nullptr),
      daResidualPtr_(nullptr)
{
    // initialize fvMesh and Time object pointer
    Info << "Initializing mesh and runtime for DASolver" << endl;
#include "setArgs.H"
#include "setRootCasePython.H"
#include "createTimePython.H"
#include "createMeshPython.H"
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DASolver> DASolver::New(
    char* argsAll,
    PyObject* pyOptions)
{
    // standard setup for runtime selectable classes

    // look up the solver name defined in pyOptions
    dictionary allOptions;
    DAUtility::pyDict2OFDict(pyOptions, allOptions);
    word modelType;
    allOptions.readEntry<word>("solverName", modelType);

    Info << "Selecting " << modelType << " for DASolver" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DASolver::New"
            "("
            "    char*,"
            "    PyObject*"
            ")")
            << "Unknown DASolver type "
            << modelType << nl << nl
            << "Valid DASolver types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DASolver>(
        cstrIter()(argsAll, pyOptions));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label DASolver::loop(Time& runTime)
{
    const scalar& endTime = runTime.endTime().value();
    const scalar& deltaT = runTime.deltaT().value();
    const scalar t = runTime.timeOutputValue();
    if (t > endTime - 0.5 * deltaT)
    {
        return 0;
    }
    else
    {
        ++runTime;
        return 1;
    }
}

void DASolver::printAllObjFuncs()
{
    /*
    Calculate the values of all objective functions and print them to screen
    NOTE: we need to call DASolver::setDAObjFuncList before calling this function!
    */

    if (daObjFuncPtrList_.size() == 0)
    {
        FatalErrorIn("printAllObjFuncs") << "daObjFuncPtrList_.size() ==0... "
                                         << "Forgot to call setDAObjFuncList?"
                                         << abort(FatalError);
    }

    forAll(daObjFuncPtrList_, idxI)
    {
        DAObjFunc& daObjFunc = daObjFuncPtrList_[idxI];
        Info << daObjFunc.getObjFuncName()
             << "-" << daObjFunc.getObjFuncPart()
             << "-" << daObjFunc.getObjFuncType()
             << ": " << daObjFunc.getObjFuncValue() << endl;
    }
}

scalar DASolver::getObjFuncValue(const word objFuncName)
{
    /*
    Return the value of the objective function.
    NOTE: we will sum up all the parts in objFuncName

    Input:
    ------
    objFuncName: the name of the objective function

    Output:
    ------
    objFuncValue: the value of the objective
    */

    if (daObjFuncPtrList_.size() == 0)
    {
        FatalErrorIn("printAllObjFuncs") << "daObjFuncPtrList_.size() ==0... "
                                         << "Forgot to call setDAObjFuncList?"
                                         << abort(FatalError);
    }

    scalar objFuncValue = 0.0;

    forAll(daObjFuncPtrList_, idxI)
    {
        DAObjFunc& daObjFunc = daObjFuncPtrList_[idxI];
        if (daObjFunc.getObjFuncName() == objFuncName)
        {
            objFuncValue += daObjFunc.getObjFuncValue();
        }
    }

    return objFuncValue;
}

void DASolver::setDAObjFuncList()
{
    /*
    NOTE: this function needs to be called before calculating any objective functions

    A typical objFunc dictionary looks like this:

    "objFunc": 
    {
        "func1": 
        {
            "part1": 
            {
                "objFuncName": "force",
                "source": "patchToFace",
                "patch": ["walls", "wallsbump"],
                "scale": 0.5,
                "addToAdjoint": False,
            },
            "part2": 
            {
                "objFuncName": "force",
                "source": "patchToFace",
                "patch": ["wallsbump", "frontandback"],
                "scale": 0.5,
                "addToAdjoint": False,
            },
        },
        "func2": 
        {
            "part1": 
            {
                "objFuncName": "force",
                "source": "patchToFace",
                "patch": ["walls", "wallsbump", "frontandback"],
                "scale": 1.0,
                "addToAdjoint": False,
            }
        },
    }
    */

    const dictionary& allOptions = daOptionPtr_->getAllOptions();

    dictionary objFuncDict = allOptions.subDict("objFunc");

    // loop over all objFuncs and parts and calc the number of
    // DAObjFunc instances we need
    label nObjFuncInstances = 0;
    forAll(objFuncDict.toc(), idxI)
    {
        word objFunI = objFuncDict.toc()[idxI];
        dictionary objFuncSubDict = objFuncDict.subDict(objFunI);
        forAll(objFuncSubDict.toc(), idxJ)
        {
            nObjFuncInstances++;
        }
    }

    daObjFuncPtrList_.setSize(nObjFuncInstances);

    // we need to repeat the loop to initialize the
    // DAObjFunc instances
    label objFuncInstanceI = 0;
    forAll(objFuncDict.toc(), idxI)
    {
        word objFunI = objFuncDict.toc()[idxI];
        dictionary objFuncSubDict = objFuncDict.subDict(objFunI);
        forAll(objFuncSubDict.toc(), idxJ)
        {

            word objPart = objFuncSubDict.toc()[idxJ];
            dictionary objFuncSubDictPart = objFuncSubDict.subDict(objPart);

            fvMesh& mesh = meshPtr_();

            daObjFuncPtrList_.set(
                objFuncInstanceI,
                DAObjFunc::New(
                    mesh,
                    daOptionPtr_(),
                    daModelPtr_(),
                    daIndexPtr_(),
                    daResidualPtr_(),
                    objFunI,
                    objPart,
                    objFuncSubDictPart)
                    .ptr());

            objFuncInstanceI++;
        }
    }
}

void DASolver::reduceStateResConLevel(
    const HashTable<label> maxResConLv4JacPCMat,
    HashTable<List<List<word>>>& stateResConInfo) const
{
    /*
    Reduce the connectivity levels for stateResConInfo
    based on maxResConLv4JacPCMat specified in DAOption

    Input:
    -----
    maxResConLv4JacPCMat: the maximal levels of connectivity for each
    state variable residual

    Output:
    ------
    stateResConInfo: reduced connectivity level.

    Example:
    -------

    If the original stateResConInfo reads:

    stateResConInfo
    {
        "URes":
        {
            {"U", "p", "phi"}, // level 0
            {"U", "p"},        // level 1
            {"U"}              // level 2
        }
    }
    And maxResConLv4JacPCMat in DAOption reads:

    maxResConLv4JacPCMat
    {
        "URes": 1
    }
    
    Then, calling reduceStateResConLevel will give:

    stateResConInfo
    {
        "URes":
        {
            {"U", "p", "phi"}, // level 0
            {"U", "p"},        // level 1
        }
    }

    Note that the level 2 of the connectivity in URes is removed becasue
    "URes"=1 in maxResConLv4JacPCMat

    */

    // if no maxResConLv4JacPCMat is specified, just return;
    if (maxResConLv4JacPCMat.size() == 0)
    {
        Info << "maxResConLv4JacPCMat is empty, just return" << endl;
        return;
    }

    // now check if maxResConLv4JacPCMat has all the maxRes level defined
    // and these max levels are <= stateResConInfo.size()
    forAll(stateResConInfo.toc(), idxJ)
    {
        word key1 = stateResConInfo.toc()[idxJ];
        bool keyFound = false;
        forAll(maxResConLv4JacPCMat.toc(), idxI)
        {
            word key = maxResConLv4JacPCMat.toc()[idxI];
            if (key == key1)
            {
                keyFound = true;
                label maxLv = maxResConLv4JacPCMat[key];
                label maxLv1 = stateResConInfo[key1].size() - 1;
                if (maxLv > maxLv1)
                {
                    FatalErrorIn("") << "maxResConLv4JacPCMat maxLevel"
                                     << maxLv << " for " << key
                                     << " larger than stateResConInfo maxLevel "
                                     << maxLv1 << " for " << key1
                                     << abort(FatalError);
                }
            }
        }
        if (!keyFound)
        {
            FatalErrorIn("") << key1 << " not found in maxResConLv4JacPCMat"
                             << abort(FatalError);
        }
    }

    Info << "Reducing max connectivity level of Jacobian PC Mat to:";
    Info << maxResConLv4JacPCMat << endl;

    // assign stateResConInfo to stateResConInfoBK
    HashTable<List<List<word>>> stateResConInfoBK;
    forAll(stateResConInfo.toc(), idxI)
    {
        word key = stateResConInfo.toc()[idxI];
        stateResConInfoBK.set(key, stateResConInfo[key]);
    }

    // now we can erase stateResConInfo
    stateResConInfo.clearStorage();

    // get the reduced stateResConInfo
    forAll(stateResConInfoBK.toc(), idxI)
    {
        word key = stateResConInfoBK.toc()[idxI];
        label maxConLevel = maxResConLv4JacPCMat[key];
        label conSize = stateResConInfoBK[key].size();
        if (conSize > maxConLevel + 1)
        {
            List<List<word>> conList;
            conList.setSize(maxConLevel + 1);
            for (label i = 0; i <= maxConLevel; i++) // NOTE: it is <=
            {
                conList[i] = stateResConInfoBK[key][i];
            }
            stateResConInfo.set(key, conList);
        }
        else
        {
            stateResConInfo.set(key, stateResConInfoBK[key]);
        }
    }
    //Info<<stateResConInfo<<endl;
}

void DASolver::calcPrimalResidualStatistics(
    const word mode,
    const label writeRes)
{
    if (mode == "print")
    {
        // print the primal residuals to screen
        Info << "Printing Primal Residual Statistics." << endl;
    }
    else
    {
        FatalErrorIn("") << "mode not valid" << abort(FatalError);
    }

    label isPC = 0;
    dictionary options;
    options.set("isPC", isPC);
    daResidualPtr_->calcResiduals(options);
    daModelPtr_->calcResiduals(options);

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        const word stateName = stateInfo_["volVectorStates"][idxI];
        const word resName = stateName + "Res";
        const volVectorField& stateRes = meshPtr_->thisDb().lookupObject<volVectorField>(resName);

        vector vecResMax(0, 0, 0);
        vector vecResNorm2(0, 0, 0);
        vector vecResMean(0, 0, 0);
        forAll(stateRes, cellI)
        {
            vecResNorm2.x() += Foam::pow(stateRes[cellI].x(), 2.0);
            vecResNorm2.y() += Foam::pow(stateRes[cellI].y(), 2.0);
            vecResNorm2.z() += Foam::pow(stateRes[cellI].z(), 2.0);
            vecResMean.x() += fabs(stateRes[cellI].x());
            vecResMean.y() += fabs(stateRes[cellI].y());
            vecResMean.z() += fabs(stateRes[cellI].z());
            if (fabs(stateRes[cellI].x()) > vecResMax.x())
            {
                vecResMax.x() = fabs(stateRes[cellI].x());
            }
            if (fabs(stateRes[cellI].y()) > vecResMax.y())
            {
                vecResMax.y() = fabs(stateRes[cellI].y());
            }
            if (fabs(stateRes[cellI].z()) > vecResMax.z())
            {
                vecResMax.z() = fabs(stateRes[cellI].z());
            }
        }
        vecResMean = vecResMean / stateRes.size();
        reduce(vecResMean, sumOp<vector>());
        vecResMean = vecResMean / Pstream::nProcs();
        reduce(vecResNorm2, sumOp<vector>());
        reduce(vecResMax, maxOp<vector>());
        vecResNorm2.x() = Foam::pow(vecResNorm2.x(), 0.5);
        vecResNorm2.y() = Foam::pow(vecResNorm2.y(), 0.5);
        vecResNorm2.z() = Foam::pow(vecResNorm2.z(), 0.5);
        if (mode == "print")
        {
            Info << stateName << " Residual Norm2: " << vecResNorm2 << endl;
            Info << stateName << " Residual Mean: " << vecResMean << endl;
            Info << stateName << " Residual Max: " << vecResMax << endl;
        }

        if (writeRes)
        {
            stateRes.write();
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        const word stateName = stateInfo_["volScalarStates"][idxI];
        const word resName = stateName + "Res";
        const volScalarField& stateRes = meshPtr_->thisDb().lookupObject<volScalarField>(resName);

        scalar scalarResMax = 0, scalarResNorm2 = 0, scalarResMean = 0;
        forAll(stateRes, cellI)
        {
            scalarResNorm2 += Foam::pow(stateRes[cellI], 2.0);
            scalarResMean += fabs(stateRes[cellI]);
            if (fabs(stateRes[cellI]) > scalarResMax)
                scalarResMax = fabs(stateRes[cellI]);
        }
        scalarResMean = scalarResMean / stateRes.size();
        reduce(scalarResMean, sumOp<scalar>());
        scalarResMean = scalarResMean / Pstream::nProcs();
        reduce(scalarResNorm2, sumOp<scalar>());
        reduce(scalarResMax, maxOp<scalar>());
        scalarResNorm2 = Foam::pow(scalarResNorm2, 0.5);
        if (mode == "print")
        {
            Info << stateName << " Residual Norm2: " << scalarResNorm2 << endl;
            Info << stateName << " Residual Mean: " << scalarResMean << endl;
            Info << stateName << " Residual Max: " << scalarResMax << endl;
        }

        if (writeRes)
        {
            stateRes.write();
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        const word stateName = stateInfo_["modelStates"][idxI];
        const word resName = stateName + "Res";
        const volScalarField& stateRes = meshPtr_->thisDb().lookupObject<volScalarField>(resName);

        scalar scalarResMax = 0, scalarResNorm2 = 0, scalarResMean = 0;
        forAll(stateRes, cellI)
        {
            scalarResNorm2 += Foam::pow(stateRes[cellI], 2.0);
            scalarResMean += fabs(stateRes[cellI]);
            if (fabs(stateRes[cellI]) > scalarResMax)
                scalarResMax = fabs(stateRes[cellI]);
        }
        scalarResMean = scalarResMean / stateRes.size();
        reduce(scalarResMean, sumOp<scalar>());
        scalarResMean = scalarResMean / Pstream::nProcs();
        reduce(scalarResNorm2, sumOp<scalar>());
        reduce(scalarResMax, maxOp<scalar>());
        scalarResNorm2 = Foam::pow(scalarResNorm2, 0.5);
        if (mode == "print")
        {
            Info << stateName << " Residual Norm2: " << scalarResNorm2 << endl;
            Info << stateName << " Residual Mean: " << scalarResMean << endl;
            Info << stateName << " Residual Max: " << scalarResMax << endl;
        }

        if (writeRes)
        {
            stateRes.write();
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        const word stateName = stateInfo_["surfaceScalarStates"][idxI];
        const word resName = stateName + "Res";
        const surfaceScalarField& stateRes = meshPtr_->thisDb().lookupObject<surfaceScalarField>(resName);

        scalar phiResMax = 0, phiResNorm2 = 0, phiResMean = 0;
        forAll(stateRes, faceI)
        {
            phiResNorm2 += Foam::pow(stateRes[faceI], 2.0);
            phiResMean += fabs(stateRes[faceI]);
            if (fabs(stateRes[faceI]) > phiResMax)
                phiResMax = fabs(stateRes[faceI]);
        }
        forAll(stateRes.boundaryField(), patchI)
        {
            forAll(stateRes.boundaryField()[patchI], faceI)
            {
                scalar bPhiRes = stateRes.boundaryField()[patchI][faceI];
                phiResNorm2 += Foam::pow(bPhiRes, 2.0);
                phiResMean += fabs(bPhiRes);
                if (fabs(bPhiRes) > phiResMax)
                    phiResMax = fabs(bPhiRes);
            }
        }
        phiResMean = phiResMean / meshPtr_->nFaces();
        reduce(phiResMean, sumOp<scalar>());
        phiResMean = phiResMean / Pstream::nProcs();
        reduce(phiResNorm2, sumOp<scalar>());
        reduce(phiResMax, maxOp<scalar>());
        phiResNorm2 = Foam::pow(phiResNorm2, 0.5);
        if (mode == "print")
        {
            Info << stateName << " Residual Norm2: " << phiResNorm2 << endl;
            Info << stateName << " Residual Mean: " << phiResMean << endl;
            Info << stateName << " Residual Max: " << phiResMax << endl;
        }

        if (writeRes)
        {
            stateRes.write();
        }
    }

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
