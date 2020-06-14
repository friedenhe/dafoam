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
                "addToAdjoint": True,
            },
            "part2": 
            {
                "objFuncName": "force",
                "source": "patchToFace",
                "patch": ["wallsbump", "frontandback"],
                "scale": 0.5,
                "addToAdjoint": True,
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

label DASolver::solveAdjoint(
    const Vec xvVec,
    const Vec wVec)
{
    // solve the adjoint linear equations

    DALinearEqn daLinearEqn(meshPtr_(), daOptionPtr_());

    // ********************** compute dRdWT **********************
    Mat dRdWT;
    {

        // initialize DAJacCon object
        word modelType = "dRdW";
        autoPtr<DAJacCon> daJacCon(DAJacCon::New(
            modelType,
            meshPtr_(),
            daOptionPtr_(),
            daModelPtr_(),
            daIndexPtr_()));

        dictionary options;
        const HashTable<List<List<word>>>& stateResConInfo = daStateInfoPtr_->getStateResConInfo();
        options.set("stateResConInfo", stateResConInfo);

        // need to first setup preallocation vectors for the dRdWCon matrix
        // because directly initializing the dRdWCon matrix will use too much memory
        daJacCon->setupJacConPreallocation(options);

        // now we can initilaize dRdWCon
        daJacCon->initializeJacCon(options);

        // setup dRdWCon
        daJacCon->setupJacCon(options);
        Info << "dRdWCon Created. " << runTimePtr_->elapsedClockTime() << " s" << endl;

        // read the coloring
        daJacCon->readJacConColoring();

        // initialize partDeriv object
        autoPtr<DAPartDeriv> daPartDeriv(DAPartDeriv::New(
            modelType,
            meshPtr_(),
            daOptionPtr_(),
            daModelPtr_(),
            daIndexPtr_(),
            daJacCon(),
            daResidualPtr_()));

        // we want transposed dRdW
        dictionary options1;
        options1.set("transposed", 1);

        // initilalize dRdWT matrix
        daPartDeriv->initializePartDerivMat(options1, &dRdWT);

        // calculate dRdWT
        daPartDeriv->calcPartDerivMat(options1, xvVec, wVec, dRdWT);

        if (daOptionPtr_->getOption<label>("debug"))
        {
            DAUtility::writeMatrixBinary(dRdWT, "dRdWT");
        }
    }

    // ********************** compute dFdW **********************
    const dictionary& allOptions = daOptionPtr_->getAllOptions();
    dictionary objFuncDict = allOptions.subDict("objFunc");

    // loop over all the objFuncName in the objFunc dictionary
    forAll(objFuncDict.toc(), idxI)
    {
        word objFuncName = objFuncDict.toc()[idxI];

        // the dFdWAllParts vector contains the sum of dFdWVec from all parts for this objFuncName
        Vec dFdWVecAllParts;
        VecDuplicate(wVec, &dFdWVecAllParts);
        VecZeroEntries(dFdWVecAllParts);

        // loop over all parts for this objFuncName
        dictionary objFuncSubDict = objFuncDict.subDict(objFuncName);
        forAll(objFuncSubDict.toc(), idxJ)
        {
            word objFuncPart = objFuncSubDict.toc()[idxJ];
            dictionary objFuncSubDictPart = objFuncSubDict.subDict(objFuncPart);

            // we only solve adjoint for objFuncs with addToAdjoint = True
            label addToAdjoint = objFuncSubDictPart.getLabel("addToAdjoint");
            if (addToAdjoint)
            {
                // append the objFuncName to objFuncNames4Adj_ such that later we can know
                // which objFuncName has addToAdjoint = True
                objFuncNames4Adj_.append(objFuncName);

                // NOTE: dFdW is a matrix here and it has nObjFuncCellSources+nObjFuncFaceSources rows
                Mat dFdW;

                // initialize DAJacCon object
                word modelType = "dFdW";
                autoPtr<DAJacCon> daJacCon(DAJacCon::New(
                    modelType,
                    meshPtr_(),
                    daOptionPtr_(),
                    daModelPtr_(),
                    daIndexPtr_()));

                // initialize objFunc to get objFuncCellSources and objFuncFaceSources
                autoPtr<DAObjFunc> daObjFunc(DAObjFunc::New(
                    meshPtr_(),
                    daOptionPtr_(),
                    daModelPtr_(),
                    daIndexPtr_(),
                    daResidualPtr_(),
                    objFuncName,
                    objFuncPart,
                    objFuncSubDictPart));

                // setup options for daJacCondFdW computation
                dictionary options;
                const List<List<word>>& objFuncConInfo = daObjFunc->getObjFuncConInfo();
                const labelList& objFuncFaceSources = daObjFunc->getObjFuncFaceSources();
                const labelList& objFuncCellSources = daObjFunc->getObjFuncCellSources();
                options.set("objFuncConInfo", objFuncConInfo);
                options.set("objFuncFaceSources", objFuncFaceSources);
                options.set("objFuncCellSources", objFuncCellSources);
                options.set("objFuncName", objFuncName);
                options.set("objFuncPart", objFuncPart);
                options.set("objFuncSubDictPart", objFuncSubDictPart);

                // now we can initilaize dFdWCon
                daJacCon->initializeJacCon(options);

                // setup dFdWCon
                daJacCon->setupJacCon(options);
                Info << "dFdWCon Created. " << meshPtr_->time().elapsedClockTime() << " s" << endl;

                // read the coloring
                word postFix = "_" + objFuncName + "_" + objFuncPart;
                daJacCon->readJacConColoring(postFix);

                // initialize DAPartDeriv to computing dFdW
                autoPtr<DAPartDeriv> daPartDeriv(DAPartDeriv::New(
                    modelType,
                    meshPtr_(),
                    daOptionPtr_(),
                    daModelPtr_(),
                    daIndexPtr_(),
                    daJacCon(),
                    daResidualPtr_()));

                // initialize dFdWMat
                daPartDeriv->initializePartDerivMat(options, &dFdW);

                // compute it
                daPartDeriv->calcPartDerivMat(options, xvVec, wVec, dFdW);

                // now we need to add all the rows of dFdW together to get dFdWVec
                // NOTE: dFdW is a matrix with nObjFuncCellSources+nObjFuncFaceSources rows
                // and nLocalAdjStates columns. So we can do dFdWVec = oneVec*dFdW
                Vec dFdWVec, oneVec;
                label objGeoSize = objFuncFaceSources.size() + objFuncCellSources.size();
                VecCreate(PETSC_COMM_WORLD, &oneVec);
                VecSetSizes(oneVec, objGeoSize, PETSC_DETERMINE);
                VecSetFromOptions(oneVec);
                // assign one to all elements
                VecSet(oneVec, 1.0);
                VecDuplicate(wVec, &dFdWVec);
                VecZeroEntries(dFdWVec);
                // dFdWVec = oneVec*dFdW
                MatMultTranspose(dFdW, oneVec, dFdWVec);

                // we need to add dFdWVec to dFdWVecAllParts because we want to sum all dFdWVec
                // for all parts of this objFuncName. When solving the adjoint equation, we use
                // dFdWVecAllParts
                VecAXPY(dFdWVecAllParts, 1.0, dFdWVec);

                if (daOptionPtr_->getOption<label>("debug"))
                {
                    word outputName = "dFdWVec_" + objFuncName + "_" + objFuncPart;
                    DAUtility::writeVectorBinary(dFdWVec, outputName);
                    DAUtility::writeVectorASCII(dFdWVec, outputName);
                }
            }
        }

        // now we should have add all dFdWVec for dFdWVecAllParts for this objFuncName
        // it will be used as the rhs for adjoint
        if (daOptionPtr_->getOption<label>("debug"))
        {
            word outputName = "dFdWVecAllParts_" + objFuncName;
            DAUtility::writeVectorBinary(dFdWVecAllParts, outputName);
            DAUtility::writeVectorASCII(dFdWVecAllParts, outputName);
        }

        // we only solve adjoint for objectives that have addToAdjoint = True
        if (DAUtility::isInList<word>(objFuncName, objFuncNames4Adj_))
        {

            // NOTE: we use dFdWVecAllParts for adjoint
            // first setup KSP options
            dictionary kspOptions;
            kspOptions.add(
                "GMRESRestart",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "gmresRestart"));
            kspOptions.add(
                "GlobalPCIters",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "globalPCIters"));
            kspOptions.add(
                "ASMOverlap",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "asmOverlap"));
            kspOptions.add(
                "LocalPCIters",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "localPCIters"));
            kspOptions.add(
                "JacMatReOrdering",
                daOptionPtr_->getSubDictOption<word>("adjEqnOption", "jacMatReOrdering"));
            kspOptions.add(
                "PCFillLevel",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "pcFillLevel"));
            kspOptions.add(
                "GMRESMaxIters",
                daOptionPtr_->getSubDictOption<label>("adjEqnOption", "gmresMaxIters"));
            kspOptions.add(
                "GMRESRelTol",
                daOptionPtr_->getSubDictOption<scalar>("adjEqnOption", "gmresRelTol"));
            kspOptions.add(
                "GMRESAbsTol",
                daOptionPtr_->getSubDictOption<scalar>("adjEqnOption", "gmresAbsTol"));
            kspOptions.add("printInfo", 1);

            KSP ksp;
            // psi is the adjoint vector; the solution
            Vec psiVec;
            VecDuplicate(wVec, &psiVec);
            VecZeroEntries(psiVec);

            // create the multi-level Richardson KSP
            daLinearEqn.createMLRKSP(kspOptions, dRdWT, dRdWT, &ksp);

            // solve the linear equation and get psiVec
            daLinearEqn.solveLinearEqn(ksp, dFdWVecAllParts, psiVec);

            // assign the psiVec to psiVecDict_ for all objFuncName such that we can
            // call getPsiVec later in DASolver::calcTotalDeriv
            this->setPsiVecDict(objFuncName, psiVec, psiVecDict_);

            if (daOptionPtr_->getOption<label>("debug"))
            {
                word outputName = "psi_" + objFuncName;
                DAUtility::writeVectorASCII(psiVec, outputName);
                DAUtility::writeVectorBinary(psiVec, outputName);
            }
        }
    }

    return 0;
}

label DASolver::calcTotalDeriv(
    const Vec xvVec,
    const Vec wVec)
{
    // compute the total derivatives

    const dictionary& allOptions = daOptionPtr_->getAllOptions();

    dictionary designVarDict = allOptions.subDict("designVar");

    // loop over all designVarNames in the designVar dict
    forAll(designVarDict.toc(), idxI)
    {
        // design variable name
        word designVarName = designVarDict.toc()[idxI];

        // get the subDict for this dvName
        dictionary dvSubDict = designVarDict.subDict(designVarName);
        // get dvType
        word designVarType = dvSubDict.getWord("designVarType");

        // *******************************************************************
        // **************************** BC dvType ****************************
        // *******************************************************************
        // boundary condition as the design variable, e.g., the inlet velocity
        if (designVarType == "BC")
        {
            // get info from dvSubDict. This needs to be defined in the pyDAFoam
            // name of the variable for changing the boundary condition
            word varName = dvSubDict.getWord("varName");
            // name of the boundary patch
            word patchName = dvSubDict.getWord("patchName");
            // the vol field type, e.g., volScalarField
            word fieldType = dvSubDict.getWord("fieldType");
            // the type of this BC, e.g., fixedValue
            word bcType = dvSubDict.getWord("bcType");
            // the compoent of a vector variable, ignore when it is a scalar
            label comp = dvSubDict.getLabel("comp");

            // no coloring is need for BC, so we create a dummy DAJacCon
            word dummyType = "dummy";
            autoPtr<DAJacCon> daJacCon(DAJacCon::New(
                dummyType,
                meshPtr_(),
                daOptionPtr_(),
                daModelPtr_(),
                daIndexPtr_()));

            // ********************** compute dRdBC **********************
            Mat dRdBC;
            {
                // create DAPartDeriv object
                word modelType = "dRdBC";
                autoPtr<DAPartDeriv> daPartDeriv(DAPartDeriv::New(
                    modelType,
                    meshPtr_(),
                    daOptionPtr_(),
                    daModelPtr_(),
                    daIndexPtr_(),
                    daJacCon(),
                    daResidualPtr_()));

                // setup options to compute dRdBC
                dictionary options;
                options.set("varName", varName);
                options.set("patchName", patchName);
                options.set("fieldType", fieldType);
                options.set("bcType", bcType);
                options.set("comp", comp);

                // initialize the dRdBC matrix
                daPartDeriv->initializePartDerivMat(options, &dRdBC);

                // compute it using brute force finite-difference
                daPartDeriv->calcPartDerivMat(options, xvVec, wVec, dRdBC);

                if (daOptionPtr_->getOption<label>("debug"))
                {
                    word outputName = "dRdBC_" + designVarName;
                    DAUtility::writeMatrixBinary(dRdBC, outputName);
                    DAUtility::writeMatrixASCII(dRdBC, outputName);
                }
            }

            // ********************** compute dFdBC **********************
            dictionary objFuncDict = allOptions.subDict("objFunc");

            // loop over all objFuncName in the objFunc dict
            forAll(objFuncDict.toc(), idxJ)
            {

                word objFuncName = objFuncDict.toc()[idxJ];

                // we only solve adjoint for objectives that have addToAdjoint = True
                if (DAUtility::isInList<word>(objFuncName, objFuncNames4Adj_))
                {
                    // the dFdBCVecAllParts vector contains the sum of dFdBCVec from all parts for this objFuncName
                    Vec dFdBCVecAllParts;
                    VecCreate(PETSC_COMM_WORLD, &dFdBCVecAllParts);
                    VecSetSizes(dFdBCVecAllParts, PETSC_DETERMINE, 1);
                    VecSetFromOptions(dFdBCVecAllParts);
                    VecZeroEntries(dFdBCVecAllParts);

                    dictionary objFuncSubDict = objFuncDict.subDict(objFuncName);
                    // loop over all parts of this objFuncName
                    forAll(objFuncSubDict.toc(), idxK)
                    {
                        word objFuncPart = objFuncSubDict.toc()[idxK];
                        dictionary objFuncSubDictPart = objFuncSubDict.subDict(objFuncPart);

                        // we only compute total derivative for objFuncs with addToAdjoint = True
                        label addToAdjoint = objFuncSubDictPart.getLabel("addToAdjoint");
                        if (addToAdjoint)
                        {

                            Mat dFdBC;

                            // initialize DAPartDeriv for dFdBC
                            word modelType = "dFdBC";
                            autoPtr<DAPartDeriv> daPartDeriv(DAPartDeriv::New(
                                modelType,
                                meshPtr_(),
                                daOptionPtr_(),
                                daModelPtr_(),
                                daIndexPtr_(),
                                daJacCon(),
                                daResidualPtr_()));

                            // initialize options
                            dictionary options;
                            options.set("objFuncName", objFuncName);
                            options.set("objFuncPart", objFuncPart);
                            options.set("objFuncSubDictPart", objFuncSubDictPart);
                            options.set("varName", varName);
                            options.set("patchName", patchName);
                            options.set("fieldType", fieldType);
                            options.set("bcType", bcType);
                            options.set("comp", comp);

                            // initialize dFdBC
                            daPartDeriv->initializePartDerivMat(options, &dFdBC);

                            // calculate it
                            daPartDeriv->calcPartDerivMat(options, xvVec, wVec, dFdBC);

                            // now we need to add all the rows of dFdBC together to get dFdBCVec
                            // NOTE: dFdBC is a 1 by 1 matrix but we just follow a similar
                            // approach for other dFd calculation
                            Vec dFdBCVec, oneVec;
                            VecDuplicate(dFdBCVecAllParts, &oneVec);
                            VecSet(oneVec, 1.0);
                            VecDuplicate(dFdBCVecAllParts, &dFdBCVec);
                            VecZeroEntries(dFdBCVec);
                            // dFdBCVec = oneVec*dFdBC
                            MatMultTranspose(dFdBC, oneVec, dFdBCVec);

                            // we need to add dFdBCVec to dFdBCVecAllParts because we want to sum
                            // all dFdBCVec for all parts of this objFuncName.
                            VecAXPY(dFdBCVecAllParts, 1.0, dFdBCVec);

                            if (daOptionPtr_->getOption<label>("debug"))
                            {
                                word outputName = "dFdBCVec_" + designVarName;
                                DAUtility::writeVectorBinary(dFdBCVec, outputName);
                                DAUtility::writeVectorASCII(dFdBCVec, outputName);
                            }
                        }
                    }

                    // now we can compute totalDeriv = dFdBCVecAllParts - psiVec * dRdBC
                    Vec psiVec, totalDerivVec;
                    VecDuplicate(dFdBCVecAllParts, &totalDerivVec);
                    VecZeroEntries(totalDerivVec);
                    VecDuplicate(wVec, &psiVec);
                    VecZeroEntries(psiVec);

                    // now we can assign DASolver::psiVecDict_ to psiVec for this objFuncName
                    // NOTE: DASolver::psiVecDict_ should be set in the DASolver::solveAdjoint
                    // function. i.e., we need to call solveAdjoint before calling calcTotalDeriv
                    this->getPsiVec(objFuncName, psiVec);

                    // totalDeriv = dFdBCVecAllParts - psiVec * dRdBC
                    MatMultTranspose(dRdBC, psiVec, totalDerivVec);
                    VecAXPY(totalDerivVec, -1.0, dFdBCVecAllParts);
                    VecScale(totalDerivVec, -1.0);

                    // assign totalDerivVec to DASolver::totalDerivDict_ such that we can
                    // get the totalDeriv in the python layer later
                    this->setTotalDerivDict(objFuncName, designVarName, totalDerivVec, totalDerivDict_);

                    if (daOptionPtr_->getOption<label>("debug"))
                    {
                        word outputName = "dFdBCTotal_" + objFuncName + "_" + designVarName;
                        DAUtility::writeVectorBinary(totalDerivVec, outputName);
                        DAUtility::writeVectorASCII(totalDerivVec, outputName);
                    }
                }
            }
        }
        // *******************************************************************
        // **************************** FFD dvType ****************************
        // *******************************************************************
        // FFD movement as the design variable
        else if (designVarType == "FFD")
        {

        }
        else
        {
            FatalErrorIn("") << "designVarType: " << designVarType << " not valid!" << abort(FatalError);
        }
    }

    //Info << totalDerivDict_ << endl;

    return 0;
}

void DASolver::setPsiVecDict(
    const word objFuncName,
    const Vec psiVec,
    dictionary& psiVecDict)
{
    // assign the psiVec to psiVecDict_
    scalarList psiList;
    const PetscScalar* psiVecArray;
    VecGetArrayRead(psiVec, &psiVecArray);
    label Istart, Iend;
    VecGetOwnershipRange(psiVec, &Istart, &Iend);
    for (label i = Istart; i < Iend; i++)
    {
        label relIdx = i - Istart;
        psiList.append(psiVecArray[relIdx]);
    }
    VecRestoreArrayRead(psiVec, &psiVecArray);

    psiVecDict.set(objFuncName, psiList);

}

void DASolver::getPsiVec(
    const word objFuncName,
    Vec psiVec)
{
    /// assign psiVecDict_ to psiVec
    VecZeroEntries(psiVec);

    scalarList psiList;
    psiVecDict_.readEntry<scalarList>(objFuncName, psiList);

    PetscScalar* psiVecArray;
    VecGetArray(psiVec, &psiVecArray);
    label Istart, Iend;
    VecGetOwnershipRange(psiVec, &Istart, &Iend);
    for (label i = Istart; i < Iend; i++)
    {
        label relIdx = i - Istart;
        psiVecArray[relIdx] = psiList[relIdx];
    }
    VecRestoreArray(psiVec, &psiVecArray);
    
}

void DASolver::setTotalDerivDict(
    const word objFuncName,
    const word designVarName,
    const Vec totalDerivVec,
    dictionary& totalDerivDict)
{
    scalarList totalDerivList;
    const PetscScalar* totalDerivVecArray;
    VecGetArrayRead(totalDerivVec, &totalDerivVecArray);
    label Istart, Iend;
    VecGetOwnershipRange(totalDerivVec, &Istart, &Iend);
    for (label i = Istart; i < Iend; i++)
    {
        label relIdx = i - Istart;
        totalDerivList.append(totalDerivVecArray[relIdx]);
    }

    VecRestoreArrayRead(totalDerivVec, &totalDerivVecArray);

    if (!totalDerivDict.found(objFuncName))
    {
        dictionary emptyDict;
        totalDerivDict.set(objFuncName, emptyDict);
    }
    dictionary& objFuncSubDict = totalDerivDict.subDict(objFuncName);
    objFuncSubDict.set(designVarName, totalDerivList);

}

scalar DASolver::getTotalDerivVal(
    const word objFuncName,
    const word designVarName,
    const label idxI) const
{
    dictionary subDict = totalDerivDict_.subDict(objFuncName);
    scalarList valList;
    subDict.readEntry<scalarList>(designVarName, valList);
    return valList[idxI];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
