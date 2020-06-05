/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAJacCon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DAJacCon, 0);
defineRunTimeSelectionTable(DAJacCon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAJacCon::DAJacCon(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DAJacCon",
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      daUtil_(),
      daOption_(mesh.thisDb().lookupObject<DAOption>("DAOption")),
      daIndex_(mesh.thisDb().lookupObject<DAIndex>("DAIndex")),
      daRegState_(mesh.thisDb().lookupObject<DARegState>("DARegState")),
      regStates_(daRegState_.getRegStates())
{
    /*
    Description:
        Construct from Foam::fvMesh
    Input:
        mesh: a fvMesh object
    */
    if (daOption_.getOption<label>("adjUseColoring"))
    {
        // Calculate the boundary connectivity
        Info << "Generating Connectivity for Boundaries:" << endl;

        this->calcNeiBFaceGlobalCompact(neiBFaceGlobalCompact_);

        this->setupStateBoundaryCon(&stateBoundaryCon_);

        this->setupStateBoundaryConID(&stateBoundaryConID_);

        if (daOption_.getOption<label>("debug"))
        {
            daUtil_.writeMatrixBinary(stateBoundaryCon_, "stateBoundaryCon");
            daUtil_.writeMatrixBinary(stateBoundaryConID_, "stateBoundaryConID");
        }

        this->initializePetscVecs();
    }
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DAJacCon> DAJacCon::New(const fvMesh& mesh)
{
    // standard setup for runtime selectable classes

    // look up the solver name defined in system/DADict
    const DAOption& daOption = mesh.thisDb().lookupObject<DAOption>("DAOption");
    word solverName = daOption.getOption<word>("solverName");

    Info << "Selecting " << solverName << " for DAJacCon" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DAJacCon::New"
            "("
            "    const fvMesh&"
            ")")
            << "Unknown DAJacCon type "
            << solverName << nl << nl
            << "Valid DAJacCon types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DAJacCon>(
        cstrIter()(mesh));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// this is a virtual function for regIOobject
bool DAJacCon::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

void DAJacCon::setupdRdWCon(
    const label isPrealloc,
    const label isPC)
{
    /*
    Calculates the state Jacobian connectivity mat DAJacCon::dRdWCon_ or 
    computing the preallocation vectors DAJacCon::dRdWPreallocOn_ and 
    DAJacCon::dRdWPreallocOff_

    Input:
    ------
    isPrealloc == 1: calculate the preallocation vectors, else, calculate dRdWCon_
    
    isPC == 1: calc dRdWConPC_, else calculate dRdWCon_

    Output:
    ------
    DAJacCon::dRdWPreallocOn_: preallocation vector that stores the number of 
    on-diagonal conectivity for each row

    DAJacCon::dRdWCon_: state Jacobian connectivity mat with dimension 
    sizeAdjStates by sizeAdjStates. dRdWCon has the same non-zero pattern as dRdW.
    The difference is that dRdWCon has values one for all non-zero values, so dRdWCon
    may look like this

                1 1 0 0 1 0 
                1 1 1 0 0 1 
                0 1 1 1 0 0 
    dRdWCon =   1 0 1 1 1 0 
                0 1 0 1 1 1
                0 0 1 0 0 1

    Example:
    -------
    The way setupdRdWCon works is that we call the DAJacCon::addStateConnections function
    to add connectivity for each row of DAJacCondRdWCon_.
    
    Here we need to loop over all cellI and add a certain number levels of connected states.
    If the connectivity list reads:

    adjStateResidualConInfo_
    {
        "URes"
        {
            {"U", "p", "phi"}, // level 0 connectivity
            {"U", "p", "phi"}, // level 1 connectivity
            {"U"},             // level 2 connectivity
        }
    }

    and the cell topology with a inter-proc boundary cen be either of the following:
    CASE 1:
                       ---------
                       | cellQ |
                -----------------------
               | cellP | cellJ | cellO |             <------ proc1
    ------------------------------------------------ <----- inter-processor boundary
       | cellT | cellK | cellI | cellL | cellU |     <------ proc0
       -----------------------------------------
               | cellN | cellM | cellR |
                ------------------------
                       | cellS |
                       ---------
    
    CASE 2:
                       ---------
                       | cellQ |                       <------ proc1
    -------------------------------------------------- <----- inter-processor boundary
               | cellP | cellJ | cellO |               <------ proc0
       ----------------------------------------- 
       | cellT | cellK | cellI | cellL | cellU |    
       -----------------------------------------
               | cellN | cellM | cellR |
                ------------------------
                       | cellS |
                       ---------
    
    Then, to add the connectivity correctly, we need to add all levels of connected
    states for cellI.
    Level 0 connectivity is straightforward becasue we don't need
    to provide connectedStatesInterProc

    To add level 1 connectivity, we need to:
    set connectedLevelLocal = 1
    set connectedStatesLocal = {U, p}
    set connectedStatesInterProc = {{U,p}, {U}}
    set addFace = 1 
    NOTE: we need set level 1 and level 2 con in connectedStatesInterProc because the 
    north face of cellI is a inter-proc boundary and there are two levels of connected
    state on the other side of the inter-proc boundary for CASE 1. This is the only chance we 
    can add all two levels of connected state across the boundary for CASE 1. For CASE 2, we won't
    add any level 1 inter-proc states because non of the faces for cellI are inter-proc
    faces so calling DAJacCon::addBoundaryFaceConnections for cellI won't add anything

    To add level 2 connectivity, we need to
    set connectedLevelLocal = 2
    set connectedStatesLocal = {U}
    set connectedStatesInterProc = {{U}}
    set addFace = 0
    NOTE 1: we need only level 2 con (U) for connectedStatesInterProc because if we are in CASE 1,
    the level 2 of inter-proc states have been added. For CASE 2, we only need to add cellQ
    by calling DAJacCon::addBoundaryFaceConnections with cellJ
    NOTE 2: If we didn't call two levels of connectedStatesInterProc in the previous call for 
    level 1 con, we can not add it for connectedLevelLocal = 2 becasue for CASE 2 there is no
    inter-proc boundary for cellI

    NOTE: how to provide connectedLevelLocal, connectedStatesLocal, and connectedStatesInterProc
    are done in DAJacCon::setupdRdWCon and DAJacCon::setupObjFuncCon

    */

    // first decide which con mat we need to set
    Mat* conMat;
    if (isPC)
    {
        conMat = &dRdWConPC_;
    }
    else
    {
        conMat = &dRdWCon_;
    }

    // if it is preconditioner mat, check if we need to reduce its
    // con level for reducing memory cost
    if (isPC)
    {
        this->reduceAdjStateResidualConLevel();
    }

    label globalIdx;
    // connectedStatesP: one row matrix that stores the actual connectivity,
    // element value 1 denotes a connected state. connectedStatesP is then used to
    // assign dRdWPreallocOn or dRdWCon
    Mat connectedStatesP;

    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    PetscInt nColsID;
    const PetscInt* colsID;
    const PetscScalar* valsID;

    if (isPrealloc)
    {
        VecZeroEntries(dRdWPreallocOn_);
        VecZeroEntries(dRdWPreallocOff_);
        VecZeroEntries(dRdWTPreallocOn_);
        VecZeroEntries(dRdWTPreallocOff_);
    }

    if (isPrealloc)
    {
        Info << "Preallocating state Jacobian connectivity mat" << endl;
    }
    else
    {
        Info << "Setup state Jacobian connectivity mat" << endl;
    }

    // loop over all cell residuals, we bascially need to compute all
    // the input parameters for the DAJacCon::addStateConnections function, then call
    // the function to get connectedStatesP (matrix that contains one row of connectivity
    // in DAJacCon::dRdWCon_). Check DAJacCon::addStateConnections for detail usages
    forAll(daIndex_.adjStateNames, idxI)
    {
        // get stateName and residual names
        word stateName = daIndex_.adjStateNames[idxI];
        word resName = stateName + "Res";

        // check if this state is a cell state, we do surfaceScalarState residuals separately
        if (daIndex_.adjStateType[stateName] == "surfaceScalarState")
        {
            continue;
        }

        // maximal connectivity level information
        // Note that adjStateResidualConInfo_ starts with level zero,
        // so the maxConLeve is its size minus one
        label maxConLevel = adjStateResidualConInfo_[resName].size() - 1;

        // if it is a vectorState, set compMax=3
        label compMax = 1;
        if (daIndex_.adjStateType[stateName] == "volVectorState")
        {
            compMax = 3;
        }

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < compMax; comp++)
            {

                // zero the connections
                this->createConnectionMat(&connectedStatesP);

                // now add the con. We loop over all the connectivity levels
                forAll(adjStateResidualConInfo_[resName], idxJ) // idxJ: con level
                {

                    // set connectedStatesLocal: the locally connected state variables for this level
                    wordList connectedStatesLocal(0);
                    forAll(adjStateResidualConInfo_[resName][idxJ], idxK)
                    {
                        word conName = adjStateResidualConInfo_[resName][idxJ][idxK];
                        // Exclude surfaceScalarState when appending connectedStatesLocal
                        // whether to add it depends on addFace parameter
                        if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                        {
                            connectedStatesLocal.append(conName);
                        }
                    }

                    // set connectedStatesInterProc: the globally connected state variables for this level
                    List<List<word>> connectedStatesInterProc;
                    if (idxJ == 0)
                    {
                        // pass a zero list, no need to add interProc connecitivity for level 0
                        connectedStatesInterProc.setSize(0);
                    }
                    else if (idxJ != maxConLevel)
                    {
                        connectedStatesInterProc.setSize(maxConLevel - idxJ + 1);
                        for (label k = 0; k < maxConLevel - idxJ + 1; k++)
                        {
                            label conSize = adjStateResidualConInfo_[resName][k + idxJ].size();
                            for (label l = 0; l < conSize; l++)
                            {
                                word conName = adjStateResidualConInfo_[resName][k + idxJ][l];
                                // Exclude surfaceScalarState when appending connectedStatesLocal
                                // whether to add it depends on addFace parameter
                                if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                                {
                                    connectedStatesInterProc[k].append(conName);
                                }
                            }
                        }
                    }
                    else
                    {
                        connectedStatesInterProc.setSize(1);
                        label conSize = adjStateResidualConInfo_[resName][maxConLevel].size();
                        for (label l = 0; l < conSize; l++)
                        {
                            word conName = adjStateResidualConInfo_[resName][maxConLevel][l];
                            // Exclude surfaceScalarState when appending connectedStatesLocal
                            // whether to add it depends on addFace parameter
                            if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                            {
                                connectedStatesInterProc[0].append(conName);
                            }
                        }
                    }

                    // check if we need to addFace for this level
                    label addFace = 0;
                    forAll(regStates_["surfaceScalarStates"], idxK)
                    {
                        word conName = regStates_["surfaceScalarStates"][idxK];
                        if (daUtil_.isInList<word>(conName, adjStateResidualConInfo_[resName][idxJ]))
                        {
                            addFace = 1;
                        }
                    }

                    // Add connectivity
                    this->addStateConnections(
                        connectedStatesP,
                        cellI,
                        idxJ,
                        connectedStatesLocal,
                        connectedStatesInterProc,
                        addFace);

                    //Info<<"lv: "<<idxJ<<" locaStates: "<<connectedStatesLocal<<" interProcStates: "
                    //    <<connectedStatesInterProc<<" addFace: "<<addFace<<endl;
                }

                // get the global index of the current state for the row index
                globalIdx = daIndex_.getGlobalAdjointStateIndex(stateName, cellI, comp);

                if (isPrealloc)
                {
                    this->allocateJacobianConnections(
                        dRdWPreallocOn_,
                        dRdWPreallocOff_,
                        dRdWTPreallocOn_,
                        dRdWTPreallocOff_,
                        connectedStatesP,
                        globalIdx);
                }
                else
                {
                    this->setupJacobianConnections(
                        *conMat,
                        connectedStatesP,
                        globalIdx);
                }
            }
        }
    }

    // loop over all face residuals
    forAll(regStates_["surfaceScalarStates"], idxI)
    {
        // get stateName and residual names
        word stateName = regStates_["surfaceScalarStates"][idxI];
        word resName = stateName + "Res";

        // maximal connectivity level information
        label maxConLevel = adjStateResidualConInfo_[resName].size() - 1;

        forAll(mesh_.faces(), faceI)
        {

            //zero the connections
            this->createConnectionMat(&connectedStatesP);

            // Get the owner and neighbour cells for this face
            label idxO = -1, idxN = -1;
            if (faceI < daIndex_.nLocalInternalFaces)
            {
                idxO = mesh_.owner()[faceI];
                idxN = mesh_.neighbour()[faceI];
            }
            else
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                label patchIdx = daIndex_.bFacePatchI[relIdx];
                label faceIdx = daIndex_.bFaceFaceI[relIdx];

                const UList<label>& pFaceCells = mesh_.boundaryMesh()[patchIdx].faceCells();
                idxN = pFaceCells[faceIdx];
            }

            // now add the con. We loop over all the connectivity levels
            forAll(adjStateResidualConInfo_[resName], idxJ) // idxJ: con level
            {

                // set connectedStatesLocal: the locally connected state variables for this level
                wordList connectedStatesLocal(0);
                forAll(adjStateResidualConInfo_[resName][idxJ], idxK)
                {
                    word conName = adjStateResidualConInfo_[resName][idxJ][idxK];
                    // Exclude surfaceScalarState when appending connectedStatesLocal
                    // whether to add it depends on addFace parameter
                    if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                    {
                        connectedStatesLocal.append(conName);
                    }
                }

                // set connectedStatesInterProc: the globally connected state variables for this level
                List<List<word>> connectedStatesInterProc;
                if (idxJ == 0)
                {
                    // pass a zero list, no need to add interProc connecitivity for level 0
                    connectedStatesInterProc.setSize(0);
                }
                else if (idxJ != maxConLevel)
                {
                    connectedStatesInterProc.setSize(maxConLevel - idxJ + 1);
                    for (label k = 0; k < maxConLevel - idxJ + 1; k++)
                    {
                        label conSize = adjStateResidualConInfo_[resName][k + idxJ].size();
                        for (label l = 0; l < conSize; l++)
                        {
                            word conName = adjStateResidualConInfo_[resName][k + idxJ][l];
                            // Exclude surfaceScalarState when appending connectedStatesLocal
                            // whether to add it depends on addFace parameter
                            if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                            {
                                connectedStatesInterProc[k].append(conName);
                            }
                        }
                    }
                }
                else
                {
                    connectedStatesInterProc.setSize(1);
                    label conSize = adjStateResidualConInfo_[resName][maxConLevel].size();
                    for (label l = 0; l < conSize; l++)
                    {
                        word conName = adjStateResidualConInfo_[resName][maxConLevel][l];
                        // Exclude surfaceScalarState when appending connectedStatesLocal
                        // whether to add it depends on addFace parameter
                        if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                        {
                            connectedStatesInterProc[0].append(conName);
                        }
                    }
                }

                // check if we need to addFace for this level
                label addFace = 0;
                forAll(regStates_["surfaceScalarStates"], idxK)
                {
                    word conName = regStates_["surfaceScalarStates"][idxK];
                    // NOTE: we need special treatment for boundary faces for level>0
                    // since addFace for boundary face should add one more extra level of faces
                    // This is because we only have idxN for a boundary face while the idxO can
                    // be on the other side of the inter-proc boundary
                    // In this case, we need to use idxJ-1 instead of idxJ information to tell whether to addFace
                    label levelCheck;
                    if (faceI < daIndex_.nLocalInternalFaces or idxJ == 0)
                    {
                        levelCheck = idxJ;
                    }
                    else
                    {
                        levelCheck = idxJ - 1;
                    }

                    if (daUtil_.isInList<word>(conName, adjStateResidualConInfo_[resName][levelCheck]))
                    {
                        addFace = 1;
                    }
                }

                // Add connectivity for idxN
                this->addStateConnections(
                    connectedStatesP,
                    idxN,
                    idxJ,
                    connectedStatesLocal,
                    connectedStatesInterProc,
                    addFace);

                if (faceI < daIndex_.nLocalInternalFaces)
                {
                    // Add connectivity for idxO
                    this->addStateConnections(
                        connectedStatesP,
                        idxO,
                        idxJ,
                        connectedStatesLocal,
                        connectedStatesInterProc,
                        addFace);
                }

                //Info<<"lv: "<<idxJ<<" locaStates: "<<connectedStatesLocal<<" interProcStates: "
                //    <<connectedStatesInterProc<<" addFace: "<<addFace<<endl;
            }

            // NOTE: if this faceI is on a coupled patch, the above connectivity is not enough to
            // cover the points on the other side of proc domain, we need to add 3 lvs of cells here
            if (faceI >= daIndex_.nLocalInternalFaces)
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                label patchIdx = daIndex_.bFacePatchI[relIdx];

                label maxLevel = adjStateResidualConInfo_[resName].size();

                if (mesh_.boundaryMesh()[patchIdx].coupled())
                {

                    label bRow = this->getLocalCoupledBFaceIndex(faceI);
                    label bRowGlobal = daIndex_.globalCoupledBFaceNumbering.toGlobal(bRow);
                    MatGetRow(stateBoundaryCon_, bRowGlobal, &nCols, &cols, &vals);
                    MatGetRow(stateBoundaryConID_, bRowGlobal, &nColsID, &colsID, &valsID);
                    for (label i = 0; i < nCols; i++)
                    {
                        PetscInt idxJ = cols[i];
                        label val = round(vals[i]);
                        // we are going to add some selective states with connectivity level <= 3
                        // first check the state
                        label stateID = round(valsID[i]);
                        word conName = daIndex_.adjStateNames[stateID];
                        label addState = 0;
                        // NOTE: we use val-1 here since phi actually has 3 levels of connectivity
                        // however, when we assign adjStateResidualConInfo_, we ignore the level 0
                        // connectivity since they are idxN and idxO
                        if (val != 10 && val < maxLevel + 1)
                        {
                            if (daUtil_.isInList<word>(conName, adjStateResidualConInfo_[resName][val - 1]))
                            {
                                addState = 1;
                            }
                        }
                        if (addState == 1 && val < maxLevel + 1 && val > 0)
                        {
                            this->setConnections(connectedStatesP, idxJ);
                        }
                    }
                    MatRestoreRow(stateBoundaryCon_, bRowGlobal, &nCols, &cols, &vals);
                    MatRestoreRow(stateBoundaryConID_, bRowGlobal, &nColsID, &colsID, &valsID);
                }
            }

            // get the global index of the current state for the row index
            globalIdx = daIndex_.getGlobalAdjointStateIndex(stateName, faceI);

            if (isPrealloc)
            {
                this->allocateJacobianConnections(
                    dRdWPreallocOn_,
                    dRdWPreallocOff_,
                    dRdWTPreallocOn_,
                    dRdWTPreallocOff_,
                    connectedStatesP,
                    globalIdx);
            }
            else
            {
                this->setupJacobianConnections(
                    *conMat,
                    connectedStatesP,
                    globalIdx);
            }
        }
    }

    if (isPC)
    {
        this->restoreAdjStateResidualConLevel();
    }

    if (isPrealloc)
    {
        VecAssemblyBegin(dRdWPreallocOn_);
        VecAssemblyEnd(dRdWPreallocOn_);
        VecAssemblyBegin(dRdWPreallocOff_);
        VecAssemblyEnd(dRdWPreallocOff_);
        VecAssemblyBegin(dRdWTPreallocOn_);
        VecAssemblyEnd(dRdWTPreallocOn_);
        VecAssemblyBegin(dRdWTPreallocOff_);
        VecAssemblyEnd(dRdWTPreallocOff_);

        //output the matrix to a file
        if (daOption_.getOption<label>("debug"))
        {
            daUtil_.writeVectorASCII(dRdWTPreallocOn_, "dRdWTPreallocOn");
            daUtil_.writeVectorASCII(dRdWTPreallocOff_, "dRdWTPreallocOff");
            daUtil_.writeVectorASCII(dRdWPreallocOn_, "dRdWPreallocOn");
            daUtil_.writeVectorASCII(dRdWPreallocOff_, "dRdWPreallocOff");
        }
    }
    else
    {
        MatAssemblyBegin(*conMat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*conMat, MAT_FINAL_ASSEMBLY);

        //output the matrix to a file
        if (daOption_.getOption<label>("debug"))
        {
            //daUtil_.writeMatRowSize(*conMat, "dRdWCon");
            daUtil_.writeMatrixBinary(*conMat, "dRdWCon");
        }
    }

    if (isPrealloc)
    {
        Info << "Preallocating state Jacobian connectivity mat: finished!" << endl;
    }
    else
    {
        Info << "Setup state Jacobian connectivity mat: finished!" << endl;
    }
}

void DAJacCon::reduceAdjStateResidualConLevel()
{
    /*
    Reduce the connectivity levels for DAJacCon::adjStateResidualConInfo_
    based on maxResConLv4JacPCMat specified in DAOption

    Input:
    -----
    maxResConLv4JacPCMat: the maximal levels of connectivity for each
    state variable residual

    Output:
    ------
    adjStateResidualConInfo_: reduced connectivity level.

    adjStateResidualConInfoBK_: Original connectivity level. Will be used
    when calling DAJacCon::restoreAdjStateResidualConLevel to restore 
    the connectivity level to DAJacCon::adjStateResidualConInfo_

    Example:
    -------

    If the original adjStateResidualConInfo_ reads:

    adjStateResidualConInfo_
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
    
    Then, calling reduceAdjStateResidualConLevel will give:

    adjStateResidualConInfo_
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
    HashTable<label> maxResConLv4JacPCMat =
        daOption_.getOption<HashTable<label>>("maxResConLv4JacPCMat");
    if (maxResConLv4JacPCMat.size() == 0)
    {
        return;
    }

    // now check if maxResConLv4JacPCMat has all the maxRes level defined
    // and these max levels are <= adjStateResidualConInfo_.size()
    forAll(adjStateResidualConInfo_.toc(), idxJ)
    {
        word key1 = adjStateResidualConInfo_.toc()[idxJ];
        bool keyFound = false;
        forAll(maxResConLv4JacPCMat.toc(), idxI)
        {
            word key = maxResConLv4JacPCMat.toc()[idxI];
            if (key == key1)
            {
                keyFound = true;
                label maxLv = maxResConLv4JacPCMat[key];
                label maxLv1 = adjStateResidualConInfo_[key1].size() - 1;
                if (maxLv > maxLv1)
                {
                    FatalErrorIn("") << "maxResConLv4JacPCMat maxLevel"
                                     << maxLv << " for " << key
                                     << " larger than adjStateResidualConInfo maxLevel "
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

    // assign adjStateResidualConInfo_ to adjStateResidualConInfoBK_
    forAll(adjStateResidualConInfo_.toc(), idxI)
    {
        word key = adjStateResidualConInfo_.toc()[idxI];
        adjStateResidualConInfoBK_.set(key, adjStateResidualConInfo_[key]);
    }

    // now we can erase adjStateResidualConInfo
    adjStateResidualConInfo_.clearStorage();

    // get the reduced adjStateResidualConInfo_
    forAll(adjStateResidualConInfoBK_.toc(), idxI)
    {
        word key = adjStateResidualConInfoBK_.toc()[idxI];
        label maxConLevel = maxResConLv4JacPCMat[key];
        label conSize = adjStateResidualConInfoBK_[key].size();
        if (conSize > maxConLevel + 1)
        {
            List<List<word>> conList;
            conList.setSize(maxConLevel + 1);
            for (label i = 0; i <= maxConLevel; i++) // NOTE: it is <=
            {
                conList[i] = adjStateResidualConInfoBK_[key][i];
            }
            adjStateResidualConInfo_.set(key, conList);
        }
        else
        {
            adjStateResidualConInfo_.set(key, adjStateResidualConInfoBK_[key]);
        }
    }
    //Info<<adjStateResidualConInfo_<<endl;
}

void DAJacCon::restoreAdjStateResidualConLevel()
{
    /*
    Assign DAJacCon::adjStateResidualConInfoBK_ to DAJacCon::adjStateResidualConInfo_
    such that the reduced connecitvity is restored.
    See DAJacCon::reduceAdjStateResidualConLevel for recuding con levels

    Input:
    -----
    adjStateResidualConInfoBK_: the back up of original adjStateResidualConInfo_

    Output:
    ------
    adjStateResidualConInfo_: original un-reduced connectivity levels

    */

    forAll(adjStateResidualConInfoBK_.toc(), idxI)
    {
        word key = adjStateResidualConInfoBK_.toc()[idxI];
        adjStateResidualConInfo_.set(key, adjStateResidualConInfoBK_[key]);
    }

    return;
}

void DAJacCon::allocateJacobianConnections(
    Vec preallocOnProc,
    Vec preallocOffProc,
    Vec preallocOnProcT,
    Vec preallocOffProcT,
    Mat connections,
    const label row)
{
    /*
    Compute the matrix allocation vector based on one row connection mat

    Input:
    -----
    connections: a one row matrix that contains all the nonzeros for one 
    row of DAJacCon::dRdWCon

    row: which row to add for the preallocation vector

    Output:
    ------
    preallocOnProc: the vector that contains the number of nonzeros for each row
    in dRdW (on-diagonal block elements)

    preallocOffProc: the vector that contains the number of nonzeros for each row
    in dRdW (off-diagonal block elements)

    preallocOnProcT: the vector that contains the number of nonzeros for each row
    in dRdWT (on-diagonal block elements)

    preallocOffProcT: the vector that contains the number of nonzeros for each row
    in dRdWT (off-diagonal block elements)

    */
    PetscScalar v = 1.0;
    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    MatAssemblyBegin(connections, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(connections, MAT_FINAL_ASSEMBLY);

    // Compute the transposed case
    // in this case connections represents a single column, so we need to
    // increment the counter in each row with a non-zero entry.

    label colMin = daIndex_.globalAdjointStateNumbering.toGlobal(0);
    label colMax = colMin + daIndex_.nLocalAdjointStates;
    // by construction rows should be limited to local rows
    MatGetRow(connections, 0, &nCols, &cols, &vals);

    // for the non-transposed case just sum up the row.
    // count up the total number of non zeros in this row
    label totalCount = 0; //2
    label localCount = 0;
    //int idx;
    for (label j = 0; j < nCols; j++)
    {
        // int idx = cols[j];
        scalar val = vals[j];
        if (daUtil_.isValueCloseToRef(val, 1.0))
        {
            // We can compute the first part of the non-transposed row here.
            totalCount++;
            label idx = cols[j];
            // Set the transposed version as well
            if (colMin <= idx && idx < colMax)
            {
                //this entry is a local entry, increment the corresponding row
                VecSetValue(preallocOnProcT, idx, v, ADD_VALUES);
                localCount++;
            }
            else
            {
                // this is an off proc entry.
                VecSetValue(preallocOffProcT, idx, v, ADD_VALUES);
            }
        }
    }

    label offProcCount = totalCount - localCount;
    VecSetValue(preallocOnProc, row, localCount, INSERT_VALUES);
    VecSetValue(preallocOffProc, row, offProcCount, INSERT_VALUES);

    // restore the row of the matrix
    MatRestoreRow(connections, 0, &nCols, &cols, &vals);
    MatDestroy(&connections);

    return;
}

void DAJacCon::setupJacobianConnections(
    Mat conMat,
    Mat connections,
    const PetscInt idxI)
{
    /*
    Assign connectivity to Jacobian conMat, e.g., dRdWCon, based on the connections input Mat
    
    Input:
    ------
    idxI: Row index to added, ad the column index to added is based on connections

    connections: the one row matrix with nonzero values to add to conMat

    Output:
    ------
    conMat: the connectivity mat to add
    */

    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    MatAssemblyBegin(connections, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(connections, MAT_FINAL_ASSEMBLY);

    MatGetRow(connections, 0, &nCols, &cols, &vals);
    MatSetValues(conMat, 1, &idxI, nCols, cols, vals, INSERT_VALUES);

    // restore the row of the matrix
    MatRestoreRow(connections, 0, &nCols, &cols, &vals);
    MatDestroy(&connections);

    return;
}

void DAJacCon::createConnectionMat(Mat* connectedStates)
{
    /*
    Initialize a serial connectivity matrix connectedStates,
    basically, it is one row of connectivity in dRdWCon
    */

    // create a local matrix to store this row's connectivity
    MatCreateSeqAIJ(
        PETSC_COMM_SELF,
        1,
        daIndex_.nGlobalAdjointStates,
        2000,
        NULL,
        connectedStates);
    //MatSetOption(*connectedStates, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*connectedStates);

    MatZeroEntries(*connectedStates);

    return;
}

void DAJacCon::addStateConnections(
    Mat connections,
    const label cellI,
    const label connectedLevelLocal,
    const wordList connectedStatesLocal,
    const List<List<word>> connectedStatesInterProc,
    const label addFace)
{
    /*
    A high level interface to add the connectivity for the row matrix connections
    Note: the connections mat is basically one row of connectivity in DAJacCon::dRdWCon_

    Input:
    ------
    cellI: cell index based on which we want to add the connectivity. We can add any level of 
    connected states to this cellI

    connectedLevelLocal: level of local connectivity, this is useually obtained from
    DAJacCon::adjStateResidualConInfo_

    connectedStatesLocal: list of connected states to add for the level: connectedLevelLocal

    connectedStatesInterProc: list of states to add for a given level of boundary connectivity 
    
    addFace: add cell faces for the current level?

    Output:
    ------
    connections: one row of connectivity in DAJacCon::dRdWCon_

    Example:
    -------
    If the connectivity list reads:

    adjStateResidualConInfo_
    {
        "URes"
        {
            {"U", "p", "phi"}, // level 0 connectivity
            {"U", "p", "phi"}, // level 1 connectivity
            {"U"},             // level 2 connectivity
        }
    }

    and the cell topology with a inter-proc boundary cen be either of the following:
    CASE 1:
                       ---------
                       | cellQ |
                -----------------------
               | cellP | cellJ | cellO |             <------ proc1
    ------------------------------------------------ <----- inter-processor boundary
       | cellT | cellK | cellI | cellL | cellU |     <------ proc0
       -----------------------------------------
               | cellN | cellM | cellR |
                ------------------------
                       | cellS |
                       ---------
    
    CASE 2:
                       ---------
                       | cellQ |                       <------ proc1
    -------------------------------------------------- <----- inter-processor boundary
               | cellP | cellJ | cellO |               <------ proc0
       ----------------------------------------- 
       | cellT | cellK | cellI | cellL | cellU |    
       -----------------------------------------
               | cellN | cellM | cellR |
                ------------------------
                       | cellS |
                       ---------
    
    Then, to add the connectivity correctly, we need to add all levels of connected
    states for cellI.
    Level 0 connectivity is straightforward becasue we don't need
    to provide connectedStatesInterProc

    To add level 1 connectivity, we need to:
    set connectedLevelLocal = 1
    set connectedStatesLocal = {U, p}
    set connectedStatesInterProc = {{U,p}, {U}}
    set addFace = 1 
    NOTE: we need set level 1 and level 2 con in connectedStatesInterProc because the 
    north face of cellI is a inter-proc boundary and there are two levels of connected
    state on the other side of the inter-proc boundary for CASE 1. This is the only chance we 
    can add all two levels of connected state across the boundary for CASE 1. For CASE 2, we won't
    add any level 1 inter-proc states because non of the faces for cellI are inter-proc
    faces so calling DAJacCon::addBoundaryFaceConnections for cellI won't add anything

    To add level 2 connectivity, we need to
    set connectedLevelLocal = 2
    set connectedStatesLocal = {U}
    set connectedStatesInterProc = {{U}}
    set addFace = 0
    NOTE 1: we need only level 2 con (U) for connectedStatesInterProc because if we are in CASE 1,
    the level 2 of inter-proc states have been added. For CASE 2, we only need to add cellQ
    by calling DAJacCon::addBoundaryFaceConnections with cellJ
    NOTE 2: If we didn't call two levels of connectedStatesInterProc in the previous call for 
    level 1 con, we can not add it for connectedLevelLocal = 2 becasue for CASE 2 there is no
    inter-proc boundary for cellI

    NOTE: how to provide connectedLevelLocal, connectedStatesLocal, and connectedStatesInterProc
    are done in DAJacCon::setupdRdWCon and DAJacCon::setupObjFuncCon

    */

    // check if the input parameters are valid
    if (connectedLevelLocal > 3 or connectedLevelLocal < 0)
    {
        FatalErrorIn("connectedLevelLocal not valid") << abort(FatalError);
    }
    if (addFace != 0 && addFace != 1)
    {
        FatalErrorIn("addFace not valid") << abort(FatalError);
    }
    if (cellI >= mesh_.nCells())
    {
        FatalErrorIn("cellI not valid") << abort(FatalError);
    }
    //if (connectedLevelLocal>=2 && addFace==1)
    //{
    //    FatalErrorIn("addFace not supported for localLevel>=2")<< abort(FatalError);
    //}

    labelList val1 = {1};
    labelList vals2 = {1, 1};
    labelList vals3 = {1, 1, 1};

    label interProcLevel = connectedStatesInterProc.size();

    if (connectedLevelLocal == 0)
    {
        // add connectedStatesLocal for level0
        forAll(connectedStatesLocal, idxI)
        {
            word stateName = connectedStatesLocal[idxI];
            label compMax = 1;
            if (daIndex_.adjStateType[stateName] == "volVectorState")
            {
                compMax = 3;
            }
            for (label i = 0; i < compMax; i++)
            {
                label idxJ = daIndex_.getGlobalAdjointStateIndex(stateName, cellI, i);
                this->setConnections(connections, idxJ);
            }
        }
        // add faces for level0
        if (addFace)
        {
            forAll(regStates_["surfaceScalarStates"], idxI)
            {
                word stateName = regStates_["surfaceScalarStates"][idxI];
                this->addConMatCellFaces(connections, 0, cellI, stateName, 1.0);
            }
        }
    }
    else if (connectedLevelLocal == 1)
    {
        // add connectedStatesLocal for level1
        forAll(connectedStatesLocal, idxI)
        {
            word stateName = connectedStatesLocal[idxI];
            this->addConMatNeighbourCells(connections, 0, cellI, stateName, 1.0);
        }

        // add faces for level1
        if (addFace)
        {
            forAll(mesh_.cellCells()[cellI], cellJ)
            {
                label localCell = mesh_.cellCells()[cellI][cellJ];
                forAll(regStates_["surfaceScalarStates"], idxI)
                {
                    word stateName = regStates_["surfaceScalarStates"][idxI];
                    this->addConMatCellFaces(connections, 0, localCell, stateName, 1.0);
                }
            }
        }
        // add inter-proc connectivity for level1
        if (interProcLevel == 0)
        {
            // pass, not adding anything
        }
        else if (interProcLevel == 1)
        {
            this->addBoundaryFaceConnections(
                connections,
                0,
                cellI,
                val1,
                connectedStatesInterProc,
                addFace);
        }
        else if (interProcLevel == 2)
        {
            this->addBoundaryFaceConnections(
                connections,
                0,
                cellI,
                vals2,
                connectedStatesInterProc,
                addFace);
        }
        else if (interProcLevel == 3)
        {
            this->addBoundaryFaceConnections(
                connections,
                0,
                cellI,
                vals3,
                connectedStatesInterProc,
                addFace);
        }
        else
        {
            FatalErrorIn("interProcLevel not valid") << abort(FatalError);
        }
    }
    else if (connectedLevelLocal == 2)
    {
        forAll(mesh_.cellCells()[cellI], cellJ)
        {
            label localCell = mesh_.cellCells()[cellI][cellJ];

            // add connectedStatesLocal for level2
            forAll(connectedStatesLocal, idxI)
            {
                word stateName = connectedStatesLocal[idxI];
                this->addConMatNeighbourCells(connections, 0, localCell, stateName, 1.0);
            }

            // add faces for level2
            if (addFace)
            {
                forAll(mesh_.cellCells()[localCell], cellK)
                {
                    label localCellK = mesh_.cellCells()[localCell][cellK];
                    forAll(regStates_["surfaceScalarStates"], idxI)
                    {
                        word stateName = regStates_["surfaceScalarStates"][idxI];
                        this->addConMatCellFaces(connections, 0, localCellK, stateName, 1.0);
                    }
                }
            }

            // add inter-proc connecitivty for level2
            if (interProcLevel == 0)
            {
                // pass, not adding anything
            }
            else if (interProcLevel == 1)
            {
                this->addBoundaryFaceConnections(
                    connections,
                    0,
                    localCell,
                    val1,
                    connectedStatesInterProc,
                    addFace);
            }
            else if (interProcLevel == 2)
            {
                this->addBoundaryFaceConnections(
                    connections,
                    0,
                    localCell,
                    vals2,
                    connectedStatesInterProc,
                    addFace);
            }
            else if (interProcLevel == 3)
            {
                this->addBoundaryFaceConnections(
                    connections,
                    0,
                    localCell,
                    vals3,
                    connectedStatesInterProc,
                    addFace);
            }
            else
            {
                FatalErrorIn("interProcLevel not valid") << abort(FatalError);
            }
        }
    }
    else if (connectedLevelLocal == 3)
    {

        forAll(mesh_.cellCells()[cellI], cellJ)
        {
            label localCell = mesh_.cellCells()[cellI][cellJ];
            forAll(mesh_.cellCells()[localCell], cellK)
            {
                label localCell2 = mesh_.cellCells()[localCell][cellK];

                // add connectedStatesLocal for level3
                forAll(connectedStatesLocal, idxI)
                {
                    word stateName = connectedStatesLocal[idxI];
                    this->addConMatNeighbourCells(connections, 0, localCell2, stateName, 1.0);
                }

                // add faces for level3
                if (addFace)
                {
                    forAll(mesh_.cellCells()[localCell2], cellL)
                    {
                        label localCellL = mesh_.cellCells()[localCell2][cellL];
                        forAll(regStates_["surfaceScalarStates"], idxI)
                        {
                            word stateName = regStates_["surfaceScalarStates"][idxI];
                            this->addConMatCellFaces(connections, 0, localCellL, stateName, 1.0);
                        }
                    }
                }

                // add inter-proc connecitivty for level3
                if (interProcLevel == 0)
                {
                    // pass, not adding anything
                }
                else if (interProcLevel == 1)
                {
                    this->addBoundaryFaceConnections(
                        connections,
                        0,
                        localCell2,
                        val1,
                        connectedStatesInterProc,
                        addFace);
                }
                else if (interProcLevel == 2)
                {
                    this->addBoundaryFaceConnections(
                        connections,
                        0,
                        localCell2,
                        vals2,
                        connectedStatesInterProc,
                        addFace);
                }
                else if (interProcLevel == 3)
                {
                    this->addBoundaryFaceConnections(
                        connections,
                        0,
                        localCell2,
                        vals3,
                        connectedStatesInterProc,
                        addFace);
                }
                else
                {
                    FatalErrorIn("interProcLevel not valid") << abort(FatalError);
                }
            }
        }
    }
    else
    {
        FatalErrorIn("connectedLevelLocal not valid") << abort(FatalError);
    }

    return;
}

void DAJacCon::setConnections(
    Mat conMat,
    const label idx) const
{

    /*
    set 1.0 for conMat, the column index is idx, the row index is
    always 1 because conMat is a row matrix
    */

    PetscInt idxI = 0;
    PetscScalar v = 1;
    MatSetValues(conMat, 1, &idxI, 1, &idx, &v, INSERT_VALUES);
    return;
}

void DAJacCon::calcNeiBFaceGlobalCompact(labelList& neiBFaceGlobalCompact)
{
    /*
    This function calculates DAJacCon::neiBFaceGlobalCompact[bFaceI]. Here neiBFaceGlobalCompact 
    stores the global coupled boundary face index for the face on the other side of the local 
    processor boundary. bFaceI is the "compact" face index. bFaceI=0 for the first boundary face
    neiBFaceGlobalCompat.size() = nLocalBoundaryFaces
    neiBFaceGlobalCompact[bFaceI] = -1 means it is not a coupled face
    NOTE: neiBFaceGlobalCompact will be used to calculate the connectivity across processors
    in DAJacCon::setupStateBoundaryCon

    Output:
    ------
    neiBFaceGlobalCompact: the global coupled boundary face index for the face on the other 
    side of the local processor boundary
   
    Example:
    -------
    On proc0, neiBFaceGlobalCompact[0] = 1024, then we have the following:
   
                         localBFaceI = 0     <--proc0
                  ---------------------------   coupled boundary face
                         globalBFaceI=1024   <--proc1   
    Taken and modified from the extended stencil code in fvMesh
    Swap the global boundary face index
    */

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    neiBFaceGlobalCompact.setSize(daIndex_.nLocalBoundaryFaces);

    // initialize the list with -1, i.e., non coupled face
    forAll(neiBFaceGlobalCompact, idx)
    {
        neiBFaceGlobalCompact[idx] = -1;
    }

    // loop over the patches and store the global indices
    label counter = 0;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            // For coupled faces set the global face index so that it can be
            // swaped across the interface.
            forAll(pp, i)
            {
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                neiBFaceGlobalCompact[bFaceI] = daIndex_.globalCoupledBFaceNumbering.toGlobal(counter);
                faceIStart++;
                counter++;
            }
        }
    }

    // Swap the cell indices, the list now contains the global index for the
    // U state for the cell on the other side of the processor boundary
    syncTools::swapBoundaryFaceList(mesh_, neiBFaceGlobalCompact);

    return;
}

label DAJacCon::getLocalCoupledBFaceIndex(const label localFaceI) const
{
    /*
    Calculate the index of the local inter-processor boundary face (bRow). 
    
    Input:
    -----
    localFaceI: The local face index. It is in a list of faces including all the
    internal and boundary faces.

    Output:
    ------
    bRow: A list of faces starts with the first inter-processor face. 
    See DAJacCon::globalBndNumbering_ for more details.
    */

    label counter = 0;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            // get the start index of this patch in the global face
            // list and the size of this patch.
            label faceStart = pp.start();
            label patchSize = pp.size();
            label faceEnd = faceStart + patchSize;
            if (localFaceI >= faceStart && localFaceI < faceEnd)
            {
                // this face is on this patch, find the exact index
                label countDelta = localFaceI - pp.start(); //-faceStart;
                PetscInt bRow = counter + countDelta;
                return bRow;
            }
            else
            {
                //increment the counter by patchSize
                counter += patchSize;
            }
        }
    }

    // no match found
    FatalErrorIn("getLocalBndFaceIndex") << abort(FatalError);
    return -1;
}

void DAJacCon::setupStateBoundaryCon(Mat* stateBoundaryCon)
{
    /*
    This function calculates DAJacCon::stateBoundaryCon_

    Output:
    -------
    stateBoundaryCon stores the level of connected states (on the other side 
    across the boundary) for a given coupled boundary face. stateBoundaryCon is 
    a matrix with sizes of nGlobalCoupledBFaces by nGlobalAdjointStates
    stateBoundaryCon is mainly used in the addBoundaryFaceConnection function
    
    Example:
    --------
    Basically, if there are 2 levels of connected states across the inter-proc boundary

                                   |<-----------proc0, globalBFaceI=1024
                           -----------------------------------  <-coupled boundary face
     globalAdjStateIdx=100 ->   | lv1 | <------ proc1
                                |_____|
     globalAdjStateIdx=200 ->   | lv2 |
                                |_____| 
                               
    
    The indices for row 1024 in the stateBoundaryCon matrix will be
    stateBoundaryCon
    rowI=1024   
    Cols: colI=0 ...... colI=100  ........ colI=200 ......... colI=nGlobalAdjointStates
    Vals (level):           1                 2           
    NOTE: globalBFaceI=1024 is owned by proc0      
           
    */

    const HashTable<wordList>& regStates = daRegState_.getRegStates();

    MatCreate(PETSC_COMM_WORLD, stateBoundaryCon);
    MatSetSizes(
        *stateBoundaryCon,
        daIndex_.nLocalCoupledBFaces,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*stateBoundaryCon);
    MatMPIAIJSetPreallocation(*stateBoundaryCon, 1000, NULL, 1000, NULL);
    MatSeqAIJSetPreallocation(*stateBoundaryCon, 1000, NULL);
    MatSetOption(*stateBoundaryCon, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*stateBoundaryCon);
    MatZeroEntries(*stateBoundaryCon);

    Mat stateBoundaryConTmp;
    MatCreate(PETSC_COMM_WORLD, &stateBoundaryConTmp);
    MatSetSizes(
        stateBoundaryConTmp,
        daIndex_.nLocalCoupledBFaces,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(stateBoundaryConTmp);
    MatMPIAIJSetPreallocation(stateBoundaryConTmp, 1000, NULL, 1000, NULL);
    MatSeqAIJSetPreallocation(stateBoundaryConTmp, 1000, NULL);
    MatSetOption(stateBoundaryConTmp, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(stateBoundaryConTmp);
    MatZeroEntries(stateBoundaryConTmp);

    // loop over the patches and set the boundary connnectivity
    // Add connectivity in reverse so that the nearer stencils take priority

    // NOTE: we need to start with level 3, then to 2, then to 1, and flush the matrix
    // for each level before going to another level This is necessary because
    // we need to make sure a proper INSERT_VALUE behavior in MatSetValues
    // i.e., we found that if you use INSERT_VALUE to insert different values (e.g., 1, 2, and 3)
    // to a same rowI and colI in MatSetValues, and call Mat_Assembly in the end. The, the actual
    // value in rowI and colI is kind of random, it does not depend on which value is
    // insert first, in this case, it can be 1, 2, or 3... This happens only in parallel and
    // only happens after Petsc-3.8.4

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    // level 3 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];

                // This cell is already a neighbour cell, so we need this plus two
                // more levels
                // Start with next to nearest neighbours
                forAll(mesh_.cellCells()[idxN], cellI)
                {
                    label localCell = mesh_.cellCells()[idxN][cellI];
                    forAll(daIndex_.adjStateNames, idxI)
                    {
                        word stateName = daIndex_.adjStateNames[idxI];
                        if (daIndex_.adjStateType[stateName] != "surfaceScalarState")
                        {
                            // Now add level 3 connectivity, add all vars except for
                            // surfaceScalarStates
                            this->addConMatNeighbourCells(
                                *stateBoundaryCon,
                                gRow,
                                localCell,
                                stateName,
                                3.0);
                            this->addConMatNeighbourCells(
                                stateBoundaryConTmp,
                                gRow,
                                localCell,
                                stateName,
                                3.0);
                        }
                    }
                }
            }
        }
    }
    // NOTE: need to flush the value before assigning the next level
    MatAssemblyBegin(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyBegin(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 2 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];

                // now add the nearest neighbour cells, add all vars for level 2 except
                // for surfaceScalarStates
                forAll(daIndex_.adjStateNames, idxI)
                {
                    word stateName = daIndex_.adjStateNames[idxI];
                    if (daIndex_.adjStateType[stateName] != "surfaceScalarState")
                    {
                        this->addConMatNeighbourCells(
                            *stateBoundaryCon,
                            gRow,
                            idxN,
                            stateName,
                            2.0);
                        this->addConMatNeighbourCells(
                            stateBoundaryConTmp,
                            gRow,
                            idxN,
                            stateName,
                            2.0);
                    }
                }
            }
        }
    }
    // NOTE: need to flush the value before assigning the next level
    MatAssemblyBegin(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyBegin(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // face con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];
                // and add the surfaceScalarStates for idxN
                forAll(regStates["surfaceScalarStates"], idxI)
                {
                    word stateName = regStates["surfaceScalarStates"][idxI];
                    this->addConMatCellFaces(
                        *stateBoundaryCon,
                        gRow,
                        idxN,
                        stateName,
                        10.0); // for faces, its connectivity level is 10
                    this->addConMatCellFaces(
                        stateBoundaryConTmp,
                        gRow,
                        idxN,
                        stateName,
                        10.0);
                }
            }
        }
    }
    // NOTE: need to flush the value before assigning the next level
    MatAssemblyBegin(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FLUSH_ASSEMBLY);
    MatAssemblyBegin(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 1 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];
                // Add all the cell states for idxN
                forAll(daIndex_.adjStateNames, idxI)
                {
                    word stateName = daIndex_.adjStateNames[idxI];
                    if (daIndex_.adjStateType[stateName] != "surfaceScalarState")
                    {
                        this->addConMatCell(
                            *stateBoundaryCon,
                            gRow,
                            idxN,
                            stateName,
                            1.0);
                        this->addConMatCell(
                            stateBoundaryConTmp,
                            gRow,
                            idxN,
                            stateName,
                            1.0);
                    }
                }
            }
        }
    }
    // Now we can do the final assembly
    MatAssemblyBegin(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);

    // Now repeat loop adding boundary connections from other procs using matrix
    // created in the first loop.
    // Add connectivity in reverse so that the nearer stencils take priority

    // level 3 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];

                // This cell is already a neighbour cell, so we need this plus two
                // more levels
                // Start with nearest neighbours
                forAll(mesh_.cellCells()[idxN], cellI)
                {
                    label localCell = mesh_.cellCells()[idxN][cellI];
                    labelList val1 = {3};
                    // pass a zero list to add all states
                    List<List<word>> connectedStates(0);
                    this->addBoundaryFaceConnections(
                        stateBoundaryConTmp,
                        gRow,
                        localCell,
                        val1,
                        connectedStates,
                        0);
                }
            }
        }
    }
    // NOTE: need to flush the value before assigning the next level
    MatAssemblyBegin(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 2 and 3 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];
                // now add the neighbour cells
                labelList vals2 = {2, 3};
                // pass a zero list to add all states
                List<List<word>> connectedStates(0);
                this->addBoundaryFaceConnections(
                    stateBoundaryConTmp,
                    gRow,
                    idxN,
                    vals2,
                    connectedStates,
                    0);
            }
        }
    }
    // NOTE: need to flush the value before assigning the next level
    MatAssemblyBegin(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 2 again, because the previous call will mess up level 2 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        // get the start index of this patch in the global face list
        label faceIStart = pp.start();

        // check whether this face is coupled (cyclic or processor?)
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                // get the necessary matrix row
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];

                // Now get the cell that borders this coupled bFace
                label idxN = pFaceCells[faceI];
                // now add the neighbour cells
                labelList vals1 = {2};
                // pass a zero list to add all states
                List<List<word>> connectedStates(0);
                this->addBoundaryFaceConnections(
                    stateBoundaryConTmp,
                    gRow,
                    idxN,
                    vals1,
                    connectedStates,
                    0);
            }
        }
    }

    MatAssemblyBegin(stateBoundaryConTmp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stateBoundaryConTmp, MAT_FINAL_ASSEMBLY);

    // the above repeat loop is not enough to cover all the stencil, we need to do more
    this->combineStateBndCon(stateBoundaryCon, &stateBoundaryConTmp);

    return;
}

void DAJacCon::combineStateBndCon(
    Mat* stateBoundaryCon,
    Mat* stateBoundaryConTmp)
{
    /*
    1. Add additional adj state connectivities if the stateBoundaryCon stencil extends through
    three or more decomposed domains, something like this:
    
    --------       ---------
           |       |       |
      Con3 |  Con2 |  Con1 |  R
           |       |       |
           ---------       --------
           
    Here R is the residual, Con1 to 3 are its connectivity, and dashed lines 
    are the inter-processor boundary
           
    2. Assign stateBoundaryConTmp to stateBoundaryCon.

    Input/Output:
    stateBoundaryCon, and stateBoundaryConTmp should come from DAJacCon::stateBoundaryCon
    */

    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    PetscInt nCols1;
    const PetscInt* cols1;
    const PetscScalar* vals1;

    // Destroy and initialize stateBoundaryCon with zeros
    MatDestroy(stateBoundaryCon);
    MatCreate(PETSC_COMM_WORLD, stateBoundaryCon);
    MatSetSizes(
        *stateBoundaryCon,
        daIndex_.nLocalCoupledBFaces,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*stateBoundaryCon);
    MatMPIAIJSetPreallocation(*stateBoundaryCon, 1000, NULL, 1000, NULL);
    MatSeqAIJSetPreallocation(*stateBoundaryCon, 1000, NULL);
    MatSetOption(*stateBoundaryCon, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*stateBoundaryCon);
    MatZeroEntries(*stateBoundaryCon); // initialize with zeros

    // assign stateBoundaryConTmp to stateBoundaryCon
    PetscInt Istart, Iend;
    MatGetOwnershipRange(*stateBoundaryConTmp, &Istart, &Iend);
    for (PetscInt i = Istart; i < Iend; i++)
    {
        MatGetRow(*stateBoundaryConTmp, i, &nCols, &cols, &vals);
        for (PetscInt j = 0; j < nCols; j++)
        {
            MatSetValue(*stateBoundaryCon, i, cols[j], vals[j], INSERT_VALUES);
        }
        MatRestoreRow(*stateBoundaryConTmp, i, &nCols, &cols, &vals);
    }
    MatAssemblyBegin(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);

    MatDestroy(stateBoundaryConTmp);

    // copy ConMat to ConMatTmp for an extract loop
    MatConvert(*stateBoundaryCon, MATSAME, MAT_INITIAL_MATRIX, stateBoundaryConTmp);

    // We need to do another loop adding boundary connections from other procs using ConMat
    // this will add missing connectivity if the stateBoundaryCon stencil extends through
    // three more more processors
    // NOTE: we need to start with level 3, then to 2, then to 1, and flush the matrix
    // for each level before going to another level This is necessary because
    // we need to make sure a proper INSERT_VALUE behavior in MatSetValues
    // i.e., we found that if you use INSERT_VALUE to insert different values (e.g., 1, 2, and 3)
    // to a same rowI and colI in MatSetValues, and call Mat_Assembly in the end. The, the actual
    // value in rowI and colI is kind of random, it does not depend on which value is
    // insert first, in this case, it can be 1, 2, or 3... This happens only in parallel and
    // only happens after Petsc-3.8.4

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    // level 3 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        label faceIStart = pp.start();
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];
                label idxN = pFaceCells[faceI];

                forAll(mesh_.cellCells()[idxN], cellI)
                {
                    label localCell = mesh_.cellCells()[idxN][cellI];
                    labelList val1 = {3};
                    // pass a zero list to add all states
                    List<List<word>> connectedStates(0);
                    this->addBoundaryFaceConnections(
                        *stateBoundaryConTmp,
                        gRow,
                        localCell,
                        val1,
                        connectedStates,
                        0);
                }
            }
        }
    }
    MatAssemblyBegin(*stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 2, 3 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        label faceIStart = pp.start();
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];
                label idxN = pFaceCells[faceI];
                // now add the neighbour cells
                labelList vals2 = {2, 3};
                // pass a zero list to add all states
                List<List<word>> connectedStates(0);
                this->addBoundaryFaceConnections(
                    *stateBoundaryConTmp,
                    gRow,
                    idxN,
                    vals2,
                    connectedStates,
                    0);
            }
        }
    }
    MatAssemblyBegin(*stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryConTmp, MAT_FLUSH_ASSEMBLY);

    // level 2 again, because the previous call will mess up level 2 con
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const UList<label>& pFaceCells = pp.faceCells();
        label faceIStart = pp.start();
        if (pp.coupled())
        {
            forAll(pp, faceI)
            {
                label bFaceI = faceIStart - daIndex_.nLocalInternalFaces;
                faceIStart++;
                label gRow = neiBFaceGlobalCompact_[bFaceI];
                label idxN = pFaceCells[faceI];
                // now add the neighbour cells
                labelList vals1 = {2};
                // pass a zero list to add all states
                List<List<word>> connectedStates(0);
                this->addBoundaryFaceConnections(
                    *stateBoundaryConTmp,
                    gRow,
                    idxN,
                    vals1,
                    connectedStates,
                    0);
            }
        }
    }
    MatAssemblyBegin(*stateBoundaryConTmp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryConTmp, MAT_FINAL_ASSEMBLY);

    // Now stateBoundaryConTmp will have all the missing stencil. However, it will also mess
    // up the existing stencil in stateBoundaryCon. So we need to do a check to make sure that
    // stateBoundaryConTmp only add stencil, not replacing any existing stencil in stateBoundaryCon.
    // If anything in stateBoundaryCon is replaced, rollback the changes.
    Mat tmpMat; // create a temp mat
    MatCreate(PETSC_COMM_WORLD, &tmpMat);
    MatSetSizes(
        tmpMat,
        daIndex_.nLocalCoupledBFaces,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(tmpMat);
    MatMPIAIJSetPreallocation(tmpMat, 1000, NULL, 1000, NULL);
    MatSeqAIJSetPreallocation(tmpMat, 1000, NULL);
    MatSetOption(tmpMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(tmpMat);
    MatZeroEntries(tmpMat); // initialize with zeros
    for (PetscInt i = Istart; i < Iend; i++)
    {
        MatGetRow(*stateBoundaryCon, i, &nCols, &cols, &vals);
        MatGetRow(*stateBoundaryConTmp, i, &nCols1, &cols1, &vals1);
        for (PetscInt j = 0; j < nCols1; j++)
        {
            // for each col in stateBoundaryConTmp, we need to check if there are any existing
            // values for the same col in stateBoundaryCon. If yes, assign the val from
            // stateBoundaryCon instead of stateBoundaryConTmp
            PetscScalar newVal = vals1[j];
            PetscInt newCol = cols1[j];
            for (PetscInt k = 0; k < nCols; k++)
            {
                if (int(cols[k]) == int(cols1[j]))
                {
                    newVal = vals[k];
                    newCol = cols[k];
                    break;
                }
            }
            MatSetValue(tmpMat, i, newCol, newVal, INSERT_VALUES);
        }
        MatRestoreRow(*stateBoundaryCon, i, &nCols, &cols, &vals);
        MatRestoreRow(*stateBoundaryConTmp, i, &nCols1, &cols1, &vals1);
    }
    MatAssemblyBegin(tmpMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tmpMat, MAT_FINAL_ASSEMBLY);

    // copy ConMat to ConMatTmp
    MatDestroy(stateBoundaryCon);
    MatConvert(tmpMat, MATSAME, MAT_INITIAL_MATRIX, stateBoundaryCon);
    MatAssemblyBegin(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryCon, MAT_FINAL_ASSEMBLY);

    MatDestroy(stateBoundaryConTmp);
    MatDestroy(&tmpMat);

    return;
}

void DAJacCon::setupStateBoundaryConID(Mat* stateBoundaryConID)
{
    /*
    This function computes DAJacCon::stateBoundaryConID_.

    Output:
    -------
    stateBoundaryConID: it has the exactly same structure as DAJacConstateBoundaryCon_ 
    except that stateBoundaryConID stores the connected stateID instead of connected 
    levels. stateBoundaryConID will be used in DAJacCon::addBoundaryFaceConnections
    */

    PetscInt nCols, colI;
    const PetscInt* cols;
    const PetscScalar* vals;
    PetscInt Istart, Iend;

    PetscScalar valIn;

    // assemble adjStateID4GlobalAdjIdx
    // adjStateID4GlobalAdjIdx stores the adjStateID for given a global adj index
    labelList adjStateID4GlobalAdjIdx;
    adjStateID4GlobalAdjIdx.setSize(daIndex_.nGlobalAdjointStates);
    daIndex_.calcAdjStateID4GlobalAdjIdx(adjStateID4GlobalAdjIdx);

    // initialize
    MatCreate(PETSC_COMM_WORLD, stateBoundaryConID);
    MatSetSizes(
        *stateBoundaryConID,
        daIndex_.nLocalCoupledBFaces,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*stateBoundaryConID);
    MatMPIAIJSetPreallocation(*stateBoundaryConID, 1000, NULL, 1000, NULL);
    MatSeqAIJSetPreallocation(*stateBoundaryConID, 1000, NULL);
    MatSetOption(*stateBoundaryConID, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*stateBoundaryConID);
    MatZeroEntries(*stateBoundaryConID);

    MatGetOwnershipRange(stateBoundaryCon_, &Istart, &Iend);

    // set stateBoundaryConID_ based on stateBoundaryCon_ and adjStateID4GlobalAdjIdx
    for (PetscInt i = Istart; i < Iend; i++)
    {
        MatGetRow(stateBoundaryCon_, i, &nCols, &cols, &vals);
        for (PetscInt j = 0; j < nCols; j++)
        {
            if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
            {
                colI = cols[j];
                valIn = adjStateID4GlobalAdjIdx[colI];
                MatSetValue(*stateBoundaryConID, i, colI, valIn, INSERT_VALUES);
            }
        }
        MatRestoreRow(stateBoundaryCon_, i, &nCols, &cols, &vals);
    }

    adjStateID4GlobalAdjIdx.clear();

    MatAssemblyBegin(*stateBoundaryConID, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*stateBoundaryConID, MAT_FINAL_ASSEMBLY);

    return;
}

void DAJacCon::addConMatCell(
    Mat conMat,
    const label gRow,
    const label cellI,
    const word stateName,
    const PetscScalar val)
{

    /* 
    Insert a value (val) to the connectivity Matrix (conMat)
    This value will be inserted at rowI=gRow
    The column index is dependent on the cellI and stateName

    Input:
    -----
    gRow: which row to insert the value for conMat
    cellI: the index of the cell to compute the column index to add
    stateName: the name of the state variable to compute the column index to add
    val: the value to add to conMat

    Output:
    ------
    conMat: the matrix to add value to

    Example:
    -------

    If we want to add value 1.0 to conMat for
    column={the U globalAdjointIndice of cellI} where cellI=5
    row = gRow = 100
    Then, call addConMatCell(conMat, 100, 5, "U", 1.0)

             -------  
            | cellI | <----------- value 1.0 will be added to 
             -------               column = {global index of U 
                                   for cellI}
    */

    PetscInt idxJ, idxI;

    idxI = gRow;

    // find the global index of this state
    label compMax = 1;
    if (daIndex_.adjStateType[stateName] == "volVectorState")
    {
        compMax = 3;
    }

    for (label i = 0; i < compMax; i++)
    {
        idxJ = daIndex_.getGlobalAdjointStateIndex(stateName, cellI, i);
        // set it in the matrix
        MatSetValues(conMat, 1, &idxI, 1, &idxJ, &val, INSERT_VALUES);
    }

    return;
}

void DAJacCon::addConMatNeighbourCells(
    Mat conMat,
    const label gRow,
    const label cellI,
    const word stateName,
    const PetscScalar val)
{

    /* 
    Insert a value (val) to the connectivity Matrix (conMat)
    This value will be inserted at rowI=gRow
    The column index is dependent on the cellI, and cellI's neibough and stateName

    Input:
    -----
    gRow: which row to insert the value for conMat
    cellI: the index of the cell to compute the column index to add
    stateName: the name of the state variable to compute the column index to add
    val: the value to add to conMat

    Output:
    ------
    conMat: the matrix to add value to

    Example:
    -------

    If we want to add value 1.0 to conMat for
    columns={the U globalAdjointIndice of all the neiboughs of cellI} where cellI=5
    row = gRow = 100
    Then, call addConMatNeighbourCells(conMat, 100, 5, "U", 1.0)

             -------  
            | cellL | <----------- value 1.0 will be added to 
     ------- ------- -------       column = {global index of U 
    | cellJ | cellI | cellK |      for cellL}, similarly for all 
     ------- ------- -------       the neiboughs of cellI
            | cellM |
             -------
    */

    label localCellJ;
    PetscInt idxJ, idxI;

    idxI = gRow;
    // Add the nearest neighbour cells for cell
    forAll(mesh_.cellCells()[cellI], cellJ)
    {
        // get the local neighbour cell
        localCellJ = mesh_.cellCells()[cellI][cellJ];

        // find the global index of this state
        label compMax = 1;
        if (daIndex_.adjStateType[stateName] == "volVectorState")
        {
            compMax = 3;
        }
        for (label i = 0; i < compMax; i++)
        {
            idxJ = daIndex_.getGlobalAdjointStateIndex(stateName, localCellJ, i);
            // set it in the matrix
            MatSetValues(conMat, 1, &idxI, 1, &idxJ, &val, INSERT_VALUES);
        }
    }

    return;
}

void DAJacCon::addConMatCellFaces(
    Mat conMat,
    const label gRow,
    const label cellI,
    const word stateName,
    const PetscScalar val)
{

    /* 
    Insert a value (val) to the connectivity Matrix (conMat)
    This value will be inserted at rowI=gRow
    The column index is dependent on the cellI's faces and stateName

    Input:
    -----
    gRow: which row to insert the value for conMat
    cellI: the index of the cell to compute the column index to add
    stateName: the name of the state variable to compute the column index to add
    val: the value to add to conMat

    Output:
    ------
    conMat: the matrix to add value to

    Example:
    -------

    If we want to add value 10.0 to conMat for
    columns={the phi globalAdjointIndice of cellI's faces} where cellI=5
    row = gRow = 100
    Then, call addConMatCell(conMat, 100, 5, "U", 1.0)

             -------  
            | cellI | <----------- value 10.0 will be added to 
             -------               column = {global adjoint index 
                                   of all cellI's faces}
    */

    PetscInt idxJ, idxI;
    idxI = gRow;

    // get the faces connected to this cell, note these are in a single
    // list that includes all internal and boundary faces
    const labelList& faces = mesh_.cells()[cellI];
    forAll(faces, idx)
    {
        //get the appropriate index for this face
        label globalState = daIndex_.getGlobalAdjointStateIndex(stateName, faces[idx]);
        idxJ = globalState;
        MatSetValues(conMat, 1, &idxI, 1, &idxJ, &val, INSERT_VALUES);
    }

    return;
}

void DAJacCon::addBoundaryFaceConnections(
    Mat conMat,
    const label gRow,
    const label cellI,
    const labelList v,
    const List<List<word>> connectedStates,
    const label addFaces)
{
    /*
    This function adds inter-proc connectivity into conMat.
    For all the inter-proc faces owned by cellI, get the global adj state indices 
    from DAJacCon::stateBoundaryCon_ and then add them into conMat
    Col index to add: the same col index for a given row (bRowGlobal) in the stateBoundaryCon 
    mat if the element value in the stateBoundaryCon mat is less than the input level, 
    i.e., v.size().
    
    Input:
    -----
    gRow: Row index to add

    cellI: the cell index for getting the faces to add inter-proc connectivity, NOTE: depending on the level
    of requested connections, we may add inter-proc face that are not belonged to cellI

    v: an array denoting the desired values to add, the size of v denotes the maximal levels to add

    connectedStates: selectively add some states into the conMat for the current level. If its size is 0,
    add all the possible states (except for surfaceStates). The dimension of connectedStates is nLevel 
    by nStates.

    addFaces: whether to add indices for face (phi) connectivity
    
    Example:
    --------
    
    labelList val2={1,2};
    PetscInt gRow=1024, idxN = 100, addFaces=1;
    wordListList connectedStates={{"U","p"},{"U"}};
    addBoundaryFaceConnections(stateBoundaryCon,gRow,idxN,vals2,connectedStates,addFaces);
    The above call will add 2 levels of connected states for all the inter-proc faces belonged to cellI=idxN
    The cols to added are: the level1 connected states (U, p) for all the inter-proc faces belonged to 
    cellI=idxN. the level2 connected states (U only) for all the inter-proc faces belonged to cellI=idxN
    The valus "1" will be added to conMat for all the level1 connected states while the value "2" will be 
    added for level2. 
    Note: this function will also add all faces belonged to level1 of the inter-proc faces, see the 
    following for reference
    
                                -------
                                | idxN|
                                |     |       proc0, idxN=100, globalBFaceI=1024 for the south face of idxN
                           -----------------  <----coupled boundary face
     add state U and p ----->   | lv1 |       proc1
     also add faces -------->   |     |
                                -------
                                | lv2 |
     add state U  ---------->   |     |
                                -------
    

    ****** NOTE: *******
    If the inter-proc boundary is like the following, calling this function will NOT add any 
    inter-proc connection for idxN because there is no inter-proc boundary for cell idxN

            -------
            | idxN|
            |     |       
            -------  <--------- there is no inter-proc boundary for idxN, not adding anything
            | lv1 |      
            |     |        proc0
      ------------------- <----coupled boundary face
            | lv2 |        proc1
            |     |
            -------
    
    */

    if (v.size() != connectedStates.size() && connectedStates.size() != 0)
    {
        FatalErrorIn("") << "size of v and connectedStates are not identical!" << abort(FatalError);
    }

    PetscInt idxJ, idxI, bRow, bRowGlobal;
    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    PetscInt nColsID;
    const PetscInt* colsID;
    const PetscScalar* valsID;

    // convert stateNames to stateIDs
    labelListList connectedStateIDs(connectedStates.size());
    forAll(connectedStates, idxI)
    {
        forAll(connectedStates[idxI], idxJ)
        {
            word stateName = connectedStates[idxI][idxJ];
            label stateID = daIndex_.adjStateID[stateName];
            connectedStateIDs[idxI].append(stateID);
        }
    }

    idxI = gRow;
    // get the faces connected to this cell, note these are in a single
    // list that includes all internal and boundary faces
    const labelList& faces = mesh_.cells()[cellI];

    //get the level
    label level = v.size();

    for (label lv = level; lv >= 1; lv--) // we need to start from the largest levels since they have higher priority
    {
        forAll(faces, faceI)
        {
            // Now deal with coupled faces
            label currFace = faces[faceI];

            if (daIndex_.isCoupledFace[currFace])
            {
                //this is a coupled face

                // use the boundary connectivity to figure out what is connected
                // to this face for this level

                // get bRow in boundaryCon for this face
                bRow = this->getLocalCoupledBFaceIndex(currFace);

                // get the global bRow index
                bRowGlobal = daIndex_.globalCoupledBFaceNumbering.toGlobal(bRow);

                // now extract the boundaryCon row
                MatGetRow(stateBoundaryCon_, bRowGlobal, &nCols, &cols, &vals);
                if (connectedStates.size() != 0)
                {
                    // check if we need to get stateID
                    MatGetRow(stateBoundaryConID_, bRowGlobal, &nColsID, &colsID, &valsID);
                }

                // now loop over the row and set any column that match this level
                // in conMat
                for (label i = 0; i < nCols; i++)
                {
                    idxJ = cols[i];
                    label val = round(vals[i]); // val is the connectivity level extracted from stateBoundaryCon_ at this col
                    // selectively add some states into conMat
                    label addState;
                    label stateID = -9999;
                    // check if we need to get stateID
                    if (connectedStates.size() != 0)
                    {
                        stateID = round(valsID[i]);
                    }

                    if (connectedStates.size() == 0)
                    {
                        addState = 1;
                    }
                    else if (daUtil_.isInList<label>(stateID, connectedStateIDs[lv - 1]))
                    {
                        addState = 1;
                    }
                    else
                    {
                        addState = 0;
                    }
                    // if the level match and the state is what you want
                    if (val == lv && addState)
                    {
                        // need to do v[lv-1] here since v is an array with starting index 0
                        PetscScalar valIn = v[lv - 1];
                        MatSetValues(conMat, 1, &idxI, 1, &idxJ, &valIn, INSERT_VALUES);
                    }
                    if (val == 10 && addFaces)
                    {
                        // this is a necessary connection
                        PetscScalar valIn = v[lv - 1];
                        MatSetValues(conMat, 1, &idxI, 1, &idxJ, &valIn, INSERT_VALUES);
                    }
                }

                // restore the row of the matrix
                MatRestoreRow(stateBoundaryCon_, bRowGlobal, &nCols, &cols, &vals);
                if (connectedStates.size() != 0)
                {
                    // check if we need to get stateID
                    MatRestoreRow(stateBoundaryConID_, bRowGlobal, &nColsID, &colsID, &valsID);
                }
            }
        }
    }

    return;
}

void DAJacCon::initializePetscVecs()
{
    // initialize the preallocation vecs
    VecCreate(PETSC_COMM_WORLD, &dRdWTPreallocOn_);
    VecSetSizes(dRdWTPreallocOn_, daIndex_.nLocalAdjointStates, PETSC_DECIDE);
    VecSetFromOptions(dRdWTPreallocOn_);

    VecDuplicate(dRdWTPreallocOn_, &dRdWTPreallocOff_);
    VecDuplicate(dRdWTPreallocOn_, &dRdWPreallocOn_);
    VecDuplicate(dRdWTPreallocOn_, &dRdWPreallocOff_);

    // initialize coloring vectors

    //dRdW Colors
    VecCreate(PETSC_COMM_WORLD, &dRdWColors_);
    VecSetSizes(dRdWColors_, daIndex_.nLocalAdjointStates, PETSC_DECIDE);
    VecSetFromOptions(dRdWColors_);
    VecDuplicate(dRdWColors_, &dRdWColoredColumns_);

    //dFdW Colors, the dFdWColoredColumns will be initialized in initializedFdWCon
    VecCreate(PETSC_COMM_WORLD, &dFdWColors_);
    VecSetSizes(dFdWColors_, daIndex_.nLocalAdjointStates, PETSC_DECIDE);
    VecSetFromOptions(dFdWColors_);

    return;
}

void DAJacCon::initializedRdWCon(const label isPC)
{

    /*
    Initialize the connectivity matrix and preallocate memory
    */

    Mat* conMat;
    if (isPC)
    {
        conMat = &dRdWConPC_;
    }
    else
    {
        conMat = &dRdWCon_;
    }

    MatCreate(PETSC_COMM_WORLD, conMat);
    MatSetSizes(
        *conMat,
        daIndex_.nLocalAdjointStates,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(*conMat);

    this->preallocateJacobianMatrix(
        *conMat,
        dRdWPreallocOn_,
        dRdWPreallocOff_);
    //MatSetOption(*conMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(*conMat);

    Info << "Connectivity matrix initialized." << endl;
}

void DAJacCon::preallocateJacobianMatrix(
    Mat dRMat,
    const Vec preallocOnProc,
    const Vec preallocOffProc)
{
    /*
    Preallocate memory for dRMat.

    Input:
    ------
    preallocOnProc, preallocOffProc: the on and off diagonal nonzeros for
    dRMat

    Output:
    dRMat: matrix to preallocate memory for
    */

    PetscScalar *onVec, *offVec;
    PetscInt onSize[daIndex_.nLocalAdjointStates], offSize[daIndex_.nLocalAdjointStates];

    VecGetArray(preallocOnProc, &onVec);
    VecGetArray(preallocOffProc, &offVec);
    for (label i = 0; i < daIndex_.nLocalAdjointStates; i++)
    {
        onSize[i] = round(onVec[i]);
        if (onSize[i] > daIndex_.nLocalAdjointStates)
        {
            onSize[i] = daIndex_.nLocalAdjointStates;
        }
        offSize[i] = round(offVec[i]) + 5; // reserve 5 more?
    }

    VecRestoreArray(preallocOnProc, &onVec);
    VecRestoreArray(preallocOffProc, &offVec);

    // MatMPIAIJSetPreallocation(dRMat,NULL,preallocOnProc,NULL,preallocOffProc);
    // MatSeqAIJSetPreallocation(dRMat,NULL,preallocOnProc);

    MatMPIAIJSetPreallocation(dRMat, NULL, onSize, NULL, offSize);
    MatSeqAIJSetPreallocation(dRMat, NULL, onSize);

    return;
}

void DAJacCon::preallocatedRdW(
    Mat dRMat,
    const label transposed)
{
    /*
    Call the DAJacCon::preallocateJacobianMatrix function with the 
    correct vectors, depending on transposed

    Input:
    -----
    transposed: whether the state Jacobian mat is transposed, i.e., it
    is for dRdW or dRdWT (transposed)

    Output:
    ------
    dRMat: the matrix to preallocate
    */
    if (transposed)
    {
        this->preallocateJacobianMatrix(dRMat, dRdWTPreallocOn_, dRdWTPreallocOff_);
    }
    else
    {
        this->preallocateJacobianMatrix(dRMat, dRdWPreallocOn_, dRdWPreallocOff_);
    }
}

void DAJacCon::initializedFdWCon(const label nLocalObjFuncGeoElements)
{
    /*
    Initialize DAJacCon::dFdWCon_

    Output:
    ------
    dFdWCon_: connectivity matrix for dFdW, here dFdWCon has a
    size of nLocalObjFuncGeoElements * nGlobalAdjointStates
    The reason that dFdWCon has nLocalObjFuncGeoElements rows is 
    because we need to divide the objective function into 
    nLocalObjFuncGeoElements discrete value such that we can
    use coloring to compute dFdW
    */

    MatCreate(PETSC_COMM_WORLD, &dFdWCon_);
    MatSetSizes(
        dFdWCon_,
        nLocalObjFuncGeoElements,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(dFdWCon_);
    MatMPIAIJSetPreallocation(dFdWCon_, 100, NULL, 100, NULL);
    MatSeqAIJSetPreallocation(dFdWCon_, 100, NULL);
    MatSetOption(dFdWCon_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dFdWCon_);
    MatZeroEntries(dFdWCon_);

    VecCreate(PETSC_COMM_WORLD, &dFdWColoredColumns_);
    VecSetSizes(
        dFdWColoredColumns_,
        nLocalObjFuncGeoElements,
        PETSC_DECIDE);
    VecSetFromOptions(dFdWColoredColumns_);
    VecZeroEntries(dFdWColoredColumns_);

    Info << "dFdWCon Created!" << endl;
}

void DAJacCon::readdRdWColoring()
{
    /*
    Read the dRdW coloring from files and 
    compute ndRdWColors. The naming convention for
    coloring vector is coloringVecName_nProcs.bin
    This is necessary because using different CPU
    cores result in different dRdWCon and therefore
    different coloring

    Output:
    ------
    dRdWColors_: read from file
    ndRdWColors: number of dRdW colors
    */

    Info << "Reading dRdW Coloring.." << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dRdWColoring_" << nProcs;
    word fileName1 = fileName.str();
    VecZeroEntries(dRdWColors_);
    daUtil_.readVectorBinary(dRdWColors_, fileName1);

    this->validateColoring(dRdWCon_, dRdWColors_);

    PetscReal maxVal;
    VecMax(dRdWColors_, NULL, &maxVal);
    ndRdWColors = maxVal + 1;
}

void DAJacCon::calcdRdWColoring()
{
    /*
    Calculate the coloring for dRdW.

    Output:
    ------
    dRdWColors_: dRdW coloring and save to files. 
    The naming convention for coloring vector is 
    coloringVecName_nProcs.bin. This is necessary because 
    using different CPU cores result in different dRdWCon 
    and therefore different coloring

    ndRdWColors: number of dRdW colors

    */

    VecZeroEntries(dRdWColors_);
    this->parallelD2Coloring(dRdWCon_, dRdWColors_, ndRdWColors);
    this->validateColoring(dRdWCon_, dRdWColors_);
    Info << " ndRdWColors: " << ndRdWColors << endl;

    // write dRdW colors
    Info << "Writing dRdW Colors.." << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dRdWColoring_" << nProcs;
    word fileName1 = fileName.str();
    daUtil_.writeVectorBinary(dRdWColors_, fileName1);

    return;
}

void DAJacCon::readdFdWColoring(const word objFunc)
{
    /*
    Read the dFdW coloring from files and 
    compute ndFdWColors. The naming convention for
    coloring vector is coloringVecName_nProcs.bin
    This is necessary because using different CPU
    cores result in different dRdWCon and therefore
    different coloring

    Output:
    ------
    dFdWColors_: read from file
    ndFdWColors: number of dFdW colors
    */

    Info << "Reading dFdW Coloring for " << objFunc << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dFdWColoring_" << objFunc << "_" << nProcs;
    word fileName1 = fileName.str();
    VecZeroEntries(dFdWColors_);
    daUtil_.readVectorBinary(dFdWColors_, fileName1);

    this->validateColoring(dFdWCon_, dFdWColors_);

    PetscReal maxVal;
    VecMax(dFdWColors_, NULL, &maxVal);
    ndFdWColors = maxVal + 1;
}

void DAJacCon::calcdFdWColoring(const word objFunc)
{
    /*
    Calculate the coloring for dFdW.

    Output:
    ------
    dFdWColors_: dFdW coloring and save to files. 
    The naming convention for coloring vector is 
    coloringVecName_nProcs.bin. This is necessary because 
    using different CPU cores result in different dRdWCon 
    and therefore different coloring

    ndFdWColors: number of dFdW colors

    */

    VecZeroEntries(dFdWColors_);
    this->parallelD2Coloring(dFdWCon_, dFdWColors_, ndFdWColors);
    this->validateColoring(dFdWCon_, dFdWColors_);
    Info << " ndFdWColors: " << ndFdWColors << endl;

    // write dFdW colors
    Info << "Writing dFdW Colors.." << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dFdWColoring_" << objFunc << "_" << nProcs;
    word fileName1 = fileName.str();
    daUtil_.writeVectorBinary(dFdWColors_, fileName1);

    return;
}

label DAJacCon::getNdRdWColors() const
{
    /*
    Return the number of colors

    Output:
    ------
    ndRdWColors: the number of colors depends on 
    whether the coloring is used
    */

    if (daOption_.getOption<label>("adjUseColoring"))
    {
        return ndRdWColors;
    }
    else
    {
        return daIndex_.nGlobalAdjointStates;
    }

    return -1;
}

label DAJacCon::getNdFdWColors() const
{
    /*
    Return the number of colors

    Output:
    ------
    ndFdWColors: the number of colors depends on 
    whether the coloring is used
    */

    if (daOption_.getOption<label>("adjUseColoring"))
    {
        return ndFdWColors;
    }
    else
    {
        return daIndex_.nGlobalAdjointStates;
    }

    return -1;
}

void DAJacCon::parallelD2Coloring(
    const Mat conMat,
    Vec colors,
    label& nColors) const
{
    /*
    A general function to compute coloring for a Jacobian matrix using a 
    paralel heuristic distance 2 algorithm

    Input:
    -----
    conMat: a Petsc matrix that have the connectivity pattern (value one for 
    all nonzero elements)

    Output:
    -------
    colors: the coloring vector to store the coloring indices, starting with 0
    
    nColors: the number of colors

    Example:
    If the conMat reads,

           color0  color1
             |     |
             1  0  0  0
    conMat = 0  1  1  0
             0  0  1  0
             0  0  0  1
                |     | 
            color0   color0

    Then, calling this function gives colors = {0, 0, 1, 0}.
    This can be done for parallel conMat
    */

    // if we end up having more than 10000 colors, something must be wrong
    label maxColors = 10000;

    PetscInt nCols, nCols2;
    const PetscInt* cols;
    const PetscScalar* vals;
    const PetscScalar* vals2;

    PetscInt colorStart, colorEnd;
    VecScatter colorScatter;

    label Istart, Iend;
    label currColor;
    label notColored = 1;
    IS globalIS;
    label maxCols = 750;
    scalar allNonZeros;
    Vec globalVec;
    PetscInt nRowG, nColG;

    Info << "Parallel Distance 2 Graph Coloring...." << endl;

    // initialize the number of colors to zero
    nColors = 0;

    // get the range of colors owned by the local prock
    VecGetOwnershipRange(colors, &colorStart, &colorEnd);

    // Set the entire color vector to -1
    VecSet(colors, -1);

    // Determine which rows are on the current processor
    MatGetOwnershipRange(conMat, &Istart, &Iend);

    //then get the global number of rows and columns
    MatGetSize(conMat, &nRowG, &nColG);
    label nRowL = Iend - Istart;
    label nColL = colorEnd - colorStart;

    /* 
    Start by looping over the rows to determine the largest
    number of non-zeros per row. This will determine maxCols
    and the minumum bound for the number of colors.
    */
    this->getMatNonZeros(conMat, maxCols, allNonZeros);
    Info << "MaxCols: " << maxCols << endl;
    Info << "AllNonZeros: " << allNonZeros << endl;

    // Create a local sparse matrix with a single row to use as a sparse vector
    Mat localCols;
    MatCreateSeqAIJ(
        PETSC_COMM_SELF,
        1,
        daIndex_.nGlobalAdjointStates,
        daIndex_.nLocalAdjointStates,
        NULL,
        &localCols);
    MatSetOption(localCols, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(localCols);
    MatZeroEntries(localCols);

    //Now loop over the owned rows and set the value in any occupied col to 1.
    PetscInt idxI = 0;
    PetscScalar v = 1;
    for (label i = Istart; i < Iend; i++)
    {
        MatGetRow(conMat, i, &nCols, &cols, &vals);
        // set any columns that have a nonzero entry into localCols
        for (label j = 0; j < nCols; j++)
        {
            if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
            {
                PetscInt idx = cols[j];
                MatSetValues(localCols, 1, &idxI, 1, &idx, &v, INSERT_VALUES);
            }
        }
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }
    MatAssemblyBegin(localCols, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localCols, MAT_FINAL_ASSEMBLY);

    // now localCols contains the unique set of local columns on each processor
    label nUniqueCols = 0;
    MatGetRow(localCols, 0, &nCols, &cols, &vals);
    for (label j = 0; j < nCols; j++)
    {
        if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
        {
            nUniqueCols++;
        }
    }
    Info << "nUniqueCols: " << nUniqueCols << endl;

    //Loop over the local vectors and set nonzero entries in a global vector
    // This lets us determine which columns are strictly local and which have
    // interproccessor overlap.
    VecCreate(PETSC_COMM_WORLD, &globalVec);
    VecSetSizes(globalVec, nColL, PETSC_DECIDE);
    VecSetFromOptions(globalVec);
    VecSet(globalVec, 0);

    for (label j = 0; j < nCols; j++)
    {
        PetscInt idx = cols[j];
        if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
        {
            VecSetValue(globalVec, idx, vals[j], ADD_VALUES);
        }
    }
    VecAssemblyBegin(globalVec);
    VecAssemblyEnd(globalVec);

    MatRestoreRow(localCols, 0, &nCols, &cols, &vals);

    // now create an index set of the strictly local columns
    label* localColumnStat = new label[nUniqueCols];
    label* globalIndexList = new label[nUniqueCols];

    PetscScalar* globalVecArray;
    VecGetArray(globalVec, &globalVecArray);

    label colCounter = 0;
    MatGetRow(localCols, 0, &nCols, &cols, &vals);
    for (label j = 0; j < nCols; j++)
    {
        label col = cols[j];
        if (daUtil_.isValueCloseToRef(vals[j], 1.0))
        {
            globalIndexList[colCounter] = col;
            localColumnStat[colCounter] = 1;
            if (col >= colorStart && col < colorEnd)
            {
                if (daUtil_.isValueCloseToRef(globalVecArray[col - colorStart], 1.0))
                {
                    localColumnStat[colCounter] = 2; // 2: strictly local
                }
            }
            colCounter++;
        }
    }
    MatRestoreRow(localCols, 0, &nCols, &cols, &vals);
    MatDestroy(&localCols);

    // create a list of the rows that have any of the strictly local columns included

    label* localRowList = new label[nRowL];
    for (label i = Istart; i < Iend; i++)
    {
        label idx = i - Istart;
        MatGetRow(conMat, i, &nCols, &cols, &vals);
        //check if this row has any strictly local columns

        /*we know that our global index lists are stored sequentially, so we
          don't need to start every index search at zero, we can start at the
          last entry found in this row */
        label kLast = -1;
        for (label j = 0; j < nCols; j++)
        {
            label matCol = cols[j];
            if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
            {
                //Info<<"j: "<<j<<vals[j]<<endl;
                kLast = this->find_index(matCol, kLast + 1, nUniqueCols, globalIndexList);

                if (kLast >= 0)
                {
                    //k was found
                    label localVal = localColumnStat[kLast];
                    if (localVal == 2)
                    {
                        localRowList[idx] = 1;
                        break;
                    }
                    else
                    {
                        localRowList[idx] = 0;
                    }
                }
                else
                {
                    localRowList[idx] = 0;
                }
            }
        }
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }
    VecRestoreArray(globalVec, &globalVecArray);

    /* Create the scatter context for the remainder of the function */
    Vec colorsLocal;
    // create a scatter context for these colors
    VecCreateSeq(PETSC_COMM_SELF, nUniqueCols, &colorsLocal);
    VecSet(colorsLocal, -1);

    // now create the Index sets
    ISCreateGeneral(PETSC_COMM_WORLD, nUniqueCols, globalIndexList, PETSC_COPY_VALUES, &globalIS);
    // Create the scatter
    VecScatterCreate(colors, globalIS, colorsLocal, NULL, &colorScatter);

    /* Create the conflict resolution scheme*/
    // create tiebreakers locally
    Vec globalTiebreaker;
    VecDuplicate(globalVec, &globalTiebreaker);
    for (label i = colorStart; i < colorEnd; i++)
    {
        srand(i);
        PetscScalar val = rand() % nColG;
        VecSetValue(globalTiebreaker, i, val, INSERT_VALUES);
    }

    // and scatter the random values
    Vec localTiebreaker;
    VecDuplicate(colorsLocal, &localTiebreaker);
    VecScatterBegin(colorScatter, globalTiebreaker, localTiebreaker, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(colorScatter, globalTiebreaker, localTiebreaker, INSERT_VALUES, SCATTER_FORWARD);

    //initialize conflict columns
    label* conflictCols = new label[maxCols];
    label* conflictLocalColIdx = new label[maxCols];
    for (label j = 0; j < maxCols; j++)
    {
        conflictCols[j] = -1;
        conflictLocalColIdx[j] = -1;
    }

    // Create a global distrbuted vector of the only the local portion of
    // localcolumnsstatus
    Vec globalColumnStat;
    VecDuplicate(colors, &globalColumnStat);
    VecSet(globalColumnStat, 0.0);
    for (label k = 0; k < nUniqueCols; k++)
    {
        label localCol = globalIndexList[k];
        PetscScalar localVal = localColumnStat[k];
        if (localCol >= colorStart && localCol < colorEnd)
        {
            VecSetValue(globalColumnStat, localCol, localVal, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(globalColumnStat);
    VecAssemblyEnd(globalColumnStat);

    /*
      create a duplicate matrix for conMat that contains its index into the
      local arrays
    */
    Mat conIndMat;
    MatDuplicate(conMat, MAT_SHARE_NONZERO_PATTERN, &conIndMat);

    // now loop over conMat locally, find the index in the local array
    // and store that value in conIndMat
    for (label i = Istart; i < Iend; i++)
    {
        MatGetRow(conMat, i, &nCols, &cols, &vals);
        label kLast = -1;
        for (label j = 0; j < nCols; j++)
        {
            label matCol = cols[j];
            if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
            {
                kLast = this->find_index(matCol, kLast + 1, nUniqueCols, globalIndexList);
                PetscScalar val = kLast;
                MatSetValue(conIndMat, i, matCol, val, INSERT_VALUES);
            }
        }
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }
    MatAssemblyBegin(conIndMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(conIndMat, MAT_FINAL_ASSEMBLY);

    /*Now color the locally independent columns using the only the rows
      that contain entries from those columns.*/

    // Retrieve the local portion of the color vector
    PetscScalar *colColor, *tbkrLocal, *tbkrGlobal, *globalStat;
    VecGetArray(colors, &colColor);
    VecGetArray(localTiebreaker, &tbkrLocal);
    VecGetArray(globalTiebreaker, &tbkrGlobal);
    VecGetArray(globalColumnStat, &globalStat);

    // Loop over the maximum number of colors
    for (label n = 0; n < maxColors; n++)
    {
        if (n % 100 == 0)
        {
            Info << "ColorSweep: " << n << "   " << mesh_.time().elapsedClockTime() << " s" << endl;
        }

        /* Set all entries for strictly local columns that are currently -1
           to the current color */
        for (label k = 0; k < nUniqueCols; k++)
        {
            label localCol = globalIndexList[k];
            label localVal = localColumnStat[k];
            if (localCol >= colorStart && localCol < colorEnd && localVal == 2)
            {
                // this is a strictly local column;
                label idx = localCol - colorStart;
                if (daUtil_.isValueCloseToRef(colColor[idx], -1.0))
                {
                    colColor[idx] = n;
                }
            }
        }

        //Now loop over the rows and resolve conflicts
        for (label i = Istart; i < Iend; i++)
        {

            // Get the row local row index
            label idx = i - Istart;

            //create the variables for later sorting
            label smallest = nColG;
            label idxKeep = -1;

            //First check if this is a row that contains strictly local columns
            if (localRowList[idx] > 0)
            {

                /* this is a row that contains strictly local columns,get the
                   row information */
                MatGetRow(conMat, i, &nCols, &cols, &vals);

                // set any columns with the current color into conflictCols
                for (label j = 0; j < nCols; j++)
                {
                    if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
                    {
                        label colIdx = cols[j];

                        // Check that this is a local column
                        if (colIdx >= colorStart && colIdx < colorEnd)
                        {
                            //now check if it is a strictly local column
                            label localVal = globalStat[colIdx - colorStart];
                            if (localVal == 2)
                            {
                                // check if the color in this column is from the
                                // current set
                                if (daUtil_.isValueCloseToRef(colColor[colIdx - colorStart], n * 1.0))
                                {
                                    /* This is a potentially conflicting column
                                       store it */
                                    conflictCols[j] = colIdx;

                                    // now check whether this is the one we keep
                                    label tbkr = tbkrGlobal[colIdx - colorStart];
                                    if (tbkr < smallest)
                                    {
                                        smallest = tbkr;
                                        idxKeep = colIdx;
                                    }
                                }
                            }
                        }
                    }
                }

                // Now reset all columns but the one that wins the tiebreak
                for (label j = 0; j < nCols; j++)
                {

                    //check if this is a conflicting column
                    label colIdx = conflictCols[j];
                    if (colIdx >= 0)
                    {
                        // Check that this is also a local column
                        if (colIdx >= colorStart && colIdx < colorEnd)
                        {
                            // and now if it is a strictly local column
                            label localVal = globalStat[colIdx - colorStart];
                            if (localVal == 2)
                            {
                                // now reset the column
                                if (colIdx >= 0 && (colIdx != idxKeep))
                                {
                                    colColor[colIdx - colorStart] = -1;
                                }
                            }
                        }
                    }
                }
                // reset the changed values in conflictCols
                for (label j = 0; j < nCols; j++)
                {
                    if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
                    {
                        //reset all values related to this row in conflictCols
                        conflictCols[j] = -1;
                    }
                }
                MatRestoreRow(conMat, i, &nCols, &cols, &vals);
            }
        }

        // now we want to check the coloring on the strictly local columns
        notColored = 0;

        //loop over the columns and check if there are any uncolored rows
        label colorCounter = 0;
        for (label k = 0; k < nUniqueCols; k++)
        {
            // get the column info
            label localVal = localColumnStat[k];
            label localCol = globalIndexList[k];
            // check if it is strictly local, if so it should be colored
            if (localVal == 2)
            {
                // confirm that it is a local column (is this redundant?
                if (localCol >= colorStart && localCol < colorEnd)
                {
                    label idx = localCol - colorStart;
                    label color = colColor[idx];
                    // now check that it has been colored
                    if (not(color >= 0))
                    {
                        // this column is not colored and coloring is not complete
                        notColored = 1;
                        colorCounter++;
                        //break;
                    }
                }
            }
        }

        // reduce the logical so that we know that all of the processors are
        // ok
        reduce(notColored, sumOp<label>());
        reduce(colorCounter, sumOp<label>());

        if (n % 100 == 0)
        {
            Info << "number of uncolored: " << colorCounter << " " << notColored << endl;
        }

        if (notColored == 0)
        {
            Info << "ColorSweep: " << n << "   " << mesh_.time().elapsedClockTime() << " s" << endl;
            Info << "number of uncolored: " << colorCounter << " " << notColored << endl;
            break;
        }
    }
    VecRestoreArray(colors, &colColor);
    /***** end of local coloring ******/

    // now redo the local row list to handle the global columns
    // create a list of the rows that have any of the global columns included

    for (label i = Istart; i < Iend; i++)
    {
        label idx = i - Istart;
        // get the row information
        MatGetRow(conMat, i, &nCols, &cols, &vals);

        //check if this row has any non-local columns

        /* We know that our localColumnStat is stored sequentially, so we don't
           need to start every index search at zero, we can start at the last
           entry found in this row.*/
        label kLast = -1;

        for (label j = 0; j < nCols; j++)
        {
            // get the column of interest
            label matCol = cols[j];
            // confirm that it has an entry
            if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
            {
                // find the index into the local arrays for this column
                kLast = this->find_index(matCol, kLast + 1, nUniqueCols, globalIndexList);
                // if this column is present (should always be true?) process the row
                if (kLast >= 0)
                {
                    // get the local column type
                    label localVal = localColumnStat[kLast];
                    /* If this is a global column, add the row and move to the next
                       row, otherwise check the next column in the row */
                    if (localVal == 1)
                    {
                        localRowList[idx] = 1;
                        break;
                    }
                    else
                    {
                        localRowList[idx] = 0;
                    }
                }
                else
                {
                    // if this column wasn't found, move to the next column
                    localRowList[idx] = 0;
                }
            }
        }
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }

    // now that we know the set of global rows, complete the coloring

    // Loop over the maximum number of colors
    for (label n = 0; n < maxColors; n++)
    {
        if (n % 100 == 0)
        {
            Info << "Global ColorSweep: " << n << "   " << mesh_.time().elapsedClockTime() << " s" << endl;
        }

        // Retrieve the local portion of the color vector
        // and set all entries that are currently -1 to the current color
        PetscScalar* colColor;
        VecGetArray(colors, &colColor);

        for (label i = colorStart; i < colorEnd; i++)
        {
            label idx = i - colorStart;
            if (daUtil_.isValueCloseToRef(colColor[idx], -1.0))
            {
                colColor[idx] = n;
            }
        }
        VecRestoreArray(colors, &colColor);

        /* We will do the confilct resolution in two passes. On the first pass
           we will keep the value with the smallest random index on the local
           column set. We will not touch the off processor columns. On the second
           pass we will keep the value with the smallest random index, regardless
           of processor. This is to prevent deadlocks in the conflict resolution.*/
        for (label conPass = 0; conPass < 2; conPass++)
        {

            // Scatter the global colors to each processor
            VecScatterBegin(colorScatter, colors, colorsLocal, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(colorScatter, colors, colorsLocal, INSERT_VALUES, SCATTER_FORWARD);

            // compute the number of local rows
            //int nRows = Iend_-Istart_;

            //Allocate a Scalar array to recieve colColors.
            PetscScalar* colColorLocal;
            VecGetArray(colorsLocal, &colColorLocal);
            VecGetArray(colors, &colColor);

            // set the iteration limits based on conPass
            label start, end;
            if (conPass == 0)
            {
                start = colorStart;
                end = colorEnd;
            }
            else
            {
                start = 0;
                end = nColG;
            }

            //Now loop over the rows and resolve conflicts
            for (label i = Istart; i < Iend; i++)
            {
                label idx = i - Istart;
                if (localRowList[idx] == 1) //this row includes at least 1 global col.
                {
                    /* Get the connectivity row as well as its index into the
                       local indices. */
                    MatGetRow(conMat, i, &nCols, &cols, &vals);
                    MatGetRow(conIndMat, i, &nCols2, NULL, &vals2);

                    // initialize the sorting variables
                    label smallest = nColG;
                    label idxKeep = -1;

                    //int localColIdx;
                    // set any columns with the current color into conflictCols
                    for (label j = 0; j < nCols; j++)
                    {
                        if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
                        {
                            label colIdx = cols[j];
                            label localColIdx = round(vals2[j]);

                            // check if the color in this column is from the
                            // current set
                            if (daUtil_.isValueCloseToRef(colColorLocal[localColIdx], n * 1.0))
                            {
                                /* This matches the current color, so set as a
                                   potential conflict */
                                conflictCols[j] = colIdx;
                                conflictLocalColIdx[j] = localColIdx;
                                /* this is one of the conflicting columns.
                                   If this is a strictly local column, keep it.
                                   Otherwise, compare its random number to the
                                   current smallest one, keep the smaller one and
                                   its index find the index of the smallest
                                   tiebreaker. On the first pass this is only
                                   for the the local columns. On pass two it is
                                   for all columns.*/
                                if (localColumnStat[localColIdx] == 2)
                                {
                                    smallest = -1;
                                    idxKeep = colIdx;
                                }
                                else if (tbkrLocal[localColIdx] < smallest and colIdx >= start and colIdx < end)
                                {
                                    smallest = tbkrLocal[localColIdx];
                                    idxKeep = cols[j];
                                }
                            }
                        }
                    }

                    // Now reset all the conflicting rows
                    for (label j = 0; j < nCols; j++)
                    {
                        label colIdx = conflictCols[j];
                        label localColIdx = conflictLocalColIdx[j];
                        // check that the column is in the range for this conPass.
                        if (colIdx >= start && colIdx < end)
                        {
                            /*this column is local. If this isn't the
                            smallest, reset it.*/
                            if (colIdx != idxKeep)
                            {
                                if (localColIdx >= 0)
                                {
                                    if (localColumnStat[localColIdx] == 2)
                                    {
                                        Pout << "local Array Index: " << colIdx << endl;
                                        Info << "Error, setting a local column!" << endl;
                                        return;
                                    }
                                    PetscScalar valIn = -1;
                                    VecSetValue(colors, colIdx, valIn, INSERT_VALUES);
                                    colColorLocal[localColIdx] = -1;
                                }
                            }
                        }
                    }

                    /* reset any columns that have been changed in conflictCols
                    and conflictLocalColIdx */
                    for (label j = 0; j < nCols; j++)
                    {
                        if (!daUtil_.isValueCloseToRef(vals[j], 0.0))
                        {
                            //reset all values related to this row in conflictCols
                            conflictCols[j] = -1;
                            conflictLocalColIdx[j] = -1;
                        }
                    }

                    // Restore the row information
                    MatRestoreRow(conIndMat, i, &nCols2, NULL, &vals2);
                    MatRestoreRow(conMat, i, &nCols, &cols, &vals);
                }
            }

            VecRestoreArray(colors, &colColor);
            VecRestoreArray(colorsLocal, &colColorLocal);
            VecAssemblyBegin(colors);
            VecAssemblyEnd(colors);
        }

        //check the coloring for completeness
        label colorCounter = 0;
        this->coloringComplete(colors, colorCounter, notColored);
        if (n % 100 == 0)
        {
            Info << "Number of Uncolored: " << colorCounter << " " << notColored << endl;
        }
        if (notColored == 0)
        {
            Info << "Global ColorSweep: " << n << "   " << mesh_.time().elapsedClockTime() << " s" << endl;
            Info << "Number of Uncolored: " << colorCounter << " " << notColored << endl;
            break;
        }
    }
    VecRestoreArray(globalTiebreaker, &tbkrGlobal);
    VecRestoreArray(globalColumnStat, &globalStat);

    // count the current colors and aggregate
    currColor = 0;
    PetscScalar color;
    for (label i = colorStart; i < colorEnd; i++)
    {
        VecGetValues(colors, 1, &i, &color);
        //Pout<<"Color: "<<i<<" "<<color<<endl;
        if (color > currColor)
        {
            currColor = color;
        }
    }

    reduce(currColor, maxOp<label>());

    nColors = currColor + 1;

    Info << "Ncolors: " << nColors << endl;

    //check the initial coloring for completeness
    //this->coloringComplete(colors, colorCounter, notColored);

    // clean up the unused memory
    VecRestoreArray(localTiebreaker, &tbkrLocal);
    delete[] conflictCols;
    delete[] conflictLocalColIdx;
    delete[] globalIndexList;
    delete[] localColumnStat;
    ISDestroy(&globalIS);
    VecScatterDestroy(&colorScatter);
    VecDestroy(&colorsLocal);
    VecDestroy(&globalTiebreaker);
    VecDestroy(&globalColumnStat);
    VecDestroy(&localTiebreaker);
    delete[] localRowList;
    VecDestroy(&globalVec);
    MatDestroy(&conIndMat);
}

void DAJacCon::getMatNonZeros(
    const Mat conMat,
    label& maxCols,
    scalar& allNonZeros) const
{
    /*
    Get the max nonzeros per row, and all the nonzeros for this matrix
    This will be used in computing coloring

    Input:
    -----
    conMat: the matrix to compute nonzeros

    Output:
    ------
    maxCols: max nonzeros per row among all rows

    allNonZeros: all non zero elements in conMat
    */

    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    label Istart, Iend;

    // set the counter
    maxCols = 0;
    allNonZeros = 0.0;

    // Determine which rows are on the current processor
    MatGetOwnershipRange(conMat, &Istart, &Iend);

    // loop over the matrix and find the largest number of cols
    for (label i = Istart; i < Iend; i++)
    {
        MatGetRow(conMat, i, &nCols, &cols, &vals);
        if (nCols < 0)
        {
            std::cout << "Warning! procI: " << Pstream::myProcNo() << " nCols <0 at rowI: " << i << std::endl;
            std::cout << "Set nCols to zero " << std::endl;
            nCols = 0;
        }
        if (nCols > maxCols) // perhaps actually check vals?
        {
            maxCols = nCols;
        }
        allNonZeros += nCols;
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }

    //reduce the maxcols value so that all procs have the same size
    reduce(maxCols, maxOp<label>());

    reduce(allNonZeros, sumOp<scalar>());

    return;
}

label DAJacCon::find_index(
    const label target,
    const label start,
    const label size,
    const label* valArray) const
{
    /*
    Find the index of a value in an array
    
    Input:
    ------
    target: the target value to find
    start: the start index in valArray
    size: the size of valArray array
    valArray: the array to check the target value

    Output:
    -------
    k: the index of the value in the array, if the value is 
    not found, return -1
    */

    // loop over the valArray from start until target is found
    for (label k = start; k < size; k++)
    {
        if (valArray[k] == target)
        {
            //Info<<"Start: "<<start<<" "<<k<<endl;
            //this is the k of interest
            return k;
        }
    }
    return -1;
}

void DAJacCon::coloringComplete(
    const Vec colors,
    label& colorCounter,
    label& notColored) const
{
    /*
    Check if the coloring process is finished and return
    the number of uncolored columns

    Input:
    -----
    colors: the current coloring vector

    Output:
    ------
    notColored: the number of uncolored columns

    */

    PetscScalar color;
    PetscInt colorStart, colorEnd;

    notColored = 0;
    // get the range of colors owned by the local prock
    VecGetOwnershipRange(colors, &colorStart, &colorEnd);
    //loop over the columns and check if there are any uncolored rows
    const PetscScalar* colColor;
    VecGetArrayRead(colors, &colColor);
    colorCounter = 0;
    for (label i = colorStart; i < colorEnd; i++)
    {
        color = colColor[i - colorStart];
        //VecGetValues(colors,1,&i,&color);
        if (not(color >= 0))
        {
            // this columns not colored and coloring is not complete
            //Pout<<"coloring incomplete...: "<<color<<" "<<i<<endl;
            //VecView(colors,PETSC_VIEWER_STDOUT_WORLD);
            notColored = 1;
            colorCounter++;
            //break;
        }
    }
    VecRestoreArrayRead(colors, &colColor);
    // reduce the logical so that we know that all of the processors are
    // ok
    //Pout<<"local number of uncolored: "<<colorCounter<<" "<<notColored<<endl;
    reduce(notColored, sumOp<label>());
    reduce(colorCounter, sumOp<label>());
}

void DAJacCon::validateColoring(
    Mat conMat,
    Vec colors) const
{
    /*
    loop over the rows and verify that no row has two columns with the same color

    Input:
    conMat: connectivity mat for check coloring

    colors: the coloring vector

    Example:
    If the conMat reads, its coloring for each column can be

           color0  color1
             |     |
             1  0  0  0
    conMat = 0  1  1  0
             0  0  1  0
             0  0  0  1
                |     | 
            color0   color0

    Then, if colors = {0, 0, 1, 0}-> no coloring conflict
    if colors = {0, 1, 0, 0}-> coloring conclict
    */

    Info << "Validating Coloring..." << endl;

    PetscInt nCols;
    const PetscInt* cols;
    const PetscScalar* vals;

    label Istart, Iend;

    // scatter colors to local array for all procs
    Vec vout;
    VecScatter ctx;
    VecScatterCreateToAll(colors, &ctx, &vout);
    VecScatterBegin(ctx, colors, vout, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, colors, vout, INSERT_VALUES, SCATTER_FORWARD);

    PetscScalar* colorsArray;
    VecGetArray(vout, &colorsArray);

    // Determine which rows are on the current processor
    MatGetOwnershipRange(conMat, &Istart, &Iend);

    // first calc the largest nCols in conMat
    label colMax = 0;
    for (label i = Istart; i < Iend; i++)
    {
        MatGetRow(conMat, i, &nCols, &cols, &vals);
        if (nCols > colMax)
        {
            colMax = nCols;
        }
        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }

    // now check if conMat has conflicting rows
    labelList rowColors(colMax);
    for (label i = Istart; i < Iend; i++)
    {
        MatGetRow(conMat, i, &nCols, &cols, &vals);

        // initialize rowColors with -1
        for (label nn = 0; nn < colMax; nn++)
        {
            rowColors[nn] = -1;
        }

        // set rowColors for this row
        for (label j = 0; j < nCols; j++)
        {
            if (daUtil_.isValueCloseToRef(vals[j], 1.0))
            {
                rowColors[j] = round(colorsArray[cols[j]]);
            }
        }

        // check if rowColors has duplicated colors
        for (label nn = 0; nn < nCols; nn++)
        {
            for (label mm = nn + 1; mm < nCols; mm++)
            {
                if (rowColors[nn] != -1 && rowColors[nn] == rowColors[mm])
                {
                    FatalErrorIn("Conflicting Colors Found!")
                        << " row: " << i << " col1: " << cols[nn] << " col2: " << cols[mm]
                        << " color: " << rowColors[nn] << abort(FatalError);
                }
            }
        }

        MatRestoreRow(conMat, i, &nCols, &cols, &vals);
    }

    VecRestoreArray(vout, &colorsArray);
    VecScatterDestroy(&ctx);
    VecDestroy(&vout);

    Info << "No Conflicting Colors Found!" << endl;

    return;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
