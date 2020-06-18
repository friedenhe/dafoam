/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAJacCondFdW.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAJacCondFdW, 0);
addToRunTimeSelectionTable(DAJacCon, DAJacCondFdW, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAJacCondFdW::DAJacCondFdW(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAJacCon(modelType, mesh, daOption, daModel, daIndex)
{
    this->initializePetscVecs();
    this->initializeStateBoundaryCon();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAJacCondFdW::initializePetscVecs()
{

    // dFdW Colors the jacConColoredColumn will be initialized in
    // DAJacCondFdW::initializeJacCon
    VecCreate(PETSC_COMM_WORLD, &jacConColors_);
    VecSetSizes(jacConColors_, daIndex_.nLocalAdjointStates, PETSC_DECIDE);
    VecSetFromOptions(jacConColors_);

    return;
}

void DAJacCondFdW::initializeJacCon(const dictionary& options)
{
    /*
    Initialize the connectivity matrix

    Input:
    ------
    options.objFuncFaceSources: a labelList that contains all the
    face indices for the objective

    options.objFuncCellSources: a labelList that contains all the
    cell indices for the objective

    Output:
    ------
    jacCon_: connectivity matrix for dFdW, here dFdWCon has a
    size of nLocalObjFuncGeoElements * nGlobalAdjointStates
    The reason that dFdWCon has nLocalObjFuncGeoElements rows is 
    because we need to divide the objective function into 
    nLocalObjFuncGeoElements discrete value such that we can
    use coloring to compute dFdW
    */

    labelList objFuncFaceSources;
    labelList objFuncCellSources;
    options.readEntry<labelList>("objFuncFaceSources", objFuncFaceSources);
    options.readEntry<labelList>("objFuncCellSources", objFuncCellSources);

    objFuncFaceSize_ = objFuncFaceSources.size();
    objFuncCellSize_ = objFuncCellSources.size();

    // nLocalObjFuncGeoElements: the number of objFunc discrete elements for local procs
    label nLocalObjFuncGeoElements = objFuncFaceSize_ + objFuncCellSize_;

    globalObjFuncGeoNumbering_ = DAUtility::genGlobalIndex(nLocalObjFuncGeoElements);

    MatCreate(PETSC_COMM_WORLD, &jacCon_);
    MatSetSizes(
        jacCon_,
        nLocalObjFuncGeoElements,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(jacCon_);
    MatMPIAIJSetPreallocation(jacCon_, 200, NULL, 200, NULL);
    MatSeqAIJSetPreallocation(jacCon_, 200, NULL);
    //MatSetOption(jacCon_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(jacCon_);
    MatZeroEntries(jacCon_);

    Info << "dFdWCon Created!" << endl;
}

void DAJacCondFdW::setupJacCon(const dictionary& options)
{
    /*
    Setup jacCon_ matrix
    
    Input:
    -----
    options.objFuncConInfo: the connectivity information for
    this objective function

    options.objFuncFaceSources: a labelList that contains all the
    face indices for the objective

    options.objFuncCellSources: a labelList that contains all the
    cell indices for the objective
    */

    Info << "Setting up dFdWCon.." << endl;

    MatZeroEntries(jacCon_);

    Mat connectedStatesP;

    List<List<word>> objFuncConInfo;
    labelList objFuncFaceSources;
    labelList objFuncCellSources;
    options.readEntry<List<List<word>>>("objFuncConInfo", objFuncConInfo);
    options.readEntry<labelList>("objFuncFaceSources", objFuncFaceSources);
    options.readEntry<labelList>("objFuncCellSources", objFuncCellSources);

    // maximal connectivity level information
    label maxConLevel = objFuncConInfo.size() - 1;

    forAll(objFuncFaceSources, idxI)
    {

        const label& objFuncFaceI = objFuncFaceSources[idxI];
        label bFaceI = objFuncFaceI - daIndex_.nLocalInternalFaces;
        const label patchI = daIndex_.bFacePatchI[bFaceI];
        const label faceI = daIndex_.bFaceFaceI[bFaceI];

        // create a shorter handle for the boundary patch
        const fvPatch& patch = mesh_.boundary()[patchI];
        // get the cells associated with this boundary patch
        const UList<label>& pFaceCells = patch.faceCells();

        // Now get the cell that borders this face
        label idxN = pFaceCells[faceI];

        //zero the connections
        this->createConnectionMat(&connectedStatesP);

        forAll(objFuncConInfo, idxJ) // idxJ: con level
        {
            // set connectedStatesLocal: the locally connected state variables for this level
            wordList connectedStatesLocal(0);
            forAll(objFuncConInfo[idxJ], idxK)
            {
                word conName = objFuncConInfo[idxJ][idxK];
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
                    label conSize = objFuncConInfo[k + idxJ].size();
                    for (label l = 0; l < conSize; l++)
                    {
                        word conName = objFuncConInfo[k + idxJ][l];
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
                label conSize = objFuncConInfo[maxConLevel].size();
                for (label l = 0; l < conSize; l++)
                {
                    word conName = objFuncConInfo[maxConLevel][l];
                    // Exclude surfaceScalarState when appending connectedStatesLocal
                    // whether to add it depends on addFace parameter
                    if (daIndex_.adjStateType[conName] != "surfaceScalarState")
                    {
                        connectedStatesInterProc[0].append(conName);
                    }
                }
            }

            // check if we need to add face
            label addFace = 0;
            forAll(stateInfo_["surfaceScalarStates"], idxK)
            {
                word conName = stateInfo_["surfaceScalarStates"][idxK];
                if (DAUtility::isInList<word>(conName, objFuncConInfo[idxJ]))
                {
                    addFace = 1;
                }
            }

            this->addStateConnections(
                connectedStatesP,
                idxN,
                idxJ,
                connectedStatesLocal,
                connectedStatesInterProc,
                addFace);
        }

        label glbRowI = this->getGlobalObjFuncGeoIndex("face", idxI);

        this->setupJacobianConnections(jacCon_, connectedStatesP, glbRowI);
    }

    forAll(objFuncCellSources, idxI)
    {
        // TODO: need to implemnt this!
        FatalErrorIn("") << "not implemented!"
                         << abort(FatalError);
    }

    MatAssemblyBegin(jacCon_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacCon_, MAT_FINAL_ASSEMBLY);
    
}

label DAJacCondFdW::getLocalObjFuncGeoIndex(
    const word idxType,
    const label idxI) const
{
    label localIdx = -9999;
    if (idxType == "face")
    {
        localIdx = idxI;
    }
    else if (idxType == "cell")
    {
        localIdx = objFuncFaceSize_ + idxI;
    }
    else
    {
        FatalErrorIn("") << "idxType: " << idxType << "not supported!"
                         << abort(FatalError);
    }
    return localIdx;
}

label DAJacCondFdW::getGlobalObjFuncGeoIndex(
    const word idxType,
    const label idxI) const
{
    label localIdx = this->getLocalObjFuncGeoIndex(idxType, idxI);

    return globalObjFuncGeoNumbering_.toGlobal(localIdx);
}

void DAJacCondFdW::setObjFuncVec(
    scalarList objFuncFaceValues,
    scalarList objFuncCellValues,
    Vec objFuncVec) const
{
    PetscScalar* objFuncVecArray;
    VecGetArray(objFuncVec, &objFuncVecArray);

    forAll(objFuncFaceValues, idxI)
    {
        scalar val = objFuncFaceValues[idxI];
        label localIdx = getLocalObjFuncGeoIndex("face", idxI);
        objFuncVecArray[localIdx] = val;
    }

    forAll(objFuncCellValues, idxI)
    {
        scalar val = objFuncCellValues[idxI];
        label localIdx = getLocalObjFuncGeoIndex("cell", idxI);
        objFuncVecArray[localIdx] = val;
    }

    VecRestoreArray(objFuncVec, &objFuncVecArray);
}

} // End namespace Foam

// ************************************************************************* //
