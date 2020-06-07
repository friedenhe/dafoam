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
    const DAModel& daModel)
    : DAJacCon(modelType, mesh, daOption, daModel)
{
    this->initializePetscVecs();
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

    // nLocalObjFuncGeoElements: the number of objFunc discrete elements for local procs
    label nLocalObjFuncGeoElements = objFuncFaceSources.size() + objFuncCellSources.size();

    DAUtility daUtil;
    globalObjFuncGeoNumbering_ = daUtil.genGlobalIndex(nLocalObjFuncGeoElements);

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

    // initilialize jacConColoredColumns_
    VecCreate(PETSC_COMM_WORLD, &jacConColoredColumns_);
    VecSetSizes(
        jacConColoredColumns_,
        nLocalObjFuncGeoElements,
        PETSC_DECIDE);
    VecSetFromOptions(jacConColoredColumns_);
    VecZeroEntries(jacConColoredColumns_);

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

    label addFace = 0;
    forAll(objFuncConInfo, idxI)
    {
        forAll(objFuncConInfo[idxI], idxJ)
        {
            const word stateName = objFuncConInfo[idxI][idxJ];
            word stateType = daIndex_.adjStateType[stateName];
            if (stateType == "surfaceScalarState")
            {
                addFace = 1;
            }
        }
    }

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
            wordList connectedStatesLocal = objFuncConInfo[idxJ];

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
                    connectedStatesInterProc[k] = objFuncConInfo[k + idxJ];
                }
            }
            else
            {
                connectedStatesInterProc.setSize(1);
                connectedStatesInterProc[0] = objFuncConInfo[maxConLevel];
            }

            this->addStateConnections(
                connectedStatesP,
                idxN,
                idxJ,
                connectedStatesLocal,
                connectedStatesInterProc,
                addFace);
        }

        label glbRowI = globalObjFuncGeoNumbering_.toGlobal(idxI);

        this->setupJacobianConnections(jacCon_, connectedStatesP, glbRowI);
    }

    MatAssemblyBegin(jacCon_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jacCon_, MAT_FINAL_ASSEMBLY);
}

label DAJacCondFdW::getNJacConColors() const
{
    /*
    Return the number of colors

    Output:
    ------
    nJacConColors_: the number of colors depends on 
    whether the coloring is used
    */

    if (daOption_.getOption<label>("adjUseColoring"))
    {
        return nJacConColors_;
    }
    else
    {
        return daIndex_.nGlobalAdjointStates;
    }

    return -1;
}

} // End namespace Foam

// ************************************************************************* //
