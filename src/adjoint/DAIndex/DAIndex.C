/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAIndex::DAIndex(const fvMesh& mesh)
    : regIOobject(
        IOobject(
            "DAIndex", // always use DAIndex for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh),
      daOption_(mesh.thisDb().lookupObject<DAOption>("DAOption")),
      daRegState_(const_cast<DARegState&>(
          mesh.thisDb().lookupObject<DARegState>("DARegState"))),
      regStates_(const_cast<HashTable<wordList>&>(
          daRegState_.getRegStates())),
      pointProcAddressing(
          IOobject(
              "pointProcAddressing",
              mesh_.facesInstance(),
              mesh_.meshSubDir,
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE),
          {})
{
    // Calculate the sizes
    // local denotes the current MPI process
    // global denotes all the MPI processes

    // Setup the adjoint state name and type lists.
    forAll(regStates_["volVectorStates"], idxI)
    {
        adjStateNames.append(regStates_["volVectorStates"][idxI]);
        adjStateType.set(regStates_["volVectorStates"][idxI], "volVectorState");
    }
    forAll(regStates_["volScalarStates"], idxI)
    {
        adjStateNames.append(regStates_["volScalarStates"][idxI]);
        adjStateType.set(regStates_["volScalarStates"][idxI], "volScalarState");
    }
    forAll(regStates_["modelStates"], idxI)
    {
        adjStateNames.append(regStates_["modelStates"][idxI]);
        adjStateType.set(regStates_["modelStates"][idxI], "modelState");
    }
    forAll(regStates_["surfaceScalarStates"], idxI)
    {
        adjStateNames.append(regStates_["surfaceScalarStates"][idxI]);
        adjStateType.set(regStates_["surfaceScalarStates"][idxI], "surfaceScalarState");
    }

    //Info<<"adjStateNames"<<adjStateNames<<endl;
    Info << "Adjoint States: " << adjStateType << endl;

    // Local mesh related sizes
    nLocalCells = mesh.nCells();
    nLocalFaces = mesh.nFaces();
    nLocalPoints = mesh.nPoints();
    nLocalXv = nLocalPoints * 3;
    nLocalInternalFaces = mesh.nInternalFaces();
    nLocalBoundaryFaces = nLocalFaces - nLocalInternalFaces;
    nLocalBoundaryPatches = mesh.boundaryMesh().size();

    // get bFacePatchI and bFaceFaceI
    // these two lists store the patchI and faceI for a given boundary mesh face index
    // the index of these lists starts from the first boundary face of the first boundary patch.
    // they will be used to quickly get the patchI and faceI from a given boundary face index
    bFacePatchI.setSize(nLocalBoundaryFaces);
    bFaceFaceI.setSize(nLocalBoundaryFaces);
    label tmpCounter = 0;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        forAll(mesh_.boundaryMesh()[patchI], faceI)
        {
            bFacePatchI[tmpCounter] = patchI;
            bFaceFaceI[tmpCounter] = faceI;
            tmpCounter++;
        }
    }

    // Initialize state local index offset, it will be used in getLocalStateIndex function
    this->calcStateLocalIndexOffset(stateLocalIndexOffset);

    // Initialize adjStateID. It stores the stateID for a given stateName
    this->calcAdjStateID(adjStateID);

    // Local adjoint state sizes
    // first get how many state variables are registered.
    // Note turbStates are treated separatedly
    nVolScalarStates = regStates_["volScalarStates"].size();
    nVolVectorStates = regStates_["volVectorStates"].size();
    nSurfaceScalarStates = regStates_["surfaceScalarStates"].size();
    nModelStates = regStates_["modelStates"].size();

    // we can now calculate adjoint state size
    label nLocalCellStates = (nVolVectorStates * 3 + nVolScalarStates + nModelStates) * nLocalCells;
    label nLocalFaceStates = nSurfaceScalarStates * nLocalFaces;
    nLocalAdjointStates = nLocalCellStates + nLocalFaceStates;

    // Setup the global numbering to convert a local index to the associated global index
    globalAdjointStateNumbering = this->genGlobalIndex(nLocalAdjointStates);
    globalCellNumbering = this->genGlobalIndex(nLocalCells);
    globalCellVectorNumbering = this->genGlobalIndex(nLocalCells * 3);
    globalFaceNumbering = this->genGlobalIndex(nLocalFaces);
    globalXvNumbering = this->genGlobalIndex(nLocalXv);

    // global Adjoint state sizes
    nGlobalAdjointStates = globalAdjointStateNumbering.size();
    nGlobalCells = globalCellNumbering.size();
    nGlobalFaces = globalFaceNumbering.size();
    nGlobalXv = globalXvNumbering.size();

    // now compute nUndecomposedPoints based on pointProcAddressing
    // this will be used in generating the total sensitivity for the undecomposed domain
    // for sensitivity map plot
    if (!Pstream::parRun())
    {
        // for serial cases, pointProcAddressing is an empty list, so manually assign it
        for (label i = 0; i < nLocalPoints; i++)
        {
            pointProcAddressing.append(i);
        }

        nUndecomposedPoints = nLocalPoints;
    }
    else
    {
        // for parallel cases, we can read pointProcAddressing
        // get the global point size
        label pointMaxIdx = max(pointProcAddressing);
        reduce(pointMaxIdx, maxOp<label>());
        // +1 since procAddressing point index starts with 0
        nUndecomposedPoints = pointMaxIdx + 1;
    }

    // Print relevant sizes to screen
    Info << "Global Cells: " << nGlobalCells << endl;
    Info << "Global Faces: " << nGlobalFaces << endl;
    Info << "Global Xv: " << nGlobalXv << endl;
    Info << "Undecomposed points: " << nUndecomposedPoints << endl;
    Info << "Global Adjoint States: " << nGlobalAdjointStates << endl;

    // initialize stuff for coloring
    label useColoring = daOption_.getOption<label>("adjUseColoring");
    if (useColoring)
    {
        // calculate nLocalCoupledBFaces and isCoupledFace
        isCoupledFace.setSize(nLocalFaces);
        for (label i = 0; i < nLocalFaces; i++)
        {
            isCoupledFace[i] = 0;
        }

        nLocalCoupledBFaces = 0;
        label faceIdx = nLocalInternalFaces;
        forAll(mesh_.boundaryMesh(), patchI)
        {
            forAll(mesh_.boundaryMesh()[patchI], faceI)
            {
                if (mesh_.boundaryMesh()[patchI].coupled())
                {
                    // this is a coupled patch
                    isCoupledFace[faceIdx] = 1;
                    nLocalCoupledBFaces++;
                }
                faceIdx++;
            }
        }

        globalCoupledBFaceNumbering = this->genGlobalIndex(nLocalCoupledBFaces);
        nGlobalCoupledBFaces = globalCoupledBFaceNumbering.size();

        // calculate nLocalCyclicAMIFaces and isCyclicAMIFace
        isCyclicAMIFace.setSize(nLocalFaces);
        for (label i = 0; i < nLocalFaces; i++)
        {
            isCyclicAMIFace[i] = 0;
        }

        nLocalCyclicAMIFaces = 0;
        faceIdx = nLocalInternalFaces;
        forAll(mesh_.boundaryMesh(), patchI)
        {
            forAll(mesh_.boundaryMesh()[patchI], faceI)
            {
                if (mesh_.boundaryMesh()[patchI].type() == "cyclicAMI")
                {
                    // this is a cyclicAMI patch
                    isCyclicAMIFace[faceIdx] = 1;
                    nLocalCyclicAMIFaces++;
                }
                faceIdx++;
            }
        }
    }

    this->calcLocalIdxLists(adjStateName4LocalAdjIdx, cellIFaceI4LocalAdjIdx);
    /*
    // check if we have user-defined patches or volumes, if yes, calculate their face and cell indices
    if (adjIO_.userDefinedPatchInfo.size() != 0)
    {
        this->calcFaceIndx4UserDefinedPatches();
    }
    if (adjIO_.userDefinedVolumeInfo.size() != 0)
    {
        this->calcCellIndx4UserDefinedVolumes();
    }
*/
}

DAIndex::~DAIndex()
{
}

// this is a virtual function for regIOobject
bool DAIndex::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

void DAIndex::calcStateLocalIndexOffset(HashTable<label>& offset)
{
    /*
    Calculate the indexing offset for all states (stateLocalIndexOffset),
    this will be used in the DAIndex::getLocalAdjointStateIndex function

    For state-by-state ordering, we set u_0, v_0, w_0, u_1, v_1, w_1,
    ...., p_0, p_1, ... nuTilda_0, nuTilda_1, ... with subscript being the
    cell index. stateLocalIndexingOffset will return how many states are
    before a specific stateName

    For cell by cell ordering, we set u_0, v_0, w_0, p_0, nuTilda_0, phi_0, .... 
    u_N, v_N, w_N, p_N, nuTilda_N, phi_N with subscript being the cell index. 
    so stateLocalIndexingOffset will return how many states are before a specific 
    stateName for a given cell index

    Output:
    ------
    offset: hash table of  local state variable index offset, This will be used in determing 
    the local indexing for adjoint states. It differs depending on whether we use state-by-state 
    or cell-by-cell ordering
    */

    word adjJacMatOrdering = daOption_.getOption<word>("adjJacMatOrdering");

    if (adjJacMatOrdering == "state")
    {

        forAll(adjStateNames, idxI)
        {
            word stateName = adjStateNames[idxI];

            label counter = 0;

            forAll(regStates_["volVectorStates"], idx)
            {
                if (regStates_["volVectorStates"][idx] == stateName)
                {
                    offset.set(stateName, counter * nLocalCells);
                }
                counter += 3;
            }

            forAll(regStates_["volScalarStates"], idx)
            {
                if (regStates_["volScalarStates"][idx] == stateName)
                {
                    offset.set(stateName, counter * nLocalCells);
                }
                counter++;
            }

            forAll(regStates_["modelStates"], idx)
            {
                if (regStates_["modelStates"][idx] == stateName)
                {
                    offset.set(stateName, counter * nLocalCells);
                }
                counter++;
            }

            forAll(regStates_["surfaceScalarStates"], idx)
            {
                if (regStates_["surfaceScalarStates"][idx] == stateName && idx == 0)
                {
                    offset.set(stateName, counter * nLocalCells);
                }
                if (regStates_["surfaceScalarStates"][idx] == stateName && idx > 0)
                {
                    offset.set(stateName, counter * nLocalFaces);
                }
                counter++;
            }
        }
    }
    else if (adjJacMatOrdering == "cell")
    {

        forAll(adjStateNames, idxI)
        {
            word stateName = adjStateNames[idxI];

            label counter = 0;

            forAll(regStates_["volVectorStates"], idx)
            {
                if (regStates_["volVectorStates"][idx] == stateName)
                {
                    offset.set(stateName, counter);
                }
                counter += 3;
            }

            forAll(regStates_["volScalarStates"], idx)
            {
                if (regStates_["volScalarStates"][idx] == stateName)
                {
                    offset.set(stateName, counter);
                }
                counter++;
            }

            forAll(regStates_["modelStates"], idx)
            {
                if (regStates_["modelStates"][idx] == stateName)
                {
                    offset.set(stateName, counter);
                }
                counter++;
            }

            forAll(regStates_["surfaceScalarStates"], idx)
            {
                if (regStates_["surfaceScalarStates"][idx] == stateName)
                {
                    offset.set(stateName, counter);
                }
                counter++;
            }
        }

        // We also need a few more offsets

        // calculate faceOwner
        faceOwner.setSize(nLocalFaces);
        const UList<label>& internalFaceOwner = mesh_.owner(); // these only include internal faces owned cellI
        forAll(faceOwner, idxI)
        {
            if (idxI < nLocalInternalFaces)
            {
                faceOwner[idxI] = internalFaceOwner[idxI];
            }
            else
            {
                label relIdx = idxI - nLocalInternalFaces;
                label patchIdx = bFacePatchI[relIdx];
                label faceIdx = bFaceFaceI[relIdx];
                const UList<label>& pFaceCells = mesh_.boundaryMesh()[patchIdx].faceCells();
                faceOwner[idxI] = pFaceCells[faceIdx];
            }
        }

        // Calculate the cell owned face index. Note: we can't use mesh.cells here since it will have
        // duplicated face indices
        List<List<label>> cellOwnedFaces;
        cellOwnedFaces.setSize(nLocalCells);
        forAll(faceOwner, idxI)
        {
            label ownedCellI = faceOwner[idxI];
            cellOwnedFaces[ownedCellI].append(idxI);
        }
        //Info<<"cellOwnedFaces "<<cellOwnedFaces<<endl;
        // check if every cells have owned faces
        // This is unnecessary since some cells may own no face.
        //forAll(cellOwnedFaces,idxI)
        //{
        //    if(cellOwnedFaces[idxI].size()==0) FatalErrorIn("")<<"cell "<<idxI<<" owns no faces"<<abort(FatalError);
        //}

        // We first calculate phiAccumulatedOffset
        phiAccumulatdOffset.setSize(nLocalCells);
        forAll(phiAccumulatdOffset, idxI) phiAccumulatdOffset[idxI] = -9999999;
        forAll(phiAccumulatdOffset, idxI)
        {
            if (idxI == 0)
                phiAccumulatdOffset[idxI] = 0;
            else
                phiAccumulatdOffset[idxI] = cellOwnedFaces[idxI - 1].size() + phiAccumulatdOffset[idxI - 1];
        }
        //Info<<"phiAccumulatdOffset "<<phiAccumulatdOffset<<endl;

        // Now calculate the phiLocalOffset
        phiLocalOffset.setSize(nLocalFaces);
        forAll(phiLocalOffset, idxI) phiLocalOffset[idxI] = -9999999;
        forAll(cellOwnedFaces, idxI) // idxI is cell Index
        {
            forAll(cellOwnedFaces[idxI], offsetI)
            {
                label ownedFace = cellOwnedFaces[idxI][offsetI];
                phiLocalOffset[ownedFace] = offsetI;
            }
        }
        //Info<<"phiLocalOffset "<<phiLocalOffset<<endl;
        //Info<<"stateLocalIndexOffset "<<stateLocalIndexOffset<<endl;
    }
    else
    {
        FatalErrorIn("") << "adjJacMatOrdering invalid" << abort(FatalError);
    }

    return;
}

void DAIndex::calcAdjStateID(HashTable<label>& adjStateID)
{
    /* 
    The stateID is an alternative for the stateNames
    stateID starts from 0 for the first volVector state
    e.g., if the state variables are U, p, nut, phi, their 
    state ID are U=0, p=1, nut=2, phi=3

    Output:
    -------
    adjStateID: the state ID list
    */

    label id = 0;
    forAll(regStates_["volVectorStates"], idx)
    {
        word stateName = regStates_["volVectorStates"][idx];
        adjStateID.set(stateName, id);
        id++;
    }

    forAll(regStates_["volScalarStates"], idx)
    {
        word stateName = regStates_["volScalarStates"][idx];
        adjStateID.set(stateName, id);
        id++;
    }

    forAll(regStates_["modelStates"], idx)
    {
        word stateName = regStates_["modelStates"][idx];
        adjStateID.set(stateName, id);
        id++;
    }

    forAll(regStates_["surfaceScalarStates"], idx)
    {
        word stateName = regStates_["surfaceScalarStates"][idx];
        adjStateID.set(stateName, id);
        id++;
    }
    return;
}

globalIndex DAIndex::genGlobalIndex(const label localIndexSize)
{
    /*
    Generate a glocal index system based on the local index size 
    such that we can use it to map a local index to a global one

    Input:
    -----
    localIndexSize: the SIZE of local index

    Output:
    ------
    globalIndex object: the global index object to map a local index
    to a global index

    Example:
    --------
    If the local index reads:

    On processor 0:
    labelList sampleList = {0, 1, 2};
    globalIndex glbSample = genGlobalIndex(sampleList.size());

    On processor 1:
    labelList sampleList = {0, 1};
    globalIndex glbSample = genGlobalIndex(sampleList.size());

    After each processor calls genGlobalIndex and get the glbSample
    object, we can use it to map a local index to a global one,
    e.g., on processor 0, if we call:

    label glxIdx = glbSample.toGlobal(1);

    it will return glbIdx = 1;
    However, on processor 1, if we call

    label glxIdx = glbSample.toGlobal(1);

    it will return glbIdx = 4;

    The date storage structure is illustrated as follows

    global index -> 0 1 2 3 4
    local index  -> 0 1 2 0 1
                    ----- ===
                    proc0 proc1

    */
    globalIndex result(localIndexSize);
    return result;
}

void DAIndex::calcLocalIdxLists(
    wordList& stateName4LocalAdjIdx,
    scalarList& cellIFaceI4LocalIdx)
{
    /*
    Initialize indexing lists:
    cellIFaceI4LocalAdjIdx
    adjStateName4LocalAdjIdx

    Output:
    -------
    cellIFaceI4LocalAdjIdx: stores the cell/face index for a local adjoint index
    For vector fields, the decima of cellIFaceI4LocalIdx denotes the vector component
    e.g., 10.1 means cellI=10, y compoent of U

    adjStateName4LocalAdjIdx: stores the state name for a local adjoint index
    */

    cellIFaceI4LocalIdx.setSize(nLocalAdjointStates);
    stateName4LocalAdjIdx.setSize(nLocalAdjointStates);

    forAll(regStates_["volVectorStates"], idx)
    {
        word stateName = regStates_["volVectorStates"][idx];
        forAll(mesh_.cells(), cellI)
        {
            for (label i = 0; i < 3; i++)
            {
                label localIdx = this->getLocalAdjointStateIndex(stateName, cellI, i);
                cellIFaceI4LocalIdx[localIdx] = cellI + i / 10.0;

                stateName4LocalAdjIdx[localIdx] = stateName;
            }
        }
    }

    forAll(regStates_["volScalarStates"], idx)
    {
        word stateName = regStates_["volScalarStates"][idx];
        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            cellIFaceI4LocalIdx[localIdx] = cellI;

            stateName4LocalAdjIdx[localIdx] = stateName;
        }
    }

    forAll(regStates_["modelStates"], idx)
    {
        word stateName = regStates_["modelStates"][idx];
        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            cellIFaceI4LocalIdx[localIdx] = cellI;

            stateName4LocalAdjIdx[localIdx] = stateName;
        }
    }

    forAll(regStates_["surfaceScalarStates"], idx)
    {
        word stateName = regStates_["surfaceScalarStates"][idx];
        forAll(mesh_.faces(), faceI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, faceI);
            cellIFaceI4LocalIdx[localIdx] = faceI;

            stateName4LocalAdjIdx[localIdx] = stateName;
        }
    }

    return;
}

label DAIndex::getLocalAdjointStateIndex(
    const word stateName,
    const label idxJ,
    const label comp) const
{
    /*
    Return the global adjoint index given a state name, a local index, 
    and vector component (optional)

    Input:
    -----
    stateName: name of the state variable for the global indexing

    idxJ: the local index for the state variable, typically it is the state's 
    local cell index or face index

    comp: if the state is a vector, give its componet for global indexing. 
    NOTE: for volVectorState, one need to set comp; while for other states, 
    comp is simply ignored in this function

    Example:
    -------
    Image we have two state variables (p, T) and we have three cells, the state
    variable vector reads (state-by-state ordering):

    w= [p0, p1, p2, T0, T1, T2]  <- p0 means p for the 0th cell
         0   1   2   3   4   5   <- adjoint local index

    Then getLocalAdjointStateIndex("p",1) returns 1 and 
    getLocalAdjointStateIndex("T",1) returns 4

    If we use cell-by-cell ordering, the state variable vector reads 
    w= [p0, T0, p1, T1, p2, T2]
         0   1   2   3   4   5   <- adjoint local index
    
    Then getLocalAdjointStateIndex("p",1) returns 2 and 
    getLocalAdjointStateIndex("T",1) returns 3

    Similarly, we can apply this functions for vector state variables, again
    we assume we have two state variables (U, p) and three cells, then the 
    state-by-state adjoint ordering gives

    w= [u0, v0, w0, u1, v1, w1, u2, v2, w2, p0, p1, p2]
         0   1   2   3   4   5   6   7   8   9  10  11 <- adjoint local index

    Then getLocalAdjointStateIndex("U", 1, 2) returns 5 and 
    getLocalAdjointStateIndex("p",1) returns 10

    NOTE: the three compoent for U are [u,v,w]

    */

    word adjJacMatOrdering = daOption_.getOption<word>("adjJacMatOrdering");

    if (adjJacMatOrdering == "state")
    {
        /*
        state by state indexing
        we set u_0, v_1, w_2, u_3, v_4, w_5, ...., p_np, p_np+1, ... nuTilda_nnu, 
        nuTilda_nnu+1, ... so getStateLocalIndexingOffset(p) will return np
        For vector, one need to provide comp, for scalar, comp is not needed.
        */
        forAll(adjStateNames, idxI)
        {
            if (adjStateNames[idxI] == stateName)
            {
                if (adjStateType[stateName] == "volVectorState")
                {
                    if (comp == -1)
                    {
                        FatalErrorIn("") << "comp needs to be set for vector states!"
                                         << abort(FatalError);
                    }
                    else
                    {
                        return stateLocalIndexOffset[stateName] + idxJ * 3 + comp;
                    }
                }
                else
                {
                    return stateLocalIndexOffset[stateName] + idxJ;
                }
            }
        }
    }
    else if (adjJacMatOrdering == "cell")
    {
        // cell by cell ordering
        // We set u_0, v_0, w_0, p_0, nuTilda_0, phi_0a,phi_0b,phi_0c.... u_N, v_N, w_N, p_N, nuTilda_N, phi_N
        // To get the local index, we need to do:
        // idxLocal =
        //   cellI*(nVectorStates*3+nScalarStates+nTurbStates)
        // + phiAccumulatedOffset (how many phis have been accumulated. Note: we have mulitple phis for a cellI)
        // + stateLocalIndexOffset+comp
        // + phiLocalOffset (only for phi idx, ie. 0a or 0b or oc. Note: we have mulitple phis for a cellI)

        label nCellStates = 3 * nVolVectorStates
            + nVolScalarStates
            + nModelStates;

        const word& stateType = adjStateType[stateName];
        label returnV = -99999999;
        if (stateType == "surfaceScalarState") // for surfaceScalarState idxJ is faceI
        {
            label idxN = faceOwner[idxJ]; // idxN is the cell who owns faceI
            returnV = idxN * nCellStates
                + stateLocalIndexOffset[stateName]
                + phiAccumulatdOffset[idxN]
                + phiLocalOffset[idxJ];
            return returnV;
        }
        else if (stateType == "volVectorState") // for other states idxJ is cellI
        {
            if (comp == -1)
            {
                FatalErrorIn("") << "comp needs to be set for vector states!"
                                 << abort(FatalError);
            }
            else
            {

                returnV = idxJ * nCellStates
                    + stateLocalIndexOffset[stateName]
                    + comp
                    + phiAccumulatdOffset[idxJ];
            }
            return returnV;
        }
        else
        {
            returnV = idxJ * nCellStates
                + stateLocalIndexOffset[stateName]
                + phiAccumulatdOffset[idxJ];
            return returnV;
        }
    }
    else
    {
        FatalErrorIn("") << "adjJacMatOrdering invalid" << abort(FatalError);
    }

    // if no stateName found, return an error
    FatalErrorIn("") << "stateName not found!" << abort(FatalError);
    return -1;
}

label DAIndex::getGlobalAdjointStateIndex(
    const word stateName,
    const label idxI,
    const label comp) const
{
    /*
    This function has the same input as DAIndex::getLocalAdjointStateIndex
    the only difference is that this function returns the global adjoint 
    state index by calling globalAdjointStateNumbering.toGlobal()

    Input:
    -----
    stateName: name of the state variable for the global indexing

    idxJ: the local index for the state variable, typically it is the state's 
    local cell index or face index

    comp: if the state is a vector, give its componet for global indexing. 
    NOTE: for volVectorState, one need to set comp; while for other states, 
    comp is simply ignored in this function


    Example:
    -------
    Image we have two state variables (p,T) and five cells, running on two CPU
    processors, the proc0 owns two cells and the proc1 owns three cells,
    then the global adjoint state variables reads (state-by-state)

    w = [p0, p1, T0, T1 | p0, p1, p2, T0, T1, T2] <- p0 means p for the 0th cell on local processor
          0   1   2   3 |  4   5   6   7   8   9  <- global adjoint index
        ---- proc0 -----|--------- proc1 ------- 
    
    Then, on proc0, getGlobalAdjointStateIndex("T", 1) returns 3 
      and on proc1, getGlobalAdjointStateIndex("T", 1) returns 8

    */
    // For vector, one need to provide comp, for scalar, comp is not needed.
    label localIdx = this->getLocalAdjointStateIndex(stateName, idxI, comp);
    label globalIdx = globalAdjointStateNumbering.toGlobal(localIdx);
    return globalIdx;
}

label DAIndex::getGlobalXvIndex(
    const label idxPoint,
    const label idxCoord) const
{
    /*
    This function has the same input as DAIndex::getLocalXvIndex except that 
    this function returns the global xv index 

    Input:
    -----
    idxPoint: local point index

    idxCoord: the compoent of the point

    Example:
    -------
    Image we have three points, running on two CPU cores, and the proc0 owns
    one point and proc1 owns two points, and the Xv vector reads

    Xv = [x0, y0, z0 | x0, y0, z0, x1, y1, z1] <- x0 means the x for the 0th point
           0   1   2    3   4   5   6   7   8  <- global Xv index
          -- proc0 --|--------- proc1 ------- 
    Then, on proc0, getGlobalXvIndex(0,1) returns 1
      and on proc1, getGlobalXvIndex(0,1) returns 4

    */

    label localXvIdx = this->getLocalXvIndex(idxPoint, idxCoord);
    label globalXvIdx = globalXvNumbering.toGlobal(localXvIdx);
    return globalXvIdx;
}

label DAIndex::getLocalXvIndex(
    const label idxPoint,
    const label idxCoord) const
{
    /*
    Returns the local xv index for a given local cell index and its component

    Input:
    -----
    idxPoint: local point index

    idxCoord: the compoent of the point

    Example:
    -------
    Image we have two points, and the Xv vector reads

    Xv = [x0, y0, z0, x1, y1, z1] <- x0 means the x for the 0th point
           0   1   2   3   4   5  <- local Xv index

    Then, getLocalXvIndex(1,1) returns 4

    */

    label localXvIdx = idxPoint * 3 + idxCoord;
    return localXvIdx;
}

void DAIndex::ofField2StateVec(Vec stateVec) const
{
    /*
    Assign values for the state variable vector based on the 
    latest OpenFOAM field values

    Input:
    ------
    OpenFOAM field variables

    Output:
    ------
    stateVec: state variable vector

    Example:
    -------
    Image we have two state variables (p,T) and five cells, running on two CPU
    processors, the proc0 owns two cells and the proc1 owns three cells,
    then calling this function gives the state vector (state-by-state ordering):

    stateVec = [p0, p1, T0, T1 | p0, p1, p2, T0, T1, T2] <- p0 means p for the 0th cell on local processor
                 0   1   2   3 |  4   5   6   7   8   9  <- global state vec index
               ---- proc0 -----|--------- proc1 ------- 
    */

    const objectRegistry& db = mesh_.thisDb();
    PetscScalar* stateVecArray;
    VecGetArray(stateVec, &stateVecArray);

    forAll(regStates_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = this->getLocalAdjointStateIndex(stateName, cellI, comp);
                stateVecArray[localIdx] = state[cellI][comp];
            }
        }
    }

    forAll(regStates_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            stateVecArray[localIdx] = state[cellI];
        }
    }

    forAll(regStates_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            stateVecArray[localIdx] = state[cellI];
        }
    }

    forAll(regStates_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < nLocalInternalFaces)
            {
                stateVecArray[localIdx] = state[faceI];
            }
            else
            {
                label relIdx = faceI - nLocalInternalFaces;
                const label& patchIdx = bFacePatchI[relIdx];
                const label& faceIdx = bFaceFaceI[relIdx];
                stateVecArray[localIdx] = state.boundaryField()[patchIdx][faceIdx];
            }
        }
    }
    VecRestoreArray(stateVec, &stateVecArray);
}

void DAIndex::stateVec2OFField(const Vec stateVec) const
{
    /*
    Assign values OpenFOAM field values based on the state variable vector

    Input:
    ------
    stateVec: state variable vector

    Output:
    ------
    OpenFoam field variables

    Example:
    -------
    Image we have two state variables (p,T) and five cells, running on two CPU
    processors, the proc0 owns two cells and the proc1 owns three cells,
    then calling this function will assign the p, and T based on the the state 
    vector (state-by-state ordering):

    stateVec = [p0, p1, T0, T1 | p0, p1, p2, T0, T1, T2] <- p0 means p for the 0th cell on local processor
                 0   1   2   3 |  4   5   6   7   8   9  <- global state vec index
               ---- proc0 -----|--------- proc1 ------- 
    */

    const objectRegistry& db = mesh_.thisDb();
    const PetscScalar* stateVecArray;
    VecGetArrayRead(stateVec, &stateVecArray);

    forAll(regStates_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = this->getLocalAdjointStateIndex(stateName, cellI, comp);
                state[cellI][comp] = stateVecArray[localIdx];
            }
        }
    }

    forAll(regStates_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            state[cellI] = stateVecArray[localIdx];
        }
    }

    forAll(regStates_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, cellI);
            state[cellI] = stateVecArray[localIdx];
        }
    }

    forAll(regStates_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(regStates_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = this->getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < nLocalInternalFaces)
            {
                state[faceI] = stateVecArray[localIdx];
            }
            else
            {
                label relIdx = faceI - nLocalInternalFaces;
                const label& patchIdx = bFacePatchI[relIdx];
                const label& faceIdx = bFaceFaceI[relIdx];
                state.boundaryFieldRef()[patchIdx][faceIdx] = stateVecArray[localIdx];
            }
        }
    }
    VecRestoreArrayRead(stateVec, &stateVecArray);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
