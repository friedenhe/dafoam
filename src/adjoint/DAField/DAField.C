/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAField::DAField(
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel),
      daIndex_(daIndex)
{
    Info << "Initialzing DAField..." << endl;
    // initialize stateInfo_
    word solverName = daOption.getOption<word>("solverName");
    autoPtr<DAStateInfo> daStateInfo(DAStateInfo::New(solverName, mesh, daOption, daModel));
    stateInfo_ = daStateInfo->getStateInfo();
}

void DAField::ofField2StateVec(Vec stateVec) const
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

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI, comp);
                stateVecArray[localIdx] = state[cellI][comp];
            }
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateVecArray[localIdx] = state[cellI];
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateVecArray[localIdx] = state[cellI];
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < daIndex_.nLocalInternalFaces)
            {
                stateVecArray[localIdx] = state[faceI];
            }
            else
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                const label& patchIdx = daIndex_.bFacePatchI[relIdx];
                const label& faceIdx = daIndex_.bFaceFaceI[relIdx];
                stateVecArray[localIdx] = state.boundaryField()[patchIdx][faceIdx];
            }
        }
    }
    VecRestoreArray(stateVec, &stateVecArray);
}

void DAField::stateVec2OFField(const Vec stateVec) const
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

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI, comp);
                state[cellI][comp] = stateVecArray[localIdx];
            }
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            state[cellI] = stateVecArray[localIdx];
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            state[cellI] = stateVecArray[localIdx];
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeState(stateInfo_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < daIndex_.nLocalInternalFaces)
            {
                state[faceI] = stateVecArray[localIdx];
            }
            else
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                const label& patchIdx = daIndex_.bFacePatchI[relIdx];
                const label& faceIdx = daIndex_.bFaceFaceI[relIdx];
                state.boundaryFieldRef()[patchIdx][faceIdx] = stateVecArray[localIdx];
            }
        }
    }
    VecRestoreArrayRead(stateVec, &stateVecArray);
}

void DAField::pointVec2OFMesh(const Vec xvVec) const
{
    /*
    Assign the points in fvMesh of OpenFOAM based on the point vector

    Input:
    -----
    xvVec: a vector that stores the x, y, and z coordinates for all
    points in the fvMesh mesh

    Output:
    ------
    New mesh metrics in fvMesh, effectively by calling mesh.movePoints

    Example:
    --------
    Image we have three points in fvMesh, running on two CPU
    processors, the proc0 owns one point and the proc1 owns two points,
    then calling this function will assign xvVec based on the the points
    coordinates in fvMesh

    xvVec = [x0, y0, z0 | x0, y0, z0, x1, y1, z1] <- x0 means x coordinate for the 0th point on local processor
             0   1   2  |  3   4   5   6   7   8  <- global point vec index
            --- proc0 --|--------- proc1 ------- 
    */

    const PetscScalar* xvVecArray;
    VecGetArrayRead(xvVec, &xvVecArray);

    pointField meshPoints(mesh_.points());

    forAll(mesh_.points(), pointI)
    {
        for (label comp = 0; comp < 3; comp++)
        {
            label localIdx = daIndex_.getLocalXvIndex(pointI, comp);
            meshPoints[pointI][comp] = xvVecArray[localIdx];
        }
    }

    VecRestoreArrayRead(xvVec, &xvVecArray);

    // movePoints update the mesh metrics such as volume, surface area and cell centers
    fvMesh& mesh = const_cast<fvMesh&>(mesh_);
    mesh.movePoints(meshPoints);
}

void DAField::ofMesh2PointVec(Vec xvVec) const
{
    /*
    Assign the point vector based on the points in fvMesh of OpenFOAM

    Input:
    ------
    Mesh coordinates in fvMesh

    Output:
    -----
    xvVec: a vector that stores the x, y, and z coordinates for all
    points in the fvMesh mesh

    Example:
    --------
    Image we have three points in fvMesh, running on two CPU
    processors, the proc0 owns one point and the proc1 owns two points,
    then calling this function will assign xvVec based on the the points
    coordinates in fvMesh

    xvVec = [x0, y0, z0 | x0, y0, z0, x1, y1, z1] <- x0 means x coordinate for the 0th point on local processor
             0   1   2  |  3   4   5   6   7   8  <- global point vec index
            --- proc0 --|--------- proc1 ------- 
    */

    PetscScalar* xvVecArray;
    VecGetArray(xvVec, &xvVecArray);

    forAll(mesh_.points(), pointI)
    {
        for (label comp = 0; comp < 3; comp++)
        {
            label localIdx = daIndex_.getLocalXvIndex(pointI, comp);
            xvVecArray[localIdx] = mesh_.points()[pointI][comp];
        }
    }

    VecRestoreArray(xvVec, &xvVecArray);
}

void DAField::ofResField2ResVec(Vec resVec) const
{
    /*
    Assign values for the residual vector based on the 
    latest OpenFOAM residual field values

    Input:
    ------
    OpenFOAM residual field variables

    Output:
    ------
    resVec: state residual  vector

    Example:
    -------
    Image we have two state residuals (pRes,TRes) and five cells, running on two CPU
    processors, the proc0 owns two cells and the proc1 owns three cells,
    then calling this function gives the residual vector (state-by-state ordering):

    resVec = [pRes0, pRes1, TRes0, TRes1 | pRes0, pRes1, pRes2, TRes0, TRes1, TRes2] 
                 0      1      2      3  |    4      5      6      7      8      9  <- global residual vec index
               ---------- proc0 ---------|------------- proc1 ----------------------
    NOTE: pRes0 means p residual for the 0th cell on local processor
    */

    const objectRegistry& db = mesh_.thisDb();
    PetscScalar* stateResVecArray;
    VecGetArray(resVec, &stateResVecArray);

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI, comp);
                stateResVecArray[localIdx] = stateRes[cellI][comp];
            }
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateResVecArray[localIdx] = stateRes[cellI];
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateResVecArray[localIdx] = stateRes[cellI];
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < daIndex_.nLocalInternalFaces)
            {
                stateResVecArray[localIdx] = stateRes[faceI];
            }
            else
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                const label& patchIdx = daIndex_.bFacePatchI[relIdx];
                const label& faceIdx = daIndex_.bFaceFaceI[relIdx];
                stateResVecArray[localIdx] = stateRes.boundaryField()[patchIdx][faceIdx];
            }
        }
    }
    VecRestoreArray(resVec, &stateResVecArray);
}

void DAField::resVec2OFResField(const Vec resVec) const
{
    /*
    Assign OpenFOAM residual values based on the residual vector

    Input:
    ------
    resVec: residual vector

    Output:
    ------
    OpenFoam field variables

    Example:
    -------
    Image we have two state residuals (pRes,TRes) and five cells, running on two CPU
    processors, the proc0 owns two cells and the proc1 owns three cells,
    then calling this function gives the residual vector (state-by-state ordering):

    resVec = [pRes0, pRes1, TRes0, TRes1 | pRes0, pRes1, pRes2, TRes0, TRes1, TRes2] 
                 0      1      2      3  |    4      5      6      7      8      9  <- global residual vec index
               ---------- proc0 ---------|------------- proc1 ----------------------
    NOTE: pRes0 means p residual for the 0th cell on local processor
    */

    const objectRegistry& db = mesh_.thisDb();
    const PetscScalar* stateResVecArray;
    VecGetArrayRead(resVec, &stateResVecArray);

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["volVectorStates"][idxI], volVectorField, db);

        forAll(mesh_.cells(), cellI)
        {
            for (label comp = 0; comp < 3; comp++)
            {
                label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI, comp);
                stateRes[cellI][comp] = stateResVecArray[localIdx];
            }
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["volScalarStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateRes[cellI] = stateResVecArray[localIdx];
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["modelStates"][idxI], volScalarField, db);

        forAll(mesh_.cells(), cellI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, cellI);
            stateRes[cellI] = stateResVecArray[localIdx];
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        // lookup state from meshDb
        makeStateRes(stateInfo_["surfaceScalarStates"][idxI], surfaceScalarField, db);

        forAll(mesh_.faces(), faceI)
        {
            label localIdx = daIndex_.getLocalAdjointStateIndex(stateName, faceI);
            if (faceI < daIndex_.nLocalInternalFaces)
            {
                stateRes[faceI] = stateResVecArray[localIdx];
            }
            else
            {
                label relIdx = faceI - daIndex_.nLocalInternalFaces;
                const label& patchIdx = daIndex_.bFacePatchI[relIdx];
                const label& faceIdx = daIndex_.bFaceFaceI[relIdx];
                stateRes.boundaryFieldRef()[patchIdx][faceIdx] = stateResVecArray[localIdx];
            }
        }
    }
    VecRestoreArrayRead(resVec, &stateResVecArray);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
