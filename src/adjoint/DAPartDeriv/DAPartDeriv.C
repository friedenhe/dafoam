/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAPartDeriv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DAPartDeriv, 0);
defineRunTimeSelectionTable(DAPartDeriv, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAPartDeriv::DAPartDeriv(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAJacCon& daJacCon,
    const DAResidual& daResidual)
    : modelType_(modelType),
      mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel),
      daIndex_(daIndex),
      daJacCon_(daJacCon),
      daResidual_(daResidual),
      allOptions_(daOption.getAllOptions())
{
    // initialize stateInfo_
    word solverName = daOption.getOption<word>("solverName");
    autoPtr<DAStateInfo> daStateInfo(DAStateInfo::New(solverName, mesh, daOption, daModel));
    stateInfo_ = daStateInfo->getStateInfo();
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DAPartDeriv> DAPartDeriv::New(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex,
    const DAJacCon& daJacCon,
    const DAResidual& daResidual)
{

    Info << "Selecting " << modelType << " for DAPartDeriv" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DAPartDeriv::New"
            "("
            "    const word,"
            "    const fvMesh&,"
            "    const DAOption&,"
            "    const DAModel&,"
            "    const DAIndex&,"
            "    const DAJacCon&,"
            "    const DAResidual&"
            ")")
            << "Unknown DAPartDeriv type "
            << modelType << nl << nl
            << "Valid DAPartDeriv types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DAPartDeriv>(
        cstrIter()(modelType,
                   mesh,
                   daOption,
                   daModel,
                   daIndex,
                   daJacCon,
                   daResidual));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAPartDeriv::perturbStates(
    const Vec jacConColors,
    const Vec normStatePerturbVec,
    const label colorI,
    const scalar delta,
    Vec wVec)
{
    /*
    Descripton:
        Perturb state variable vec such that it can be used to compute 
        perturbed residual vector later
    
    Input:
        jacConColors: the coloring vector for this Jacobian, obtained from Foam::DAJacCon
    
        normStatePerturbVec: state normalization vector, 1.0 means no normalization, the 
        actually perturbation added on the wVec is normStatePerturbVec * delta

        colorI: we perturb the rows that associated with the coloring index colorI

        delta: the delta value for perturbation the actually perturbation added on 
        the wVec is normStatePerturbVec * delta

    Input/Output:
        wVec: the perturbed state variable vector

    Example:
        Assuming we have five element in wVec {0.1, 0.5, 1.0, 1.5, 2.0}
        If jacConColors reads {0, 0, 1, 0, 1},
        normStatePerturbVec reads {1.0, 1.0, 0.5, 1.0, 0.1}
        colorI = 1, delta = 0.01

        Then after calling this function wVec reads
        {0.1, 0.5, 1.005, 1.5, 2.001}
    */
    PetscInt Istart, Iend;
    VecGetOwnershipRange(jacConColors, &Istart, &Iend);

    const PetscScalar* colorArray;
    VecGetArrayRead(jacConColors, &colorArray);

    const PetscScalar* normStateArray;
    VecGetArrayRead(normStatePerturbVec, &normStateArray);

    PetscScalar* wVecArray;
    VecGetArray(wVec, &wVecArray);

    for (label i = Istart; i < Iend; i++)
    {
        label relIdx = i - Istart;
        label colorJ = colorArray[relIdx];
        if (colorI == colorJ)
        {
            wVecArray[relIdx] += delta * normStateArray[relIdx];
        }
    }

    VecRestoreArrayRead(jacConColors, &colorArray);
    VecRestoreArray(wVec, &wVecArray);
    VecRestoreArrayRead(normStatePerturbVec, &normStateArray);

    return;
}

void DAPartDeriv::setPartDerivMat(
    const Vec resVec,
    const Vec coloredColumn,
    const label transposed,
    Mat jacMat) const
{
    /*
    Description:
        Set the values from resVec to jacMat
    
    Input:
        resVec: residual vector, obtained after calling the DAResidual::masterFunction

        coloredColumn: a vector to determine the element in resVec is associated with 
        which column in the jacMat. coloredColumn is computed in DAJacCon::calcColoredColumns
        -1 in coloredColumn means don't set values to jacMat for this column

        transposed: whether the jacMat is transposed or not, for dRdWT transpoed = 1, for 
        all the other cases, set it to 0

    Output:
        jacMat: the jacobian matrix to set

    Example:
        Condering a 5 by 5 jacMat with all zeros, and

        resVec
        {
            0.1,
            0.2,
            0.3,
            0.4,
            0.5
        }

        coloredColumn
        {
            -1
            1
            3
            -1
            0
        }

        transposed = 0

        Then, after calling this function, jacMat reads

        0.0  0.0  0.0  0.0  0.0  
        0.0  0.2  0.0  0.0  0.0  
        0.0  0.0  0.0  0.3  0.0  
        0.0  0.0  0.0  0.0  0.0  
        0.5  0.0  0.0  0.0  0.0  
    */

    label rowI, colI;
    scalar val;
    PetscInt Istart, Iend;
    const PetscScalar* resVecArray;
    const PetscScalar* coloredColumnArray;

    VecGetArrayRead(resVec, &resVecArray);
    VecGetArrayRead(coloredColumn, &coloredColumnArray);

    // get the local ownership range
    VecGetOwnershipRange(resVec, &Istart, &Iend);

    // Loop over the owned values of this row and set the corresponding
    // Jacobian entries
    for (PetscInt i = Istart; i < Iend; i++)
    {
        label relIdx = i - Istart;
        colI = coloredColumnArray[relIdx];
        if (colI >= 0)
        {
            rowI = i;
            val = resVecArray[relIdx];

            if (transposed)
            {
                MatSetValue(jacMat, colI, rowI, val, INSERT_VALUES);
            }
            else
            {
                MatSetValue(jacMat, rowI, colI, val, INSERT_VALUES);
            }
        }
    }

    VecRestoreArrayRead(resVec, &resVecArray);
    VecRestoreArrayRead(coloredColumn, &coloredColumnArray);
}

void DAPartDeriv::perturbBC(
    const dictionary options,
    const scalar delta)
{
    /*
    Descripton:
        Perturb values in boundary conditions of the OpenFOAM fields
    
    Input:
        options.varName: the name of the variable to perturb

        options.patchName: the name of the boundary patches to perturb, do not support
        multiple patches yet
    
        options.fieldType: the type of the OpenFOAM field for this variable

        options.bcType: the boundary type, fixedValue?

        optoins.comp: the component to perturb

        delta: the delta value to perturb
    */

    word varName, patchName, fieldType, bcType;
    label comp;
    options.readEntry<word>("varName", varName);
    options.readEntry<word>("patchName", patchName);
    options.readEntry<word>("fieldType", fieldType);
    options.readEntry<word>("bcType", bcType);
    options.readEntry<label>("comp", comp);

    label patchI = mesh_.boundaryMesh().findPatchID(patchName);

    if (bcType == "fixedValue")
    {
        if (fieldType == "volVectorField")
        {
            volVectorField& var(
                const_cast<volVectorField&>(
                    mesh_.thisDb().lookupObject<volVectorField>(varName)));

            if (mesh_.boundaryMesh()[patchI].size() > 0)
            {
                forAll(var.boundaryField()[patchI], faceI)
                {
                    var.boundaryFieldRef()[patchI][faceI][comp] += delta;
                }
            }
        }
        else if (fieldType == "volScalarField")
        {
            volScalarField& var(
                const_cast<volScalarField&>(
                    mesh_.thisDb().lookupObject<volScalarField>(varName)));

            if (mesh_.boundaryMesh()[patchI].size() > 0)
            {
                forAll(var.boundaryField()[patchI], faceI)
                {
                    var.boundaryFieldRef()[patchI][faceI] += delta;
                }
            }
        }
        else
        {
            FatalErrorIn("") << "fieldType: "
                             << fieldType << " not supported"
                             << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("") << "bcType: "
                         << bcType << " not supported"
                         << abort(FatalError);
    }
}

void DAPartDeriv::setdXvdFFDMat(const Mat dXvdFFDMat)
{
    /*
    Description:
        Assign value to dXvdFFDMat_ basically we do a MatConvert
    */

    MatConvert(dXvdFFDMat, MATSAME, MAT_INITIAL_MATRIX, &dXvdFFDMat_);
    //MatDuplicate(dXvdFFDMat, MAT_COPY_VALUES, &dXvdFFDMat_);
    MatAssemblyBegin(dXvdFFDMat_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dXvdFFDMat_, MAT_FINAL_ASSEMBLY);
}

void DAPartDeriv::setNormStatePerturbVec(Vec* normStatePerturbVec)
{
    /*
    Description:
        Set the values for the normStatePerturbVec

    Input/Output:
        normStatePerturbVec: the vector to store the state normalization
        values, this will be used in DAPartDeriv::perturbStates

    The normalization referene values are set in "normalizeStates" in DAOption
    */
    label localSize = daIndex_.nLocalAdjointStates;
    VecCreate(PETSC_COMM_WORLD, normStatePerturbVec);
    VecSetSizes(*normStatePerturbVec, localSize, PETSC_DETERMINE);
    VecSetFromOptions(*normStatePerturbVec);
    VecSet(*normStatePerturbVec, 1.0);

    dictionary normStateDict = allOptions_.subDict("normalizeStates");

    wordList normStateNames = normStateDict.toc();

    forAll(stateInfo_["volVectorStates"], idxI)
    {
        word stateName = stateInfo_["volVectorStates"][idxI];
        if (DAUtility::isInList<word>(stateName, normStateNames))
        {
            scalar scale = normStateDict.getScalar(stateName);
            forAll(mesh_.cells(), idxI)
            {
                for (label comp = 0; comp < 3; comp++)
                {
                    label glbIdx = daIndex_.getGlobalAdjointStateIndex(stateName, idxI, comp);
                    VecSetValue(*normStatePerturbVec, glbIdx, scale, INSERT_VALUES);
                }
            }
        }
    }

    forAll(stateInfo_["volScalarStates"], idxI)
    {
        const word stateName = stateInfo_["volScalarStates"][idxI];
        if (DAUtility::isInList<word>(stateName, normStateNames))
        {
            scalar scale = normStateDict.getScalar(stateName);
            forAll(mesh_.cells(), idxI)
            {
                label glbIdx = daIndex_.getGlobalAdjointStateIndex(stateName, idxI);
                VecSetValue(*normStatePerturbVec, glbIdx, scale, INSERT_VALUES);
            }
        }
    }

    forAll(stateInfo_["modelStates"], idxI)
    {
        const word stateName = stateInfo_["modelStates"][idxI];
        if (DAUtility::isInList<word>(stateName, normStateNames))
        {
            scalar scale = normStateDict.getScalar(stateName);
            forAll(mesh_.cells(), idxI)
            {
                label glbIdx = daIndex_.getGlobalAdjointStateIndex(stateName, idxI);
                VecSetValue(*normStatePerturbVec, glbIdx, scale, INSERT_VALUES);
            }
        }
    }

    forAll(stateInfo_["surfaceScalarStates"], idxI)
    {
        const word stateName = stateInfo_["surfaceScalarStates"][idxI];
        if (DAUtility::isInList<word>(stateName, normStateNames))
        {
            scalar scale = normStateDict.getScalar(stateName);

            forAll(mesh_.faces(), faceI)
            {
                if (faceI < daIndex_.nLocalInternalFaces)
                {
                    scale = mesh_.magSf()[faceI];
                }
                else
                {
                    label relIdx = faceI - daIndex_.nLocalInternalFaces;
                    label patchIdx = daIndex_.bFacePatchI[relIdx];
                    label faceIdx = daIndex_.bFaceFaceI[relIdx];
                    scale = mesh_.magSf().boundaryField()[patchIdx][faceIdx];
                }

                label glbIdx = daIndex_.getGlobalAdjointStateIndex(stateName, faceI);
                VecSetValue(*normStatePerturbVec, glbIdx, scale, INSERT_VALUES);
            }
        }
    }

    VecAssemblyBegin(*normStatePerturbVec);
    VecAssemblyEnd(*normStatePerturbVec);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
