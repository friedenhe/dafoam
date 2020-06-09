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
    const DAJacCon& daJacCon)
    : modelType_(modelType),
      mesh_(mesh),
      daOption_(daOption),
      daModel_(daModel),
      daJacCon_(daJacCon),
      daIndex_(daIndex),
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
    const DAJacCon& daJacCon)
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
            "    const DAJacCon&"
            ")")
            << "Unknown DAPartDeriv type "
            << modelType << nl << nl
            << "Valid DAPartDeriv types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DAPartDeriv>(
        cstrIter()(modelType, mesh, daOption, daModel, daIndex, daJacCon));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void AdjointDerivative::perturbStates(
    const label colorI,
    const word mode)
{

    // perturb volVectorStates
    forAll(adjReg_.volVectorStates, idxI)
    {
        // create state and stateRef
        makeState(volVectorStates, volVectorField, adjReg_);
        makeStateRef(volVectorStates, volVectorField, adjReg_);

        forAll(mesh_.cells(), cellI)
        {
            // check if this state's color = coloI
            for (label i = 0; i < 3; i++)
            {
                label color = adjCon_.getStateColor(mode, stateName, cellI, i);
                if (color == colorI)
                {
                    scalar normScaling = this->getStateScaling(stateName);
                    state[cellI][i] = stateRef[cellI][i] + adjIO_.epsDeriv * normScaling;
                }
            }
        }
        // correct BC
        state.correctBoundaryConditions();
    }

    // perturb volScalarStates
    forAll(adjReg_.volScalarStates, idxI)
    {
        // create state and stateRef
        makeState(volScalarStates, volScalarField, adjReg_);
        makeStateRef(volScalarStates, volScalarField, adjReg_);

        forAll(mesh_.cells(), cellI)
        {
            // check if this state's color = coloI
            label color = adjCon_.getStateColor(mode, stateName, cellI);
            if (color == colorI)
            {
                scalar normScaling = this->getStateScaling(stateName);
                state[cellI] = stateRef[cellI] + adjIO_.epsDeriv * normScaling;
            }
        }

        // correct BC
        state.correctBoundaryConditions();
    }

    // perturb turbStates
    forAll(adjRAS_.turbStates, idxI)
    {
        // create state and stateRef
        makeState(turbStates, volScalarField, adjRAS_);
        makeStateRef(turbStates, volScalarField, adjRAS_);

        forAll(mesh_.cells(), cellI)
        {
            // check if this state's color = coloI
            label color = adjCon_.getStateColor(mode, stateName, cellI);
            if (color == colorI)
            {
                scalar normScaling = this->getStateScaling(stateName);
                state[cellI] = stateRef[cellI] + adjIO_.epsDeriv * normScaling;
            }
        }
    }
    // BC for turbStates are implemented in the AdjRAS class
    adjRAS_.correctTurbBoundaryConditions();

    // perturb surfaceScalarStates
    forAll(adjReg_.surfaceScalarStates, idxI)
    {
        // create state and stateRef
        makeState(surfaceScalarStates, surfaceScalarField, adjReg_);
        makeStateRef(surfaceScalarStates, surfaceScalarField, adjReg_);

        forAll(mesh_.faces(), faceI)
        {
            // check if this state's color = coloI
            label color = adjCon_.getStateColor(mode, stateName, faceI);
            if (color == colorI)
            {
                scalar normScaling = this->getStateScaling(stateName, faceI);
                if (faceI < adjIdx_.nLocalInternalFaces)
                {
                    state[faceI] = stateRef[faceI] + adjIO_.epsDeriv * normScaling;
                }
                else
                {
                    label relIdx = faceI - adjIdx_.nLocalInternalFaces;
                    label patchIdx = adjIdx_.bFacePatchI[relIdx];
                    label faceIdx = adjIdx_.bFaceFaceI[relIdx];
                    state.boundaryFieldRef()[patchIdx][faceIdx] =
                        stateRef.boundaryField()[patchIdx][faceIdx] + adjIO_.epsDeriv * normScaling;
                }
            }
        }
    }

    // NOTE: we also need to update states that are related to the adjoint states but not
    // perturbed here. For example, in buoyantBoussinesqFoam, p is related to p_rgh;
    // however, we only perturb p_rgh in this function. To calculate the perturbed
    // p due to the p_rgh perturbation for calculating force, we need to do p=p_rgh+rhok*gh
    // Similar treatment is needed for rhok and alphat. Basically, any variables apprear in flow residual
    // calculation or objection function calculation that are not state variables need to be updated.
    // This function is implemented in child class
    this->updateIntermediateVariables();

    // for some special boundary conditions, e.g. pressureInletVelocity, the inlet U depends on
    // rho and phi, so we need to call U.correctBoundaryConditions again
    this->checkSpecialBCs();

    return;
}
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
