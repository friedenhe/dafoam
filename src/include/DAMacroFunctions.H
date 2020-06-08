/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Macro functions for DAFoam

\*---------------------------------------------------------------------------*/

#ifndef DAMacroFunctions_H
#define DAMacroFunctions_H

// make state
#define makeState(regState, fieldType, db) \
    const word stateName = regState;       \
    fieldType& state(                      \
        const_cast<fieldType&>(            \
            db.lookupObject<fieldType>(stateName)));

#define setResidualClassMemberScalar(stateName, stateUnit)            \
    stateName##_(                                                     \
        const_cast<volScalarField&>(                                  \
            mesh.thisDb().lookupObject<volScalarField>(#stateName))), \
        stateName##Ref_(                                              \
            IOobject(                                                 \
                #stateName "Ref",                                     \
                mesh.time().timeName(),                               \
                mesh,                                                 \
                IOobject::NO_READ,                                    \
                IOobject::NO_WRITE),                                  \
            stateName##_),                                            \
        stateName##Res_(                                              \
            IOobject(                                                 \
                #stateName "Res",                                     \
                mesh.time().timeName(),                               \
                mesh,                                                 \
                IOobject::NO_READ,                                    \
                IOobject::NO_WRITE),                                  \
            mesh,                                                     \
            dimensionedScalar(#stateName "Res", stateUnit, 0.0),      \
            zeroGradientFvPatchScalarField::typeName),                \
        stateName##ResRef_(                                           \
            IOobject(                                                 \
                #stateName "ResRef",                                  \
                mesh.time().timeName(),                               \
                mesh,                                                 \
                IOobject::NO_READ,                                    \
                IOobject::NO_WRITE),                                  \
            stateName##Res_),                                         \
        stateName##ResPartDeriv_(                                     \
            IOobject(                                                 \
                #stateName "ResPartDeriv",                            \
                mesh.time().timeName(),                               \
                mesh,                                                 \
                IOobject::NO_READ,                                    \
                IOobject::NO_WRITE),                                  \
            stateName##Res_)

#define setResidualClassMemberVector(stateName, stateUnit)                \
    stateName##_(                                                         \
        const_cast<volVectorField&>(                                      \
            mesh.thisDb().lookupObject<volVectorField>(#stateName))),     \
        stateName##Ref_(                                                  \
            IOobject(                                                     \
                #stateName "Ref",                                         \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            stateName##_),                                                \
        stateName##Res_(                                                  \
            IOobject(                                                     \
                #stateName "Res",                                         \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            mesh,                                                         \
            dimensionedVector(#stateName "Res", stateUnit, vector::zero), \
            zeroGradientFvPatchVectorField::typeName),                    \
        stateName##ResRef_(                                               \
            IOobject(                                                     \
                #stateName "ResRef",                                      \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            stateName##Res_),                                             \
        stateName##ResPartDeriv_(                                         \
            IOobject(                                                     \
                #stateName "ResPartDeriv",                                \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            stateName##Res_)

#define setResidualClassMemberPhi(stateName)                              \
    stateName##_(                                                         \
        const_cast<surfaceScalarField&>(                                  \
            mesh.thisDb().lookupObject<surfaceScalarField>(#stateName))), \
        stateName##Ref_(                                                  \
            IOobject(                                                     \
                #stateName "Ref",                                         \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            stateName##_),                                                \
        stateName##Res_(                                                  \
            IOobject(                                                     \
                #stateName "Res",                                         \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            mesh.thisDb().lookupObject<surfaceScalarField>(#stateName)),  \
        stateName##ResRef_(                                               \
            IOobject(                                                     \
                #stateName "ResRef",                                      \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            mesh.thisDb().lookupObject<surfaceScalarField>(#stateName)),  \
        stateName##ResPartDeriv_(                                         \
            IOobject(                                                     \
                #stateName "ResPartDeriv",                                \
                mesh.time().timeName(),                                   \
                mesh,                                                     \
                IOobject::NO_READ,                                        \
                IOobject::NO_WRITE),                                      \
            mesh.thisDb().lookupObject<surfaceScalarField>(#stateName))

#endif
