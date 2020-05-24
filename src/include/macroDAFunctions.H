/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    Macro functions used in the AdjointDerivative class

\*---------------------------------------------------------------------------*/

// make state
#define makeState(regState, fieldType, db) \
    const word stateName = regState;       \
    fieldType& state(                      \
        const_cast<fieldType&>(            \
            db.lookupObject<fieldType>(stateName)));
