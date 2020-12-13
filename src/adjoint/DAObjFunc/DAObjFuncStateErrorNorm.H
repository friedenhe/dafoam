/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for force 

\*---------------------------------------------------------------------------*/

#ifndef DAObjFuncStateErrorNorm_H
#define DAObjFuncStateErrorNorm_H

#include "DAObjFunc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAObjFuncStateErrorNorm Declaration
\*---------------------------------------------------------------------------*/

class DAObjFuncStateErrorNorm
    : public DAObjFunc
{

protected:

    /// name of the state variable to compute the error norm
    word stateName_;

    /// name of the reference state variable to compute the error norm
    word stateRefName_;

    /// state variable type: can be either scalar or vector
    word stateType_;

    /// the name of the argument variable, e.g., betaSA
    word argVarName_;

    /// reference value of the argument variable to compute the error norm
    scalar argVarRefVal_;

    /// the lambda coefficient of the argument variable
    scalar argVarCoeff_;

    /// scaling factor for the whole objective function
    scalar scale_;

public:
    TypeName("stateErrorNorm");
    // Constructors

    //- Construct from components
    DAObjFuncStateErrorNorm(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAResidual& daResidual,
        const word objFuncName,
        const word objFuncPart,
        const dictionary& objFuncDict);

    //- Destructor
    virtual ~DAObjFuncStateErrorNorm()
    {
    }

    /// calculate the value of objective function
    virtual void calcObjFunc(
        const labelList& objFuncFaceSources,
        const labelList& objFuncCellSources,
        scalarList& objFuncFaceValues,
        scalarList& objFuncCellValues,
        scalar& objFuncValue);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
