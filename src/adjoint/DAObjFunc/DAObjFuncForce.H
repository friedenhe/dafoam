/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for force 

\*---------------------------------------------------------------------------*/

#ifndef DAObjFuncForce_H
#define DAObjFuncForce_H

#include "DAObjFunc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAObjFuncForce Declaration
\*---------------------------------------------------------------------------*/

class DAObjFuncForce
    : public DAObjFunc
{

protected:
    /// the direction of the force
    vector forceDir_;

    /// if dynamically adjusting the angle what mode to use
    word dirMode_;

    /// the patch name from the alpha design variable dict
    word inoutPatch_;

    /// xAxisIndex from the alpha design variable dict with aoa=tan(U_y/U_x)
    label xAxisIndex_;

    /// yAxisIndex from the alpha design variable dict with aoa=tan(U_y/U_x)
    label yAxisIndex_;

    /// dynamically adjust the force direction
    void updateForceDir(vector& forceDir);

    /// DATurbulenceModel object
    const DATurbulenceModel& daTurb_;

public:
    TypeName("force");
    // Constructors

    //- Construct from components
    DAObjFuncForce(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAResidual& daResidual,
        const word objFuncName,
        const word objFuncPart,
        const dictionary& objFuncDict);

    //- Destructor
    virtual ~DAObjFuncForce()
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
