/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for drag coefficient CD

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

public:
    TypeName("force");
    // Constructors

    //- Construct from components
    DAObjFuncForce(
        const fvMesh& mesh,
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
