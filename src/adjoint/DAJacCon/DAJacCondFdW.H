/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for dFdW

\*---------------------------------------------------------------------------*/

#ifndef DAJacCondFdW_H
#define DAJacCondFdW_H

#include "DAJacCon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAJacCondFdW Declaration
\*---------------------------------------------------------------------------*/

class DAJacCondFdW
    : public DAJacCon
{

protected:

    /// initialize petsc vectors
    void initializePetscVecs();

    /// the global numbering for the discrete source of objective function
    globalIndex globalObjFuncGeoNumbering_;

public:
    TypeName("dFdW");
    // Constructors

    //- Construct from components
    DAJacCondFdW(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAJacCondFdW()
    {
    }

    // Member functions

    /// calculate the
    /// initialize the state Jacobian connectivity matrix
    virtual void initializeJacCon(const dictionary& options);

    /// assign 1 to all non-zero elements for the Jacobian connecitivyt matrix
    virtual void setupJacCon(const dictionary& options);

    /// get the number of JacCon colors
    virtual label getNJacConColors() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
