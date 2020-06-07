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

    /// dFdW connectivity mat
    Mat dFdWCon_;

    /// dFdWPC connectivity mat
    Mat dFdWConPC_;

    /// dFdW colors
    Vec dFdWColors_;

    /// dFdW colored columns. It contains the global column index of perturbed state for the current color
    Vec dFdWColoredColumns_;

    /// number of dFdW colors
    label ndFdWColors_;

public:
    TypeName("DAJacCondFdW");
    // Constructors

    //- Construct from components
    DAJacCondFdW(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAJacCondFdW()
    {
    }

    // Member functions

    /// initialize petsc vectors
    void initializePetscVecs();

    /// initialize the state Jacobian connectivity matrix
    virtual void initializeJacCon(const label isPC = 0);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
