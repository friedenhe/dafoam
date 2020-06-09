/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DARhoSimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DAResidualRhoSimpleFoam_H
#define DAResidualRhoSimpleFoam_H

#include "DAResidual.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAResidualRhoSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAResidualRhoSimpleFoam
    : public DAResidual
{

public:
    TypeName("DARhoSimpleFoam");
    // Constructors

    //- Construct from components
    DAResidualRhoSimpleFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAResidualRhoSimpleFoam()
    {
    }

    // Members

    /// compute residual
    virtual void calcResiduals(const dictionary& options);

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the boundary condition for all the states in the selected solver
    virtual void correctBoundaryConditions();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
