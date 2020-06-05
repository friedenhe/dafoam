/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DARegStateRhoSimpleFoam_H
#define DARegStateRhoSimpleFoam_H

#include "DARegState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DARegStateRhoSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DARegStateRhoSimpleFoam
    : public DARegState
{

public:
    TypeName("DARhoSimpleFoam");
    // Constructors

    //- Construct from components
    DARegStateRhoSimpleFoam(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DARegStateRhoSimpleFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
