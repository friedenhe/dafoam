/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DARegStateSimpleFoam_H
#define DARegStateSimpleFoam_H

#include "DARegState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DARegStateSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DARegStateSimpleFoam
    : public DARegState
{

public:
    TypeName("DASimpleFoam");
    // Constructors

    //- Construct from components
    DARegStateSimpleFoam(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DARegStateSimpleFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
