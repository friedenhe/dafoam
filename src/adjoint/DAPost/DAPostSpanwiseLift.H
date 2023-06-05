/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
        Output spanwise lift

\*---------------------------------------------------------------------------*/

#ifndef DAPostSpanwiseLift_H
#define DAPostSpanwiseLift_H

#include "DAPost.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAPostSpanwiseLift Declaration
\*---------------------------------------------------------------------------*/

class DAPostSpanwiseLift
    : public DAPost
{

public:
    TypeName("spanwiseLift");
    // Constructors

    //- Construct from components
    DAPostSpanwiseLift(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAPostSpanwiseLift()
    {
    }

    virtual void run();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
