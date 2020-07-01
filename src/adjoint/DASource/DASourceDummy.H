/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for Dummy disk source

\*---------------------------------------------------------------------------*/

#ifndef DASourceDummy_H
#define DASourceDummy_H

#include "DASource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DASourceDummy Declaration
\*---------------------------------------------------------------------------*/

class DASourceDummy
    : public DASource
{

public:
    TypeName("dummy");
    // Constructors

    //- Construct from components
    DASourceDummy(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DASourceDummy()
    {
    }

    /// compute the source term
    virtual void calcSource(volVectorField& source);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
