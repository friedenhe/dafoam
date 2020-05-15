/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for simpleFoam

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

private:
public:
    TypeName("DASimpleFoam");
    // Constructors

    //- Construct from components
    DARegStateSimpleFoam(
        const fvMesh& mesh);

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
