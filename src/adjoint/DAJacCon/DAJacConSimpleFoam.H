/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DAJacConSimpleFoam_H
#define DAJacConSimpleFoam_H

#include "DAJacCon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAJacConSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAJacConSimpleFoam
    : public DAJacCon
{

public:
    TypeName("DASimpleFoam");
    // Constructors

    //- Construct from components
    DAJacConSimpleFoam(
        const fvMesh& mesh);

    //- Destructor
    virtual ~DAJacConSimpleFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
