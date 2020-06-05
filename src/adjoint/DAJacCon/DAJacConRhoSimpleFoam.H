/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for DASimpleFoam

\*---------------------------------------------------------------------------*/

#ifndef DAJacConRhoSimpleFoam_H
#define DAJacConRhoSimpleFoam_H

#include "DAJacCon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAJacConRhoSimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAJacConRhoSimpleFoam
    : public DAJacCon
{

public:
    TypeName("DARhoSimpleFoam");
    // Constructors

    //- Construct from components
    DAJacConRhoSimpleFoam(
        const fvMesh& mesh);

    //- Destructor
    virtual ~DAJacConRhoSimpleFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
