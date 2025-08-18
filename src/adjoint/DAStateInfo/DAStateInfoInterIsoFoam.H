/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DAInterIsoFoam

\*---------------------------------------------------------------------------*/

#ifndef DAStateInfoInterIsoFoam_H
#define DAStateInfoInterIsoFoam_H

#include "DAStateInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAStateInfoInterIsoFoam Declaration
\*---------------------------------------------------------------------------*/

class DAStateInfoInterIsoFoam
    : public DAStateInfo
{

public:
    TypeName("DAInterIsoFoam");
    // Constructors

    //- Construct from components
    DAStateInfoInterIsoFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAStateInfoInterIsoFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
