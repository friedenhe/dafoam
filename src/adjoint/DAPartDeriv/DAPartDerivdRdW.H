/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for dRdW

\*---------------------------------------------------------------------------*/

#ifndef DAPartDerivdRdW_H
#define DAPartDerivdRdW_H

#include "DAPartDeriv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAPartDerivdRdW Declaration
\*---------------------------------------------------------------------------*/

class DAPartDerivdRdW
    : public DAPartDeriv
{

protected:
public:
    TypeName("dRdW");
    // Constructors

    //- Construct from components
    DAPartDerivdRdW(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAJacCon& daJacCon);

    //- Destructor
    virtual ~DAPartDerivdRdW()
    {
    }

    // Member functions

    /// initialize partial derivative matrix/vector
    virtual void initializePartDeriv(
        Mat* jacMat,
        const dictionary& options);
    
    /// compute the partial derivative matrix/vector
    virtual void calcPartDeriv(const dictionary& options);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
