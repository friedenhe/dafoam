/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Child class for the SpalartAllmaras model

\*---------------------------------------------------------------------------*/

#ifndef DASpalartAllmaras_H
#define DASpalartAllmaras_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class DASpalartAllmaras
    : public DATurbulenceModel
{

protected:

public:
    TypeName("SpalartAllmaras");
    // Constructors

    //- Construct from components
    DASpalartAllmaras(const fvMesh& mesh);

    //- Destructor
    virtual ~DASpalartAllmaras()
    {
    }

    // Member functions


    /// update nut
    virtual void correctNut();

    /// update the turbulence state for DARegState::regStates_
    virtual void correctModelStates(wordList& modelStates);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
