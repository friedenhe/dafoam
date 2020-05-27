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

    /// update the turbulence state for DARegState::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void updateNut();

    /// update turbulence variable boundary values
    virtual void correctTurbBoundaryConditions();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
