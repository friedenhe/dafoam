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
    DASpalartAllmaras(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    //- Destructor
    virtual ~DASpalartAllmaras()
    {
    }

    // Member functions

    /// update the turbulence state for DAStateInfo::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void updateNut();

    /// update turbulence variable boundary values
    virtual void correctTurbBoundaryConditions();

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
