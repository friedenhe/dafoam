/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        DAModel is a general template for handling physical models in 
        DAFoam. Here we can simultaneously have multiple models for adjoint,
        e.g., turbulence model + radiation model

\*---------------------------------------------------------------------------*/

#ifndef DAModel_H
#define DAModel_H

#include "fvOptions.H"
#include "DATurbulenceModel.H"
#include "DARadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAModel Declaration
\*---------------------------------------------------------------------------*/

class DAModel
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DAModel(const DAModel&);

    /// Disallow default bitwise assignment
    void operator=(const DAModel&);

    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// have turbulence models for the flow and adjoint?
    label hasTurbulenceModel_ = 0;

    /// have radiation models for the flow and adjoint?
    label hasRadiationModel_ = 0;

public:
    /// Constructors
    DAModel(const fvMesh& mesh);

    /// Destructor
    virtual ~DAModel();

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;

    /// update the name in modelStates based on the selected physical model at runtime
    void correctModelStates(wordList& modelStates);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
