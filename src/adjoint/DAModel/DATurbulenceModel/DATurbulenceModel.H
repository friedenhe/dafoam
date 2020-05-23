/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Augmented turbulence model for the adjoint method, including residual 
    calculation functions, etc

    NOTE:
    Instead of inheriting from the OpenFOAM turbulence implementation, in 
    RASModel's children, we re-write all the correspndong functions for each 
    turbulence model. This is to avoid using template classes and template 
    functions for all the other classes in DAFoam. The downside is that we 
    need to update all the RASModel's children when upgrading to a new version 
    of OpenFOAM. Hopefully, the turbulence model part does not change too much 
    from version to version so the modification will be minimal. 

\*---------------------------------------------------------------------------*/

#ifndef DATurbulenceModel_H
#define DATurbulenceModel_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "geometricOneField.H"
#include "wallDist.H"
#include "turbulenceModel.H"
#include "nearWallDist.H"
#include "bound.H"
#include "DAOption.H"
#include "DARegState.H"

#ifdef IncompressibleFlow
#include "singlePhaseTransportModel.H"
#endif

#ifdef CompressibleFlow
#include "fluidThermo.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DATurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class DATurbulenceModel
    : public regIOobject
{

private:
    //- Disallow default bitwise copy construct
    DATurbulenceModel(const DATurbulenceModel&);

    //- Disallow default bitwise assignment
    void operator=(const DATurbulenceModel&);

protected:
    const fvMesh& mesh_;

public:
    //- Runtime type information
    TypeName("DATurbulenceModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DATurbulenceModel,
        dictionary,
        (const fvMesh& mesh),
        (mesh));

    // Constructors

    //- Construct from components
    DATurbulenceModel(const fvMesh& mesh);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DATurbulenceModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~DATurbulenceModel()
    {
    }

    // Members

    /// update nut
    virtual void correctNut()=0;

    /// update the turbulence state for DARegState::regStates_
    virtual void correctModelStates(wordList& modelStates)=0;

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
