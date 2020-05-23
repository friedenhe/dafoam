/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Augmented radiation model for the adjoint method, including residual 
    calculation functions, etc

\*---------------------------------------------------------------------------*/

#ifndef DARadiationModel_H
#define DARadiationModel_H

#include "runTimeSelectionTables.H"
#include "DAOption.H"
#include "DARegState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DARadiationModel Declaration
\*---------------------------------------------------------------------------*/

class DARadiationModel
    : public regIOobject
{

private:
    //- Disallow default bitwise copy construct
    DARadiationModel(const DARadiationModel&);

    //- Disallow default bitwise assignment
    void operator=(const DARadiationModel&);

protected:
    const fvMesh& mesh_;

public:
    //- Runtime type information
    TypeName("DARadiationModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DARadiationModel,
        dictionary,
        (const fvMesh& mesh),
        (mesh));

    // Constructors

    //- Construct from components
    DARadiationModel(const fvMesh& mesh);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DARadiationModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~DARadiationModel()
    {
    }

    // Members

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
