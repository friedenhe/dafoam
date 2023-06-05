/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
        Post processing

\*---------------------------------------------------------------------------*/

#ifndef DAPost_H
#define DAPost_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAModel.H"
#include "DATurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAPost Declaration
\*---------------------------------------------------------------------------*/

class DAPost
{

private:
    /// Disallow default bitwise copy construct
    DAPost(const DAPost&);

    /// Disallow default bitwise assignment
    void operator=(const DAPost&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

public:
    /// Runtime type information
    TypeName("DAPost");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAPost,
        dictionary,
        (const word modelType,
         const fvMesh& mesh,
         const DAOption& daOption,
         const DAModel& daModel),
        (modelType, mesh, daOption, daModel));

    // Constructors

    //- Construct from components
    DAPost(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAPost> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAPost()
    {
    }

    virtual void run() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
