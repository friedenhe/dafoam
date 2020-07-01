/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Compute the source term (e.g., actuator disk) for the primal 
        governing equations

\*---------------------------------------------------------------------------*/

#ifndef DASource_H
#define DASource_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAModel.H"
#include "DAIndex.H"
#include "topoSetSource.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DASource Declaration
\*---------------------------------------------------------------------------*/

class DASource
{

private:
    /// Disallow default bitwise copy construct
    DASource(const DASource&);

    /// Disallow default bitwise assignment
    void operator=(const DASource&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// DAIndex object
    const DAIndex& daIndex_;

public:
    /// Runtime type information
    TypeName("DASource");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DASource,
        dictionary,
        (const word modelType,
         const fvMesh& mesh,
         const DAOption& daOption,
         const DAModel& daModel,
         const DAIndex& daIndex),
        (modelType, mesh, daOption, daModel, daIndex));

    // Constructors

    //- Construct from components
    DASource(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DASource> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DASource()
    {
    }

    /// compute the source term
    virtual void calcSource(volVectorField& source) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
