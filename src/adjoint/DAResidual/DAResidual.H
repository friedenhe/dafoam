/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Compute the residual for a given solver

\*---------------------------------------------------------------------------*/

#ifndef DAResidual_H
#define DAResidual_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAModel.H"
#include "DAMacroFunctions.H"
#include "DAUtility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAResidual Declaration
\*---------------------------------------------------------------------------*/

class DAResidual
{

private:
    /// Disallow default bitwise copy construct
    DAResidual(const DAResidual&);

    /// Disallow default bitwise assignment
    void operator=(const DAResidual&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// DAUtility
    DAUtility daUtil_;

public:
    /// Runtime type information
    TypeName("DAResidual");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAResidual,
        dictionary,
        (const word modelType,
         const fvMesh& mesh,
         const DAOption& daOption,
         const DAModel& daModel),
        (modelType, mesh, daOption, daModel));

    // Constructors

    //- Construct from components
    DAResidual(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAResidual> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAResidual()
    {
    }

    // Members

    virtual void calcResiduals(const dictionary& options) = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
