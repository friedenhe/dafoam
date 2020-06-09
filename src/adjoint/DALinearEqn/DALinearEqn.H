/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Solve the linear euqation using Petsc

\*---------------------------------------------------------------------------*/

#ifndef DALinearEqn_H
#define DALinearEqn_H

#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAUtility.H"
#include "DAStateInfo.H"
#include "DAModel.H"
#include "DAIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DALinearEqn Declaration
\*---------------------------------------------------------------------------*/

class DALinearEqn
{

private:
    /// Disallow default bitwise copy construct
    DALinearEqn(const DALinearEqn&);

    /// Disallow default bitwise assignment
    void operator=(const DALinearEqn&);

protected:
    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// Foam::DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// DAIndex object
    const DAIndex& daIndex_;

public:
    /// Constructors
    DALinearEqn(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    /// Destructor
    virtual ~DALinearEqn()
    {
    }

    // Members
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
