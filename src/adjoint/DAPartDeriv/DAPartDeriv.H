/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Compute partial derivatives using the finite-difference method

\*---------------------------------------------------------------------------*/

#ifndef DAPartDeriv_H
#define DAPartDeriv_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "DAUtility.H"
#include "DAOption.H"
#include "DAIndex.H"
#include "DAModel.H"
#include "DAStateInfo.H"
#include "syncTools.H"
#include "DAObjFunc.H"
#include "DAJacCon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAPartDeriv Declaration
\*---------------------------------------------------------------------------*/

class DAPartDeriv
{

private:
    /// Disallow default bitwise copy construct
    DAPartDeriv(const DAPartDeriv&);

    /// Disallow default bitwise assignment
    void operator=(const DAPartDeriv&);

protected:
    /// the name of the jacCon matrix
    const word modelType_;

    /// fvMesh
    const fvMesh& mesh_;

    /// DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// DAJacCon object
    const DAJacCon& daJacCon_;

    /// DAIndex object
    const DAIndex& daIndex_;

    /// all the DAFoam option
    const dictionary& allOptions_;

    /// the regState_ list from DAStateInfo object
    HashTable<wordList> stateInfo_;

    /// perturb state variables given a color index
    void perturbStates(
        const label colorI,
        const word mode);

public:
    /// Runtime type information
    TypeName("DAPartDeriv");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAPartDeriv,
        dictionary,
        (const word modelType,
         const fvMesh& mesh,
         const DAOption& daOption,
         const DAModel& daModel,
         const DAIndex& daIndex,
         const DAJacCon& daJacCon),
        (modelType,
         mesh,
         daOption,
         daModel,
         daIndex,
         daJacCon));

    // Constructors

    //- Construct from components
    DAPartDeriv(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAJacCon& daJacCon);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAPartDeriv> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAJacCon& daJacCon);

    //- Destructor
    virtual ~DAPartDeriv()
    {
    }

    // Member functions

    label getRunTime()
    {
        return mesh_.time().elapsedClockTime();
    }

    /// initialize partial derivative matrix
    virtual void initializePartDeriv(
        Mat* jacMat,
        const dictionary& options) = 0;

    /// compute the partial derivative matrix
    virtual void calcPartDeriv(const dictionary& options) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
