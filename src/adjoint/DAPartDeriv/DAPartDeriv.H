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
#include "DAResidual.H"

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

    /// DAIndex object
    const DAIndex& daIndex_;

    /// DAJacCon object
    const DAJacCon& daJacCon_;

    /// DAResidual object
    const DAResidual& daResidual_;

    /// all the DAFoam option
    const dictionary& allOptions_;

    /// the regState_ list from DAStateInfo object
    HashTable<wordList> stateInfo_;

    /// perturb state variables given a color index
    void perturbStates(
        const Vec jacConColors,
        const label colorI,
        const scalar delta,
        Vec wVec);

    void perturbBC(
        const dictionary options,
        const scalar delta);

    void setPartDerivMat(
        const Vec resVec,
        const Vec coloredColumn,
        const label transposed,
        Mat jacMat) const;

    Mat dXvdFFDMat_;

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
         const DAJacCon& daJacCon,
         const DAResidual& daResidual),
        (modelType,
         mesh,
         daOption,
         daModel,
         daIndex,
         daJacCon,
         daResidual));

    // Constructors

    //- Construct from components
    DAPartDeriv(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAJacCon& daJacCon,
        const DAResidual& daResidual);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAPartDeriv> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const DAJacCon& daJacCon,
        const DAResidual& daResidual);

    //- Destructor
    virtual ~DAPartDeriv()
    {
    }

    // Member functions

    /// initialize partial derivative matrix
    virtual void initializePartDerivMat(
        const dictionary& options,
        Mat* jacMat) = 0;

    /// compute the partial derivative matrix
    virtual void calcPartDerivMat(
        const dictionary& options,
        const Vec xvVec,
        const Vec wVec,
        Mat jacMat) = 0;

    void setdXvdFFDMat(const Mat dXvdFFDMat);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
