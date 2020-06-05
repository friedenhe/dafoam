/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register the names of state variables and how do they connect
        with their residuals

\*---------------------------------------------------------------------------*/

#ifndef DARegState_H
#define DARegState_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DARegState Declaration
\*---------------------------------------------------------------------------*/

class DARegState
{

private:
    /// Disallow default bitwise copy construct
    DARegState(const DARegState&);

    /// Disallow default bitwise assignment
    void operator=(const DARegState&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    const DAOption& daOption_;

    const DAModel& daModel_;

    /// registered states 1st hash: solverName, 2nd hash: fieldType, 3nd list, stateNames
    HashTable<wordList> regStates_;

    /// table to specify how the states are connected to the residuals for a given solver
    HashTable<List<List<word>>> adjStateResidualConInfo_;

    /// a backup for DAJacCon::adjStateResidualConInfo_ when reducing the connectivity levels
    HashTable<List<List<word>>> adjStateResidualConInfoBK_;

public:
    /// Runtime type information
    TypeName("DARegState");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DARegState,
        dictionary,
        (const fvMesh& mesh, const DAOption& daOption, const DAModel& daModel),
        (mesh, daOption, daModel));

    // Constructors

    //- Construct from components
    DARegState(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DARegState> New(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DARegState()
    {
    }

    /// return the reference of regStates_
    const HashTable<wordList>& getRegStates() const
    {
        return regStates_;
    }

    /// return the reference of adjStateResidualConInfo_
    const HashTable<List<List<word>>>& getAdjStateResidualConInfo() const
    {
        return adjStateResidualConInfo_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
