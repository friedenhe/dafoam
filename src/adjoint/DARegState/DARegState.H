/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register states for the adjoint, compute the derived state information
        such as reference state and state residuals

\*---------------------------------------------------------------------------*/

#ifndef DARegState_H
#define DARegState_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DARegState Declaration
\*---------------------------------------------------------------------------*/

class DARegState
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DARegState(const DARegState&);

    /// Disallow default bitwise assignment
    void operator=(const DARegState&);

protected:
    /// registered states 1st hash: solverName, 2nd hash: fieldType, 3nd list, stateNames
    HashTable<wordList> regStates_;

public:
    /// Runtime type information
    TypeName("DARegState");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DARegState,
        dictionary,
        (
            const fvMesh& mesh),
        (mesh));

    // Constructors

    //- Construct from components
    DARegState(const fvMesh& mesh);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DARegState> New(const fvMesh& mesh);

    //- Destructor
    virtual ~DARegState()
    {
    }

    /// return the reference of regStates
    const HashTable<wordList>& getRegStates()
    {
        return regStates_;
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
