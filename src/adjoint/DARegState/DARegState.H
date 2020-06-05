/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register state variables for the adjoint

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
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DARegState(const DARegState&);

    /// Disallow default bitwise assignment
    void operator=(const DARegState&);

protected:

    /// fvMesh
    const fvMesh& mesh_;
    
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
        (const fvMesh& mesh),
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
    const HashTable<wordList>& getRegStates() const
    {
        return regStates_;
    }

    void correctModelStates(wordList& modelStates);

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;

    /// return the name of pressure field, it can be either p or p_rgh
    word getPName() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
