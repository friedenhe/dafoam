/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Compute object functions and their derivatives

\*---------------------------------------------------------------------------*/

#ifndef DAObjFunc_H
#define DAObjFunc_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "DAOption.H"
#include "DATurbulenceModel.H"
#include "DAIndex.H"
#include "topoSetSource.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAObjFunc Declaration
\*---------------------------------------------------------------------------*/

class DAObjFunc
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DAObjFunc(const DAObjFunc&);

    /// Disallow default bitwise assignment
    void operator=(const DAObjFunc&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// dictionary containing the information for the objective function
    const dictionary& objFuncDict_;

    /// DAOption object
    const DAOption& daOption_;

    /// DATurbuleneModel object
    const DATurbulenceModel& daTurb_;

    /// DAIndex object
    const DAIndex& daIndex_;

    /// a sorted list of all face sources for the objective function
    labelList objFuncFaceSources_;

    /// a sorted list of all cell sources for the objective function
    labelList objFuncCellSources_;

    /// value of the obj function on the discrete faces
    scalarList objFuncFaceValues_;

    /// value of the obj function in the discrete cells
    scalarList objFuncCellValues_;

    /// scale of the objective function
    scalar scale_;

    /// the value of the objective function
    scalar objFuncValue_;

    /// the name of the objective function
    word objFuncName_;

public:
    /// Runtime type information
    TypeName("DAObjFunc");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAObjFunc,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& objFuncDict),
        (
            mesh,
            objFuncDict));

    // Constructors

    //- Construct from components
    DAObjFunc(
        const fvMesh& mesh,
        const dictionary& objFuncDict);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAObjFunc> New(
        const fvMesh& mesh,
        const dictionary& objFuncDict);

    //- Destructor
    virtual ~DAObjFunc()
    {
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;

    /// calculate DAObjFunc::objFuncFaceSources_ and DAObjFunc::objFuncCellSources_
    void calcObjFuncSources(
        labelList& faceSources,
        labelList& cellSources);
    
    /// calculate the value of objective function
    virtual void calcObjFunc(
        const labelList& objFuncFaceSources,
        const labelList& objFuncCellSources,
        scalarList& objFuncFaceValues,
        scalarList& objFuncCellValues,
        scalar& objFuncValue) = 0;

    /// calcluate the value of objective function
    virtual scalar getObjFuncValue()
    {
        // calculate 
        this->calcObjFunc(
            objFuncFaceSources_,
            objFuncCellSources_,
            objFuncFaceValues_,
            objFuncCellValues_,
            objFuncValue_);
        
        // return 
        return objFuncValue_;
    }

    /// return the name of objective function
    word getObjFuncName()
    {
        return objFuncName_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
