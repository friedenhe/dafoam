/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Compute Jacobian connectivity and coloring

\*---------------------------------------------------------------------------*/

#ifndef DAJacCon_H
#define DAJacCon_H

#include "runTimeSelectionTables.H"
#include "fvOptions.H"
#include "DAOption.H"
#include "DAIndex.H"
#include "DAModel.H"
#include "DAStateInfo.H"
#include "syncTools.H"
#include "DAObjFunc.H"
#include "DAColoring.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAJacCon Declaration
\*---------------------------------------------------------------------------*/

class DAJacCon
{

private:
    /// Disallow default bitwise copy construct
    DAJacCon(const DAJacCon&);

    /// Disallow default bitwise assignment
    void operator=(const DAJacCon&);

protected:
    /// fvMesh
    const fvMesh& mesh_;

    /// DAUtility object
    DAUtility daUtil_;

    /// DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// DAIndex object
    DAIndex daIndex_;

    /// DAColoring object
    DAColoring daColoring_;

    /// the regState_ list from DAStateInfo object
    HashTable<wordList> stateInfo_;

    /// matrix to store boundary connectivity levels for state Jacobians
    Mat stateBoundaryCon_;

    /// matrix to store boundary connectivity ID for state Jacobians
    Mat stateBoundaryConID_;

    /// neibough face global index for a given local boundary face
    labelList neiBFaceGlobalCompact_;

public:
    /// Runtime type information
    TypeName("DAJacCon");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAJacCon,
        dictionary,
        (const word modelType,
         const fvMesh& mesh,
         const DAOption& daOption,
         const DAModel& daModel),
        (modelType,
         mesh,
         daOption,
         daModel));

    // Constructors

    //- Construct from components
    DAJacCon(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAJacCon> New(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAJacCon()
    {
    }

    // Member functions

    /// initialize state boundary connection
    void initializeStateBoundaryCon();

    /// calculate DAJacCon::neiBFaceGlobalCompact_
    void calcNeiBFaceGlobalCompact(labelList& neiBFaceGlobalCompact);

    /// given a local face index, return the local index of the coupled boundary face
    label getLocalCoupledBFaceIndex(const label localFaceI) const;

    /// calculate DAJacCon::stateBoundaryCon_
    void setupStateBoundaryCon(Mat* stateBoundaryCon);

    /// calculate DAJacCon::stateBoundaryConID_
    void setupStateBoundaryConID(Mat* stateBoundaryConID);

    /// allocate connectedState matrix
    void createConnectionMat(Mat* connectedStates);

    /// a high-level function to add connected state column indices to the connectivity matrix
    void addStateConnections(
        Mat connections,
        const label cellI,
        const label connectedLevelLocal,
        const wordList connectedStatesLocal,
        const List<List<word>> connectedStateInterProc,
        const label addFace);

    /// add value 1 for the colume idx to conMat
    void setConnections(
        Mat conMat,
        const label idx) const;

    /// assign values in connections to a specific row idxI in conMat
    void setupJacobianConnections(
        Mat conMat,
        Mat connections,
        const label idxI);

    /// combine stateBoundaryCon and stateBoundaryConTmp, this is to ensure including all connected states for parallel cases
    void combineStateBndCon(
        Mat* stateBoundaryCon,
        Mat* stateBoundaryConTmp);

    /// add val to the gRow row in conMat, the column indice are the state (given by stateName) for a given cellI
    void addConMatCell(
        Mat conMat,
        const label gRow,
        const label cellI,
        const word stateName,
        const scalar val);

    /// add val to gRow row in conMat, column indice are the neighbouring states (given by stateName) for a given cellI
    void addConMatNeighbourCells(
        Mat conMat,
        const label gRow,
        const label cellI,
        const word stateName,
        const scalar val);

    /// add val to gRow row in conMat, column indice are the faces (given by stateName) for a given cellI
    void addConMatCellFaces(
        Mat conMat,
        const label gRow,
        const label cellI,
        const word stateName,
        const scalar val);

    /// add the column index of the (iner-proc) connected states and faces to conMat, given a local face index
    void addBoundaryFaceConnections(
        Mat conMat,
        const label gRow,
        const label cellI,
        const labelList v,
        const List<List<word>> connectedStates,
        const label addFaces);

    // ****************** Pure virutal function ********************

    /// assign 1 to all non-zero elements for the Jacobian connecitivyt matrix
    virtual void setupJacCon(const dictionary& options) = 0;

    /// initialize the state Jacobian connectivity matrix
    virtual void initializeJacCon(const dictionary& options) = 0;

    virtual void deleteJacCon(const dictionary& options) = 0;

    /// get the number of dRdW colors
    virtual label getNJacConColors() const = 0;

    /// return DAJacCon::dRdWColoredColumns_
    virtual Vec getJacConColoredColumns() const = 0;

    /// compute graph coloring for Jacobian connectivity matrix
    virtual void calcJacConColoring() = 0;

    /// read colors for dRdW
    virtual void readJacConColoring() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
