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
#include "DARegState.H"
#include "syncTools.H"
#include "DAObjFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DAJacCon Declaration
\*---------------------------------------------------------------------------*/

class DAJacCon
    : public regIOobject
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

    /// DAIndex object
    const DAIndex& daIndex_;

    /// DARegState object
    const DARegState& daRegState_;

    /// the regState_ list from DARegState object
    const HashTable<wordList>& regStates_;

    /// table to specify how the states are connected to the residuals for a given solver
    HashTable<List<List<word>>> adjStateResidualConInfo_;

    /// a backup for DAJacCon::adjStateResidualConInfo_ when reducing the connectivity levels
    HashTable<List<List<word>>> adjStateResidualConInfoBK_;

    /// matrix to store boundary connectivity levels for state Jacobians
    Mat stateBoundaryCon_;

    /// matrix to store boundary connectivity ID for state Jacobians
    Mat stateBoundaryConID_;

    /// neibough face global index for a given local boundary face
    labelList neiBFaceGlobalCompact_;

    /// dRdW connectivity mat
    Mat dRdWCon_;

    /// dRdWPC connectivity mat
    Mat dRdWConPC_;

    /// \name on- and off-diagonal preallocation sizes for each row of dRdW, dRdWT, and dRdXv
    //@{
    Vec dRdWTPreallocOn_;
    Vec dRdWTPreallocOff_;
    Vec dRdWPreallocOn_;
    Vec dRdWPreallocOff_;
    //@}

    /// dRdW colors
    Vec dRdWColors_;

    /// dRdW colored columns. It contains the global column index of perturbed state for the current color
    Vec dRdWColoredColumns_;

    /// dFdW colors
    Vec dFdWColors_;

    /// dFdW colored columns. It contains the global column index of perturbed state for the current color
    Vec dFdWColoredColumns_;

    /// dFdW connectivity
    Mat dFdWCon_;

public:
    /// Runtime type information
    TypeName("DAJacCon");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DAJacCon,
        dictionary,
        (const fvMesh& mesh),
        (mesh));

    // Constructors

    //- Construct from components
    DAJacCon(const fvMesh& mesh);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DAJacCon> New(const fvMesh& mesh);

    //- Destructor
    virtual ~DAJacCon()
    {
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;

    // public variables
    /// number of dRdW colors
    label ndRdWColors;

    /// number of dFdW colors
    label ndFdWColors;

    // Member functions

    /// calculate DAJacCon::neiBFaceGlobalCompact_
    void calcNeiBFaceGlobalCompact(labelList& neiBFaceGlobalCompact);

    /// given a local face index, return the local index of the coupled boundary face
    label getLocalCoupledBFaceIndex(const label localFaceI) const;

    /// calculate DAJacCon::stateBoundaryCon_
    void setupStateBoundaryCon(Mat* stateBoundaryCon);

    /// calculate DAJacCon::stateBoundaryConID_
    void setupStateBoundaryConID(Mat* stateBoundaryConID);

    /// assign 1 to all non-zero elements for DAJacCon::dRdWCon_
    void setupdRdWCon(
        const label isPrealloc,
        const label isPC);

    /// reduce the connectivity level for DAJacCon::dRdWCon_
    void reduceAdjStateResidualConLevel();

    /// restore the connectivity level for DAJacCon::dRdWCon_
    void restoreAdjStateResidualConLevel();

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

    /// compute preallocation vectors
    void allocateJacobianConnections(
        Vec preallocOnProc,
        Vec preallocOffProc,
        Vec preallocOnProcT,
        Vec preallocOffProcT,
        Mat connections,
        const label row);

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

    /// initialize petsc vectors
    void initializePetscVecs();

    /// initialize the state Jacobian connectivity matrix
    void initializedRdWCon(const label isPC);

    /// compute preallocation vectors
    void preallocateJacobianMatrix(
        Mat dRMat,
        const Vec preallocOnProc,
        const Vec preallocOffProc);

    /// preallocate memory for dRdW using the computed preallocation vectors
    void preallocatedRdW(
        Mat dRMat,
        const label transposed);

    /// initialize DAJacCon::dFdWCon_
    void initializedFdWCon(const label nLocalObjFuncGeoElements);

    /// assign 1 to all non-zero elements for DAJacCon::dFdWCon_
    void setupObjFuncCon(
        const word objFunc,
        const word mode);

    /// delete DAJacCon::dRdWCon_
    void deletedRdWCon()
    {
        MatDestroy(&dRdWCon_);
        //VecDestroy(&dRdWColoredColumns_);
    }

    /// delete DAJacCon::dRdWConPC_
    void deletedRdWConPC()
    {
        MatDestroy(&dRdWConPC_);
    }

    /// delete DAJacCon::dFdWCon_
    void deletedFdWCon()
    {
        MatDestroy(&dFdWCon_);
        //VecDestroy(&dFdWColoredColumns_);
    }

    /// get the number of dRdW colors
    label getNdRdWColors() const;

    /// get the number of dFdW colors
    label getNdFdWColors() const;

    /// return DAJacCon::dRdWColoredColumns_
    Vec getdRdWColoredColumns() const
    {
        return dRdWColoredColumns_;
    }

    /// return DAJacCon::dFdWColoredColumns_
    Vec getdFdWColoredColumns() const
    {
        return dFdWColoredColumns_;
    }

    /// compute graph coloring for dRdW: DAJacCon::dRdWColors_
    void calcdRdWColoring();

    /// read colors for dRdW
    void readdRdWColoring();

    /// calculate DAJacCon::dFdWColors_
    void calcdFdWColoring(const word objFunc);

    /// read DAJacCon::dFdWColors_
    void readdFdWColoring(const word objFunc);

    /// validate if there is coloring conflict
    void validateColoring(
        Mat conMat,
        Vec colors) const;

    /// number of non-zero elements in a matrix
    void getMatNonZeros(
        const Mat conMat,
        label& maxCols,
        scalar& allNonZeros) const;

    /// find the index of a prescribed value in an array
    label find_index(
        const label target,
        const label start,
        const label size,
        const label* valArray) const;

    /// check if there is non-colored columns
    void coloringComplete(
        const Vec colors,
        label& colorCounter,
        label& notColored) const;

    /// a parallel distance-2 graph coloring function
    void parallelD2Coloring(
        const Mat conMat,
        Vec colors,
        label& nColors) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
