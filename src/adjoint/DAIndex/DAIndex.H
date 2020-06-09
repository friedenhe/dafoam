/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Handel state indexing in serial and in parallel

\*---------------------------------------------------------------------------*/

#ifndef DAIndex_H
#define DAIndex_H

#include "fvOptions.H"
#include "surfaceFields.H"
#include "DAOption.H"
#include "DAUtility.H"
#include "DAStateInfo.H"
#include "DAModel.H"
#include "globalIndex.H"
#include "DAMacroFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAIndex Declaration
\*---------------------------------------------------------------------------*/

class DAIndex
{

private:
    /// Disallow default bitwise copy construct
    DAIndex(const DAIndex&);

    /// Disallow default bitwise assignment
    void operator=(const DAIndex&);

protected:
    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// Foam::DAOption object
    const DAOption& daOption_;

    /// DAModel object
    const DAModel& daModel_;

    /// the StateInfo_ list from DAStateInfo object
    HashTable<wordList> stateInfo_;

public:
    /// Constructors
    DAIndex(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    /// Destructor
    virtual ~DAIndex()
    {
    }

    // Members

    // adjoint indexing info
    /// list of adjoint state names for a specific solver
    List<word> adjStateNames;

    /// hash table of adjoint state types, e.g., volVectorState for a given state name
    HashTable<word> adjStateType;

    // mesh sizes
    // local sizes

    /// local cell size
    label nLocalCells;

    /// local face size
    label nLocalFaces;

    /// local point size
    label nLocalPoints;

    /// local Xv size (point size*3)
    label nLocalXv;

    /// local boundary face size
    label nLocalBoundaryFaces;

    /// local internal face size
    label nLocalInternalFaces;

    /// local boundary patch size
    label nLocalBoundaryPatches;

    /// local coupled boundary patch size
    label nLocalCoupledBFaces;

    /** given a local boundary face index (starting with zero for the first boundary face, 
        and ends with the last boundary faces), return its boundary patchI, 
        e.g., for face index = 100, its patchI = 2
    */
    labelList bFacePatchI;

    /** given a local boundary face index (starting with zero for the first boundary face, 
        and ends with the last boundary faces), return its boundary faceI associated with a boundary patchI,
        e.g., for face index = 100, its patchI = 2, faceI=20
    */
    labelList bFaceFaceI;

    /** hash table of  local state variable index offset, see initializeStateLocalIndexOffset for definition. 
        This will be used in determing the local indexing  for adjoint states. It differs depending on whether
        we use state-by-state or cell-by-cell ordering
    */
    HashTable<label> stateLocalIndexOffset;

    /// a unique number ID for adjoint states, it depends on the sequence of adjStateNames
    HashTable<label> adjStateID;

    // glocal sizes
    /// global cell size
    label nGlobalCells;

    /// global face size
    label nGlobalFaces;

    /// global Xv size (global face size*3)
    label nGlobalXv;

    /// global objective function face size
    label nGlobalObjFuncFaces;

    /// global coupled boundary face size
    label nGlobalCoupledBFaces;

    // adjoint sizes
    /// number of local adjoint states (including all cells and faces)
    label nLocalAdjointStates;

    /// number of global adjoint states (including all cells and faces)
    label nGlobalAdjointStates;

    /// number of state variables for volScalarField
    label nVolScalarStates;

    /// number of state variables for volVectorField
    label nVolVectorStates;

    /// number of state variables for surfaceScalarField
    label nSurfaceScalarStates;

    /// number of model states, NOTE: they are counted separately
    label nModelStates;

    /// \name global indexing variables, they will be used to transfer local and global indices
    //@{
    globalIndex globalAdjointStateNumbering;
    globalIndex globalCellNumbering;
    globalIndex globalCellVectorNumbering; // similar to globalCellNumbering but has 3 components per cell
    globalIndex globalFaceNumbering;
    globalIndex globalCoupledBFaceNumbering;
    globalIndex globalObjFuncGeoNumbering;
    globalIndex globalXvNumbering;
    //@}

    /// a list to map the point index from decomposed domains to the original un-decomposed domain
    labelIOList pointProcAddressing;

    /// number of points for the un-decomposed domain
    label nUndecomposedPoints;

    /// for a given face index, return whether this face is a coupled boundary face
    labelList isCoupledFace;

    /// given a local adjoint state index, return its state name
    wordList adjStateName4LocalAdjIdx;

    /// given a local adjoint state index, return its cell/face index
    scalarList cellIFaceI4LocalAdjIdx;

    /// owner cell of a given face
    labelList faceOwner;

    /// the accumulated phi indexing offset for cell-by-cell indexing
    labelList phiAccumulatdOffset;

    /// phi local indexing offset for cell-by-cell indexing
    labelList phiLocalOffset;

    /*
    // Member functions

    /// get local segregated adjoint index for a given state name, cell/face indxI and its component (optional, only for vector states)
    label getLocalSegregatedAdjointStateIndex
    (
        const word stateName,
        const label idxI,
        label comp=-1
    );
    
    /// get global segregated adjoint index for a given state name, cell/face indxI and its component (optional, only for vector states)
    label getGlobalSegregatedAdjointStateIndex
    (
        const word stateName,
        const label idxI,
        label comp=-1
    );
        
    /// get the objective function geometry information from the adjIO object
    List<word> getObjFuncGeoInfo(const word objFunc);

    /// get the number of objective function geometry elements, this could be the total number of discrete surfaces for all the design surfaces, or total cells for all the user-defined volumes. 
    label getNLocalObjFuncGeoElements(const word objFunc);
    
    void calcCellIFaceI4GlobalAdjIdx(scalarList& cellIFaceI4GlobalAdjIdx);
    
    /// given a boundary patchI and faceI, return its local boundary face index
    label BFacePatchIFaceI2LocalIndex(const label patchI, const label faceI);
    
    /// given a local boundary face index, return its boundary patchI and faceI
    void BFaceLocalIndex2PatchIFaceI
    (
        const label idxI,
        label& patchIdx,
        label& faceIdx
    );
    
    /// write adjoint indexing information to disk for debugging
    void writeAdjointIndexing();
    
    /// how many non-zeros in the petsc matrix
    void getMatNonZeros
    (
        Mat matIn,
        label& maxCols, 
        scalar& allNonZeros
    );
    
    /// print matrix non-zero information (e.g., diag/off-diag ratio)
    void printMatChars(Mat matIn);

    /// initialize globalObjFuncGeoNumbering for a given objective function   
    void initializeObjFuncGeoNumbering(const word objFunc);

    /// delete globalObjFuncGeoNumbering
    void deleteObjFuncGeoNumbering();
    
    /// check if a patch is a user-defined patch
    label isUserDefinedPatch(word geoInfo);

    /// check if a volume is user-defined volume
    label isUserDefinedVolume(word geoInfo);

    /// calculate the face indices for all the user-defined patches: faceIdx4UserDefinedPatches
    void calcFaceIndx4UserDefinedPatches();

    /// calculate the cell indices for all the user-defined volumes: cellIdx4UserDefinedVolumes
    void calcCellIndx4UserDefinedVolumes();
    
    ///  face indices for all the user-defined patches
    HashTable< labelList > faceIdx4UserDefinedPatches;

    ///  cell indices for all the user-defined volumes
    HashTable< labelList > cellIdx4UserDefinedVolumes;
    */

    // Member functions

    /// calculate stateLocalIndexOffset
    void calcStateLocalIndexOffset(HashTable<label>& offset);

    /// set adjoint state unique ID: adjStateID
    void calcAdjStateID(HashTable<label>& adjStateID);

    /// compute local lists such as adjStateName4LocalAdjIdx and  cellIFaceI4LocalAdjIdx;
    void calcLocalIdxLists(
        wordList& adjStateName4LocalAdjIdx,
        scalarList& cellIFaceI4LocalAdjIdx);

    /// get local adjoint index for a given state name, cell/face indxI and its component (optional, only for vector states)
    label getLocalAdjointStateIndex(
        const word stateName,
        const label idxI,
        const label comp = -1) const;

    /// get global adjoint index for a given state name, cell/face indxI and its component (optional, only for vector states)
    label getGlobalAdjointStateIndex(
        const word stateName,
        const label idxI,
        const label comp = -1) const;

    /// get global Xv index for a given point index and coordinate component (x, y, or z)
    label getGlobalXvIndex(
        const label idxPoint,
        const label idxCoord) const;

    /// get local Xv index for a given point index and coordinate component (x, y, or z)
    label getLocalXvIndex(
        const label idxPoint,
        const label idxCoord) const;

    /// compute global list adjStateID4GlobalAdjIdx that stores the stateID for a given globalAdjIndx
    void calcAdjStateID4GlobalAdjIdx(labelList& adjStateID4GlobalAdjIdx) const;

    /// print all the sizes and indices
    void printIndices() const
    {
        //Info<<"adjStateNames"<<adjStateNames<<endl;
        Info << "Adjoint States: " << adjStateType << endl;

        // Print relevant sizes to screen
        Info << "Global Cells: " << nGlobalCells << endl;
        Info << "Global Faces: " << nGlobalFaces << endl;
        Info << "Global Xv: " << nGlobalXv << endl;
        Info << "Undecomposed points: " << nUndecomposedPoints << endl;
        Info << "Global Adjoint States: " << nGlobalAdjointStates << endl;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
