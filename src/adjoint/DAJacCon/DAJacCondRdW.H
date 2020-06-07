/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for dRdW

\*---------------------------------------------------------------------------*/

#ifndef DAJacCondRdW_H
#define DAJacCondRdW_H

#include "DAJacCon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAJacCondRdW Declaration
\*---------------------------------------------------------------------------*/

class DAJacCondRdW
    : public DAJacCon
{

protected:
    /// table to specify how the states are connected to the residuals for a given solver
    HashTable<List<List<word>>> stateResConInfo_;

    /// a backup for DAJacCondRdW::adjstateResConInfo_ when reducing the connectivity levels
    HashTable<List<List<word>>> stateResConInfoBK_;

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

    /// number of dRdW colors
    label ndRdWColors_;

public:
    TypeName("DAJacCondRdW");
    // Constructors

    //- Construct from components
    DAJacCondRdW(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAJacCondRdW()
    {
    }

    // Member functions

    /// initialize petsc vectors
    void initializePetscVecs();

    /// reduce the connectivity level for DAJacCondRdW::dRdWCon_
    void reduceStateResConLevel();

    /// restore the connectivity level for DAJacCondRdW::dRdWCon_
    void restoreStateResConLevel();

    /// compute preallocation vectors
    void allocateJacobianConnections(
        Vec preallocOnProc,
        Vec preallocOffProc,
        Vec preallocOnProcT,
        Vec preallocOffProcT,
        Mat connections,
        const label row);

    /// compute preallocation vectors
    void preallocateJacobianMatrix(
        Mat dRMat,
        const Vec preallocOnProc,
        const Vec preallocOffProc);

    /// preallocate memory for dRdW using the computed preallocation vectors
    void preallocatedRdW(
        Mat dRMat,
        const label transposed);
    
     /// initialize the state Jacobian connectivity matrix
    virtual void initializeJacCon(const dictionary& options);

    /// assign 1 to all non-zero elements for the Jacobian connecitivyt matrix
    virtual void setupJacCon(const dictionary& options);
    
    /// get the number of dRdW colors
    virtual label getNJacConColors() const;

    /// return DAJacCondRdW::dRdWColoredColumns_
    virtual Vec getJacConColoredColumns() const
    {
        return dRdWColoredColumns_;
    }

    /// delete DAJacCondRdW::dRdWCon_
    virtual void deleteJacCon(const dictionary& options);

    /// compute graph coloring for Jacobian connectivity matrix
    virtual void calcJacConColoring();

    /// read colors for dRdW
    virtual void readJacConColoring();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
