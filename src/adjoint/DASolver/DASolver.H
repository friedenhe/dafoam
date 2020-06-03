/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Solver class that calls primal and adjoint solvers,
        and compute the total derivatives

\*---------------------------------------------------------------------------*/

#ifndef DASolver_H
#define DASolver_H

#include <petscksp.h>
#include "Python.h"
#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "DAUtility.H"
#include "DAOption.H"
#include "DARegState.H"
#include "DAModel.H"
#include "DATurbulenceModel.H"
#include "DAIndex.H"
#include "DAObjFunc.H"
#include "DACheckMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DASolver Declaration
\*---------------------------------------------------------------------------*/

class DASolver
{

private:
    /// Disallow default bitwise copy construct
    DASolver(const DASolver&);

    /// Disallow default bitwise assignment
    void operator=(const DASolver&);

protected:
    /// all the arguments
    char* argsAll_;

    /// all options in DAFoam
    PyObject* pyOptions_;

    /// args pointer
    autoPtr<argList> argsPtr_;

    /// runTime pointer
    autoPtr<Time> runTimePtr_;

    /// fvMesh pointer
    autoPtr<fvMesh> meshPtr_;

    /// DAUtility pointer
    autoPtr<DAUtility> daUtilPtr_;

    /// DAOption pointer
    autoPtr<DAOption> daOptionPtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// DAModel pointer
    autoPtr<DAModel> daModelPtr_;

    /// DARegState pointer
    autoPtr<DARegState> daRegStatePtr_;

    /// DAIndex pointer
    autoPtr<DAIndex> daIndexPtr_;

    /// a list of DAObjFunc pointers
    UPtrList<DAObjFunc> daObjFuncPtrList_;

    /// DACheckMesh object pointer
    autoPtr<DACheckMesh> daCheckMeshPtr_;

public:
    /// Runtime type information
    TypeName("DASolver");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DASolver,
        dictionary,
        (
            char* argsAll,
            PyObject* pyOptions),
        (argsAll, pyOptions));

    // Constructors

    //- Construct from components
    DASolver(
        char* argsAll,
        PyObject* pyOptions);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<DASolver> New(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DASolver()
    {
    }

    // Member functions

    /// initialize fields and variables
    virtual void initSolver() = 0;

    /// solve the primal equations
    virtual label solvePrimal(
        const Vec xvVec,
        Vec wVec) = 0;

    /// solve the adjoint linear equations
    virtual void solveAdjoint() = 0;

    /// compute the total derivatives
    virtual void calcTotalDerivs() = 0;

    /// return whether to loop the primal solution, similar to runTime::loop() except we don't do file IO
    label loop(Time& runTime);

    /// basically, we call DAIndex::getGlobalXvIndex
    label getGlobalXvIndex(
        const label idxPoint,
        const label idxCoord) const
    {
        return daIndexPtr_->getGlobalXvIndex(idxPoint, idxCoord);
    }

    /// set the state vector based on the latest fields in OpenFOAM
    void ofField2StateVec(Vec stateVec) const
    {
        return daIndexPtr_->ofField2StateVec(stateVec);
    }

    /// assign the fields in OpenFOAM based on the state vector
    void stateVec2OFField(const Vec stateVec) const
    {
        return daIndexPtr_->stateVec2OFField(stateVec);
    }

    /// assign the points in fvMesh of OpenFOAM based on the point vector
    void pointVec2OFMesh(const Vec xvVec) const
    {
        return daIndexPtr_->pointVec2OFMesh(xvVec);
    }

    /// assign the point vector based on the points in fvMesh of OpenFOAM
    void ofMesh2PointVec(Vec xvVec) const
    {
        return daIndexPtr_->ofMesh2PointVec(xvVec);
    }

    /// return the number of local adjoint states
    label getNLocalAdjointStates() const
    {
        return daIndexPtr_->nLocalAdjointStates;
    }

    /// initialize DASolver::daObjFuncPtrList_ one needs to call this before calling printAllObjFuncs
    void setDAObjFuncList();

    /// calculate the values of all objective functions and print them to screen
    void printAllObjFuncs();

    /// check the mesh quality and return meshOK
    label checkMesh() const
    {
        return daCheckMeshPtr_->run();
    }

    /// return the value of the objective function
    scalar getObjFuncValue(const word objFuncName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
