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
#include "DACheckMesh.H"
#include "DAOption.H"
#include "DAStateInfo.H"
#include "DAModel.H"
#include "DATurbulenceModel.H"
#include "DAIndex.H"
#include "DAObjFunc.H"
#include "DAJacCon.H"
#include "DAColoring.H"
#include "DAResidual.H"
#include "DAField.H"
#include "DAPartDeriv.H"
#include "DALinearEqn.H"

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

    /// DAOption pointer
    autoPtr<DAOption> daOptionPtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// DAModel pointer
    autoPtr<DAModel> daModelPtr_;

    /// DAIndex pointer
    autoPtr<DAIndex> daIndexPtr_;

    /// DAField pointer
    autoPtr<DAField> daFieldPtr_;

    /// a list of DAObjFunc pointers
    UPtrList<DAObjFunc> daObjFuncPtrList_;

    /// DACheckMesh object pointer
    autoPtr<DACheckMesh> daCheckMeshPtr_;

    /// DAResidual pointer
    autoPtr<DAResidual> daResidualPtr_;

    /// DAStateInfo pointer
    autoPtr<DAStateInfo> daStateInfoPtr_;

    /// the stateInfo_ list from DAStateInfo object
    HashTable<wordList> stateInfo_;

    /// the list of objective function names that requires adjoint solution
    wordList objFuncNames4Adj_;

    /// the dictionary of adjoint vector (psi) values for all objectives
    dictionary psiVecDict_;

    dictionary totalDerivDict_;

    /// assign the psiVec to psiVecArray_
    void setPsiVecDict(
        const word objFuncName,
        const Vec psiVec,
        dictionary& psiVecDict);

    void setTotalDerivDict(
        const word objFuncName,
        const word designVarName,
        const Vec totalDerivVec,
        dictionary& totalDerivDict);

    Mat dXvdFFDMat_;

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
    label solveAdjoint(
        const Vec xvVec,
        const Vec wVec);

    /// compute the total derivatives
    label calcTotalDeriv(
        const Vec xvVec,
        const Vec wVec,
        const word designVarName);

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
        return daFieldPtr_->ofField2StateVec(stateVec);
    }

    /// assign the fields in OpenFOAM based on the state vector
    void stateVec2OFField(const Vec stateVec) const
    {
        return daFieldPtr_->stateVec2OFField(stateVec);
    }

    /// assign the points in fvMesh of OpenFOAM based on the point vector
    void pointVec2OFMesh(const Vec xvVec) const
    {
        return daFieldPtr_->pointVec2OFMesh(xvVec);
    }

    /// assign the point vector based on the points in fvMesh of OpenFOAM
    void ofMesh2PointVec(Vec xvVec) const
    {
        return daFieldPtr_->ofMesh2PointVec(xvVec);
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

    /// print all DAOption
    void printAllOptions()
    {
        Info << "DAFoam option dictionary: ";
        Info << daOptionPtr_->getAllOptions() << endl;
    }

    /// reduce the connectivity level for Jacobian connectivity mat
    void reduceStateResConLevel(
        const HashTable<label> maxResConLv4JacPCMat,
        HashTable<List<List<word>>>& stateResConInfo) const;

    void calcPrimalResidualStatistics(
        const word mode,
        const label writeRes = 0);

    /// assign the psiVecDict_ to psiVec
    void getPsiVec(
        const word objFuncName,
        Vec psiVec);

    scalar getTotalDerivVal(
        const word objFuncName,
        const word designVarName,
        const label idxI) const;

    void setdXvdFFDMat(const Mat dXvdFFDMat);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
