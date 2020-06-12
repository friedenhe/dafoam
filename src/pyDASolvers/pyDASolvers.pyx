
# distutils: language = c++
# distutils: sources = DASolvers.C

"""
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Cython wrapper functions that call OpenFOAM libraries defined
        in the *.C and *.H files. The python naming convention is to
        add "py" before the C++ class name
"""

# for using Petsc
from petsc4py.PETSc cimport Vec, PetscVec

# declear cpp functions
cdef extern from "DASolvers.H" namespace "Foam":
    cppclass DASolvers:
        DASolvers(char *, object) except +
        void initSolver()
        int solvePrimal(PetscVec, PetscVec)
        int solveAdjoint(PetscVec, PetscVec, PetscVec)
        int calcTotalDerivs(PetscVec, PetscVec, PetscVec, PetscVec)
        int getGlobalXvIndex(int, int)
        void ofField2StateVec(PetscVec)
        void stateVec2OFField(PetscVec)
        int getNLocalAdjointStates()
        int checkMesh()
        double getObjFuncValue(char *)
        void printAllOptions()

# create python wrappers that call cpp functions
cdef class pyDASolvers:

    # define a class pointer for cpp functions
    cdef:
        DASolvers * _thisptr

    # initialize this class pointer with NULL
    def __cinit__(self):
        self._thisptr = NULL

    # deallocate the class pointer, and
    # make sure we don't have memory leak
    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    # point the class pointer to the cpp class constructor
    def __init__(self, argsAll, pyOptions):
        """
        Parameters
        ----------

        argsAll: char
            Chars that contains all the arguments
            for running OpenFOAM solvers, including
            the name of the solver.

        pyOptions: dict
            Dictionary that defines all the options
            in DAFoam

        Examples
        --------
        solver = pyDASolvers("DASolvers -parallel -python", aeroOptions)
        """
        self._thisptr = new DASolvers(argsAll, pyOptions)

    # wrap all the other memeber functions in the cpp class
    def initSolver(self):
        self._thisptr.initSolver()

    def solvePrimal(self, Vec xvVec, Vec wVec):
        return self._thisptr.solvePrimal(xvVec.vec, wVec.vec)
    
    def solveAdjoint(self, Vec xvVec, Vec wVec, Vec psiVec):
        self._thisptr.solveAdjoint(xvVec.vec, wVec.vec, psiVec.vec)
    
    def calcTotalDerivs(self, Vec xvVec, Vec wVec, Vec psiVec, Vec totalDerivVec):
        return self._thisptr.calcTotalDerivs(xvVec.vec, wVec.vec, psiVec.vec, totalDerivVec.vec)
    
    def getGlobalXvIndex(self, pointI, coordI):
        return self._thisptr.getGlobalXvIndex(pointI, coordI)
    
    def ofField2StateVec(self, Vec stateVec):
        self._thisptr.ofField2StateVec(stateVec.vec)
    
    def stateVec2OFField(self, Vec stateVec):
        self._thisptr.stateVec2OFField(stateVec.vec)
    
    def getNLocalAdjointStates(self):
        return self._thisptr.getNLocalAdjointStates()
    
    def checkMesh(self):
        return self._thisptr.checkMesh()
    
    def getObjFuncValue(self, objFuncName):
        return self._thisptr.getObjFuncValue(objFuncName)

    def printAllOptions(self):
        self._thisptr.printAllOptions()
