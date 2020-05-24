
# distutils: language = c++
# distutils: sources = DASolverCompressible.C

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
cdef extern from "DASolverCompressible.H" namespace "Foam":
    cppclass DASolverCompressible:
        DASolverCompressible(char *, object) except +
        void initSolver()
        void solvePrimal()
        void solveAdjoint()
        void calcTotalDerivs()
        int getGlobalXvIndex(int, int)
        void ofField2StateVec(PetscVec)
        void stateVec2OFField(PetscVec)
        int getNLocalAdjointStates()

# create python wrappers that call cpp functions
cdef class pyDASolverCompressible:

    # define a class pointer for cpp functions
    cdef:
        DASolverCompressible * _thisptr

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
        solver = pyDASolverCompressible("DASolverCompressible -parallel -python", aeroOptions)
        """
        self._thisptr = new DASolverCompressible(argsAll, pyOptions)

    # wrap all the other memeber functions in the cpp class
    def initSolver(self):
        self._thisptr.initSolver()

    def solvePrimal(self):
        self._thisptr.solvePrimal()
    
    def solveAdjoint(self):
        self._thisptr.solveAdjoint()
    
    def calcTotalDerivs(self):
        self._thisptr.calcTotalDerivs()
    
    def getGlobalXvIndex(self, pointI, coordI):
        return self._thisptr.getGlobalXvIndex(pointI, coordI)
    
    def ofField2StateVec(self, Vec stateVec):
        self._thisptr.ofField2StateVec(stateVec.vec)
    
    def stateVec2OFField(self, Vec stateVec):
        self._thisptr.stateVec2OFField(stateVec.vec)
    
    def getNLocalAdjointStates(self):
        return self._thisptr.getNLocalAdjointStates()
