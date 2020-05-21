
# distutils: language = c++
# distutils: sources = DASolverIncompressible.C

"""
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Cython wrapper functions that call OpenFOAM libraries defined
        in the *.C and *.H files. The python naming convention is to
        add "py" before the C++ class name
"""

# for using Petsc
# from petsc4py.PETSc cimport Vec, PetscVec

# declear cpp functions
cdef extern from "DASolverIncompressible.H" namespace "Foam":
    cppclass DASolverIncompressible:
        DASolverIncompressible(char *, object) except +
        void initSolver()
        void solvePrimal()

# create python wrappers that call cpp functions
cdef class pyDASolverIncompressible:

    # define a class pointer for cpp functions
    cdef:
        DASolverIncompressible * _thisptr

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
        solver = pyDASolverIncompressible("DASolverIncompressible -parallel -python", aeroOptions)
        """
        self._thisptr = new DASolverIncompressible(argsAll, pyOptions)

    # wrap all the other memeber functions in the cpp class
    def initSolver(self):
        self._thisptr.initSolver()

    def solvePrimal(self):
        self._thisptr.solvePrimal()
