
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
from petsc4py.PETSc cimport Vec, PetscVec, Mat, PetscMat, KSP, PetscKSP

# declear cpp functions
cdef extern from "DASolvers.H" namespace "Foam":
    cppclass DASolvers:
        DASolvers(char *, object) except +
        void initSolver()
        int solvePrimal(PetscVec, PetscVec)
        int solveAdjoint(PetscVec, PetscVec)
        int calcTotalDeriv(PetscVec, PetscVec, char *)
        void calcdRdWT(PetscVec, PetscVec, int, PetscMat)
        void calcdFdW(PetscVec, PetscVec, char *, PetscVec)
        void createMLRKSP(PetscMat, PetscMat, PetscKSP)
        void solveLinearEqn(PetscKSP, PetscVec, PetscVec)
        void calcdRdBC(PetscVec, PetscVec, char *, PetscMat)
        void calcdFdBC(PetscVec, PetscVec, char *, char *, PetscVec)
        void calcdRdAOA(PetscVec, PetscVec, char *, PetscMat)
        void calcdFdAOA(PetscVec, PetscVec, char *, char *, PetscVec)
        void calcdRdFFD(PetscVec, PetscVec, char *, PetscMat)
        void calcdFdFFD(PetscVec, PetscVec, char *, char *, PetscVec)
        void calcdRdACT(PetscVec, PetscVec, char *, char *, PetscMat)
        void multiPointTreatment(PetscVec)
        void setdXvdFFDMat(PetscMat)
        int getGlobalXvIndex(int, int)
        void ofField2StateVec(PetscVec)
        void stateVec2OFField(PetscVec)
        int getNLocalAdjointStates()
        int checkMesh()
        double getObjFuncValue(char *)
        void printAllOptions()
        double getTotalDerivVal(char *, char *, int)
        void updateDAOption(object)
        double getPrevPrimalSolTime()
        # functions for unit tests
        void pointVec2OFMesh(PetscVec)
        void ofMesh2PointVec(PetscVec)
        void resVec2OFResField(PetscVec)
        void ofResField2ResVec(PetscVec)
        void writeMatrixBinary(PetscMat, char *)
        void writeMatrixASCII(PetscMat, char *)
        void readMatrixBinary(PetscMat, char *)
        void writeVectorASCII(PetscVec, char *)
        void readVectorBinary(PetscVec, char *)
        void writeVectorBinary(PetscVec, char *)
        void setTimeInstanceField(int)
        double getTimeInstanceObjFunc(int, char *)
    
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
    
    def solveAdjoint(self, Vec xvVec, Vec wVec):
        return self._thisptr.solveAdjoint(xvVec.vec, wVec.vec)
    
    def calcTotalDeriv(self, Vec xvVec, Vec wVec, designVarName):
        return self._thisptr.calcTotalDeriv(xvVec.vec, wVec.vec, designVarName)
    
    def calcdRdWT(self, Vec xvVec, Vec wVec, isPC, Mat dRdWT):
        self._thisptr.calcdRdWT(xvVec.vec, wVec.vec, isPC, dRdWT.mat)
    
    def calcdFdW(self, Vec xvVec, Vec wVec, objFuncName, Vec dFdW):
        self._thisptr.calcdFdW(xvVec.vec, wVec.vec, objFuncName, dFdW.vec)
    
    def createMLRKSP(self, Mat jacMat, Mat jacPCMat, KSP myKSP):
        self._thisptr.createMLRKSP(jacMat.mat, jacPCMat.mat, myKSP.ksp)
    
    def solveLinearEqn(self, KSP myKSP, Vec rhsVec, Vec solVec):
        self._thisptr.solveLinearEqn(myKSP.ksp, rhsVec.vec, solVec.vec)

    def calcdRdBC(self, Vec xvVec, Vec wVec, designVarName, Mat dRdBC):
        self._thisptr.calcdRdBC(xvVec.vec, wVec.vec, designVarName, dRdBC.mat)
    
    def calcdFdBC(self, Vec xvVec, Vec wVec, objFuncName, designVarName, Vec dFdBC):
        self._thisptr.calcdFdBC(xvVec.vec, wVec.vec, objFuncName, designVarName, dFdBC.vec)

    def calcdRdAOA(self, Vec xvVec, Vec wVec, designVarName, Mat dRdAOA):
        self._thisptr.calcdRdAOA(xvVec.vec, wVec.vec, designVarName, dRdAOA.mat)

    def calcdFdAOA(self, Vec xvVec, Vec wVec, objFuncName, designVarName, Vec dFdAOA):
        self._thisptr.calcdFdAOA(xvVec.vec, wVec.vec, objFuncName, designVarName, dFdAOA.vec)

    def calcdRdFFD(self, Vec xvVec, Vec wVec, designVarName, Mat dRdFFD):
        self._thisptr.calcdRdFFD(xvVec.vec, wVec.vec, designVarName, dRdFFD.mat)

    def calcdFdFFD(self, Vec xvVec, Vec wVec, objFuncName, designVarName, Vec dFdFFD):
        self._thisptr.calcdFdFFD(xvVec.vec, wVec.vec, objFuncName, designVarName, dFdFFD.vec)

    def calcdRdACT(self, Vec xvVec, Vec wVec, designVarName, designVarType, Mat dRdACT):
        self._thisptr.calcdRdACT(xvVec.vec, wVec.vec, designVarName, designVarType, dRdACT.mat)
    
    def multiPointTreatment(self, Vec wVec):
        self._thisptr.multiPointTreatment(wVec.vec)

    def setdXvdFFDMat(self, Mat dXvdFFDMat):
        self._thisptr.setdXvdFFDMat(dXvdFFDMat.mat)
    
    def getGlobalXvIndex(self, pointI, coordI):
        return self._thisptr.getGlobalXvIndex(pointI, coordI)
    
    def ofField2StateVec(self, Vec stateVec):
        self._thisptr.ofField2StateVec(stateVec.vec)
    
    def stateVec2OFField(self, Vec stateVec):
        self._thisptr.stateVec2OFField(stateVec.vec)
    
    def pointVec2OFMesh(self, Vec xvVec):
        self._thisptr.pointVec2OFMesh(xvVec.vec)
    
    def ofMesh2PointVec(self, Vec xvVec):
        self._thisptr.ofMesh2PointVec(xvVec.vec)
    
    def resVec2OFResField(self, Vec rVec):
        self._thisptr.resVec2OFResField(rVec.vec)
    
    def ofResField2ResVec(self, Vec rVec):
        self._thisptr.ofResField2ResVec(rVec.vec)
    
    def getNLocalAdjointStates(self):
        return self._thisptr.getNLocalAdjointStates()
    
    def checkMesh(self):
        return self._thisptr.checkMesh()
    
    def getObjFuncValue(self, objFuncName):
        return self._thisptr.getObjFuncValue(objFuncName)
    
    def getTotalDerivVal(self, objFuncName, designVarName, idxI):
        return self._thisptr.getTotalDerivVal(objFuncName, designVarName, idxI)

    def printAllOptions(self):
        self._thisptr.printAllOptions()

    def updateDAOption(self, pyOptions):
        self._thisptr.updateDAOption(pyOptions)
    
    def getPrevPrimalSolTime(self):
        return self._thisptr.getPrevPrimalSolTime()
    
    def writeMatrixBinary(self, Mat magIn, prefix):
        self._thisptr.writeMatrixBinary(magIn.mat, prefix)
    
    def writeMatrixASCII(self, Mat magIn, prefix):
        self._thisptr.writeMatrixASCII(magIn.mat, prefix)
    
    def readMatrixBinary(self, Mat magIn, prefix):
        self._thisptr.readMatrixBinary(magIn.mat, prefix)
    
    def writeVectorASCII(self, Vec vecIn, prefix):
        self._thisptr.writeVectorASCII(vecIn.vec, prefix)
    
    def readVectorBinary(self, Vec vecIn, prefix):
        self._thisptr.readVectorBinary(vecIn.vec, prefix)
    
    def writeVectorBinary(self, Vec vecIn, prefix):
        self._thisptr.writeVectorBinary(vecIn.vec, prefix)
    
    def setTimeInstanceField(self, instanceI):
        self._thisptr.setTimeInstanceField(instanceI)
    
    def getTimeInstanceObjFunc(self, instanceI, objFuncName):
        return self._thisptr.getTimeInstanceObjFunc(instanceI, objFuncName)
