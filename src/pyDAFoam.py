#!/usr/bin/env python

"""

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    The Python interface to DAFoam. It controls the adjoint
    solvers and external modules for design optimization

"""

import subprocess
import os
import sys
import copy
import shutil
import numpy as np
from pprint import pprint as pp
from mpi4py import MPI
from collections import OrderedDict
import petsc4py
from petsc4py import PETSc

petsc4py.init(sys.argv)


class PYDAFOAM(object):

    """
    Main class for pyDAFoam

    Parameters
    ----------

    comm : mpi4py communicator
        An optional argument to pass in an external communicator.

    options : dictionary
        The list of options to use with pyDAFoam.

    """

    def _getDefOptions(self):
        """
        Setup default options

        Returns
        -------

        defOpts : dict
            All the DAFoam options.
        """
        defOpts = {
            # primal options
            "primalEndTime": [float, 1.0],
            "primalDeltaT": [float, 1.0],
            "primalVarBounds": [dict, {}],
            "flowCondition": [str, "Incompressible"],
            "turbulenceModel": [str, "SpalartAllmaras"],
            "primalBC": [dict, {}],
            "fvSource": [dict, {}],
            "printInterval": [int, 100],
            "primalMinResTol": [float, 1.0e-8],
            "primalMinResTolDiff": [float, 1.0e2],
            # adjoint options
            "adjUseColoring": [bool, True],
            "adjEpsDerivState": [float, 1.0e-5],
            "adjEpsDerivFFD": [float, 1.0e-3],
            "adjEpsDerivBC": [float, 1.0e-2],
            "adjEpsDerivAOA": [float, 1.0e-3],
            "adjStateOrdering": [str, "state"],
            "adjEqnOption": [
                dict,
                {
                    "globalPCIters": 0,
                    "asmOverlap": 1,
                    "localPCIters": 1,
                    "jacMatReOrdering": "rcm",
                    "pcFillLevel": 1,
                    "gmresMaxIters": 1000,
                    "gmresRestart": 1000,
                    "gmresRelTol": 1.0e-6,
                    "gmresAbsTol": 1.0e-14,
                    "gmresTolDiff": 1.0e2,
                },
            ],
            "normalizeStates": [dict, {}],
            "normalizeResiduals": [list, ["URes", "pRes", "nuTildaRes", "phiRes", "TRes"]],
            "maxResConLv4JacPCMat": [
                dict,
                {
                    "pRes": 2,
                    "phiRes": 1,
                    "URes": 2,
                    "TRes": 2,
                    "nuTildaRes": 2,
                    "kRes": 2,
                    "epsilonRes": 2,
                    "omegaRes": 2,
                    "p_rghRes": 2,
                },
            ],
            "transonicPCOption": [int, -1],
            # optimization options
            "designVar": [dict, {}],
            # system options
            "rootDir": [str, "./"],
            "solverName": [str, "DASimpleFoam"],
            "printAllOptions": [bool, True],
            "objFunc": [dict, {}],
            "debug": [bool, False],
            # surface definition
            "meshSurfaceFamily": [str, "None"],
            "designSurfaceFamily": [str, "None"],
            "designSurfaces": [list, ["body"]],
        }

        return defOpts

    def __init__(self, comm=None, options=None):
        """
        Initialize class members
        """

        assert not os.getenv("WM_PROJECT") is None, "$WM_PROJECT not found. Please source OpenFOAM-v1812/etc/bashrc"

        # name
        self.name = "PYDAFOAM"

        # initialize options for adjoints
        self._initializeOptions(options)

        # initialize comm for parallel communication
        self._initializeComm(comm)

        # Initialize families
        self.families = OrderedDict()

        # Default it to fault, after calling setSurfaceCoordinates, set it to true
        self._updateGeomInfo = False

        # Use double data type: 'd'
        self.dtype = "d"

        # write all the setup files
        self._writeOFCaseFiles()

        # Remind the user of all the DAFoam options:
        if self.getOption("printAllOptions"):
            self._printCurrentOptions()

        # run decomposePar for parallel runs
        self.runDecomposePar()

        # initialize the pySolvers
        self.solverInitialized = 0
        self._initSolver()

        # initialize mesh information and read grids
        self._readMeshInfo()

        # initialize the mesh point vector xvVec
        self._initializeMeshPointVec()

        # initialize state variable vector self.wVec
        self._initializeStateVec()

        # get the reduced point connectivities for the base patches in the mesh
        self._computeBasicFamilyInfo()

        # Add a couple of special families.
        self.allFamilies = "allSurfaces"
        self.addFamilyGroup(self.allFamilies, self.basicFamilies)

        self.allWallsGroup = "allWalls"
        self.addFamilyGroup(self.allWallsGroup, self.wallList)

        # Set the design families if given, otherwise default to all
        # walls
        self.designFamilyGroup = self.getOption("designSurfaceFamily")
        if self.designFamilyGroup == "None":
            self.designFamilyGroup = self.allWallsGroup

        # Set the mesh families if given, otherwise default to all
        # walls
        self.meshFamilyGroup = self.getOption("meshSurfaceFamily")
        if self.meshFamilyGroup == "None":
            self.meshFamilyGroup = self.allWallsGroup

        # get the surface coordinate of allFamilies
        self.xs0 = self.getSurfaceCoordinates(self.allFamilies)

        # By Default we don't have an external mesh object or a
        # geometric manipulation object
        self.mesh = None
        self.DVGeo = None

        # initialize the number of primal and adjoint calls
        self.nSolvePrimals = 0
        self.nSolveAdjoints = 0

        self.primalFail = 0
        self.adjointFail = 0

        if self.comm.rank == 0:
            print("pyDAFoam initialization done!")

        return

    def __del__(self):

        self.primalSolver = None

        return

    def __call__(self):
        """
        Solve the primal
        """

        # update the mesh coordinates if DVGeo is set
        # add point set and update the mesh based on the DV values
        self.ptSetName = self.getPointSetName("dummy")
        ptSetName = self.ptSetName
        if self.DVGeo is not None:

            # if the point set is not in DVGeo add it first
            if ptSetName not in self.DVGeo.points:

                xs0 = self.mapVector(self.xs0, self.allFamilies, self.designFamilyGroup)

                self.DVGeo.addPointSet(xs0, self.ptSetName)
                self.pointsSet = True

            # set the surface coords xs
            if self.comm.rank == 0:
                print("DVGeo PointSet UpToDate: " + str(self.DVGeo.pointSetUpToDate(ptSetName)))
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                if self.comm.rank == 0:
                    print("Updating DVGeo PointSet....")
                xs = self.DVGeo.update(ptSetName, config=None)
                self.setSurfaceCoordinates(xs, self.designFamilyGroup)
                if self.comm.rank == 0:
                    print("DVGeo PointSet UpToDate: " + str(self.DVGeo.pointSetUpToDate(ptSetName)))

                # warp the mesh to get the new volume coordinates
                if self.comm.rank == 0:
                    print("Warping the volume mesh....")
                self.mesh.warpMesh()

                xvNew = self.mesh.getSolverGrid()
                self.xvFlatten2XvVec(xvNew, self.xvVec)

        # solve the primal to get new state variables
        self.solvePrimal()

        return

    def evalFunctions(self, funcs, evalFuncs=None, ignoreMissing=False):
        """
        Evaluate the desired functions given in iterable object,
        'evalFuncs' and add them to the dictionary 'funcs'. The keys
        in the funcs dictioary will be have an _<ap.name> appended to
        them. Additionally, information regarding whether or not the
        last analysis with the solvePrimal was sucessful is
        included. This information is included as "funcs['fail']". If
        the 'fail' entry already exits in the dictionary the following
        operation is performed:

        funcs['fail'] = funcs['fail'] or <did this problem fail>

        In other words, if any one problem fails, the funcs['fail']
        entry will be False. This information can then be used
        directly in the pyOptSparse.

        Parameters
        ----------
        funcs : dict
            Dictionary into which the functions are saved.

        evalFuncs : iterable object containing strings
          If not None, use these functions to evaluate.

        ignoreMissing : bool
            Flag to supress checking for a valid function. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> CFDsolver()
        >>> CFDsolver.evalFunctions(funcs, ['CD', 'CL'])
        >>> funcs
        >>> # Result will look like:
        >>> # {'CD':0.501, 'CL':0.02750}
        """

        for funcName in evalFuncs:
            objFuncValue = self.solver.getObjFuncValue(funcName.encode())
            funcs[funcName] = objFuncValue

        if self.primalFail:
            funcs["fail"] = True
        else:
            funcs["fail"] = False

        return

    def evalFunctionsSens(self, funcsSens, evalFuncs=None):
        """
        Evaluate the sensitivity of the desired functions given in
        iterable object,'evalFuncs' and add them to the dictionary
        'funcSens'.

        Parameters
        ----------
        funcSens : dict
        Dictionary into which the function derivatives are saved.

        evalFuncs : iterable object containing strings
            The functions the user wants the derivatives of

        Examples
        --------
        >>> funcSens = {}
        >>> CFDsolver.evalFunctionsSens(funcSens, ['CD', 'CL'])
        """

        if self.DVGeo is None:
            raise Error("DVGeo not set!")

        dvs = self.DVGeo.getValues()

        for funcName in evalFuncs:
            funcsSens[funcName] = {}
            for dvName in dvs:
                nDVs = len(dvs[dvName])
                funcsSens[funcName][dvName] = np.zeros(nDVs, self.dtype)
                for i in range(nDVs):
                    sensVal = self.solver.getTotalDerivVal(funcName.encode(), dvName.encode(), i)
                    funcsSens[funcName][dvName][i] = sensVal

        if self.adjointFail:
            funcsSens["fail"] = True
        else:
            funcsSens["fail"] = False

        return

    def setDVGeo(self, DVGeo):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that <SOLVER> does not **strictly** need a
        DVGeometry object, but if optimization with geometric
        changes is desired, then it is required.
        Parameters
        ----------
        dvGeo : A DVGeometry object.
            Object responsible for manipulating the constraints that
            this object is responsible for.
        Examples
        --------
        >>> CFDsolver = <SOLVER>(comm=comm, options=CFDoptions)
        >>> CFDsolver.setDVGeo(DVGeo)
        """

        self.DVGeo = DVGeo

    def addFamilyGroup(self, groupName, families):
        """
        Add a custom grouping of families called groupName. The groupName
        must be distinct from the existing families. All families must
        in the 'families' list must be present in the CGNS file.
        Parameters
        ----------
        groupName : str
            User-supplied custom name for the family groupings
        families : list
            List of string. Family names to combine into the family group
        """

        # Do some error checking
        if groupName in self.families:
            raise Error(
                "The specified groupName '%s' already exists in the " "cgns file or has already been added." % groupName
            )

        # We can actually allow for nested groups. That is, an entry
        # in families may already be a group added in a previous call.
        indices = []
        for fam in families:
            if fam not in self.families:
                raise Error(
                    "The specified family '%s' for group '%s', does "
                    "not exist in the cgns file or has "
                    "not already been added. The current list of "
                    "families (original and grouped) is: %s" % (fam, groupName, repr(self.families.keys()))
                )

            indices.extend(self.families[fam])

        # It is very important that the list of families is sorted
        # becuase in fortran we always use a binary search to check if
        # a famID is in the list.
        self.families[groupName] = sorted(np.unique(indices))

    def setMesh(self, mesh):
        """
        Set the mesh object to the aero_solver to do geometric deformations
        Parameters
        ----------
        mesh : MBMesh or USMesh object
            The mesh object for doing the warping
        """

        # Store a reference to the mesh
        self.mesh = mesh

        # Setup External Warping with volume indices
        meshInd = self.getSolverMeshIndices()
        self.mesh.setExternalMeshIndices(meshInd)

        # Set the surface the user has supplied:
        conn, faceSizes = self.getSurfaceConnectivity(self.meshFamilyGroup)
        pts = self.getSurfaceCoordinates(self.meshFamilyGroup)
        self.mesh.setSurfaceDefinition(pts, conn, faceSizes)

    def getSurfaceConnectivity(self, groupName=None):
        """
        Return the connectivity of the coordinates at which the forces (or tractions) are
        defined. This is the complement of getForces() which returns
        the forces at the locations returned in this routine.

        Parameters
        ----------
        groupName : str
            Group identifier to get only forces cooresponding to the
            desired group. The group must be a family or a user-supplied
            group of families. The default is None which corresponds to
            all wall-type surfaces.
        """

        if groupName is None:
            groupName = self.allWallsGroup

        # loop over the families in this group and populate the connectivity
        famInd = self.families[groupName]
        conn = []
        faceSizes = []

        pointOffset = 0
        for Ind in famInd:
            # select the face from the basic families
            name = self.basicFamilies[Ind]

            # get the size of this
            bc = self.boundaries[name]
            nPts = len(bc["indicesRed"])

            # get the number of reduced faces associated with this boundary
            nFace = len(bc["facesRed"])

            # check that this isn't an empty boundary
            if nFace > 0:
                # loop over the faces and add them to the connectivity and faceSizes array
                for iFace in range(nFace):
                    face = copy.copy(bc["facesRed"][iFace])
                    for i in range(len(face)):
                        face[i] += pointOffset
                    conn.extend(face)
                    faceSizes.append(len(face))

                pointOffset += nPts

        return conn, faceSizes

    def getTriangulatedMeshSurface(self, groupName=None, **kwargs):
        """
        This function returns a trianguled verision of the surface
        mesh on all processors. The intent is to use this for doing
        constraints in DVConstraints.
        Returns
        -------
        surf : list
           List of points and vectors describing the surface. This may
           be passed directly to DVConstraint setSurface() function.
        """

        if groupName is None:
            groupName = self.allWallsGroup

        # Obtain the points and connectivity for the specified
        # groupName
        pts = self.comm.allgather(self.getSurfaceCoordinates(groupName, **kwargs))
        conn, faceSizes = self.getSurfaceConnectivity(groupName)
        conn = np.array(conn).flatten()
        conn = self.comm.allgather(conn)
        faceSizes = self.comm.allgather(faceSizes)

        # Triangle info...point and two vectors
        p0 = []
        v1 = []
        v2 = []

        # loop over the faces
        for iProc in range(len(faceSizes)):

            connCounter = 0
            for iFace in range(len(faceSizes[iProc])):
                # Get the number of nodes on this face
                faceSize = faceSizes[iProc][iFace]
                faceNodes = conn[iProc][connCounter : connCounter + faceSize]

                # Start by getting the centerpoint
                ptSum = [0, 0, 0]
                for i in range(faceSize):
                    # idx = ptCounter+i
                    idx = faceNodes[i]
                    ptSum += pts[iProc][idx]

                avgPt = ptSum / faceSize

                # Now go around the face and add a triangle for each adjacent pair
                # of points. This assumes an ordered connectivity from the
                # meshwarping
                for i in range(faceSize):
                    idx = faceNodes[i]
                    p0.append(avgPt)
                    v1.append(pts[iProc][idx] - avgPt)
                    if i < (faceSize - 1):
                        idxp1 = faceNodes[i + 1]
                        v2.append(pts[iProc][idxp1] - avgPt)
                    else:
                        # wrap back to the first point for the last element
                        idx0 = faceNodes[0]
                        v2.append(pts[iProc][idx0] - avgPt)

                # Now increment the connectivity
                connCounter += faceSize

        return [p0, v1, v2]

    def printFamilyList(self):
        """
        Print a nicely formatted dictionary of the family names
        """
        pp(self.families)

    def setDesignVars(self, x):
        """
        Set the internal design variables.
        At the moment we don't have any internal DVs to set.
        """
        pass

        return

    def _initializeOptions(self, options):
        """
        Initialize the options passed into pyDAFoam

        Parameters
        ----------

        options : dictionary
            The list of options to use with pyDAFoam.
        """

        # If 'options' is None raise an error
        if options is None:
            raise Error("The 'options' keyword argument must be passed pyDAFoam.")

        # set immutable options that users should not change during the optimization
        self.imOptions = self._getImmutableOptions()

        # Load all the option information:
        self.defaultOptions = self._getDefOptions()

        # Set options based on defaultOptions
        # we basically overwrite defaultOptions with the given options
        # first assign self.defaultOptions to self.options
        self.options = OrderedDict()
        for key in self.defaultOptions:
            if len(self.defaultOptions[key]) != 2:
                raise Error(
                    "key %s has wrong format! \
                    Example: {'iters' : [int, 1]}"
                    % key
                )
            self.options[key] = self.defaultOptions[key]
        # now set options to self.options
        for key in options:
            self.setOption(key, options[key])

        return

    def _initializeComm(self, comm):
        """
        Initialize MPI COMM and setup parallel flags
        """

        # Set the MPI Communicators and associated info
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Check whether we are running in parallel
        nProc = self.comm.size
        self.parallel = False
        if nProc > 1:
            self.parallel = True

        # Save the rank and number of processors
        self.rank = self.comm.rank
        self.nProcs = self.comm.size

        # Setup the parallel flag for OpenFOAM executives
        self.parallelFlag = ""
        if self.parallel:
            self.parallelFlag = "-parallel"

        return

    def _writeOFCaseFiles(self):

        return

    def solvePrimal(self):
        """
        Run primal solver to compute state variables and objectives

        Input:
        ------
        xvVec: vector that contains all the mesh point coordinates

        Output:
        -------
        wVec: vector that contains all the state variables

        self.primalFail: if the primal solution fails, assigns 1, otherwise 0
        """

        if self.comm.rank == 0:
            print("Running Primal Solver %03d" % self.nSolvePrimals)
        
        self.deletePrevPrimalSolTime()

        self.primalFail = 0
        self.primalFail = self.solver.solvePrimal(self.xvVec, self.wVec)

        self.nSolvePrimals += 1

        return

    def solveAdjoint(self):
        """
        Run adjoint solver to compute the adjoint vector psiVec

        Input:
        ------
        xvVec: vector that contains all the mesh point coordinates

        wVec: vector that contains all the state variables

        Output:
        -------
        psiVec: the adjoint vector

        self.adjointFail: if the primal solution fails, assigns 1, otherwise 0
        """

        # save the point vector and state vector to disk
        """
        if self.comm.rank == 0:
            print("Saving the xvVec and wVec vectors to disk....")
        self.comm.Barrier()
        viewerXv = PETSc.Viewer().createBinary("xvVec_%03d.bin" % self.nSolveAdjoints, mode="w", comm=PETSc.COMM_WORLD)
        viewerXv(self.xvVec)
        viewerW = PETSc.Viewer().createBinary("wVec_%03d.bin" % self.nSolveAdjoints, mode="w", comm=PETSc.COMM_WORLD)
        viewerW(self.wVec)
        """

        self.renameSolution(self.nSolveAdjoints)

        if self.comm.rank == 0:
            print("Running adjoint Solver %03d" % self.nSolveAdjoints)

        self.adjointFail = 0
        self.adjointFail = self.solver.solveAdjoint(self.xvVec, self.wVec)

        self.nSolveAdjoints += 1

        return

    def calcTotalDeriv(self):
        """
        Compute total derivative

        Input:
        ------
        xvVec: vector that contains all the mesh point coordinates

        wVec: vector that contains all the state variables

        psiVec: the adjoint vector

        Output:
        -------
        totalDerivVec: the total derivative vector

        self.adjointFail: if the total derivative computation fails, assigns 1, otherwise 0
        """

        if self.comm.rank == 0:
            print("Computing total derivatives....")

        designVarDict = self.getOption("designVar")
        for key in designVarDict:
            if designVarDict[key]["designVarType"] == "FFD":
                self.setdXvdFFDMat(key)
            self.solver.calcTotalDeriv(self.xvVec, self.wVec, key.encode())

        return

    def _initSolver(self):
        """
        Initialize the solvers. This needs to be called before calling any runs
        """

        if self.solverInitialized == 1:
            raise Error("pyDAFoam: self._initSolver has been called! One shouldn't initialize solvers twice!")

        solverName = self.getOption("solverName")
        solverArg = solverName + " -python " + self.parallelFlag
        if self.getOption("flowCondition") == "Incompressible":

            from .pyDASolverIncompressible import pyDASolvers

            self.solver = pyDASolvers(solverArg.encode(), self.options)
        elif self.getOption("flowCondition") == "Compressible":

            from .pyDASolverCompressible import pyDASolvers

            self.solver = pyDASolvers(solverArg.encode(), self.options)
        else:
            raise Error("pyDAFoam: flowCondition %s: not valid!" % self.getOption("flowCondition"))

        self.solver.initSolver()

        self.solver.printAllOptions()

        self.solverInitialized = 1

        return

    def runCheckMesh(self):
        """
        Run checkMesh for mesh quality

        Returns
        -------

        meshOK : int
            meshOK=1 means the mesh quality check passes
        """

        if self.comm.rank == 0:
            print("\n")
            print("+--------------------------------------------------------------------------+")
            print("|                        Checking Mesh Quality                             |")
            print("+--------------------------------------------------------------------------+")

        meshOK = self.solver.checkMesh()

        return meshOK

    def runColoring(self):
        """
        Run coloring solver
        """

        if self.comm.rank == 0:
            print("\n")
            print("+--------------------------------------------------------------------------+")
            print("|                       Running Coloring Solver                            |")
            print("+--------------------------------------------------------------------------+")

        if self.getOption("flowCondition") == "Incompressible":

            from .pyColoringIncompressible import pyColoringIncompressible

            solverArg = "ColoringIncompressible -python " + self.parallelFlag
            solver = pyColoringIncompressible(solverArg.encode(), self.options)
        elif self.getOption("flowCondition") == "Compressible":

            from .pyColoringCompressible import pyColoringCompressible

            solverArg = "ColoringCompressible -python " + self.parallelFlag
            solver = pyColoringCompressible(solverArg.encode(), self.options)
        else:
            raise Error("pyDAFoam: flowCondition %s: not valid!" % self.getOption("flowCondition"))
        solver.run()

        solver = None

        return

    def runDecomposePar(self):
        """
        Run decomposePar to parallel run
        """

        # don't run it if it is a serial case
        if self.comm.size == 1:
            return

        if self.comm.rank == 0:
            status = subprocess.call("decomposePar", stdout=sys.stdout, stderr=subprocess.STDOUT, shell=False)
            if status != 0:
                # raise Error('pyDAFoam: status %d: Unable to run decomposePar'%status)
                print("\nUnable to run decomposePar, the domain has been already decomposed?\n")
        self.comm.Barrier()

        return

    def deletePrevPrimalSolTime(self):
        """
        Delete the previous primal solution time folder
        """

        solTime = self.solver.getPrevPrimalSolTime()

        rootDir = os.getcwd()
        if self.parallel:
            checkPath = os.path.join(rootDir, "processor%d/%g" % (self.comm.rank, solTime))
        else:
            checkPath = os.path.join(rootDir, "%g" % solTime)

        if os.path.isdir(checkPath):
            try:
                shutil.rmtree(checkPath)
            except Exception:
                raise Error("Can not delete %s" % checkPath)

            if self.comm.rank == 0:
                print("Previous solution time %g found and deleted." % solTime)
        else:
            if self.comm.rank == 0:
                print("Previous solution time %g not found and nothing deleted." % solTime)

        return

    def renameSolution(self, solIndex):
        """
        Rename the primal solution folder to specific format for post-processing. The renamed time has the
        format like 0.00000001, 0.00000002, etc. One can load these intermediate shapes and fields and
        plot them in paraview.
        The way it is implemented is that we sort the solution folder and consider the largest time folder
        as the solution folder and rename it

        Parameters
        ----------
        solIndex: int
            The major interation index
        """

        allSolutions = []
        rootDir = os.getcwd()
        if self.parallel:
            checkPath = os.path.join(rootDir, "processor%d" % self.comm.rank)
        else:
            checkPath = rootDir

        folderNames = os.listdir(checkPath)
        for folderName in folderNames:
            try:
                float(folderName)
                allSolutions.append(folderName)
            except ValueError:
                continue
        allSolutions.sort(reverse=True)
        # choose the latst solution to rename
        solutionTime = allSolutions[0]

        distTime = "%.8f" % ((solIndex + 1) / 1e8)

        src = os.path.join(checkPath, solutionTime)
        dst = os.path.join(checkPath, distTime)

        if self.comm.rank == 0:
            print("Moving time %s to %s" % (solutionTime, distTime))

        try:
            shutil.move(src, dst)
        except Exception:
            raise Error("Can not move %s to %s" % (src, dst))

        return

    def setdXvdFFDMat(self, designVarName, deltaVPointThreshold=1.0e-16):
        """
        Perturb each design variable and save the delta volume point coordinates
        to a mat, this will be used to calculate dRdFFD and dFdFFD in DAFoam

        Parameters
        ----------
        deltaVPointThreshold: float
            A threshold, any delta volume coordinates smaller than this value will be ignored

        """

        if self.DVGeo is None:
            raise Error("DVGeo not set!")

        # Get the FFD size
        nDVs = -9999
        xDV = self.DVGeo.getValues()
        nDVs = len(xDV[designVarName])

        # get the unperturbed point coordinates
        oldVolPoints = self.mesh.getSolverGrid()
        # get the size of xv, it is the number of points * 3
        nXvs = len(oldVolPoints)
        # get eps
        epsFFD = self.getOption("adjEpsDerivFFD")
        if self.comm.rank == 0:
            print("Caclculating the dXvdFFD matrix with epsFFD: " + str(epsFFD))

        dXvdFFDMat = PETSc.Mat().create(PETSc.COMM_WORLD)
        dXvdFFDMat.setSizes(((nXvs, None), (None, nDVs)))
        dXvdFFDMat.setFromOptions()
        dXvdFFDMat.setPreallocationNNZ((nDVs, nDVs))
        dXvdFFDMat.setUp()
        Istart, Iend = dXvdFFDMat.getOwnershipRange()

        # for each DV, perturb epsFFD and save the delta vol point coordinates
        for i in range(nDVs):
            # perturb
            xDV[designVarName][i] += epsFFD
            # set the dv to DVGeo
            self.DVGeo.setDesignVars(xDV)
            # update the vol points according to the new DV values
            self.updateVolumePoints()
            # get the new vol points
            newVolPoints = self.mesh.getSolverGrid()
            # assign the delta vol coords to the mat
            for idx in range(Istart, Iend):
                idxRel = idx - Istart
                deltaVal = newVolPoints[idxRel] - oldVolPoints[idxRel]
                if abs(deltaVal) > deltaVPointThreshold:  # a threshold
                    dXvdFFDMat[idx, i] = deltaVal
            # reset the perturbation of the dv
            xDV[designVarName][i] -= epsFFD

        # reset the volume mesh coordinates
        self.DVGeo.setDesignVars(xDV)
        self.updateVolumePoints()

        # assemble
        dXvdFFDMat.assemblyBegin()
        dXvdFFDMat.assemblyEnd()

        # viewer = PETSc.Viewer().createASCII("dXvdFFDMat_%s_%s.dat" % (designVarName, self.comm.size), "w")
        # viewer(dXvdFFDMat)

        self.solver.setdXvdFFDMat(dXvdFFDMat)

        return

    def updateVolumePoints(self):
        """
        Update the vol mesh point coordinates based on the current values of design variables
        """

        # update the CFD Coordinates
        self.ptSetName = self.getPointSetName("dummy")
        ptSetName = self.ptSetName
        if self.DVGeo is not None:
            if ptSetName not in self.DVGeo.points:
                coords0 = self.mapVector(self.coords0, self.allFamilies, self.designFamilyGroup)
                self.DVGeo.addPointSet(coords0, self.ptSetName)
                self.pointsSet = True

            # set the surface coords
            # if self.comm.rank == 0:
            #     print ('DVGeo PointSet UpToDate: '+str(self.DVGeo.pointSetUpToDate(ptSetName)))
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                # if self.comm.rank == 0:
                #     print 'Updating DVGeo PointSet....'
                coords = self.DVGeo.update(ptSetName, config=None)
                self.setSurfaceCoordinates(coords, self.designFamilyGroup)
                # if self.comm.rank == 0:
                #     print ('DVGeo PointSet UpToDate: '+str(self.DVGeo.pointSetUpToDate(ptSetName)))

            # warp the mesh
            self.mesh.warpMesh()

        return

    def _readMeshInfo(self):
        """
        Initialize mesh information and read mesh information
        """

        dirName = os.getcwd()

        self.fileNames, self.xv0, self.faces, self.boundaries, self.owners, self.neighbours = self._readOFGrid(dirName)
        self.xv = copy.copy(self.xv0)

        return

    def setSurfaceCoordinates(self, coordinates, groupName=None):
        """
        Set the updated surface coordinates for a particular group.
        Parameters
        ----------
        coordinates : numpy array
            Numpy array of size Nx3, where N is the number of coordinates on this processor.
            This array must have the same shape as the array obtained with getSurfaceCoordinates()
        groupName : str
            Name of family or group of families for which to return coordinates for.
        """
        if self.mesh is None:
            return

        if groupName is None:
            groupName = self.allWallsGroup

        self._updateGeomInfo = True
        if self.mesh is None:
            raise Error("Cannot set new surface coordinate locations without a mesh" "warping object present.")

        # First get the surface coordinates of the meshFamily in case
        # the groupName is a subset, those values will remain unchanged.

        meshSurfCoords = self.getSurfaceCoordinates(self.meshFamilyGroup)
        meshSurfCoords = self.mapVector(coordinates, groupName, self.meshFamilyGroup, meshSurfCoords)

        self.mesh.setSurfaceCoordinates(meshSurfCoords)

    def getSurfaceCoordinates(self, groupName=None):
        """
        Return the coordinates for the surfaces defined by groupName.

        Parameters
        ----------
        groupName : str
            Group identifier to get only coordinates cooresponding to
            the desired group. The group must be a family or a
            user-supplied group of families. The default is None which
            corresponds to all wall-type surfaces.

        Output
        ------
        xs: numpy array of size nPoints * 3 for surface points
        """

        if groupName is None:
            groupName = self.allWallsGroup

        # Get the required size
        npts, ncell = self._getSurfaceSize(groupName)
        xs = np.zeros((npts, 3), self.dtype)

        # loop over the families in this group and populate the surface
        famInd = self.families[groupName]
        counter = 0
        for Ind in famInd:
            name = self.basicFamilies[Ind]
            bc = self.boundaries[name]
            for ptInd in bc["indicesRed"]:
                xs[counter, :] = self.xv[ptInd]
                counter += 1

        return xs

    def _getSurfaceSize(self, groupName):
        """
        Internal routine to return the size of a particular surface. This
        does *NOT* set the actual family group
        """
        if groupName is None:
            groupName = self.allFamilies

        if groupName not in self.families:
            raise Error(
                "'%s' is not a family in the OpenFoam Case or has not been added"
                " as a combination of families" % groupName
            )

        # loop over the basic surfaces in the family group and sum up the number of
        # faces and nodes

        famInd = self.families[groupName]
        nPts = 0
        nCells = 0
        for Ind in famInd:
            name = self.basicFamilies[Ind]
            bc = self.boundaries[name]
            nCells += len(bc["facesRed"])
            nPts += len(bc["indicesRed"])

        return nPts, nCells

    def _computeBasicFamilyInfo(self):
        """
        Loop over the boundary data and compute necessary family
        information for the basic patches

        """
        # get the list of basic families
        self.basicFamilies = sorted(self.boundaries.keys())

        # save and return a list of the wall boundaries
        self.wallList = []
        counter = 0
        # for each boundary, figure out the unique list of volume node indices it uses
        for name in self.basicFamilies:
            # setup the basic families dictionary
            self.families[name] = [counter]
            counter += 1

            # Create a handle for this boundary
            bc = self.boundaries[name]

            # get the number of faces associated with this boundary
            nFace = len(bc["faces"])

            # create the index list
            indices = []

            # check that this isn't an empty boundary
            if nFace > 0:
                for iFace in bc["faces"]:
                    # get the node information for the current face
                    face = self.faces[iFace]
                    indices.extend(face)

            # Get the unique entries
            indices = np.unique(indices)

            # now create the reverse dictionary to connect the reduced set with the original
            inverseInd = {}
            for i in range(len(indices)):
                inverseInd[indices[i]] = i

            # Now loop back over the faces and store the connectivity in terms of the reduces index set
            facesRed = []
            for iFace in bc["faces"]:
                # get the node information for the current face
                face = self.faces[iFace]
                nNodes = len(face)
                # Generate the reduced connectivity.
                faceReduced = []
                for j in range(nNodes):
                    indOrig = face[j]
                    indRed = inverseInd[indOrig]
                    faceReduced.append(indRed)
                facesRed.append(faceReduced)

            # Check that the length of faces and facesRed are equal
            if not (len(bc["faces"]) == len(facesRed)):
                raise Error("Connectivity for faces on reduced index set is not the same length as original.")

            # put the reduced faces and index list in the boundary dict
            bc["facesRed"] = facesRed
            bc["indicesRed"] = list(indices)

            # now check for walls
            if bc["type"] == "wall" or bc["type"] == "slip" or bc["type"] == "cyclic":
                self.wallList.append(name)

        return

    def getPointSetName(self, apName):
        """
        Take the apName and return the mangled point set name.
        """
        return "openFoamCoords"

    def getSolverMeshIndices(self):
        """
        Get the list of indices to pass to the mesh object for the
        volume mesh mapping
        """

        # Setup External Warping
        nCoords = len(self.xv0.flatten())

        nCoords = self.comm.allgather(nCoords)
        offset = 0
        for i in range(self.comm.rank):
            offset += nCoords[i]

        meshInd = np.arange(nCoords[self.comm.rank]) + offset

        return meshInd

    def mapVector(self, vec1, groupName1, groupName2, vec2=None):
        """This is the main workhorse routine of everything that deals with
        families in pyDAFoam. The purpose of this routine is to convert a
        vector 'vec1' (of size Nx3) that was evaluated with
        'groupName1' and expand or contract it (and adjust the
        ordering) to produce 'vec2' evaluated on groupName2.

        A little ascii art might help. Consider the following "mesh"
        . Family 'fam1' has 9 points, 'fam2' has 10 pts and 'fam3' has
        5 points.  Consider that we have also also added two
        additional groups: 'f12' containing 'fam1' and 'fma2' and a
        group 'f23' that contains families 'fam2' and 'fam3'. The vector
        we want to map is 'vec1'. It is length 9+10. All the 'x's are
        significant values.

        The call: mapVector(vec1, 'f12', 'f23')

        will produce the "returned vec" array, containing the
        significant values from 'fam2', where the two groups overlap,
        and the new values from 'fam3' set to zero. The values from
        fam1 are lost. The returned vec has size 15.

            fam1     fam2      fam3
        |---------+----------+------|

        |xxxxxxxxx xxxxxxxxxx|        <- vec1
                  |xxxxxxxxxx 000000| <- returned vec (vec2)

        Parameters
        ----------
        vec1 : Numpy array
            Array of size Nx3 that will be mapped to a different family set.

        groupName1 : str
            The family group where the vector vec1 is currently defined

        groupName2 : str
            The family group where we want to the vector to mapped into

        vec2 : Numpy array or None
            Array containing existing values in the output vector we want to keep.
            If this vector is not given, the values will be filled with zeros.

        Returns
        -------
        vec2 : Numpy array
            The input vector maped to the families defined in groupName2.
        """
        if groupName1 not in self.families or groupName2 not in self.families:
            raise Error(
                "'%s' or '%s' is not a family in the CGNS file or has not been added"
                " as a combination of families" % (groupName1, groupName2)
            )

        # Shortcut:
        if groupName1 == groupName2:
            return vec1

        if vec2 is None:
            npts, ncell = self._getSurfaceSize(groupName2)
            vec2 = np.zeros((npts, 3), self.dtype)

        famList1 = self.families[groupName1]
        famList2 = self.families[groupName2]

        """
        This functionality is predicated on the surfaces being traversed in the
        same order every time. Loop over the allfamilies list, keeping track of sizes
        as we go and if the family is in both famLists, copy the values from vec1 to vec2.

        """

        vec1counter = 0
        vec2counter = 0

        for ind in self.families[self.allFamilies]:
            npts, ncell = self._getSurfaceSize(self.basicFamilies[ind])

            if ind in famList1 and ind in famList2:
                vec2[vec2counter : npts + vec2counter] = vec1[vec1counter : npts + vec1counter]

            if ind in famList1:
                vec1counter += npts

            if ind in famList2:
                vec2counter += npts

        return vec2

    def _initializeMeshPointVec(self):
        """
        Initialize the mesh point vec: xvVec
        """

        xvSize = len(self.xv) * 3
        self.xvVec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        self.xvVec.setSizes((xvSize, PETSc.DECIDE), bsize=1)
        self.xvVec.setFromOptions()

        self.xv2XvVec(self.xv, self.xvVec)

        # viewer = PETSc.Viewer().createASCII("xvVec", comm=PETSc.COMM_WORLD)
        # viewer(self.xvVec)

        return

    def _initializeStateVec(self):
        """
        Initialize state variable vector
        """

        wSize = self.solver.getNLocalAdjointStates()
        self.wVec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        self.wVec.setSizes((wSize, PETSc.DECIDE), bsize=1)
        self.wVec.setFromOptions()

        self.solver.ofField2StateVec(self.wVec)

        # viewer = PETSc.Viewer().createASCII("wVec", comm=PETSc.COMM_WORLD)
        # viewer(self.wVec)

        # adjoint vector
        self.psiVec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        self.psiVec.setSizes((wSize, PETSc.DECIDE), bsize=1)
        self.psiVec.setFromOptions()

        # total deriv vector
        # NOTE ***** need to change this!
        self.totalDerivVec = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
        self.totalDerivVec.setSizes((1, PETSc.DECIDE), bsize=1)
        self.totalDerivVec.setFromOptions()

        return

    def xv2XvVec(self, xv, xvVec):
        """
        Convert a Nx3 mesh point numpy array to a Petsc xvVec
        """

        xSize = len(xv)

        for i in range(xSize):
            for j in range(3):
                globalIdx = self.solver.getGlobalXvIndex(i, j)
                xvVec[globalIdx] = xv[i][j]

        xvVec.assemblyBegin()
        xvVec.assemblyEnd()

        return

    def xvFlatten2XvVec(self, xv, xvVec):
        """
        Convert a 3Nx1 mesh point numpy array to a Petsc xvVec
        """

        xSize = len(xv)
        xSize = int(xSize // 3)

        counterI = 0
        for i in range(xSize):
            for j in range(3):
                globalIdx = self.solver.getGlobalXvIndex(i, j)
                xvVec[globalIdx] = xv[counterI]
                counterI += 1

        xvVec.assemblyBegin()
        xvVec.assemblyEnd()

        return

    def xvVec2Xv(self, xvVec, xv):
        """
        Convert a Petsc xvVec to a Nx3 mesh point numpy array
        """

        xSize = len(xv)

        for i in range(xSize):
            for j in range(3):
                globalIdx = self.solver.getGlobalXvIndex(i, j)
                xv[i][j] = xvVec[globalIdx]

        return

    def _coloringComputationRequired(self):
        """
        check whether any of the required colorings are missing, if so
        recompute.
        """
        missingColorings = False

        if self.getOption("adjUseColoring"):
            # We need colorings, check if they exist
            requiredColorings = []

            requiredColorings.append("dRdWColoring_%d.bin" % self.nProcs)
            for objFunc in self.getOption("objFunc"):
                requiredColorings.append("dFdWColoring_%s_%d.bin" % (objFunc, self.nProcs))

            # now check for the require coloring
            for coloring in requiredColorings:
                if not os.path.exists(coloring):
                    missingColorings = True
                    break

        return missingColorings

    # base case files
    def _readOFGrid(self, caseDir):
        """
        Read in the mesh information we need to run the case using pyofm

        Parameters
        ----------
        caseDir : str
            The directory containing the openFOAM Mesh files
        """

        if self.comm.rank == 0:
            print("Reading OpenFOAM mesh information...")

        from pyofm import PYOFM

        # Initialize pyOFM
        ofm = PYOFM(comm=self.comm)

        # generate the file names
        fileNames = ofm.getFileNames(caseDir, comm=self.comm)

        # Read in the volume points
        x0 = ofm.readVolumeMeshPoints()

        # Read the face info for the mesh
        faces = ofm.readFaceInfo()

        # Read the boundary info
        boundaries = ofm.readBoundaryInfo(faces)

        # Read the cell info for the mesh
        owners, neighbours = ofm.readCellInfo()

        return fileNames, x0, faces, boundaries, owners, neighbours

    def setOption(self, name, value):
        """
        Set a value to options.
        NOTE: calling this function will only change the values in self.options
        It will NOT change values for allOptions_ in DAOption. To make the options changes
        from pyDAFoam to DASolvers, call self.updateDAOption()

        Parameters
        ----------
        name : str
           Name of option to set. Not case sensitive
        value : varies
           Value to set. Type is checked for consistency.

        Examples
        --------
        If self.options reads:
        self.options =
        {
            'solverName': [str, 'DASimpleFoam'],
            'flowEndTime': [float, 1.0]
        }
        Then, calling self.options('solverName', 'DARhoSimpleFoam') will give:
        self.options =
        {
            'solverName': [str, 'DARhoSimpleFoam'],
            'flowEndTime': [float, 1.0]
        }


        NOTE: if 'value' is of dict type, we will set all the subKey values in
        'value' dict to self.options, instead of overiding it

        For example, if self.options reads
        self.options =
        {
            'objFunc': [dict, {
                'name': 'CD',
                'direction': [1.0, 0.0, 0.0],
                'scale': 1.0}]
        }

        Then, calling self.setOption('objFunc', {'name': 'CL'}) will give:

        self.options =
        {
            'objFunc': [dict, {
                'name': 'CL',
                'direction': [1.0, 0.0, 0.0],
                'scale': 1.0}]
        }

        INSTEAD OF

        self.options =
        {
            'objFunc': [dict, {'name': 'CL'}]
        }
        """

        try:
            self.defaultOptions[name]
        except KeyError:
            Error("Option '%-30s' is not a valid %s option." % (name, self.name))

        # Make sure we are not trying to change an immutable option if
        # we are not allowed to.
        if name in self.imOptions:
            raise Error("Option '%-35s' cannot be modified after the solver " "is created." % name)

        # Now we know the option exists, lets check if the type is ok:
        if isinstance(value, self.defaultOptions[name][0]):
            # the type matches, now we need to check if the 'value' is of dict type, if yes, we only
            # replace the subKey values of 'value', instead of overiding all the subKey values
            if isinstance(value, dict):
                for subKey in value:
                    # no need to set self.options[name][0] since it has the right type
                    self.options[name][1][subKey] = value[subKey]
            else:
                # It is not dict, just set
                # no need to set self.options[name][0] since it has the right type
                self.options[name][1] = value
        else:
            raise Error(
                "Datatype for Option %-35s was not valid \n "
                "Expected data type is %-47s \n "
                "Received data type is %-47s" % (name, self.defaultOptions[name][0], type(value))
            )

    def getOption(self, name):
        """
        Get a value from options

        Parameters
        ----------
        name : str
           Name of option to get. Not case sensitive

        Returns
        -------
        value : varies
           Return the value of the option.
        """

        if name in self.defaultOptions:
            return self.options[name][1]
        else:
            raise Error("%s is not a valid option name." % name)

    def updateDAOption(self):
        """
        Update the allOptions_ in DAOption based on the latest self.options in
        pyDAFoam. This will pass the changes of self.options from pyDAFoam
        to DASolvers. NOTE: need to call this function after calling
        self.initSolver
        """

        if self.solverInitialized == 0:
            raise Error("self._initSolver not called!")

        self.solver.updateDAOption(self.options)

    def _printCurrentOptions(self):
        """
        Prints a nicely formatted dictionary of all the current solver
        options to the stdout on the root processor
        """

        if self.comm.rank == 0:
            print("+---------------------------------------+")
            print("|         All %s Options:         |" % self.name)
            print("+---------------------------------------+")
            # Need to assemble a temporary dictionary
            tmpDict = {}
            for key in self.options:
                tmpDict[key] = self.getOption(key)
            pp(tmpDict)

    def _getImmutableOptions(self):
        """
        We define the list of options that *cannot* be changed after the
        object is created. pyDAFoam will raise an error if a user tries to
        change these. The strings for these options are placed in a set
        """

        # return ("meshSurfaceFamily", "designSurfaceFamily")
        return ()


class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| pyDAFoam Error: "
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)
        Exception.__init__(self)

        return
