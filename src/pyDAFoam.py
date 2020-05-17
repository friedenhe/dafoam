#!/usr/bin/env python

"""
pyDAFoam: Python wrapper for DAFoam.
"""

import subprocess
import os
import sys
from pprint import pprint as pp
from mpi4py import MPI
from collections import OrderedDict
import petsc4py

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

        return

    def __del__(self):

        self.primalSolver = None

        return

    def _getDefOptions(self):
        """
        Setup default options

        Returns
        -------

        defOpts : dict
            All the DAFoam options.
            NOTE: We support only one level of sub dictionary
        """
        defOpts = {
            # flow options
            "flowEndTime": [float, 1.0],
            "flowDeltaT": [float, 1.0],
            # adjoint options
            "adjUseColoring": [bool, True],
            "adjEpsDerivFFD": [float, 1.0e-6],
            # system options
            "rootDir": [str, "./"],
            "solverName": [str, "DASimpleFoam"],
            "printAllOptions": [bool, False],

        }

        return defOpts

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
            raise Error(
                "The 'options' keyword argument must be passed \
            pyDAFoam. The options dictionary must contain (at least) the gridFile \
            entry for the grid"
            )

        # set immutable options that users should not change during the optimization
        self.imOptions = self._getImmutableOptions()

        # Load all the option information:
        self.defaultOptions = self._getDefOptions()

        # Set options based on defaultOptions
        # we basically overwrite defaultOptions with the given options
        self.options = OrderedDict()
        for key in self.defaultOptions:
            if len(self.defaultOptions[key]) != 2:
                raise Error(
                    "key %s has wrong format! \
                    Example: {'iters' : [int, 1]}"
                    % key
                )

            self.setOption(key, self.defaultOptions[key][1])
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

    def readMesh(self):
        """
        Read the OpenFOAM mesh information using the pyofm repo
        """

        from pyofm import PYOFM

        # Initialize pyOFM
        self.ofm = PYOFM(comm=self.comm)

        return

    def runPrimalSolver(self):
        """
        Run primal solver to compute state variables and objectives
        """

        if self.comm.rank == 0:
            print("\n")
            print("+--------------------------------------------------------------------------+")
            print("|                        Running Primal Solver                             |")
            print("+--------------------------------------------------------------------------+")

        self.solver.solvePrimal()

        return

    def initSolver(self):
        """
        Initialize solver
        """

        solverName = self.getOption("solverName")
        solverArg = solverName + " -python " + self.parallelFlag
        if solverName == "DASimpleFoam":
            from .pyDASimpleFoam import pyDASimpleFoam

            self.solver = pyDASimpleFoam(solverArg.encode(), self.options)
            self.solver.init()

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

        solverName = "CheckMesh"
        solverArg = solverName + " -python " + self.parallelFlag
        from .pyCheckMesh import pyCheckMesh

        checkMesh = pyCheckMesh(solverArg.encode())
        meshOK = checkMesh.run()
        checkMesh = None

        return meshOK

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

    def setOption(self, name, value):
        """
        Set a value to options

        Parameters
        ----------
        name : str
           Name of option to set. Not case sensitive
        value : varies
           Value to set. Type is checked for consistency.
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
            # Just set:
            self.options[name] = [type(value), value]
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

        return ("meshSurfaceFamily", "designSurfaceFamily")


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
