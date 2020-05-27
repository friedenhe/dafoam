#!/usr/bin/env python
"""
Run C++ tests
"""

from mpi4py import MPI
from pyTestDAFoamIncompressible import pyTestDAFoamIncompressible
import sys
import os
import petsc4py

petsc4py.init(sys.argv)

comm = MPI.COMM_WORLD


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


def checkErrors(testName, errorCode):
    if errorCode != 0:
        raise Error("Tests Failed for %s! Rank %d " % (testName, comm.rank))
    else:
        print("Tests Passed for %s! Rank %d" % (testName, comm.rank))


pyDict = {
    "key1": [int, 15],
    "key2": [float, 5.5],
    "key3": [str, "solver1"],
    "key4": [bool, False],
    "key5": [list, [1, 2, 3]],
    "key6": [list, [1.5, 2.3, 3.4]],
    "key7": [list, ["ele1", "ele2", "ele3"]],
    "key8": [list, [False, True, False]],
    "key9": [
        dict,
        {
            "subkey1": [int, 30],
            "subkey2": [float, 3.5],
            "subkey3": [str, "solver2"],
            "subkey4": [bool, True],
            "subkey5": [list, [4, 5, 6]],
            "subkey6": [list, [2.5, 7.7, 8.9]],
            "subkey7": [list, ["ele4", "ele5", "ele6"]],
            "subkey8": [list, [True, False, True]],
        },
    ],
}

parallelFlag = ""
if comm.size > 1:
    parallelFlag = "-parallel"
solverArg = "TestDAFoamIncompressible -python " + parallelFlag
tests = pyTestDAFoamIncompressible(solverArg.encode())

# Test1: DAUtility
os.chdir("../input/CurvedCubeHexMesh")
testErrors = tests.testDAUtility(pyDict)
checkErrors("DAUtility", testErrors)
os.chdir("../../DAFoamIncompressible")

# Test2: DAOption
os.chdir("../input/CurvedCubeHexMesh")
testErrors = tests.testDAOption(pyDict)
checkErrors("DAOption", testErrors)
os.chdir("../../DAFoamIncompressible")

# Test2: DARegState
os.chdir("../input/CurvedCubeHexMesh")
testDict = {"solverName": [str, "DASimpleFoam"], "turbulenceModel": [str, "SpalartAllmaras"]}
testErrors = tests.testDARegState(testDict)
checkErrors("DARegState", testErrors)
os.chdir("../../DAFoamIncompressible")

# Test2: DARegState
os.chdir("../input/CurvedCubeHexMesh")
testDict = {
    "solverName": [str, "DASimpleFoam"],
    "turbulenceModel": [str, "SpalartAllmaras"],
    "adjUseColoring": [bool, True],
    "adjJacMatOrdering": [str, "state"],
    "objFunc": [
        dict,
        {
            "func1": {
                "part1": {
                    "objFuncName": "force",
                    "source": "patchToFace",
                    "patch": ["walls", "wallsbump"],
                    "scale": 0.5,
                    "addToAdjoint": False,
                },
                "part2": {
                    "objFuncName": "force",
                    "source": "patchToFace",
                    "patch": ["wallsbump", "frontandback"],
                    "scale": 0.5,
                    "addToAdjoint": False,
                },
            },
            "func2": {
                "part1": {
                    "objFuncName": "force",
                    "source": "patchToFace",
                    "patch": ["walls", "wallsbump", "frontandback"],
                    "scale": 1.0,
                    "addToAdjoint": False,
                }
            },
        },
    ],
}
testErrors = tests.testDAObjFunc(testDict)
checkErrors("DAObjFunc", testErrors)
os.chdir("../../DAFoamIncompressible")
