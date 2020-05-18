#!/usr/bin/env python
"""
Run C++ tests
"""

from mpi4py import MPI
from pyTestDAFoamIncompressible import pyTestDAFoamIncompressible
import sys
import petsc4py

petsc4py.init(sys.argv)

comm = MPI.COMM_WORLD


def checkErrors(testName, errorCode):
    if errorCode != 0:
        print("%s Failed! Rank %d" % (testName, comm.rank))
        exit(1)
    else:
        print("%s Passed! Rank %d" % (testName, comm.rank))
        exit(0)


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


solverArg = "tests -python"
tests = pyTestDAFoamIncompressible(solverArg.encode())

# Test1: DAUtility
testErrors = tests.testDAUtility(pyDict)
checkErrors("Test DAUtility", testErrors)
