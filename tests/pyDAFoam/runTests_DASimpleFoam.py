#!/usr/bin/env python
"""
Run Python tests
"""

from mpi4py import MPI
from dafoam.src import pyDAFoam
import sys
import os

comm = MPI.COMM_WORLD

os.chdir("../input/CurvedCubeHexMesh")

# test incompressible solvers
aeroOptions = {
    "solverName": "DASimpleFoam",
    "flowCondition": "Incompressible",
    "printAllOptions": False,
    "designSurfaceFamily": "designSurface",
    "designSurfaces": ["walls"],
    "objFunc": {
        "func1": {
            "part1": {
                "objFuncName": "force",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 0.5,
                "addToAdjoint": False,
            },
            "part2": {
                "objFuncName": "force",
                "source": "patchToFace",
                "patches": ["wallsbump", "frontandback"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 0.5,
                "addToAdjoint": False,
            },
        },
        "func2": {
            "part1": {
                "objFuncName": "force",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump", "frontandback"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": False,
            }
        },
    },
}

if comm.rank == 0:
    os.system("rm -rf  0/* processor*")
    os.system("cp 0.incompressible/* 0/")

dafoam = pyDAFoam.PYDAFOAM(options=aeroOptions)

meshOK = dafoam.runCheckMesh()
dafoam.runPrimalSolver()
dafoam.runPrimalSolver()
dafoam = None
