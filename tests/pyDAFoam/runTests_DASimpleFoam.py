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
}

if comm.rank == 0:
    os.system("rm -rf  0/* processor*")
    os.system("cp 0.incompressible/* 0/")

dafoam = pyDAFoam.PYDAFOAM(options=aeroOptions)

meshOK = dafoam.runCheckMesh()
dafoam.runPrimalSolver()
dafoam.runPrimalSolver()
dafoam = None
