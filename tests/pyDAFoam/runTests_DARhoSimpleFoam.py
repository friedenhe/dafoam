#!/usr/bin/env python
"""
Run Python tests
"""

from mpi4py import MPI
from dafoam.src import pyDAFoam
import sys
import os

comm = MPI.COMM_WORLD

os.chdir('../input/CurvedCubeHexMesh')

# test compressible solvers
aeroOptions = {
    "solverName": "DARhoSimpleFoam",
    "flowCondition": "Compressible",
    "printAllOptions": False,
}
if comm.rank == 0:
    os.system('rm -rf  0/* processor*')
    os.system('cp 0.compressible/* 0/')

dafoam = pyDAFoam.PYDAFOAM(options=aeroOptions)
dafoam.runPrimalSolver()
dafoam.runPrimalSolver()
dafoam = None
