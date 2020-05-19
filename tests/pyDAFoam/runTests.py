#!/usr/bin/env python
"""
Run Python tests
"""

from mpi4py import MPI
from dafoam.src import pyDAFoam
import sys
import os

comm = MPI.COMM_WORLD

aeroOptions = {
    "solverName": "DASimpleFoam",
    "printAllOptions": False,
}

os.chdir('../input/CurvedCubeHexMesh')
if comm.rank == 0:
    os.system('cp 0.incompressible/* 0/')

dafoam = pyDAFoam.PYDAFOAM(options=aeroOptions)

meshOK = dafoam.runCheckMesh()
dafoam.runPrimalSolver()
dafoam.runPrimalSolver()

if comm.rank == 0:
    os.system('rm -rf  0/* processor*')
    os.system('cp 0.compressible/* 0/')

dafoam.runDecomposePar()
dafoam.setOption("solverName", "DARhoSimpleFoam")
dafoam.runPrimalSolver()

os.chdir('../../pyDAFoam')
