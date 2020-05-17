#!/usr/bin/env python
"""
Run Python tests
"""

from mpi4py import MPI
from dafoam import *
import sys
import os

os.chdir('../input/CurvedCubeHexMesh')

aeroOptions = {
    "solverName": "DASimpleFoam",
    "printAllOptions": False,
}

pyDAFoam = PYDAFOAM(options=aeroOptions)

meshOK = pyDAFoam.runCheckMesh()
pyDAFoam.initSolver()
pyDAFoam.runPrimalSolver()
pyDAFoam.runPrimalSolver()

os.chdir('../../python')