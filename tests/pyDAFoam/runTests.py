#!/usr/bin/env python
"""
Run Python tests
"""

from mpi4py import MPI
from dafoam.src import pyDAFoam
import sys
import os

os.chdir('../input/CurvedCubeHexMesh')

aeroOptions = {
    "solverName": "DASimpleFoam",
    "printAllOptions": False,
}

dafoam = pyDAFoam.PYDAFOAM(options=aeroOptions)

meshOK = dafoam.runCheckMesh()
dafoam.initSolver()
dafoam.runPrimalSolver()
dafoam.runPrimalSolver()

os.chdir('../../pyDAFoam')
