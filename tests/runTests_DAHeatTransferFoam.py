#!/usr/bin/env python
"""
Run Python tests for DASolidDisplacementFoam
"""

from mpi4py import MPI
from dafoam import PYDAFOAM
import sys
import os
from pygeo import *
from pyspline import *
from idwarp import *
import numpy as np
from testFuncs import *

gcomm = MPI.COMM_WORLD

os.chdir("./input/ChannelConjugateHeat/thermal")

aeroOptions = {
    "designSurfaces": ["channel_outer", "channel_inner", "channel_sides"],
    "solverName": "DAHeatTransferFoam",
    "objFunc": {
        "HF_INNER": {
            "part1": {
                "type": "wallHeatFlux",
                "source": "patchToFace",
                "patches": ["channel_inner"],
                "scale": 1,
                "addToAdjoint": True,
            }
        },
        "HF_OUTER": {
            "part1": {
                "type": "wallHeatFlux",
                "source": "patchToFace",
                "patches": ["channel_outer"],
                "scale": 1,
                "addToAdjoint": True,
            }
        },
    },
    "debug": False,
    "primalMinResTol": 1e-12,
}
DASolver = PYDAFOAM(options=aeroOptions, comm=MPI.COMM_WORLD)

DASolver()
funcs = {}
evalFuncs = ["HF_INNER", "HF_OUTER"]
DASolver.evalFunctions(funcs, evalFuncs)

if gcomm.rank == 0:
    reg_write_dict(funcs, 1e-8, 1e-10)