#!/usr/bin/env python
"""
DAFoam run script for the Onera M4 wing at subsonic speed
"""

# =================================================================================================
# Imports
# =================================================================================================
import os
import time
import argparse
import sys
import numpy as np
from mpi4py import MPI
from dafoam import PYDAFOAM, optFuncs
from pygeo import *
from pyspline import *
from idwarp import *
from pyoptsparse import Optimization, OPT
import numpy as np

# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help="Output directory", type=str, default="./")
parser.add_argument("--opt", help="optimizer to use", type=str, default="slsqp")
parser.add_argument("--task", help="type of run to do", type=str, default="opt")
parser.add_argument("--optVars", type=str, help="Vars for the optimizer", default="['shape']")
args = parser.parse_args()
exec("optVars=%s" % args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

UmagIn = 285.0
pIn = 101325.0
nuTildaIn = 4.5e-5
TIn = 300.0
CL_target = 0.270
alpha0 = 3.0
LRef = 0.64607
ARef = 0.7575
rhoRef = 1.0

# Set the parameters for optimization
aeroOptions = {
    # output options
    # design surfaces and cost functions
    "designSurfaceFamily": "designSurfaces",
    "designSurfaces": ["wing"],
    # flow setup
    "solverName": "DARhoSimpleCFoam",
    "turbulenceModel": "SpalartAllmaras",
    "flowCondition": "Compressible",
    "primalMinResTol": 1.0e-8,
    # adjoint setup
    "primalBC": {
        "UIn": {"variable": "U", "patch": "inout", "value": [UmagIn, 0.0, 0.0]},
        "pIn": {"variable": "p", "patch": "inout", "value": [pIn]},
        "TIn": {"variable": "T", "patch": "inout", "value": [TIn]},
        "nuTildaIn": {"variable": "nuTilda", "patch": "inout", "value": [nuTildaIn], "useWallFunction": True},
    },
    "primalVarBounds": {
        "UUpperBound": 1000.0,
        "ULowerBound": -1000.0,
        "pUpperBound": 500000.0,
        "pLowerBound": 20000.0,
        "eUpperBound": 500000.0,
        "eLowerBound": 100000.0,
        "rhoUpperBound": 5.0,
        "rhoLowerBound": 0.2,
    },
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "parallelToFlow",
                "alphaName": "alpha",
                "scale": 1.0 / (0.5 * rhoRef * UmagIn * UmagIn * ARef),
                "addToAdjoint": True,
            }
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "normalToFlow",
                "alphaName": "alpha",
                "scale": 1.0 / (0.5 * rhoRef * UmagIn * UmagIn * ARef),
                "addToAdjoint": True,
            }
        },
    },
    "adjEqnOption": {"pcFillLevel": 1, "jacMatReOrdering": "rcm"},
    "transonicPCOption": 1,
    "normalizeStates": {"U": UmagIn, "p": pIn, "nuTilda": nuTildaIn * 10.0, "phi": 1.0, "T": TIn},
    "adjPartDerivFDStep": {"State": 1e-6, "FFD": 1e-3},
    # Design variable setup
    "designVar": {
        "shapey": {"designVarType": "FFD"},
        "twist": {"designVarType": "FFD"},
        "alpha": {"designVarType": "AOA", "patch": "inout", "xAxisIndex": 0, "yAxisIndex": 1},
    },
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    "userotations": False,
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]],
}

# options for optimizers
outPrefix = outputDirectory + task + optVars[0]
if args.opt == "snopt":
    optOptions = {
        "Major feasibility tolerance": 1.0e-6,  # tolerance for constraint
        "Major optimality tolerance": 1.0e-6,  # tolerance for gradient
        "Minor feasibility tolerance": 1.0e-6,  # tolerance for constraint
        "Verify level": -1,
        "Function precision": 1.0e-6,
        "Major iterations limit": 40,
        "Nonderivative linesearch": None,
        "Print file": os.path.join(outPrefix + "_SNOPT_print.out"),
        "Summary file": os.path.join(outPrefix + "_SNOPT_summary.out"),
    }
elif args.opt == "slsqp":
    optOptions = {
        "ACC": 1.0e-6,  # convergence accuracy
        "MAXIT": 40,  # max optimization iterations
        "IFILE": os.path.join(outPrefix + "_SLSQP.out"),
    }
elif args.opt == "ipopt":
    optOptions = {
        "tol": 1.0e-6,  # convergence accuracy
        "max_iter": 40,  # max optimization iterations
        "output_file": os.path.join(outPrefix + "_IPOPT.out"),
    }
else:
    print("opt arg not valid!")
    exit(0)


# =================================================================================================
# DVGeo
# =================================================================================================
FFDFile = "./FFD/wingFFD.xyz"
DVGeo = DVGeometry(FFDFile)

# nTwists is the number of FFD points in the spanwise direction
nTwists = DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")

# twist function, we keep the root twist constant so the first
# element in the twist design variable is the twist at the 2nd
# spanwise location
def twist(val, geo):
    for i in range(1, nTwists):
        geo.rot_z["bodyAxis"].coef[i] = val[i - 1]

def alpha(val, geo):
    aoa = val[0] * np.pi / 180.0
    inletU = [float(UmagIn * np.cos(aoa)), float(UmagIn * np.sin(aoa)), 0]
    DASolver.setOption("primalBC", {"UIn": {"variable": "U", "patch": "inout", "value": inletU}})
    DASolver.updateDAOption()

# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
DVGeo.addGeoDVGlobal("twist", np.zeros(nTwists-1), twist, lower=-10.0, upper=10.0, scale=1.0)
DVGeo.addGeoDVGlobal("alpha", [alpha0], alpha, lower=0.0, upper=10.0, scale=1.0)

# =================================================================================================
# DAFoam
# =================================================================================================
DASolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
DASolver.printFamilyList()
DASolver.setMesh(mesh)
# set evalFuncs
evalFuncs = []
objFuncs = DASolver.getOption("objFunc")
for funcName in objFuncs:
    for funcPart in objFuncs[funcName]:
        if objFuncs[funcName][funcPart]["addToAdjoint"] == True:
            if not funcName in evalFuncs:
                evalFuncs.append(funcName)

# =================================================================================================
# DVCon
# =================================================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)

leList = [[0.01, 0.0, 1e-3], [0.7, 0.0, 1.19]]
teList = [[0.79, 0.0, 1e-3], [1.135, 0.0, 1.19]]

DVCon.addVolumeConstraint(leList, teList, nSpan=10, nChord=10, lower=1.0, upper=3, scaled=True)
DVCon.addThicknessConstraints2D(leList, teList, nSpan=10, nChord=10, lower=0.8, upper=3.0, scaled=True)

# Le/Te constraints
DVCon.addLeTeConstraints(0, "iLow")
DVCon.addLeTeConstraints(0, "iHigh")

# ================================================================================================
# optFuncs
# =================================================================================================
optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# =================================================================================================
# Task
# =================================================================================================
if task == "opt":

    DASolver.runColoring()
    optProb = Optimization("opt", optFuncs.calcObjFuncValues, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    # Add objective
    optProb.addObj("CD", scale=1)
    # Add physical constraints
    optProb.addCon("CL", lower=CL_target, upper=CL_target, scale=1)

    if gcomm.rank == 0:
        print(optProb)

    opt = OPT(args.opt, options=optOptions)
    histFile = os.path.join(outputDirectory, "%s_hist.hst" % args.opt)
    sol = opt(optProb, sens=optFuncs.calcObjFuncSens, storeHistory=histFile)
    if gcomm.rank == 0:
        print(sol)

elif task == "run":

    optFuncs.run()

elif task == "solveCL":

    optFuncs.solveCL(CL_star, "pitch", "CL")

elif task == "testSensShape":

    optFuncs.testSensShape()

else:
    print("task arg not found!")
    exit(0)
