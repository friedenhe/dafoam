#!/usr/bin/env python
"""
DAFoam run script for the NACA0012 airfoil at subsonic speed
"""

# =================================================================================================
# Imports
# =================================================================================================
import os
import argparse
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

UmagIn = 238.0
pIn = 101325.0
nuTildaIn = 4.5e-5
TIn = 300.0
ARef = 0.1
rhoRef = 1.0
CL_target = 0.5
alpha0 = 2.764442

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
        "rhoLowerBound": 0.2},
    # adjoint setup
    "adjUseColoring": True,
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
    "adjEpsDerivState": 1e-7,
    "adjEpsDerivFFD": 1e-3,
    # Design variable setup
    "designVar": {
        "shapey": {"designVarType": "FFD"},
        "alpha": {"designVarType": "AOA", "patch": "inout", "xAxisIndex": 0, "yAxisIndex": 1},
    }
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    "userotations": False,
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
}

# options for optimizers
outPrefix = outputDirectory + task + optVars[0]
if args.opt == "snopt":
    optOptions = {
        "Major feasibility tolerance": 1.0e-7,  # tolerance for constraint
        "Major optimality tolerance": 1.0e-7,  # tolerance for gradient
        "Minor feasibility tolerance": 1.0e-7,  # tolerance for constraint
        "Verify level": -1,
        "Function precision": 1.0e-7,
        "Major iterations limit": 50,
        "Nonderivative linesearch": None,
        "Print file": os.path.join(outPrefix + "_SNOPT_print.out"),
        "Summary file": os.path.join(outPrefix + "_SNOPT_summary.out"),
    }
elif args.opt == "slsqp":
    optOptions = {
        "ACC": 1.0e-7,  # convergence accuracy
        "MAXIT": 50,  # max optimization iterations
        "IFILE": os.path.join(outPrefix + "_SLSQP.out"),
    }
elif args.opt == "ipopt":
    optOptions = {
        "tol": 1.0e-7,  # convergence accuracy
        "max_iter": 50,  # max optimization iterations
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
DVGeo.addGeoDVGlobal("alpha", [alpha0], alpha, lower=-10.0, upper=10.0, scale=1.0)

# =================================================================================================
# DAFoam
# =================================================================================================
DASolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
if MPI.COMM_WORLD.rank == 0:
    DASolver.printFamilyList()
DASolver.setMesh(mesh)
# set evalFuncs
evalFuncs = []
objFuncs = DASolver.getOption("objFunc")
for funcName in objFuncs:
    for funcPart in objFuncs[funcName]:
        if objFuncs[funcName][funcPart]["addToAdjoint"] is True:
            if funcName not in evalFuncs:
                evalFuncs.append(funcName)

# =================================================================================================
# DVCon
# =================================================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)

leList = [[1e-4, 0.0, 1e-4], [1e-4, 0.0, 0.1 - 1e-4]]
teList = [[0.998 - 1e-4, 0.0, 1e-4], [0.998 - 1e-4, 0.0, 0.1 - 1e-4]]

DVCon.addVolumeConstraint(leList, teList, nSpan=2, nChord=10, lower=1.0, upper=3, scaled=True)
DVCon.addThicknessConstraints2D(leList, teList, nSpan=2, nChord=10, lower=0.8, upper=3.0, scaled=True)

# Create a linear constraint so that the curvature at the symmetry plane is zero
nFFDs_x = 5
pts1 = DVGeo.getLocalIndex(0)
indSetA = []
indSetB = []
for i in range(nFFDs_x):
    for j in [0, 1]:
        indSetA.append(pts1[i, j, 1])
        indSetB.append(pts1[i, j, 0])
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=-1.0, lower=0.0, upper=0.0)

# Create a linear constraint so that the leading and trailing edges do not change
pts1 = DVGeo.getLocalIndex(0)
indSetA = []
indSetB = []
for i in [0, nFFDs_x - 1]:
    for k in [0]:  # do not constrain k=1 because it is linked in the above symmetry constraint
        indSetA.append(pts1[i, 0, k])
        indSetB.append(pts1[i, 1, k])
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=1.0, lower=0.0, upper=0.0)


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

    optFuncs.solveCL(CL_target, "alpha", "CL")

elif task == "testSensShape":

    optFuncs.testSensShape()

else:
    print("task arg not found!")
    exit(0)
