#!/usr/bin/env python
"""
DAFoam run script for the NACA0012 airfoil at low-speed
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

UmagIn = 10.0
CL_star = 0.5
twist0 = 5.113547

# Set the parameters for optimization
aeroOptions = {
    # output options
    # design surfaces and cost functions
    "designSurfaceFamily": "designSurfaces",
    "designSurfaces": ["wing"],
    # flow setup
    "solverName": "DASimpleFoam",
    "turbulenceModel": "SpalartAllmaras",
    "flowCondition": "Incompressible",
    # adjoint setup
    "adjUseColoring": True,
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 0.2,
                "addToAdjoint": True,
            }
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "direction": [0.0, 1.0, 0.0],
                "scale": 0.2,
                "addToAdjoint": True,
            }
        },
    },
    "designVar": {"shapey": {"designVarType": "FFD"}, "twist": {"designVarType": "FFD"},},
    "normalizeStates": {"U": UmagIn, "p": UmagIn*UmagIn/2.0, "nuTilda": 0.001, "phi": 1.0},
    "adjEpsDerivState": 1e-6,
    "adjEpsDerivFFD": 1e-3,
    "maxResConLv4JacPCMat": {"pRes": 2, "phiRes": 1, "URes": 2, "nuTildaRes": 2},
    ########## misc setup ##########
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
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
        "Major iterations limit": 50,
        "Nonderivative linesearch": None,
        "Major step limit": 2.0,
        "Penalty parameter": 0.0,  # initial penalty parameter
        "Print file": os.path.join(outPrefix + "_SNOPT_print.out"),
        "Summary file": os.path.join(outPrefix + "_SNOPT_summary.out"),
    }
elif args.opt == "psqp":
    optOptions = {
        "TOLG": 1.0e-6,  # tolerance for gradient
        "TOLC": 1.0e-6,  # tolerance for constraint
        "MIT": 50,  # max optimization iterations
        "IFILE": os.path.join(outPrefix + "_PSQP.out"),
    }
elif args.opt == "slsqp":
    optOptions = {
        "ACC": 1.0e-6,  # convergence accuracy
        "MAXIT": 50,  # max optimization iterations
        "IFILE": os.path.join(outPrefix + "_SLSQP.out"),
    }
elif args.opt == "ipopt":
    optOptions = {
        "tol": 1.0e-6,  # convergence accuracy
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

# ref axis
x = [0.25, 0.25]
y = [0.00, 0.00]
z = [0.00, 0.10]
c1 = pySpline.Curve(x=x, y=y, z=z, k=2)
DVGeo.addRefAxis("bodyAxis", curve=c1, axis="z")


def twist(val, geo):
    # Set all the twist values
    for i in range(2):
        geo.rot_z["bodyAxis"].coef[i] = -val[0]


# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
DVGeo.addGeoDVGlobal("twist", twist0, twist, lower=-10.0, upper=10.0, scale=1.0)

# =================================================================================================
# DAFoam
# =================================================================================================
CFDSolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
CFDSolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
CFDSolver.addFamilyGroup(CFDSolver.getOption("designSurfaceFamily"), CFDSolver.getOption("designSurfaces"))
if MPI.COMM_WORLD.rank == 0:
    CFDSolver.printFamilyList()
CFDSolver.setMesh(mesh)
# set evalFuncs
evalFuncs = []
objFuncs = CFDSolver.getOption("objFunc")
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
[p0, v1, v2] = CFDSolver.getTriangulatedMeshSurface(groupName=CFDSolver.getOption("designSurfaceFamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)

leList = [[1e-4, 0.0, 1e-4], [1e-4, 0.0, 0.1 - 1e-4]]
teList = [[0.998 - 1e-4, 0.0, 1e-4], [0.998 - 1e-4, 0.0, 0.1 - 1e-4]]

DVCon.addVolumeConstraint(leList, teList, nSpan=2, nChord=20, lower=1.0, upper=3, scaled=True)
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
for i in [0, nFFDs_x-1]:
    for k in [0]:  # do not constrain k=1 because it is linked in the above symmetry constraint
        indSetA.append(pts1[i, 0, k])
        indSetB.append(pts1[i, 1, k])
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=1.0, lower=0.0, upper=0.0)

# ================================================================================================
# optFuncs
# =================================================================================================
optFuncs.CFDSolver = CFDSolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# =================================================================================================
# Task
# =================================================================================================
if task == "opt":

    CFDSolver.runColoring()
    optProb = Optimization("opt", optFuncs.getObjFuncValues, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    # Add objective
    optProb.addObj("CD", scale=1)
    # Add physical constraints
    optProb.addCon("CL", lower=CL_star, upper=CL_star, scale=1)

    if gcomm.rank == 0:
        print(optProb)

    opt = OPT(args.opt, options=optOptions)
    histFile = os.path.join(outputDirectory, "%s_hist.hst" % args.opt)
    sol = opt(optProb, sens=optFuncs.getObjFuncSens, storeHistory=histFile)
    if gcomm.rank == 0:
        print(sol)

elif task == "run":

    optFuncs.run()

elif task == "solveCL":

    optFuncs.solveCL(CL_star, "twist", "CL")

elif task == "testSensShape":

    optFuncs.testSensShape()

else:
    print("task arg not found!")
    exit(0)
