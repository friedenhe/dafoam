#!/usr/bin/env python
"""
Run Python tests for DASimpleFoam
"""

from mpi4py import MPI
from dafoam import PYDAFOAM, optFuncs
import sys
import os
from pygeo import *
from pyspline import *
from idwarp import *
from pyoptsparse import Optimization, OPT
import numpy as np

checkRegVal = 1
if len(sys.argv) == 1:
    checkRegVal = 1
elif sys.argv[1] == "noCheckVal":
    checkRegVal = 0
else:
    print("sys.argv %s not valid!" % sys.argv[1])
    exit(1)

gcomm = MPI.COMM_WORLD

CL_target = 12.6
CMZ_target = 5.82
CMY_target = 1.71


os.chdir("../input/CurvedCubeHexMesh")

if gcomm.rank == 0:
    os.system("rm -rf  0/* processor*")
    os.system("cp 0.incompressible/* 0/")

# test incompressible solvers
aeroOptions = {
    "solverName": "DASimpleFoam",
    "flowCondition": "Incompressible",
    "designSurfaceFamily": "designSurface",
    "designSurfaces": ["wallsbump"],
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["walls"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": True,
            },
            "part2": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wallsbump", "frontandback"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": True,
            },
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump", "frontandback"],
                "direction": [0.0, 1.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
        "CMZ": {
            "part1": {
                "type": "moment",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump", "frontandback"],
                "direction": [0.0, 0.0, 1.0],
                "center": [0.5, 0.5, 0.5],
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
        "CMY": {
            "part1": {
                "type": "moment",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump", "frontandback"],
                "direction": [0.0, 1.0, 0.0],
                "center": [0.0, 0.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
    },
    "designVar": {"shapex": {"designVarType": "FFD"}, "shapey": {"designVarType": "FFD"}},
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
}

optOptions = {
    "ACC": 1.0e-5,  # convergence accuracy
    "MAXIT": 2,  # max optimization iterations
    "IFILE": "opt_SLSQP.out",
}

# DVGeo
FFDFile = "./FFD/bumpFFD.xyz"
DVGeo = DVGeometry(FFDFile)

# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[1:3, 1, 1:3].flatten()
PS = geo_utils.PointSelect("list", indexList)
DVGeo.addGeoDVLocal("shapey", lower=-0.1, upper=0.1, axis="y", scale=1.0, pointSelect=PS)
DVGeo.addGeoDVLocal("shapex", lower=-0.1, upper=0.1, axis="x", scale=1.0, pointSelect=PS)

# DAFoam
CFDSolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
CFDSolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
CFDSolver.addFamilyGroup(CFDSolver.getOption("designSurfaceFamily"), CFDSolver.getOption("designSurfaces"))
if gcomm.rank == 0:
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

# DVCon
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = CFDSolver.getTriangulatedMeshSurface(groupName=CFDSolver.getOption("designSurfaceFamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)

# Test a linear constraint
pts1 = DVGeo.getLocalIndex(0)
indSetA = []
indSetB = []
for i in range(3):
    for k in range(3):
        indSetA.append(pts1[i, 0, k])
        indSetB.append(pts1[i, 1, k])
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=-1.0, lower=0.0, upper=0.0)

# optFuncs
optFuncs.CFDSolver = CFDSolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# Opt
CFDSolver.runColoring()
optProb = Optimization("opt", optFuncs.calcObjFuncValues, comm=gcomm)
DVGeo.addVariablesPyOpt(optProb)
DVCon.addConstraintsPyOpt(optProb)

# Add objective
optProb.addObj("CD", scale=1)
# Add physical constraints
optProb.addCon("CL", lower=CL_target, upper=CL_target, scale=1)
optProb.addCon("CMZ", lower=CMZ_target, upper=CMZ_target, scale=1)
optProb.addCon("CMY", lower=CMY_target, upper=CMY_target, scale=1)

if gcomm.rank == 0:
    print(optProb)

opt = OPT("slsqp", options=optOptions)
histFile = "slsqp_hist.hst"
sol = opt(optProb, sens=optFuncs.calcObjFuncSens, storeHistory=histFile)
if gcomm.rank == 0:
    print(sol)

if checkRegVal:
    xDVs = DVGeo.getValues()

    l2_shapex = np.linalg.norm(xDVs["shapex"])
    l2_shapey = np.linalg.norm(xDVs["shapey"])

    if gcomm.rank == 0:
        print("l2_shapex: ", l2_shapex)
        print("l2_shapey: ", l2_shapey)

    ref_shapex = 0.18792724386851453
    ref_shapey = 0.1646036795038631

    diff_shapey = abs(l2_shapey - ref_shapey)
    diff_shapex = abs(l2_shapex - ref_shapex)

    dvTol = 1.0e-8
    if diff_shapey > dvTol or diff_shapex > dvTol:
        print("Failed!")
        exit(1)
    else:
        print("Succes!")

