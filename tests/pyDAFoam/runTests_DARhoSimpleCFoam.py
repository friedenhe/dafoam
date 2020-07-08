#!/usr/bin/env python
"""
Run Python tests for DARhoSimpleCFoam
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

os.chdir("../input/CurvedCubeHexMesh")

if gcomm.rank == 0:
    os.system("rm -rf  0/* processor*")
    os.system("cp 0.compressible/* 0/")

# test incompressible solvers
aeroOptions = {
    "solverName": "DARhoSimpleFoam",
    "flowCondition": "Compressible",
    "designSurfaceFamily": "designSurface",
    "designSurfaces": ["wallsbump"],
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["walls"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 0.1,
                "addToAdjoint": True,
            },
            "part2": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wallsbump", "frontandback"],
                "direction": [1.0, 0.0, 0.0],
                "scale": 0.1,
                "addToAdjoint": True,
            },
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["walls", "wallsbump", "frontandback"],
                "direction": [0.0, 1.0, 0.0],
                "scale": 0.1,
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
DASolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
if gcomm.rank == 0:
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

# DVCon
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)

# optFuncs
optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# Opt
DASolver.runColoring()
xDV = DVGeo.getValues()
funcs = {}
funcs, fail = optFuncs.calcObjFuncValues(xDV)
funcsSens = {}
funcsSens, fail = optFuncs.calcObjFuncSens(xDV, funcs)

if checkRegVal:

    CD = funcs["CD"]
    CL = funcs["CL"]
    l2_CD_shapex = np.linalg.norm(funcsSens["CD"]["shapex"])
    l2_CD_shapey = np.linalg.norm(funcsSens["CD"]["shapey"])
    l2_CL_shapex = np.linalg.norm(funcsSens["CL"]["shapex"])
    l2_CL_shapey = np.linalg.norm(funcsSens["CL"]["shapey"])

    if gcomm.rank == 0:
        print(CD, CL, l2_CD_shapex, l2_CD_shapey, l2_CL_shapex, l2_CL_shapey)

    CD_ref = 10.050080671067752
    CL_ref = 37.87716335434652

    l2_CD_shapex_ref = 3.2663055157149183
    l2_CD_shapey_ref = 25.162005608326258
    l2_CL_shapex_ref = 8.056355533817573
    l2_CL_shapey_ref = 64.8422882564767

    diff_CD = abs(CD - CD_ref) / CD_ref
    diff_CL = abs(CL - CL_ref) / CL_ref
    diff_CD_shapex = abs(l2_CD_shapex - l2_CD_shapex_ref) / l2_CD_shapex_ref
    diff_CD_shapey = abs(l2_CD_shapey - l2_CD_shapey_ref) / l2_CD_shapey_ref
    diff_CL_shapex = abs(l2_CL_shapex - l2_CL_shapex_ref) / l2_CL_shapex_ref
    diff_CL_shapey = abs(l2_CL_shapey - l2_CL_shapey_ref) / l2_CL_shapey_ref

    checkFail = 0
    funcTol = 1.0e-10
    funcSensTol = 1.0e-8
    if diff_CD > funcTol or diff_CL > funcTol:
        checkFail += 1
    if diff_CD_shapex > funcSensTol or diff_CD_shapey > funcSensTol or diff_CL_shapex > funcSensTol or diff_CL_shapey > funcSensTol:
        checkFail += 1

    if checkFail > 0:
        print("Failed!")
        exit(1)
    else:
        print("Succes!")

