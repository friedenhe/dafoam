#!/usr/bin/env python
"""
DAFoam run script for the NACA0012 airfoil at low-speed (multi point)
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
from idwarp import USMesh
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

nMultiPoints = 2
MPWeights = [0.25, 0.25, 0.5]
UmagIn = [10.0, 10.0, 10.0]
URef = 10.0
CL_target = [0.2, 0.8, 0.5]
alpha0 = [1.992939, 5.139186]
pIn = 0.0
nuTildaIn = 4.5e-5
ARef = 0.1

# Set the parameters for optimization
daOptions = {
    # output options
    # design surfaces and cost functions
    "designSurfaceFamily": "designSurfaces",
    "designSurfaces": ["wing"],
    # primal setup
    "multiPoint": True,
    "nMultiPoints": nMultiPoints,
    "solverName": "DASimpleFoam",
    "turbulenceModel": "SpalartAllmaras",
    "flowCondition": "Incompressible",
    "primalMinResTol": 1.0e-8,
    "primalBC": {
        "UIn": {"variable": "U", "patch": "inout", "value": [URef, 0.0, 0.0]},
        "pIn": {"variable": "p", "patch": "inout", "value": [pIn]},
        "nuTildaIn": {"variable": "nuTilda", "patch": "inout", "value": [nuTildaIn], "useWallFunction": True},
    },
    # adjoint setup
    "adjUseColoring": True,
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "parallelToFlow",
                "alphaName": "mp0_alpha",
                "scale": 1.0 / (0.5 * URef * URef * ARef),
                "addToAdjoint": True,
            }
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "normalToFlow",
                "alphaName": "mp0_alpha",
                "scale": 1.0 / (0.5 * URef * URef * ARef),
                "addToAdjoint": True,
            }
        },
    },
    "adjEqnOption": {"pcFillLevel": 1, "jacMatReOrdering": "rcm"},
    "normalizeStates": {"U": URef, "p": URef * URef / 2.0, "nuTilda": nuTildaIn * 10.0, "phi": 1.0},
    "adjPartDerivFDStep": {"State": 1e-7, "FFD": 1e-3},
    # Design variable setup
    "designVar": {
        "shapey": {"designVarType": "FFD"},
        # will add alpha in addGeoDVGlobal
    }
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


def dummyFunc(val, geo):
    pass


# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
for i in range(nMultiPoints):
    # NOTE: here we don't need to implement the alpha function because the alpha values will be changed
    # in setMultiPointCondition. So we provide a dummyFunc
    DVGeo.addGeoDVGlobal("mp%d_alpha" % i, alpha0[i], dummyFunc, lower=-10.0, upper=10.0, scale=1.0)
    # add alpha for designVar
    daOptions["designVar"]["mp%d_alpha" % i] = {
        "designVarType": "AOA",
        "patch": "inout",
        "xAxisIndex": 0,
        "yAxisIndex": 1,
    }

# =================================================================================================
# DAFoam
# =================================================================================================
DASolver = PYDAFOAM(options=daOptions, comm=gcomm)
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
# provide a function to set primal conditions
def setMultiPointCondition(xDV, index):
    aoa = xDV["mp%d_alpha" % index].real * np.pi / 180.0
    inletU = [float(UmagIn[index] * np.cos(aoa)), float(UmagIn[index] * np.sin(aoa)), 0]
    DASolver.setOption("primalBC", {"UIn": {"variable": "U", "patch": "inout", "value": inletU}})
    DASolver.updateDAOption()
    return


# provide a function to assemble the funcs from MP
def setMultiPointObjFuncs(funcs, funcsMP, index):
    for key in funcs:
        if "fail" in key:
            pass
        elif "DVCon" in key:
            funcsMP[key] = funcs[key]
        elif "CD" in key:
            # funcsMP["mp%d_CD" % index] = funcs[key]
            try:
                funcsMP["obj"] += funcs[key] * MPWeights[index]
            except Exception:
                funcsMP["obj"] = 0.0
                funcsMP["obj"] += funcs[key] * MPWeights[index]
        elif "CL" in key:
            funcsMP["mp%d_CL" % index] = funcs[key]
    return


# provide a function to assemble the funcs from MP
def setMultiPointObjFuncsSens(xDVs, funcsMP, funcsSens, funcsSensMP, index):
    for key in funcsMP:
        try:
            keySize = len(funcsMP[key])
        except Exception:
            keySize = 1
        try:
            funcsSensMP[key]
        except Exception:
            funcsSensMP[key] = {}

        if "fail" in key:
            pass
        elif "DVCon" in key:
            funcsSensMP[key]["mp%d_alpha" % index] = np.zeros((keySize, 1), "d")
            funcsSensMP[key]["shapey"] = funcsSens[key]["shapey"]
        elif "obj" in key:
            funcsSensMP[key]["mp%d_alpha" % index] = funcsSens["CD"]["mp%d_alpha" % index] * MPWeights[index]
            try:
                funcsSensMP[key]["shapey"] += funcsSens["CD"]["shapey"] * MPWeights[index]
            except Exception:
                funcsSensMP[key]["shapey"] = np.zeros(len(xDVs["shapey"]), "d")
                funcsSensMP[key]["shapey"] += funcsSens["CD"]["shapey"] * MPWeights[index]
        elif "mp%d_CL" % index in key:
            for alphaI in range(nMultiPoints):
                if alphaI == index:
                    funcsSensMP[key]["mp%d_alpha" % alphaI] = funcsSens["CL"]["mp%d_alpha" % index]
                else:
                    funcsSensMP[key]["mp%d_alpha" % alphaI] = np.zeros((keySize, 1), "d")
            funcsSensMP[key]["shapey"] = funcsSens["CL"]["shapey"]

    return


optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm
optFuncs.setMultiPointCondition = setMultiPointCondition
optFuncs.setMultiPointObjFuncs = setMultiPointObjFuncs
optFuncs.setMultiPointObjFuncsSens = setMultiPointObjFuncsSens

# =================================================================================================
# Task
# =================================================================================================
if task == "opt":

    DASolver.runColoring()
    optProb = Optimization("opt", optFuncs.calcObjFuncValuesMP, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    # Add objective
    optProb.addObj("obj", scale=1)
    # Add physical constraints
    for i in range(nMultiPoints):
        optProb.addCon("mp%d_CL" % i, lower=CL_target[i], upper=CL_target[i], scale=1)

    if gcomm.rank == 0:
        print(optProb)

    opt = OPT(args.opt, options=optOptions)
    histFile = os.path.join(outputDirectory, "%s_hist.hst" % args.opt)
    sol = opt(optProb, sens=optFuncs.calcObjFuncSensMP, storeHistory=histFile)
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
