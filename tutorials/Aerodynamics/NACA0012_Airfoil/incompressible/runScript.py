#!/usr/bin/env python
"""
DAFoam run script for the NACA0012 airfoil at low-speed
"""

# =============================================================================
# Imports
# =============================================================================
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
# which optimizer to use. Options are: slsqp (default), snopt, or ipopt
parser.add_argument("--opt", help="optimizer to use", type=str, default="slsqp")
# which task to run. Options are: opt (default), run, testSensShape, or solveCL
parser.add_argument("--task", help="type of run to do", type=str, default="opt")
args = parser.parse_args()
gcomm = MPI.COMM_WORLD

# Define the global parameters here
# far field velocity
U0 = 10.0
# far field pressure, for incompressible p is relative to reference pressure
p0 = 0.0
# far field turbulence variable
nuTilda0 = 4.5e-5
# target lift coefficient to maintain during optimization
CL_target = 0.5
# initial angle of attack
alpha0 = 5.139186
# reference wing area for computing drag and lift coefficients
A0 = 0.1

# Set the parameters for optimization
daOptions = {
    # Primal solver setup
    # list of patch names for the design surface, need to be wall type in constant/polyMesh/boundary.
    # Here "wing" is a wall boundary patch defined in constant/polyMesh/boundary, we will change
    # the shape of the "wing" patch during optimization
    "designSurfaces": ["wing"],
    # name of the DASolver to run analysis and adjoint
    "solverName": "DASimpleFoam",
    # flow condition of the selected DASolver
    "flowCondition": "Incompressible",
    # turbulence model, need to be consistant with constant/turbulenceProperties
    "turbulenceModel": "SpalartAllmaras",
    # the residual convergence tolerance for the primal solver (DASimpleFoam)
    "primalMinResTol": 1.0e-8,
    # boundary condition for primal solution, if primalBC is defined, it will overwrite any defined
    # values in the "0" folder in OpenFOAM
    "primalBC": {
        # define the far field velocity for the "inout" patch, it will overwrite the value in 0/U
        "U0": {"variable": "U", "patch": "inout", "value": [U0, 0.0, 0.0]},
        # define the far field pressure for the "inout" patch, it will overwrite the value in 0/p
        "p0": {"variable": "p", "patch": "inout", "value": [p0]},
        # define the far field pressure for the "inout" patch, it will overwrite the value in 0/nuTilda
        # here we use the wallFunction
        "nuTilda0": {"variable": "nuTilda", "patch": "inout", "value": [nuTilda0], "useWallFunction": True},
    },
    # define the objective functions
    "objFunc": {
        # name of the objective function
        "CD": {
            # the first part of the objective function, most of the time the objective has only one part,
            # but one can also combine two parts of objectives, e.g., combining force and moment
            "part1": {
                # type of the objective
                "type": "force",
                # how to select the discrete mesh faces to compute the objective
                # we select them from the name of a patch
                "source": "patchToFace",
                # name of the patch for "patchToFace"
                "patches": ["wing"],
                # this is only for force objective because we need to project the force vector
                # to a specific direction. Here we say CD is the force that is parallel to the flow direction
                # Alternative, we can also use "fixedDirection" and provide a "direction" key for force:
                # "directionMode": "fixedDirection", "direction": [1.0, 0.0, 0.0]
                "directionMode": "parallelToFlow",
                # since we select "parallelToFlow" we need to prescribe the name of alpha (angle of attack)
                # variable to determine the flow direction. Here "alpha" is the name of design variable
                # that will be define later in:
                # DVGeo.addGeoDVGlobal("alpha", [alpha0], alpha, lower=-10.0, upper=10.0, scale=1.0)
                # NOTE: if no alpha is added in DVGeo.addGeoDVGlobal, we cant use "parallelToFlow",
                # In this case, use "directionMode": "fixedDirection" instead.
                "alphaName": "alpha",
                # the scaling factor for the force: CD = force / (0.5*U0*U0*A0)
                "scale": 1.0 / (0.5 * U0 * U0 * A0),
                # if addToAdjoint is True, the adjoint solver will compute the derivative for this objective
                # otherwise, it will only calculate the objective value for primal solver, no adjoint will
                # be computed for this objective
                "addToAdjoint": True,
            }
        },
        # the definition of CL is similar to CD except that we use "normalToFlow" for "directionMode"
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "normalToFlow",
                "alphaName": "alpha",
                "scale": 1.0 / (0.5 * U0 * U0 * A0),
                "addToAdjoint": True,
            }
        },
    },
    # Adjoint solver setup
    # adjoint linear equation solution options. If the adjoint does not converge, increase
    # pcFillLevel to 2. Or try "jacMatReOrdering" : "nd"
    "adjEqnOption": {"pcFillLevel": 1, "jacMatReOrdering": "rcm"},
    # state normalization value, we use far field value as reference
    # NOTE: since p is relative, we use the dynamic pressure U0*U0/2.
    # "phi" = 1.0 works for most of the case
    "normalizeStates": {"U": U0, "p": U0 * U0 / 2.0, "nuTilda": nuTilda0 * 10.0, "phi": 1.0},
    # the finite difference step size for computing partial derivatives in adjoint
    # the delta for state is typically 1e-8 to 1e-6 and the delta for the FFD point displacement
    # is typically LRef/1000
    "adjPartDerivFDStep": {"State": 1e-7, "FFD": 1e-3},
    # reserve an empty designVar key here and we will add values for it later in DVGeo.addGeoDVLocal
    "designVar": {},
}

# mesh warping parameters, users need to manually specify the symmetry plane and their normals
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane, for the airfoil case we have two symmetry planes
    # at z=0 and z=0.1 and the normal direction is [0.0, 0.0, 1.0]
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
}

# options for optimizers
if args.opt == "snopt":
    optOptions = {
        # tolerance for constraint
        "Major feasibility tolerance": 1.0e-7,
        # tolerance for gradient
        "Major optimality tolerance": 1.0e-7,
        # precision of objective function
        "Function precision": 1.0e-7,
        # do not verify gradients
        "Verify level": -1,
        # max optimization iteration to run
        "Major iterations limit": 50,
        # we do gradient-free line search
        "Nonderivative linesearch": None,
        # names of the optimization output file
        "Print file": "opt_SNOPT_print.out",
        "Summary file": "opt_SNOPT_summary.out",
    }
elif args.opt == "slsqp":
    optOptions = {
        # convergence accuracy
        "ACC": 1.0e-7,
        # max optimization iteration to run
        "MAXIT": 50,
        # name of the optimization output file
        "IFILE": "opt_SLSQP.out",
    }
elif args.opt == "ipopt":
    optOptions = {
        # convergence accuracy
        "tol": 1.0e-7,
        # max optimization iteration to run
        "max_iter": 50,
        # name of the optimization output file
        "output_file": "opt_IPOPT.out",
    }
else:
    print("opt arg not valid!")
    exit(0)


# =============================================================================
# Design variable setup
# =============================================================================
# define alpha function for angle of attack
# this function takes input "val" (can be an array) as angle of attack
# and change the "geo" object (output). Here geo is an class object to
# change the displacement of FFD points. For alpha function, we don't change
# the geometry or FFD so the geo object is not used.
def alpha(val, geo):
    # compute angle of attack in radian
    aoa = val[0] * np.pi / 180.0
    # compute the component of velocities based on U0 (global variable defined in the begining)
    # and the aoa
    inletU = [float(U0 * np.cos(aoa)), float(U0 * np.sin(aoa)), 0]
    # update the primalBC key in daOptions
    DASolver.setOption("primalBC", {"UIn": {"variable": "U", "patch": "inout", "value": inletU}})
    # apply the update to DAOption class in DAFoam
    DASolver.updateDAOption()


# create a DVGeo object to manipulate the design surface geometry using the free-form deformation (FFD)
# NOTE: the FFD volume should completely contain the design surface. The FFD file needs to be in Plot3D
# format (*.xyz), wingFFD.xyz is generated by running `python genFFD.py` in the FFD folder.
DVGeo = DVGeometry("./FFD/wingFFD.xyz")

# add reference axis for twist, in this case no twist is defined so the refAxis is not used
DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")

# select FFD point to move, the FFD file supports multi block meshes but in this case we have
# only one block in the FFD so iVol = 0
iVol = 0
pts = DVGeo.getLocalIndex(iVol)
# we allow all the points to move so we set pts[:, :, :]. Alternatively, you can select a subset
# of indices to move by setting a range for pts
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
# add local shape change variables "shapey" that moves in y direction with displacement bounds [-1.0:1.0]
# the points we select to move is provided by pointSelect=PS
DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
# set the design variable type for "shapey" (FFD)
daOptions["designVar"]["shapey"] = {"designVarType": "FFD"}
# add angle of attack variable "alpha" with bounds [0.0:10.0], the initial angle of attack is [alpha0]
# the the function to change angle of attack is provided as func=alpha (alpha is defined above)
DVGeo.addGeoDVGlobal("alpha", value=[alpha0], func=alpha, lower=0.0, upper=10.0, scale=1.0)
# set the design variable type for "alpha", which is "AOA". Here we also need to set at what patch
# the aoa is applied to (typically the name of far field patch) and the "flowAxis" and "normalAxis"
# for computing angle of attack: aoa = atan(U_normal/U_flow)
daOptions["designVar"]["alpha"] = {"designVarType": "AOA", "patch": "inout", "flowAxis": "x", "normalAxis": "y"}

# =============================================================================
# DAFoam initialization
# =============================================================================
# Users don't need to change this section
DASolver = PYDAFOAM(options=daOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
DASolver.printFamilyList()
DASolver.setMesh(mesh)
evalFuncs = []
DASolver.setEvalFuncs(evalFuncs)

# =============================================================================
# Constraint setup
# =============================================================================
# create DVCon object and set surface for constraints, these lines are common for most of applications
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily")))

# define leading (leList) and trailing (teList) edge lists for thickness and volume constraints
# NOTE: the leList and teList should be completely within the wing surface mesh
leList = [[1e-4, 0.0, 1e-4], [1e-4, 0.0, 0.1 - 1e-4]]
teList = [[0.998 - 1e-4, 0.0, 1e-4], [0.998 - 1e-4, 0.0, 0.1 - 1e-4]]

# volume constraint with bounds [1.0:3.0], here we use relative value with respect to the initial volume
# by setting scaled=True. The volume is computed based on the 2D mesh (2 x 10) constructed from the leList and teList
# nSpan = 2 and nChord = 10 means we use two points in the spanwise (z) and 10 points in the chordwise (x) to
# construct the 2D mesh, then we will project the 2D mesh upward and downward to the wing surface mesh
# and form a trapezoid volume to appoximate the wing volume. The more the leList and teList are close to the
# actual leading and trailing edges of the aifoil mesh, the better the volume approximation will be.
# Also, increasing the nSpan and nChord gives a better approximate. We recommend nSpan and nChord be similar
# to the number of FFD points in the spanwise and chordwise
DVCon.addVolumeConstraint(leList, teList, nSpan=2, nChord=10, lower=1.0, upper=3, scaled=True)

# thickness constraint with bounds [0.8:3.0], here we use relative value with respect to the initial thickness
# by setting scaled=True. Similar to the volume constraint, we sample 2 by 10 points and project them upward
# and downward to the wing surface to compute thickness at these 20 locations. We will have 20 thickness constraints
DVCon.addThicknessConstraints2D(leList, teList, nSpan=2, nChord=10, lower=0.8, upper=3.0, scaled=True)

# Create linear constraints to link the shape change between k=0 and k=1
# so that the shape changes are same in the spanwise direction, this is needed only for the airfoil case
# where we have two symmetry planes
nFFDs_x = pts.shape[0]  # number of FFD points in the x direction
indSetA = []
indSetB = []
# link shape change between k=0 and k =1 for all i and j indices
for i in range(nFFDs_x):
    for j in [0, 1]:
        indSetA.append(pts[i, j, 1])
        indSetB.append(pts[i, j, 0])
# here we impose: lower <= factorA * dy_{k=0} + factorB * dy_{k=1} <= upper
# Here dy is the displacement of FFD point in the y direction, and is defined in the DVGeo.addGeoDVLocal function
# Subsituting the parameters into the above equation, we have:
# 0 <= dy_{k=0} - dy_{k=1} <= 0
# In other words: dy_{k=0} = dy_{k=1}
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=-1.0, lower=0.0, upper=0.0)

# Create a linear constraint to fix the leading and trailing point. This is done by requiring the upper
# and lower FFD point on the leading and trailing edges to move in opposite directions.
indSetA = []
indSetB = []
for i in [0, nFFDs_x - 1]:
    for k in [0]:  # do not constrain k=1 because it is linked in the above symmetry constraint
        indSetA.append(pts[i, 0, k])
        indSetB.append(pts[i, 1, k])
# this imposes dy_{j=0} = - dy_{j=1}
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=1.0, lower=0.0, upper=0.0)

# =============================================================================
# Initialize optFuncs for optimization
# =============================================================================
# we need to assign the objects created in this script to optFuncs
optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# =============================================================================
# Task
# =============================================================================
# Run optimization
if args.task == "opt":

    # create an optimization problem and provide a function to compute objective (objFun=...)
    # here optFuncs.calcObjFuncValues is defined in optFuncs.py
    optProb = Optimization("opt", objFun=optFuncs.calcObjFuncValues, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    # add objective
    optProb.addObj("CD", scale=1)
    # add physical constraints with a lower and upper bounds
    # here we need to keep CL = CL_target during the optimization
    optProb.addCon("CL", lower=CL_target, upper=CL_target, scale=1)

    if gcomm.rank == 0:
        print(optProb)

    # run the coloring solver before running the optimization
    DASolver.runColoring()

    # create the optimization and provide a function to compute derivatives (sens=..)
    # here optFuncs.calcObjFuncSens is defined in optFuncs.py
    opt = OPT(args.opt, options=optOptions)
    histFile = "./%s_hist.hst" % args.opt
    sol = opt(optProb, sens=optFuncs.calcObjFuncSens, storeHistory=histFile)
    if gcomm.rank == 0:
        print(sol)

elif args.task == "run":

    # just run the primal and adjoint once
    optFuncs.run()

elif args.task == "solveCL":

    # change the angle of attack to meet the CL_target. This function is useful
    # to determine alpha0 for CL_target
    optFuncs.solveCL(CL_target, "alpha", "CL")

elif args.task == "testSensShape":

    # verify the sensitivity computed by adjoint
    optFuncs.testSensShape()

else:
    print("task arg not found!")
    exit(0)
