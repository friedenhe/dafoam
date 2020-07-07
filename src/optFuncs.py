#!/usr/bin/env python
"""

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
    Common functions for DAFoam optimization setup.

"""


# =============================================================================
# Imports
# =============================================================================
import time
import sys
import numpy as np

np.set_printoptions(precision=16)


def calcObjFuncValues(xDV):
    """
    Update the design surface and run the primal solver to get objective function values.
    """

    if gcomm.rank == 0:
        print("\n")
        print("+--------------------------------------------------------------------------+")
        print("|                  Evaluating Objective Functions %03d                      |" % DASolver.nSolvePrimals)
        print("+--------------------------------------------------------------------------+")
        print("Design Variables: ", xDV)

    a = time.time()

    # Setup an empty dictionary for the evaluated function values
    funcs = {}

    # Set the current design variables in the DV object
    DVGeo.setDesignVars(xDV)
    DASolver.setDesignVars(xDV)

    # Evaluate the geometric constraints and add them to the funcs dictionary
    DVCon.evalFunctions(funcs)

    # Solve the CFD problem
    DASolver()

    # Populate the required values from the CFD problem
    DASolver.evalFunctions(funcs, evalFuncs=evalFuncs)

    b = time.time()

    # Print the current solution to the screen
    if gcomm.rank == 0:
        print("Objective Functions: ", funcs)
        print("Flow Runtime: ", b - a)

    fail = funcs["fail"]

    # flush the output to the screen/file
    sys.stdout.flush()

    return funcs, fail


def calcObjFuncSens(xDV, funcs):
    """
    Run the adjoint solver and get objective function sensitivities.
    """

    if gcomm.rank == 0:
        print("\n")
        print("+--------------------------------------------------------------------------+")
        print(
            "|              Evaluating Objective Function Sensitivities %03d             |" % DASolver.nSolveAdjoints
        )
        print("+--------------------------------------------------------------------------+")

    a = time.time()

    # Setup an empty dictionary for the evaluated derivative values
    funcsSens = {}

    # Evaluate the geometric constraint derivatives
    DVCon.evalFunctionsSens(funcsSens)

    # Solve the adjoint
    DASolver.solveAdjoint()
    DASolver.calcTotalDeriv()

    # Evaluate the CFD derivatives
    DASolver.evalFunctionsSens(funcsSens, evalFuncs=evalFuncs)

    b = time.time()

    # Print the current solution to the screen
    if gcomm.rank == 0:
        print("Objective Function Sensitivity: ", funcsSens)
        print("Adjoint Runtime: ", b - a)

    fail = funcsSens["fail"]

    # flush the output to the screen/file
    sys.stdout.flush()

    return funcsSens, fail


def run():
    """
    Just run the primal and adjoint
    """

    DASolver.runColoring()
    xDV = DVGeo.getValues()
    funcs = {}
    funcs, fail = calcObjFuncValues(xDV)
    funcsSens = {}
    funcsSens, fail = calcObjFuncSens(xDV, funcs)


def solveCL(CL_star, alphaName, liftName):

    DASolver.setOption("adjUseColoring", False)

    xDVs = DVGeo.getValues()
    alpha = xDVs[alphaName]

    for i in range(10):
        # Solve the CFD problem
        xDVs[alphaName] = alpha
        funcs = {}
        funcs, fail = calcObjFuncValues(xDVs)
        CL0 = funcs[liftName]
        if gcomm.rank == 0:
            print("alpha: %f, CL: %f" % (alpha.real, CL0))
        if abs(CL0 - CL_star) / CL_star < 1e-5:
            if gcomm.rank == 0:
                print("Completed! alpha = %f" % alpha.real)
            break
        # compute sens
        eps = 1e-2
        alphaVal = alpha + eps
        xDVs[alphaName] = alphaVal
        funcsP = {}
        funcsP, fail = calcObjFuncValues(xDVs)
        CLP = funcsP[liftName]
        deltaAlpha = (CL_star - CL0) * eps / (CLP - CL0)
        alpha += deltaAlpha


def testSensShape():
    """
    Verify the FFD sensitivity against finite-difference references
    """

    DASolver.runColoring()

    xDV = DVGeo.getValues()

    if gcomm.rank == 0:
        fOut = open("./testFFDSens.txt", "w")

    # gradAdj

    funcs = {}
    funcsSens = {}
    funcs, fail = calcObjFuncValues(xDV)
    funcsSens, fail = calcObjFuncSens(xDV, funcs)
    if gcomm.rank == 0:
        for funcName in evalFuncs:
            for shapeVar in xDV:
                fOut.write(funcName + " " + shapeVar + "\n")
                try:
                    nDVs = len(funcsSens[funcName][shapeVar])
                except Exception:
                    nDVs = 1
                for n in range(nDVs):
                    line = str(funcsSens[funcName][shapeVar][n]) + "\n"
                    fOut.write(line)
                    fOut.flush()

    # gradFD
    deltaX = DASolver.getOption("adjEpsDerivFFD")
    # initialize gradFD
    gradFD = {}
    for funcName in evalFuncs:
        gradFD[funcName] = {}
        for shapeVar in xDV:
            gradFD[funcName][shapeVar] = np.zeros(len(xDV[shapeVar]))
    if gcomm.rank == 0:
        print("-------FD----------", deltaX)
        fOut.write("DeltaX: " + str(deltaX) + "\n")
    for shapeVar in xDV:
        try:
            nDVs = len(xDV[shapeVar])
        except Exception:
            nDVs = 1
        for i in range(nDVs):
            funcp = {}
            funcm = {}
            xDV[shapeVar][i] += deltaX
            funcp, fail = calcObjFuncValues(xDV)
            xDV[shapeVar][i] -= 2.0 * deltaX
            funcm, fail = calcObjFuncValues(xDV)
            xDV[shapeVar][i] += deltaX
            for funcName in evalFuncs:
                gradFD[funcName][shapeVar][i] = (funcp[funcName] - funcm[funcName]) / (2.0 * deltaX)
            if gcomm.rank == 0:
                print(gradFD)
    # write FD results
    if gcomm.rank == 0:
        for funcName in evalFuncs:
            for shapeVar in xDV:
                fOut.write(funcName + " " + shapeVar + "\n")
                try:
                    nDVs = len(gradFD[funcName][shapeVar])
                except Exception:
                    nDVs = 1
                for n in range(nDVs):
                    line = str(gradFD[funcName][shapeVar][n]) + "\n"
                    fOut.write(line)
                    fOut.flush()

    if gcomm.rank == 0:
        fOut.close()
