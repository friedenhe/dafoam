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

    Info("\n")
    Info("+--------------------------------------------------------------------------+")
    Info("|                  Evaluating Objective Functions %03d                      |" % DASolver.nSolvePrimals)
    Info("+--------------------------------------------------------------------------+")
    Info("Design Variables: ")
    Info(xDV)

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
    Info("Objective Functions: ")
    Info(funcs)
    Info("Flow Runtime: %g" % (b - a))

    fail = funcs["fail"]

    return funcs, fail


def calcObjFuncValuesMP(xDV):
    """
    Update the design surface and run the primal solver to get objective function values.
    This is the multipoint version of calcObjFuncValues
    """

    Info("\n")
    Info("+--------------------------------------------------------------------------+")
    Info("|                  Evaluating Objective Functions %03d                      |" % DASolver.nSolvePrimals)
    Info("+--------------------------------------------------------------------------+")
    Info("Design Variables: ")
    Info(xDV)

    a = time.time()

    fail = False

    # Setup an empty dictionary for the evaluated function values
    funcsMP = {}

    # Set the current design variables in the DV object
    DVGeo.setDesignVars(xDV)
    DASolver.setDesignVars(xDV)

    nMultiPoints = DASolver.getOption("nMultiPoints")

    for i in range(nMultiPoints):

        Info("--Solving Primal for Configuration %d--" % i)

        funcs = {}

        # Evaluate the geometric constraints and add them to the funcs dictionary
        DVCon.evalFunctions(funcs)

        # set the multi point condition provided by users in the
        # runScript.py script. This function should define what
        # conditions to change for each case i
        setMultiPointCondition(xDV, i)

        # Solve the CFD problem
        DASolver()

        # Populate the required values from the CFD problem
        DASolver.evalFunctions(funcs, evalFuncs=evalFuncs)

        # save the state vector for case i and will be used in solveAdjoint
        DASolver.saveMultiPointField(i)

        # if any of the multipoint primal fails, return fail=True
        if funcs["fail"] is True:
            fail = True

        if DASolver.getOption("debug"):
            Info("Objective Functions for Configuration %d: " % i)
            Info(funcs)

        # assign funcs to funcsMP
        setMultiPointObjFuncs(funcs, funcsMP, i)

    funcsMP["fail"] = fail

    Info("Objective Functions MultiPoint: ")
    Info(funcsMP)

    b = time.time()
    Info("Flow Runtime: %g" % (b - a))

    return funcsMP, fail


def calcObjFuncSens(xDV, funcs):
    """
    Run the adjoint solver and get objective function sensitivities.
    """

    Info("\n")
    Info("+--------------------------------------------------------------------------+")
    Info("|              Evaluating Objective Function Sensitivities %03d             |" % DASolver.nSolveAdjoints)
    Info("+--------------------------------------------------------------------------+")

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
    Info("Objective Function Sensitivity: ")
    Info(funcsSens)
    Info("Adjoint Runtime: %g s" % (b - a))

    fail = funcsSens["fail"]

    return funcsSens, fail


def calcObjFuncSensMP(xDV, funcs):
    """
    Run the adjoint solver and get objective function sensitivities.
    This is the multipoint version of calcObjFuncSens
    """

    Info("\n")
    Info("+--------------------------------------------------------------------------+")
    Info("|              Evaluating Objective Function Sensitivities %03d             |" % DASolver.nSolveAdjoints)
    Info("+--------------------------------------------------------------------------+")

    fail = False

    a = time.time()

    # Setup an empty dictionary for the evaluated derivative values
    funcsSensMP = {}

    nMultiPoints = DASolver.getOption("nMultiPoints")

    for i in range(nMultiPoints):

        Info("--Solving Adjoint for Configuration %d--" % i)

        funcsSens = {}

        # Evaluate the geometric constraint derivatives
        DVCon.evalFunctionsSens(funcsSens)

        # set the state vector for case i
        DASolver.setMultiPointField(i)

        # set the multi point condition provided by users in the
        # runScript.py script. This function should define what
        # conditions to change for each case i
        setMultiPointCondition(xDV, i)

        # Solve the adjoint
        DASolver.solveAdjoint()
        DASolver.calcTotalDeriv()

        # Evaluate the CFD derivatives
        DASolver.evalFunctionsSens(funcsSens, evalFuncs=evalFuncs)

        if funcsSens["fail"] is True:
            fail = True

        if DASolver.getOption("debug"):
            Info("Objective Function Sensitivity: ")
            Info(funcsSens)

        # assign funcs to funcsMP
        setMultiPointObjFuncsSens(xDV, funcs, funcsSens, funcsSensMP, i)

    funcsSensMP["fail"] = fail

    # Print the current solution to the screen
    Info("Objective Function Sensitivity MultiPoiint: ")
    Info(funcsSensMP)

    b = time.time()
    Info("Adjoint Runtime: %g s" % (b - a))

    return funcsSensMP, fail


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
        Info("alpha: %f, CL: %f" % (alpha.real, CL0))
        if abs(CL0 - CL_star) / CL_star < 1e-5:
            Info("Completed! alpha = %f" % alpha.real)
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
    deltaX = DASolver.getOption("adjPartDerivFDStep")["FFD"]
    # initialize gradFD
    gradFD = {}
    for funcName in evalFuncs:
        gradFD[funcName] = {}
        for shapeVar in xDV:
            gradFD[funcName][shapeVar] = np.zeros(len(xDV[shapeVar]))
    if gcomm.rank == 0:
        print("-------FD----------", deltaX, flush=True)
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
            Info(gradFD)
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


class Info(object):
    """
    Print information and flush to screen for parallel cases
    """

    def __init__(self, message):
        if gcomm.rank == 0:
            print(message, flush=True)
        gcomm.Barrier()
