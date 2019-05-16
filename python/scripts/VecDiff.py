#!/usr/bin/env python
'''
Compare values in two different vectors
Example:
python VecDiff.py dFdW1.bin dFdW2.bin
'''
import os,sys
import argparse
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc


def evalVecDiff(vecName1,vecName2,diffTol=1e-16):
    # read  vec 1
    vec1 = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    viewer = PETSc.Viewer().createBinary(vecName1,comm=PETSc.COMM_WORLD)
    vec1.load(viewer)
    
    # read  vec 2
    vec2 = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    viewer = PETSc.Viewer().createBinary(vecName2,comm=PETSc.COMM_WORLD)
    vec2.load(viewer)
    
    # read  vec Diff
    vecDiff = PETSc.Vec().create(comm=PETSc.COMM_WORLD)
    viewer = PETSc.Viewer().createBinary(vecName1,comm=PETSc.COMM_WORLD)
    vecDiff.load(viewer)
    
    Istart, Iend = vec1.getOwnershipRange()
    
    vecDiff.axpy(-1.0, vec2)
    
    maxDiff = -1.0e12
    maxVal1 = -1.0e12
    maxVal2 = -1.0e12
    maxRowI = -1e12
    l2norm=0.0
    foundDiff=0
    for i in xrange(Istart, Iend):
        valDiff=abs(vecDiff[i])
        l2norm = l2norm + valDiff**2
        if valDiff>diffTol:
            if valDiff>maxDiff:
                maxDiff=valDiff
                maxRowI=i
                maxVal1 = vec1.getValue(i)
                maxVal2 = vec2.getValue(i)
                foundDiff=1
    l2norm = l2norm**0.5
    
    if foundDiff==1:
        maxDiffRel = maxDiff/abs(maxVal1+1e-16) # relative value for the maxDiff
        print('L2Norm: %20.16e'%l2norm)
        print('MaxDiff: %20.16e' %maxDiff)
        print('MaxRelD: %20.16e' %maxDiffRel)
        print('MaxVal1: %20.16e' %maxVal1)
        print('MaxVal2: %20.16e' %maxVal2)
        print('MaxrowI: %d' %maxRowI)
        return maxDiff,maxDiffRel
    else:
        print('Two vectors are exactly same with tolerance: %e'%diffTol)
        return 0.0,0.0

if __name__ == '__main__':
    print("\nUsage: python VecDiff.py vecName1 vecName2")
    print("Example python VecDiff.py dFdW1.bin dFdW2.bin\n")
    evalVecDiff(sys.argv[1],sys.argv[2])
