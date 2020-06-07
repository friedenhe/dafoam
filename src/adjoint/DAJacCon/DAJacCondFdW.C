/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAJacCondFdW.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAJacCondFdW, 0);
addToRunTimeSelectionTable(DAJacCon, DAJacCondFdW, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAJacCondFdW::DAJacCondFdW(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAJacCon(modelType, mesh, daOption, daModel)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DAJacCondRdW::initializePetscVecs()
{

    //dFdW Colors, the dFdWColoredColumns will be initialized in initializedFdWCon
    VecCreate(PETSC_COMM_WORLD, &dFdWColors_);
    VecSetSizes(dFdWColors_, daIndex_.nLocalAdjointStates, PETSC_DECIDE);
    VecSetFromOptions(dFdWColors_);

    return;
}

void DAJacCon::initializeJacCon(const label isPC)
{
    /*
    Initialize DAJacCon::dFdWCon_

    Output:
    ------
    dFdWCon_: connectivity matrix for dFdW, here dFdWCon has a
    size of nLocalObjFuncGeoElements * nGlobalAdjointStates
    The reason that dFdWCon has nLocalObjFuncGeoElements rows is 
    because we need to divide the objective function into 
    nLocalObjFuncGeoElements discrete value such that we can
    use coloring to compute dFdW
    */

    MatCreate(PETSC_COMM_WORLD, &dFdWCon_);
    MatSetSizes(
        dFdWCon_,
        nLocalObjFuncGeoElements,
        daIndex_.nLocalAdjointStates,
        PETSC_DETERMINE,
        PETSC_DETERMINE);
    MatSetFromOptions(dFdWCon_);
    MatMPIAIJSetPreallocation(dFdWCon_, 100, NULL, 100, NULL);
    MatSeqAIJSetPreallocation(dFdWCon_, 100, NULL);
    MatSetOption(dFdWCon_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dFdWCon_);
    MatZeroEntries(dFdWCon_);

    VecCreate(PETSC_COMM_WORLD, &dFdWColoredColumns_);
    VecSetSizes(
        dFdWColoredColumns_,
        nLocalObjFuncGeoElements,
        PETSC_DECIDE);
    VecSetFromOptions(dFdWColoredColumns_);
    VecZeroEntries(dFdWColoredColumns_);

    Info << "dFdWCon Created!" << endl;
}

void DAJacCon::readdFdWColoring(const word objFunc)
{
    /*
    Read the dFdW coloring from files and 
    compute ndFdWColors. The naming convention for
    coloring vector is coloringVecName_nProcs.bin
    This is necessary because using different CPU
    cores result in different dRdWCon and therefore
    different coloring

    Output:
    ------
    dFdWColors_: read from file
    ndFdWColors: number of dFdW colors
    */

    Info << "Reading dFdW Coloring for " << objFunc << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dFdWColoring_" << objFunc << "_" << nProcs;
    word fileName1 = fileName.str();
    VecZeroEntries(dFdWColors_);
    daUtil_.readVectorBinary(dFdWColors_, fileName1);

    this->validateColoring(dFdWCon_, dFdWColors_);

    PetscReal maxVal;
    VecMax(dFdWColors_, NULL, &maxVal);
    ndFdWColors = maxVal + 1;
}

void DAJacCon::calcdFdWColoring(const word objFunc)
{
    /*
    Calculate the coloring for dFdW.

    Output:
    ------
    dFdWColors_: dFdW coloring and save to files. 
    The naming convention for coloring vector is 
    coloringVecName_nProcs.bin. This is necessary because 
    using different CPU cores result in different dRdWCon 
    and therefore different coloring

    ndFdWColors: number of dFdW colors

    */

    VecZeroEntries(dFdWColors_);
    this->parallelD2Coloring(dFdWCon_, dFdWColors_, ndFdWColors);
    this->validateColoring(dFdWCon_, dFdWColors_);
    Info << " ndFdWColors: " << ndFdWColors << endl;

    // write dFdW colors
    Info << "Writing dFdW Colors.." << endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream fileName("");
    fileName << "dFdWColoring_" << objFunc << "_" << nProcs;
    word fileName1 = fileName.str();
    daUtil_.writeVectorBinary(dFdWColors_, fileName1);

    return;
}

label DAJacCon::getNdFdWColors() const
{
    /*
    Return the number of colors

    Output:
    ------
    ndFdWColors: the number of colors depends on 
    whether the coloring is used
    */

    if (daOption_.getOption<label>("adjUseColoring"))
    {
        return ndFdWColors;
    }
    else
    {
        return daIndex_.nGlobalAdjointStates;
    }

    return -1;
}

} // End namespace Foam

// ************************************************************************* //
