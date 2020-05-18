/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/
#include "TestDAFoamIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
TestDAFoamIncompressible::TestDAFoamIncompressible(char* argsAll)
    : argsAll_(argsAll)
{
}

TestDAFoamIncompressible::~TestDAFoamIncompressible()
{
}

label TestDAFoamIncompressible::testDAUtility(PyObject* pyDict)
{
    label testErrors = 0;

    DAUtility daUtil;

    // ********************* pyDict2OFDict *********************
    dictionary ofDict;
    daUtil.pyDict2OFDict(pyDict, ofDict);

    // This is how ofDict should look like, it contains
    // all the supported types
    /*
    {
        key1            15;
        key2            5.5;
        key3            solver1;
        key4            0;
        key5            3 ( 1 2 3 );
        key6            3 ( 1.5 2.3 3.4 );
        key7            3 ( ele1 ele2 ele3 );
        key8            3 ( 0 1 0 );
        key9
        {
            subkey1         30;
            subkey2         3.5;
            subkey3         solver2;
            subkey4         1;
            subkey5         3 ( 4 5 6 );
            subkey6         3 ( 2.5 7.7 8.9 );
            subkey7         3 ( ele4 ele5 ele6 );
            subkey8         3 ( 1 0 1 );
        }
    }
    */

    // Now check if ofDict is properly set
    if (readLabel(ofDict.lookup("key1")) != 15)
    {
        Pout << "key1 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readScalar(ofDict.lookup("key2")) != 5.5)
    {
        Pout << "key2 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (word(ofDict.lookup("key3")) != "solver1")
    {
        Pout << "key3 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readLabel(ofDict.lookup("key4")) != 0)
    {
        Pout << "key4 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList list1;
    ofDict.readEntry("key5", list1);
    if (list1[0] != 1 || list1[1] != 2 || list1[2] != 3)
    {
        Pout << "key5 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    scalarList list2;
    ofDict.readEntry("key6", list2);
    if (list2[0] != 1.5 || list2[1] != 2.3 || list2[2] != 3.4)
    {
        Pout << "key6 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    wordList list3;
    ofDict.readEntry("key7", list3);
    if (list3[0] != "ele1" || list3[1] != "ele2" || list3[2] != "ele3")
    {
        Pout << "key7 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList list4;
    ofDict.readEntry("key8", list4);
    if (list4[0] != 0 || list4[1] != 1 || list4[2] != 0)
    {
        Pout << "key8 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    // check subDict
    dictionary subDict1;
    subDict1 = ofDict.subDict("key9");

    if (readLabel(subDict1.lookup("subkey1")) != 30)
    {
        Pout << "subkey1 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readScalar(subDict1.lookup("subkey2")) != 3.5)
    {
        Pout << "subkey2 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (word(subDict1.lookup("subkey3")) != "solver2")
    {
        Pout << "subkey3 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readLabel(subDict1.lookup("subkey4")) != 1)
    {
        Pout << "subkey4 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList subList1;
    subDict1.readEntry("subkey5", subList1);
    if (subList1[0] != 4 || subList1[1] != 5 || subList1[2] != 6)
    {
        Pout << "subkey5 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    scalarList subList2;
    subDict1.readEntry("subkey6", subList2);
    if (subList2[0] != 2.5 || subList2[1] != 7.7 || subList2[2] != 8.9)
    {
        Pout << "subkey6 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    wordList subList3;
    subDict1.readEntry("subkey7", subList3);
    if (subList3[0] != "ele4" || subList3[1] != "ele5" || subList3[2] != "ele6")
    {
        Pout << "subkey7 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList subList4;
    subDict1.readEntry("subkey8", subList4);
    if (subList4[0] != 1 || subList4[1] != 0 || subList4[2] != 1)
    {
        Pout << "key8 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (testErrors != 0)
    {
        Pout << "ofDict" << ofDict << endl;
    }

    // ********************* List functions *********************
    // isInList
    scalarList sList = {1.0, 2.3, 3.5};
    labelList lList = {1, 2, 3};
    wordList wList = {"e1", "e2", "e3"};
    if (!daUtil.isInList<scalar>(1.0, sList) || daUtil.isInList<scalar>(5.0, sList))
    {
        Pout << "scalar error in isInList!" << endl;
        testErrors += 1;
    }
    if (!daUtil.isInList<label>(1, lList) || daUtil.isInList<label>(4, lList))
    {
        Pout << "label error in isInList!" << endl;
        testErrors += 1;
    }
    if (!daUtil.isInList<word>("e1", wList) || daUtil.isInList<word>("e5", wList))
    {
        Pout << "word error in isInList!" << endl;
        testErrors += 1;
    }

    // listReplaceVal
    scalarList sListRe = {1.0, 2.3, 3.5};
    labelList lListRe = {1, 2, 3};
    wordList wListRe = {"e1", "e2", "e3"};

    scalarList sListReNew = {5.5, 2.3, 3.5};
    labelList lListReNew = {1, 5, 3};
    wordList wListReNew = {"e1", "e2", "e10"};

    daUtil.listReplaceVal<scalar>(sListRe, 1.0, 5.5);
    daUtil.listReplaceVal<label>(lListRe, 2, 5);
    daUtil.listReplaceVal<word>(wListRe, "e3", "e10");

    if (sListRe != sListReNew)
    {
        Pout << "scalar error in listReplaceVal!" << endl;
        testErrors += 1;
    }
    if (lListRe != lListReNew)
    {
        Pout << "label error in listReplaceVal!" << endl;
        testErrors += 1;
    }
    if (wListRe != wListReNew)
    {
        Pout << "word error in listReplaceVal!" << endl;
        testErrors += 1;
    }

    // listDeleteVal
    scalarList sListDel = {1.0, 2.3, 3.5, 1.0}; // delete multiple
    labelList lListDel = {1, 2, 3};
    wordList wListDel = {"e1", "e2", "e3"};

    scalarList sListDelNew = {2.3, 3.5};
    labelList lListDelNew = {1, 3};
    wordList wListDelNew = {"e1", "e2"};

    daUtil.listDeleteVal<scalar>(sListDel, 1.0);
    daUtil.listDeleteVal<label>(lListDel, 2);
    daUtil.listDeleteVal<word>(wListDel, "e3");

    if (sListDel != sListDelNew)
    {
        Pout << "scalar error in listDeleteVal!" << endl;
        testErrors += 1;
    }
    if (lListDel != lListDelNew)
    {
        Pout << "label error in listDeleteVal!" << endl;
        testErrors += 1;
    }
    if (wListDel != wListDelNew)
    {
        Pout << "word error in listDeleteVal!" << endl;
        testErrors += 1;
    }

    // ********************* MAT VEC IO *********************
    Vec tVec;
    label myProcNo = Pstream::myProcNo();
    VecCreate(PETSC_COMM_WORLD, &tVec);
    VecSetSizes(tVec, myProcNo, PETSC_DETERMINE);
    VecSetFromOptions(tVec);
    VecZeroEntries(tVec);

    PetscScalar* tVecArray;
    VecGetArray(tVec, &tVecArray);
    for (label i = 0; i < myProcNo; i++)
    {
        tVecArray[i] = i * 2.0;
    }

    VecRestoreArray(tVec, &tVecArray);

    VecAssemblyBegin(tVec);
    VecAssemblyEnd(tVec);

    daUtil.writeVectorBinary(tVec, "tVec");
    daUtil.writeVectorASCII(tVec, "tVec");

    Vec tVecRead;
    VecCreate(PETSC_COMM_WORLD, &tVecRead);
    daUtil.readVectorBinary(tVecRead, "tVec");

    PetscBool equalVec;
    VecEqual(tVec, tVecRead, &equalVec);
    if (!equalVec)
    {
        Pout << "error in read/writeVectorBinary!" << endl;
        testErrors += 1;
    }

    Mat tMat;
    MatCreate(PETSC_COMM_WORLD, &tMat);
    MatSetSizes(tMat, myProcNo, PETSC_DECIDE, PETSC_DETERMINE, 2);
    MatSetFromOptions(tMat);
    MatMPIAIJSetPreallocation(tMat, 2, NULL, 2, NULL);
    MatSeqAIJSetPreallocation(tMat, 2, NULL);
    MatSetUp(tMat);

    PetscInt Istart, Iend;
    MatGetOwnershipRange(tMat, &Istart, &Iend);
    for (PetscInt i = Istart; i < Iend; i++)
    {
        for (label j = 0; j < 2; j++)
        {
            scalar val = i*2.0+j;
            MatSetValue(tMat, i, j, val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(tMat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tMat,MAT_FINAL_ASSEMBLY);

    daUtil.writeMatrixBinary(tMat, "tMat");
    daUtil.writeMatrixASCII(tMat, "tMat");

    Mat tMatRead;
    MatCreate(PETSC_COMM_WORLD, &tMatRead);
    daUtil.readMatrixBinary(tMatRead, "tMat");

    PetscBool equalMat;
    MatEqual(tMat, tMatRead, &equalMat);
    if (!equalMat)
    {
        Pout << "error in read/writeMatrixBinary!" << endl;
        testErrors += 1;
    }

    return testErrors;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
