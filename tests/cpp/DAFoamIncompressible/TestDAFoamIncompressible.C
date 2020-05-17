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

label TestDAFoamIncompressible::test1(PyObject* pyDict)
{
    label testErrors = 0;

    DAUtility daUtil;

    // **************************** pyDict2OFDict ****************************
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
    Info << ofDict << endl;
    // Now check if ofDict is properly set
    if (readLabel(ofDict.lookup("key1")) != 15)
    {
        Info << "key1 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readScalar(ofDict.lookup("key2")) != 5.5)
    {
        Info << "key2 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (word(ofDict.lookup("key3")) != "solver1")
    {
        Info << "key3 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readLabel(ofDict.lookup("key4")) != 0)
    {
        Info << "key4 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList list1;
    ofDict.readEntry("key5", list1);
    if (list1[0] != 1 || list1[1] != 2 || list1[2] != 3)
    {
        Info << "key5 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    scalarList list2;
    ofDict.readEntry("key6", list2);
    if (list2[0] != 1.5 || list2[1] != 2.3 || list2[2] != 3.4)
    {
        Info << "key6 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    wordList list3;
    ofDict.readEntry("key7", list3);
    if (list3[0] != "ele1" || list3[1] != "ele2" || list3[2] != "ele3")
    {
        Info << "key7 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList list4;
    ofDict.readEntry("key8", list4);
    if (list4[0] != 0 || list4[1] != 1 || list4[2] != 0)
    {
        Info << "key8 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    // check subDict
    dictionary subDict1;
    subDict1 = ofDict.subDict("key9");

    if (readLabel(subDict1.lookup("subkey1")) != 30)
    {
        Info << "subkey1 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readScalar(subDict1.lookup("subkey2")) != 3.5)
    {
        Info << "subkey2 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (word(subDict1.lookup("subkey3")) != "solver2")
    {
        Info << "subkey3 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    if (readLabel(subDict1.lookup("subkey4")) != 1)
    {
        Info << "subkey4 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList subList1;
    subDict1.readEntry("subkey5", subList1);
    if (subList1[0] != 4 || subList1[1] != 5 || subList1[2] != 6)
    {
        Info << "subkey5 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    scalarList subList2;
    subDict1.readEntry("subkey6", subList2);
    if (subList2[0] != 2.5 || subList2[1] != 7.7 || subList2[2] != 8.9)
    {
        Info << "subkey6 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    wordList subList3;
    subDict1.readEntry("subkey7", subList3);
    if (subList3[0] != "ele4" || subList3[1] != "ele5" || subList3[2] != "ele6")
    {
        Info << "subkey7 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    labelList subList4;
    subDict1.readEntry("subkey8", subList4);
    if (subList4[0] != 1 || subList4[1] != 0 || subList4[2] != 1)
    {
        Info << "key8 error in pyDict2OFDict!" << endl;
        testErrors += 1;
    }

    // **************************** pyDict2OFDict ****************************

    return testErrors;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
