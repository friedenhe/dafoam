/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAUtility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAUtility::DAUtility()
{
}

DAUtility::~DAUtility()
{
}

void DAUtility::pyDict2OFDict(
    PyObject* pyDict,
    dictionary& ofDict)
{
    /*
    Description:
        Parse a Python dictionary to an OpenFOAM dictionary
        Only support a certain number of types and one nested dictionary
    Input:
        pyDict: Pytion dictionary
    Output:
        ofDict: OpenFoam dictionary
    */

    //PyObject_Print(pyDict,stdout,0);Info<<endl;

    // size of the pyOpions keys
    Py_ssize_t dictSize = PyDict_Size(pyDict);
    // all the keys
    PyObject* keys = PyDict_Keys(pyDict);
    // loop over all the keys in pyDict and assign their values
    // to ofDict
    for (label i = 0; i < dictSize; i++)
    {
        // the ith key
        PyObject* keyI = PyList_GetItem(keys, i);
        // convert it to UTF8 such that we can use it in C++
        char* keyUTF8 = PyUnicode_AsUTF8(keyI);

        //std::cout << "Key is "<<keyUTF8<<std::endl;
        // the actual value of this key, NOTE: it is a list
        // the 1st one is its type and the 2nd one is its value
        PyObject* value = PyDict_GetItem(pyDict, keyI);
        const char* valueTypeTmp = Py_TYPE(value)->tp_name;
        if (word(valueTypeTmp) != "list")
        {
            FatalErrorIn("pyDict2OFDict") << keyUTF8 << " needs to be in a list format with "
                                          << "the 1st element being its type and the 2nd element being its value."
                                          << "Example: " << keyUTF8 << ":[str,\"solverName\"]" << abort(FatalError);
        }

        // the second element of value is the actual value
        PyObject* value1 = PyList_GetItem(value, 1);

        //PyObject_Print(value1,stdout,0);Info<<endl;
        // get the type of value1
        const char* valueType = Py_TYPE(value1)->tp_name;
        if (word(valueType) == "str")
        {
            char* valSet = PyUnicode_AsUTF8(value1);
            ofDict.add(keyUTF8, word(valSet));
        }
        else if (word(valueType) == "int")
        {
            long valSet = PyLong_AsLong(value1);
            ofDict.add(keyUTF8, label(valSet));
        }
        else if (word(valueType) == "bool")
        {
            label valSet = PyObject_IsTrue(value1);
            ofDict.add(keyUTF8, valSet);
        }
        else if (word(valueType) == "float")
        {
            scalar valSet = PyFloat_AS_DOUBLE(value1);
            ofDict.add(keyUTF8, valSet);
        }
        else if (word(valueType) == "list")
        {
            // size of the list
            Py_ssize_t listSize = PyList_Size(value1);

            //Info<<listSize<<endl;
            // create OpenFOAM lists to hold this list
            // we create all the possible types
            scalarList valSetScalar;
            valSetScalar.setSize(label(listSize));
            labelList valSetLabel;
            valSetLabel.setSize(label(listSize));
            List<word> valSetWord;
            valSetWord.setSize(label(listSize));

            // need to check what type of list this is
            // by checking its first element
            PyObject* tmp = PyList_GetItem(value1, 0);
            const char* tmpType = Py_TYPE(tmp)->tp_name;
            word tmpTypeWord = word(tmpType);
            // assign value to the OpenFOAM list
            for (label j = 0; j < listSize; j++)
            {
                PyObject* valueListJ = PyList_GetItem(value1, j);
                if (tmpTypeWord == "str")
                {
                    char* valSet = PyUnicode_AsUTF8(valueListJ);
                    valSetWord[j] = word(valSet);
                }
                else if (tmpTypeWord == "int")
                {
                    long valSet = PyLong_AsLong(valueListJ);
                    valSetLabel[j] = label(valSet);
                }
                else if (tmpTypeWord == "float")
                {
                    scalar valSet = PyFloat_AS_DOUBLE(valueListJ);
                    valSetScalar[j] = valSet;
                }
                else if (tmpTypeWord == "bool")
                {
                    label valSet = PyObject_IsTrue(valueListJ);
                    valSetLabel[j] = valSet;
                }
                else
                {
                    FatalErrorIn("pyDict2OFDict") << "Type: " << tmpTypeWord << " for " << keyUTF8
                                                  << " list is not supported! Options are: str, int, bool, and float!"
                                                  << abort(FatalError);
                }
            }

            // add the list to the ofDict dict
            if (tmpTypeWord == "str")
            {
                ofDict.add(keyUTF8, valSetWord);
            }
            else if (tmpTypeWord == "int")
            {
                ofDict.add(keyUTF8, valSetLabel);
            }
            else if (tmpTypeWord == "float")
            {
                ofDict.add(keyUTF8, valSetScalar);
            }
            else if (tmpTypeWord == "bool")
            {
                ofDict.add(keyUTF8, valSetLabel);
            }
        }
        else if (word(valueType) == "dict")
        {
            dictionary subDict;
            this->pyDict2OFDict(value1, subDict);
            ofDict.add(keyUTF8, subDict);
        }
        else
        {
            FatalErrorIn("pyDict2OFDict") << "Type: " << valueType << " for " << keyUTF8
                                          << " is not supported! Options are: str, int, float, bool, list, and dict!"
                                          << abort(FatalError);
        }
        //std::cout << "My type is " << valueType << std::endl;
    }
}

void DAUtility::readVectorBinary(
    Vec vecIn,
    const word prefix) const
{
    /*
    Description:
        Read a vector in binary form
    Input:
        vecIn: a Petsc vector to read values into (also output)
        prefix: Name of the Petsc vector from disk
    Example:
        If the vector storing in the disk reads: dFdWVector.bin
        Then read the vector using:
        Vec dFdW;
        codes to initialize vector dFdW....
        readVectorBinary(dFdW,"dFdwVector");
        NOTE: the prefix does not include ".bin"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".bin";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(vecIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

void DAUtility::writeVectorBinary(
    const Vec vecIn,
    const word prefix) const
{
    /*
    Description:
        Write a vector in binary form
    Input:
        vecIn: a Petsc vector to write to disk
        prefix: Name of the Petsc vector to write to disk
    Example:
        Vec dFdW;
        codes to initialize vector dFdW....
        writeVectorBinary(dFdW,"dFdwVector");
        This will write the dFdW vector to disk with name "dFdWVector.bin"
        NOTE: the prefix does not include ".bin"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".bin";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(vecIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

void DAUtility::writeVectorASCII(
    const Vec vecIn,
    const word prefix) const
{
    /*
    Description:
        Write a vector in ASCII form
    Input:
        vecIn: a Petsc vector to write to disk
        prefix: Name of the Petsc vector to write to disk
    Example:
        Vec dFdW;
        codes to initialize vector dFdW....
        writeVectorASCII(dFdW,"dFdwVector");
        This will write the dFdW vector to disk with name "dFdWVector.dat"
        NOTE: the prefix does not include ".dat"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".dat";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, fileName.c_str(), &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); // write all the digits
    VecView(vecIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

void DAUtility::readMatrixBinary(
    Mat matIn,
    const word prefix) const
{
    /*
    Description:
        Read a matrix in binary form
    Input:
        matIn: a Petsc matrix to read values into (also output)
        prefix: Name of the Petsc matrix from disk
    Example:
        If the matrix storing in the disk reads: dRdWMat.bin
        Then read the matrix using:
        Mat dRdW;
        codes to initialize matrix dRdW....
        readMatrixBinary(dRdW,"dRdwVector");
        NOTE: the prefix does not include ".bin"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".bin";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer);
    MatLoad(matIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

void DAUtility::writeMatrixBinary(
    const Mat matIn,
    const word prefix) const
{
    /*
    Description:
        Write a matrix in binary form
    Input:
        matIn: a Petsc matrix to write to disk
        prefix: Name of the Petsc matrix to write to disk
    Example:
        Mat dRdW;
        codes to initialize matrix dRdW....
        writeMatrixBinary(dRdW,"dRdWMat");
        This will write the dRdW matrix to disk with name "dRdWMat.bin"
        NOTE: the prefix does not include ".bin"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".bin";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer);
    MatView(matIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

void DAUtility::writeMatrixASCII(
    const Mat matIn,
    const word prefix) const
{
    /*
    Description:
        Write a matrix in ASCII form
    Input:
        matIn: a Petsc matrix to write to disk
        prefix: Name of the Petsc matrix to write to disk
    Example:
        Mat dRdW;
        codes to initialize matrix dRdW....
        writeMatrixASCII(dRdW,"dRdWMat");
        This will write the dRdW matrix to disk with name "dRdWMat.dat"
        NOTE: the prefix does not include ".dat"
    */

    std::ostringstream fileNameStream("");
    fileNameStream << prefix << ".dat";
    word fileName = fileNameStream.str();

    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, fileName.c_str(), &viewer);
    MatView(matIn, viewer);
    PetscViewerDestroy(&viewer);

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
