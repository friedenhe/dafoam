/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DAOption::DAOption(
    const fvMesh& mesh,
    PyObject* pyOptions)
    : regIOobject(
        IOobject(
            "DAOption", // always use DAOption for the db name
            mesh.time().timeName(),
            mesh, // register to mesh
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true // always register object
            )),
      mesh_(mesh)
{
    /*
    Description:
        Construct from fvMesh and a dict object from Python
    Input:
        mesh: Foam::fvMesh object
        pyOptions: dictionary from Python, which contains all
        options for DAFoam 
    */

    // now we need to convert the pyOptions<PyObject*> to allOptions_<dictionary> in OpenFOAM
    this->setAllOptions(pyOptions, allOptions_);

    Info << "All DAFoam Options:";
    Info << allOptions_ << endl;
}

DAOption::~DAOption()
{
}

void DAOption::setAllOptions(
    PyObject* pyOptions,
    dictionary& allOptions)
{
    /*
    Description:
        Assign the OpenFOAM dictionary allOptions based on a
        Python dictionary pyOptoins. Here we parse a Python dictionary
        to an OpenFOAM dictionary
    Input:
        pyOptions: Pytion dictionary
    Output:
        allOptions: OpenFoam dictionary
    */

    //PyObject_Print(pyOptions,stdout,0);Info<<endl;

    // size of the pyOpions keys
    Py_ssize_t dictSize = PyDict_Size(pyOptions);
    // all the keys
    PyObject* keys = PyDict_Keys(pyOptions);
    // loop over all the keys in pyOptions and assign their values
    // to allOptions
    for (label i = 0; i < dictSize; i++)
    {
        // the ith key
        PyObject* keyI = PyList_GetItem(keys, i);
        // convert it to UTF8 such that we can use it in C++
        char* keyUTF8 = PyUnicode_AsUTF8(keyI);

        //std::cout << "Key is "<<keyUTF8<<std::endl;
        // the actual value of this key, NOTE: it is a list
        // the 1st one is its type and the 2nd one is its value
        PyObject* value = PyDict_GetItem(pyOptions, keyI);
        const char* valueTypeTmp = Py_TYPE(value)->tp_name;
        if (word(valueTypeTmp) != "list")
        {
            FatalErrorIn("setAllOptions") << keyUTF8 << " needs to be in a list format with "
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
            allOptions.add(keyUTF8, word(valSet));
        }
        else if (word(valueType) == "int")
        {
            long valSet = PyLong_AsLong(value1);
            allOptions.add(keyUTF8, label(valSet));
        }
        else if (word(valueType) == "bool")
        {
            label valSet = PyObject_IsTrue(value1);
            allOptions.add(keyUTF8, valSet);
        }
        else if (word(valueType) == "float")
        {
            scalar valSet = PyFloat_AS_DOUBLE(value1);
            allOptions.add(keyUTF8, valSet);
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
                    FatalErrorIn("setAllOptions") << "Type: " << tmpTypeWord << " for " << keyUTF8
                                                  << " list is not supported! Options are: str, int, bool, and float!"
                                                  << abort(FatalError);
                }
            }

            // add the list to the allOptions dict
            if (tmpTypeWord == "str")
            {
                allOptions.add(keyUTF8, valSetWord);
            }
            else if (tmpTypeWord == "int")
            {
                allOptions.add(keyUTF8, valSetLabel);
            }
            else if (tmpTypeWord == "float")
            {
                allOptions.add(keyUTF8, valSetScalar);
            }
            else if (tmpTypeWord == "bool")
            {
                allOptions.add(keyUTF8, valSetLabel);
            }
        }
        else if (word(valueType) == "dict")
        {
            dictionary subDict;
            this->setAllOptions(value1, subDict);
            allOptions.add(keyUTF8, subDict);
        }
        else
        {
            FatalErrorIn("setAllOptions") << "Type: " << valueType << " for " << keyUTF8
                                          << " is not supported! Options are: str, int, float, bool, list, and dict!"
                                          << abort(FatalError);
        }
        //std::cout << "My type is " << valueType << std::endl;
    }
}

// this is a virtual function for regIOobject
bool DAOption::writeData(Ostream& os) const
{
    // do nothing
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
