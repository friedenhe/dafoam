/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        DAUtility contains basic functions such as matrix file IO, and
        it is independent of fvMesh.

\*---------------------------------------------------------------------------*/

#ifndef DAUtility_H
#define DAUtility_H

#include <petscksp.h>
#include "Python.h"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAUtility Declaration
\*---------------------------------------------------------------------------*/
class DAUtility
{

private:
    /// Disallow default bitwise copy construct
    DAUtility(const DAUtility&);

    /// Disallow default bitwise assignment
    void operator=(const DAUtility&);

public:
    /// Constructors
    DAUtility();

    /// Destructor
    virtual ~DAUtility();

    /// convert a python dictionary object to OpenFoam dictionary
    void pyDict2OFDict(
        PyObject* pyDict,
        dictionary& ofDict);

    /// check whether a value is in the list
    template<class classType>
    label isInList(
        const classType value,
        const List<classType>& list) const;

    /// replace a value in the list
    template<class classType>
    label listReplaceVal(
        List<classType>& listIn,
        const classType valOrig,
        const classType valNew) const;

    /// delete a value in the list
    template<class classType>
    label listDeleteVal(
        List<classType>& listIn,
        const classType valDel) const;

    /// write petsc matrix in binary format
    void writeMatrixBinary(
        const Mat matIn,
        const word prefix) const;

    /// write petsc matrix in ascii format
    void writeMatrixASCII(
        const Mat matIn,
        const word prefix) const;

    /// read petsc matrix in binary format
    void readMatrixBinary(
        Mat matIn,
        const word prefix) const;

    /// write petsc vector in ascii format
    void writeVectorASCII(
        const Vec vecIn,
        const word prefix) const;

    /// read petsc vector in binary format
    void readVectorBinary(
        Vec vecIn,
        const word prefix) const;

    /// write petsc vector in binary format
    void writeVectorBinary(
        const Vec vecIn,
        const word prefix) const;

    /// bound a volScalar variable based on parametes defined in DAOption::allOptions_
    void boundVar(
        const dictionary& allOptions,
        volScalarField& var);

    /// bound a volVector variable based on parametes defined in DAOption::allOptions_
    void boundVar(
        const dictionary& allOptions,
        volVectorField& var);
};

template<class classType>
label DAUtility::isInList(
    const classType value,
    const List<classType>& list) const
{
    /*
    Description:
        Check whether a value is in the list
    Input:
        value: a value to check with type: classType
        list: list to check the value
    Example:
        If the list reads
        list1={"test","run","compare"};
        label val=isInList<word>("run",list1);
        Here the return value: val==1 (value found in the list)
    */

    forAll(list, idxI)
    {
        if (list[idxI] == value)
        {
            return 1;
        }
    }
    return 0;
}

template<class classType>
label DAUtility::listReplaceVal(
    List<classType>& listIn,
    const classType valOrig,
    const classType valNew) const
{
    /*
    Description:
        Replace a value in hte list, the list sequence is preserved
    Input:
        listIn: the list to replace
        valOrig: original value in the list
        valNew: new value to replace
    Example:
        List<word> listIn={"apple","orange"};
        listReplaceVal<word>(listIn,"apple","banana");
        Now listIn will be {"banana","orange"}
        NOTE: if list has multiple valOrig, they will be all replaced
    */

    label foundVal = 0;
    forAll(listIn, idxI)
    {
        const classType& val = listIn[idxI];
        if (val == valOrig)
        {
            listIn[idxI] = valNew;
            foundVal = 1;
        }
    }

    return foundVal;
}

template<class classType>
label DAUtility::listDeleteVal(
    List<classType>& listIn,
    const classType valDel) const
{
    /*
    Description:
        Delete a value in the list, the sequence of other
        elements will be preserved
    Input:
        listIn: the list to delete value
        valDel: value to delete
    Example:
        List<word> listIn={"apple","orange"};
        listDeleteVal<word>(listIn,"apple");
        Now listIn will be {orange"}
        NOTE: if list has multiple val, they will be all deleted
    */

    // we will create a new list and delete the value
    List<classType> listNew;

    label foundVal = 0;
    forAll(listIn, idxI)
    {
        const classType& val = listIn[idxI];
        if (val == valDel)
        {
            foundVal = 1;
            // do not append
        }
        else
        {
            listNew.append(val);
        }
    }

    listIn.clear();
    listIn = listNew;

    return foundVal;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
