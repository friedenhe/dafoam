/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Child class for ActuatorDisk disk source

\*---------------------------------------------------------------------------*/

#ifndef DASourceActuatorDisk_H
#define DASourceActuatorDisk_H

#include "DASource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DASourceActuatorDisk Declaration
\*---------------------------------------------------------------------------*/

class DASourceActuatorDisk
    : public DASource
{

protected:

    /// HashTable that contains lists of cell indices that are within the actuator disk space
    HashTable<labelList> sourceCellIndices_;

    /// calculate DASourceActuatorDisk::sourceCellIndices_
    void calcSourceCellIndices(HashTable<labelList>& sourceCellIndices);

    /// the subDict of fvSource in DAOption
    dictionary fvSourceSubDict_;

public:
    TypeName("actuatorDisk");
    // Constructors

    //- Construct from components
    DASourceActuatorDisk(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DASourceActuatorDisk()
    {
    }

    /// compute the source term
    virtual void calcSource(volVectorField& source);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
