/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/
#include "TestDAFoamCompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
TestDAFoamCompressible::TestDAFoamCompressible(char* argsAll)
    : argsAll_(argsAll)
{
}

TestDAFoamCompressible::~TestDAFoamCompressible()
{
}

label TestDAFoamCompressible::testDARegState(PyObject* pyDict)
{
#include "setArgs.H"
#include "setRootCasePython.H"
#include "createTime.H"
#include "createMesh.H"

    label testErrors = 0;

    DAOption daOption(mesh, pyDict);

    autoPtr<DARegState> daRegState(DARegState::New(mesh));

    const HashTable<wordList>& regStates = daRegState->getRegStates();

    HashTable<wordList> regStatesRef;

    regStatesRef.set("volScalarField", {});
    regStatesRef.set("volVectorField", {});
    regStatesRef.set("surfaceScalarField", {});
    regStatesRef.set("surfaceVectorField", {});
    regStatesRef["volScalarField"].append("p");
    regStatesRef["volScalarField"].append("T");
    regStatesRef["volScalarField"].append("nut");
    regStatesRef["volVectorField"].append("U");
    regStatesRef["surfaceScalarField"].append("phi");

    if( regStates != regStatesRef)
    {
        Pout << "compressible error in DARegState!" << endl;
        testErrors += 1;
    }

    return testErrors;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
