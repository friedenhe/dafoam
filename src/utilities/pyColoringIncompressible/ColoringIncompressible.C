/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "ColoringIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
ColoringIncompressible::ColoringIncompressible(
    char* argsAll,
    PyObject* pyOptions)
    : argsAll_(argsAll),
      pyOptions_(pyOptions),
      argsPtr_(nullptr),
      runTimePtr_(nullptr),
      meshPtr_(nullptr)
{
}

ColoringIncompressible::~ColoringIncompressible()
{
}

void ColoringIncompressible::run()
{
#include "setArgs.H"
#include "setRootCasePython.H"
#include "createTimePython.H"
#include "createMeshPython.H"
#include "createFields.H"
#include "createAdjoint.H"

    label isPrealloc = 1;
    label isPC = 0;
    daJacCon->setupdRdWCon(isPrealloc, isPC);
    daJacCon->initializedRdWCon(isPC);
    isPrealloc = 0;
    daJacCon->setupdRdWCon(isPrealloc, isPC);
    Info << "dRdWCon Created. " << mesh.time().elapsedClockTime() << " s" << endl;
    Info << "Calculating dRdW Coloring... " << mesh.time().elapsedClockTime() << " s" << endl;
    daJacCon->calcdRdWColoring();
    Info << "Calculating dRdW Coloring... Completed! " << mesh.time().elapsedClockTime() << " s" << endl;
    daJacCon->deletedRdWCon();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //