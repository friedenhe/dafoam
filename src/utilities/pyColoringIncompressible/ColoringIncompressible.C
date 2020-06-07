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

    // a dictionary to pass the parameters int daJacCon
    dictionary options;
    label isPC = 0;
    label isPrealloc = 1;

    // need to first set isPrealloc to true to calcualte the preallocation for the dRdWCon matrix
    // because directly initializing the dRdWCon matrix will use too much memory
    options.set("isPC", isPC);
    options.set("isPrealloc", isPrealloc);
    daJacCon->setupJacCon(options);

    // now we can initilaize dRdWCon
    daJacCon->initializeJacCon(options);

    // setup dRdWCon
    isPrealloc = 0;
    options.set("isPrealloc", isPrealloc);
    daJacCon->setupJacCon(options);
    Info << "dRdWCon Created. " << mesh.time().elapsedClockTime() << " s" << endl;

    // compute the coloring
    Info << "Calculating dRdW Coloring... " << mesh.time().elapsedClockTime() << " s" << endl;
    daJacCon->calcJacConColoring();
    Info << "Calculating dRdW Coloring... Completed! " << mesh.time().elapsedClockTime() << " s" << endl;

    // clean up
    daJacCon->deleteJacCon(options);

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //