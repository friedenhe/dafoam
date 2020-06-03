/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DACheckMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
DACheckMesh::DACheckMesh(
    const Time& runTime1,
    const fvMesh& mesh1)
    : runTime(runTime1),
      mesh(mesh1)
{
}

DACheckMesh::~DACheckMesh()
{
}

label DACheckMesh::run() const
{
    /*
    Run checkMesh and return meshOK
    */

    label meshOK=1;

    Info << "Checking mesh quality for time = " << runTime.timeName() << nl << endl;

    label nFailedChecks = checkGeometry(mesh);

    if (nFailedChecks)
    {
        Info << "\nFailed " << nFailedChecks << " mesh checks.\n"
             << endl;
        meshOK = 0;
    }
    else
    {
        Info << "\nMesh OK.\n"
             << endl;
    }

    return meshOK;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //