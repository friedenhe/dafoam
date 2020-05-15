/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

Info << "Create mesh for time = " << runTime.timeName() << endl;

word regionName = fvMesh::defaultRegion;
meshPtr_.reset(
    new fvMesh(
        IOobject(
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ)));

fvMesh& mesh = meshPtr_();
