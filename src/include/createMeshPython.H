/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    This file is modified from OpenFOAM's source code
    src/OpenFOAM/include/createMesh.H

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

        OpenFOAM is free software: you can redistribute it and/or modify it
        under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
        for more details.
    
        You should have received a copy of the GNU General Public License
        along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


\*---------------------------------------------------------------------------*/

Info << "Create mesh for time = " << runTimePtr_().timeName() << endl;

word regionName = fvMesh::defaultRegion;
meshPtr_.reset(
    new fvMesh(
        IOobject(
            regionName,
            runTimePtr_().timeName(),
            runTimePtr_(),
            IOobject::MUST_READ)));

daGlobalVarPtr_.reset(new DAGlobalVar(meshPtr_()));
