#!/bin/bash

if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit
fi

# pre-processing

# generate mesh
echo "Generating mesh.."

if [ -f "m6_surfaceMesh_fine.cgns.tar.gz" ]; then
  echo "Surface mesh m6_surfaceMesh_fine.cgns.tar.gz already exists."
else
  echo "Downloading surface mesh m6_surfaceMesh_fine.cgns.tar.gz"
  wget https://github.com/mdolab/dafoam/releases/download/v1.1.2/m6_surfaceMesh_fine.cgns.tar.gz --no-check-certificate
fi
tar -xvf m6_surfaceMesh_fine.cgns.tar.gz
# coarsen the surface mesh three times
cgns_utils coarsen m6_surfaceMesh_fine.cgns surfaceMesh.cgns
cgns_utils coarsen surfaceMesh.cgns
cgns_utils coarsen surfaceMesh.cgns
python genWingMesh.py > log.meshGeneration
plot3dToFoam -noBlank volumeMesh.xyz >> log.meshGeneration
autoPatch 60 -overwrite >> log.meshGeneration
createPatch -overwrite >> log.meshGeneration
renumberMesh -overwrite >> log.meshGeneration
echo "Generating mesh.. Done!"

# copy initial and boundary condition files
cp -r 0.orig 0
