#!/bin/bash

if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit
fi

# pre-processing

# generate mesh
echo "Generating mesh.."
python genAirFoilMesh.py > log.meshGeneration
plot3dToFoam -noBlank volumeMesh.xyz >> log.meshGeneration
autoPatch 30 -overwrite >> log.meshGeneration
createPatch -overwrite >> log.meshGeneration
renumberMesh -overwrite >> log.meshGeneration
echo "Generating mesh.. Done!"

# copy initial and boundary condition files
cp -r 0.orig 0
