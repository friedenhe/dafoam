#!/bin/bash

# Check if the OpenFOAM enviroments are loaded
if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit
fi

# pre-processing

# generate mesh
echo "Generating mesh.."

# Run the python script to read the airfoil coordinates from 
# the profiles folder and generate structured mesh using pyHyp
# and output the mesh as volumeMesh.xyz.  
python genAirFoilMesh.py &> logMeshGeneration.txt

# Use the plot3dToFoam utility to convert the plot3D mesh 
# volumeMesh.xyz to OpenFOAM format and store the mesh in 
# constant/polyMesh
plot3dToFoam -noBlank volumeMesh.xyz >> logMeshGeneration.txt

# The above generated mesh has one boundary patch, we 
# need to use the autoPatch utility to split boundaries
# here the argument 30 means split the mesh based on 30
# degree feature angle. Essentially, this call will modify
# the constant/polyMesh/boundary file
autoPatch 30 -overwrite >> logMeshGeneration.txt

# The above generated boundary file has boundary names such 
# as auto0, auto1, etc.Now rename the above generated patches 
# to meaningful names, e.g. wing, symmetry, etc, 
# see system/createPathDict
createPatch -overwrite >> logMeshGeneration.txt

# Renumber the mesh points to reduce the bandwith of Jacobians
# and reduce the memory usage
renumberMesh -overwrite >> logMeshGeneration.txt

echo "Generating mesh.. Done!"

# copy initial and boundary condition files
cp -r 0.orig 0
