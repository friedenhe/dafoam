/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        DAModel is a general template for handling physical models in DAFoam. 
        The constructor for DAModel does not explicitly contain objects of 
        physical models. Instead, it will look up any model class, e.g., 
        DATurbulenceModel, DARaditionModel, that has been registred in fvMesh. 
        This enables the flexibility such that we  can simultaneously have 
        multiple models for adjoint, e.g., turbulence model + radiation model.
        
        NOTE 1: in adjoint functions, we should call functions in DAModel to handle
        general model operation, e.g., call daModel.correctBoundaryCondition(), 
        instead of calling daTurbulenceModel.correctBoundaryCondition(). Because
        daModel.correctBoundaryCondition() will call the correctBoundaryCondition()
        in daTurbulenceModel and daRadiationModel, if these two models are used.

        NOTE 2: Because of the above feature, we need to initialize any physical
        model BEFORE initializing the DAModel, e.g.,

        // initialize DATurbulence and automatically register it to fvMesh
        autoPtr<DATurbulenceModel> daTurb(DATurbulenceModel::New(mesh));
        
        ............ // codes to initialize any other physical models

        // now we can initialize DAModel and lookup DATurbulenceModel in mesh
        DAModel daModel(mesh);

        Otherwise, the physical models will NOT be picked up by the DAModel

\*---------------------------------------------------------------------------*/

#ifndef DAModel_H
#define DAModel_H

#include "fvOptions.H"
#include "DATurbulenceModel.H"
#include "DARadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAModel Declaration
\*---------------------------------------------------------------------------*/

class DAModel
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DAModel(const DAModel&);

    /// Disallow default bitwise assignment
    void operator=(const DAModel&);

    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// have turbulence models for the flow and adjoint?
    label hasTurbulenceModel_ = 0;

    /// have radiation models for the flow and adjoint?
    label hasRadiationModel_ = 0;

public:
    /// Constructors
    DAModel(const fvMesh& mesh);

    /// Destructor
    virtual ~DAModel();

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;

    /// update the name in modelStates based on the selected physical model at runtime
    void correctModelStates(wordList& modelStates) const;

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
