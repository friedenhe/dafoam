/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register singlePhaseTransportModel to Foam::fvMesh

        Example:

        To register a singlePhaseTransportModel object to fvMesh:

        singlePhaseTransportModel laminarTransport(U, phi);
        DASinglePhaseTransportModel sg(mesh,laminarTransport);

        Then, to get the registered nearWallDist object from fvMesh:

        DASinglePhaseTransportModel& sg1= const_cast<DASinglePhaseTransportModel&>(
            mesh.thisDb().lookupObject<DASinglePhaseTransportModel>(
            "DASinglePhaseTransportModel"));
        singlePhaseTransportModel& lt1 = const_cast<singlePhaseTransportModel&>(
            sg1.getSinglePhaseTransportModel());

\*---------------------------------------------------------------------------*/

#ifndef DASinglePhaseTransportModel_H
#define DASinglePhaseTransportModel_H

#include "fvOptions.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DASinglePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class DASinglePhaseTransportModel
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DASinglePhaseTransportModel(const DASinglePhaseTransportModel&);

    /// Disallow default bitwise assignment
    void operator=(const DASinglePhaseTransportModel&);

    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// fluid thermo
    singlePhaseTransportModel& singlePhaseTransportModel_;

public:
    /// Constructors
    DASinglePhaseTransportModel(
        const fvMesh& mesh,
        singlePhaseTransportModel& singlePhaseTransportModel);

    /// Destructor
    virtual ~DASinglePhaseTransportModel();

    /// return fluidThermo
    const singlePhaseTransportModel& getSinglePhaseTransportModel()
    {
        return singlePhaseTransportModel_;
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
