/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register wall distance to Foam::fvMesh

        Example:

        To register a nearWallDist object to fvMesh:

        nearWallDist& d = const_cast<nearWallDist&>(turbulence->y());
        DANearWallDist daD(mesh,d);

        Then, to get the registered nearWallDist object from fvMesh:

        DANearWallDist& daD1= const_cast<DANearWallDist&>(
            mesh.thisDb().lookupObject<DANearWallDist>("DANearWallDist"));
        nearWallDist& d1 = const_cast<nearWallDist&>(daD1.getNearWallDist());

\*---------------------------------------------------------------------------*/

#ifndef DANearWallDist_H
#define DANearWallDist_H

#include "fvOptions.H"
#include "nearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DANearWallDist Declaration
\*---------------------------------------------------------------------------*/

class DANearWallDist
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DANearWallDist(const DANearWallDist&);

    /// Disallow default bitwise assignment
    void operator=(const DANearWallDist&);

    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// nearWallDist
    nearWallDist& nearWallDist_;

public:
    /// Constructors
    DANearWallDist(
        const fvMesh& mesh,
        nearWallDist& nearWallDist);

    /// Destructor
    virtual ~DANearWallDist();

    /// return nearWallDist_
    const nearWallDist& getNearWallDist()
    {
        return nearWallDist_;
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
