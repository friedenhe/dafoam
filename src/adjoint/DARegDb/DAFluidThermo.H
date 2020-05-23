/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Register fluidThermo to Foam::fvMesh

        Example:

        To register a fluidThermo object to fvMesh:

        autoPtr<fluidThermo> fluidThermo = fluidThermo::New(mesh);
        fluidThermo& thermo = fluidThermo();
        DAFluidThermo daThermo(mesh,thermo);

        Then, to get the registered fluidThermo object from fvMesh:

        DAFluidThermo& daThermo1 = const_cast<DAFluidThermo&>(
            mesh.thisDb().lookupObject<DAFluidThermo>("DAFluidThermo"));
        fluidThermo& thermo1 = const_cast<fluidThermo&>(daThermo1.getFuildThermo());

\*---------------------------------------------------------------------------*/

#ifndef DAFluidThermo_H
#define DAFluidThermo_H

#include "fvOptions.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class DAFluidThermo
    : public regIOobject
{

private:
    /// Disallow default bitwise copy construct
    DAFluidThermo(const DAFluidThermo&);

    /// Disallow default bitwise assignment
    void operator=(const DAFluidThermo&);

    /// Foam::fvMesh object
    const fvMesh& mesh_;

    /// fluid thermo
    fluidThermo& fluidThermo_;

public:
    /// Constructors
    DAFluidThermo(
        const fvMesh& mesh,
        fluidThermo& fluidThermo);

    /// Destructor
    virtual ~DAFluidThermo();

    /// return fluidThermo
    const fluidThermo& getFluidThermo()
    {
        return fluidThermo_;
    }

    /// this is a virtual function for regIOobject
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
