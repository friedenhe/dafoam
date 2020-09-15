/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAFvSourceActuatorPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAFvSourceActuatorPoint, 0);
addToRunTimeSelectionTable(DAFvSource, DAFvSourceActuatorPoint, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAFvSourceActuatorPoint::DAFvSourceActuatorPoint(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAFvSource(modelType, mesh, daOption, daModel, daIndex)
{
    const dictionary& allOptions = daOption_.getAllOptions();
    fvSourceSubDict_ = allOptions.subDict("fvSource");

    printIntervalUnsteady_ = daOption.getOption<label>("printIntervalUnsteady");
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAFvSourceActuatorPoint::calcFvSource(volVectorField& fvSource)
{
    /*
    Description:
        Compute the actuator point source term. We support two types of smooth functions:

        1. Hyperbolic:
            thrust = xT*yT*zT where xT is
            xT = tanh( eps*(meshCX-actCX+0.5*actSX) ) -  tanh( eps*(meshCX-actCX-0.5*actSX) )
            if the mesh cell center is within the actuator box, xT~2 other wise xT~0
            eps is a parameter to smooth the distribution

        2. Gaussian: 
            We use a 2D Gaussian function for thrust source distribution:
            thrust = 1/(2*pi*eps^2)*e^{-d^2/(2*eps^2)} where d is the
            distance between the cell center and actuator source center

            NOTE: this option has some issues for U field... it has many spikes
    
    Example:
        An example of the fvSource in pyOptions in pyDAFoam can be
        defOptions = 
        {
            "fvSource"
            {
                "line1"
                {
                    "type": "actuatorPoint",
                    "smoothFunction": "hyperbolic",
                    "center": [0.0, 0.0, 0.0], # center and size define a rectangular
                    "size": [0.5, 0.3, 0.5],
                    "amplitude": [0.0, 1.0, 0.0],
                    "thrustDirIdx": 0,
                    "periodicity": 1.0,
                    "eps": 1.0,
                    "scale": 1.0  # scale the source such the integral equals desired thrust
                },
                "line2"
                {
                    "type": "actuatorDisk",
                    "smoothFunction": "gaussian",
                    "center": [0.0, 0.0, 0.0],
                    "amplitude": [0.0, 1.0, 0.0],
                    "thrustDirIdx": 0,
                    "periodicity": 1.0,
                    "eps": 1.0,
                    "scale": 1.0  # scale the source such the integral equals desired thrust
                }
            }
        }
    */

    forAll(fvSource, idxI)
    {
        fvSource[idxI] = vector::zero;
    }

    // loop over all the cell indices for all actuator disks
    forAll(fvSourceSubDict_.toc(), idxI)
    {

        // name of this point model
        word pointName = fvSourceSubDict_.toc()[idxI];

        // sub dictionary with all parameters for this disk
        dictionary pointSubDict = fvSourceSubDict_.subDict(pointName);

        word smoothFunction = pointSubDict.get<word>("smoothFunction");

        if (smoothFunction == "hyperbolic")
        {
            // now read in all parameters for this actuator point
            scalarList center1;
            scalarList size1;
            scalarList amp1;
            pointSubDict.readEntry<scalarList>("center", center1);
            pointSubDict.readEntry<scalarList>("size", size1);
            pointSubDict.readEntry<scalarList>("amplitude", amp1);
            vector center = {center1[0], center1[1], center1[2]};
            vector size = {size1[0], size1[1], size1[2]};
            vector amp = {amp1[0], amp1[1], amp1[2]};
            scalar period = pointSubDict.get<scalar>("periodicity");
            scalar eps = pointSubDict.get<scalar>("eps");
            scalar scale = pointSubDict.get<scalar>("scale");
            scalar thrustDirIdx = pointSubDict.get<label>("thrustDirIdx");

            scalar t = mesh_.time().timeOutputValue();
            center += amp * sin(constant::mathematical::twoPi * t / period);

            scalar xTerm, yTerm, zTerm, s;
            scalar thrustTotal = 0.0;
            forAll(mesh_.C(), cellI)
            {
                const vector& meshC = mesh_.C()[cellI];
                xTerm = (tanh(eps * (meshC[0] + 0.5 * size[0] - center[0])) - tanh(eps * (meshC[0] - 0.5 * size[0] - center[0])));
                yTerm = (tanh(eps * (meshC[1] + 0.5 * size[1] - center[1])) - tanh(eps * (meshC[1] - 0.5 * size[1] - center[1])));
                zTerm = (tanh(eps * (meshC[2] + 0.5 * size[2] - center[2])) - tanh(eps * (meshC[2] - 0.5 * size[2] - center[2])));

                s = xTerm * yTerm * zTerm;
                fvSource[cellI][thrustDirIdx] = s * scale;

                thrustTotal += s * scale * mesh_.V()[cellI];
            }
            reduce(thrustTotal, sumOp<scalar>());

            if (daOption_.getOption<word>("runStatus") == "solvePrimal")
            {
                if (mesh_.time().timeIndex() % printIntervalUnsteady_ == 0 || mesh_.time().timeIndex() == 0)
                {
                    Info << "Actuator point source: " << pointName << endl;
                    Info << "Total thrust source: " << thrustTotal << endl;
                }
            }
        }
        else if (smoothFunction == "gaussian")
        {
            // now read in all parameters for this actuator point
            scalarList center1;
            scalarList amp1;
            pointSubDict.readEntry<scalarList>("center", center1);
            pointSubDict.readEntry<scalarList>("amplitude", amp1);
            vector center = {center1[0], center1[1], center1[2]};
            vector amp = {amp1[0], amp1[1], amp1[2]};
            scalar period = pointSubDict.get<scalar>("periodicity");
            scalar eps = pointSubDict.get<scalar>("eps");
            scalar scale = pointSubDict.get<scalar>("scale");
            scalar thrustDirIdx = pointSubDict.get<label>("thrustDirIdx");
            scalar t = mesh_.time().timeOutputValue();
            center += amp * sin(constant::mathematical::twoPi * t / period);

            scalar thrustTotal = 0.0;
            scalar coeff = 1.0 / constant::mathematical::twoPi / eps / eps;
            forAll(mesh_.C(), cellI)
            {
                const vector& meshC = mesh_.C()[cellI];
                scalar d = mag(meshC - center);
                scalar s = coeff * Foam::exp(-d * d / 2.0 / eps / eps);
                fvSource[cellI][thrustDirIdx] = s * scale;
                thrustTotal += s * scale * mesh_.V()[cellI];
            }
            reduce(thrustTotal, sumOp<scalar>());

            if (daOption_.getOption<word>("runStatus") == "solvePrimal")
            {
                if (mesh_.time().timeIndex() % printIntervalUnsteady_ == 0 || mesh_.time().timeIndex() == 0)
                {
                    Info << "Actuator point source: " << pointName << endl;
                    Info << "Total thrust source: " << thrustTotal << endl;
                }
            }
        }
        else
        {
            FatalErrorIn("") << "smoothFunction should be either hyperbolic or gaussian" << abort(FatalError);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //