/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DASolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(DASolver, 0);
defineRunTimeSelectionTable(DASolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASolver::DASolver(
    char* argsAll,
    PyObject* pyOptions)
    : argsAll_(argsAll),
      pyOptions_(pyOptions)
{
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<DASolver> DASolver::New(
    char* argsAll,
    PyObject* pyOptions)
{
    // standard setup for runtime selectable classes

    // look up the solver name defined in pyOptions
    DAUtility daUtil;
    dictionary allOptions;
    daUtil.pyDict2OFDict(pyOptions, allOptions);
    word solverName;
    allOptions.readEntry<word>("solverName", solverName);

    Info << "Selecting " << solverName << " for DASolver" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverName);

    // if the solver name is not found in any child class, print an error
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "DASolver::New"
            "("
            "    char*,"
            "    PyObject*"
            ")")
            << "Unknown DASolver type "
            << solverName << nl << nl
            << "Valid DASolver types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // child class found
    return autoPtr<DASolver>(
        cstrIter()(argsAll, pyOptions));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label DASolver::loop(Time& runTime)
{
    const scalar& endTime = runTime.endTime().value();
    const scalar& deltaT = runTime.deltaT().value();
    const scalar t = runTime.timeOutputValue();
    if (t > endTime - 0.5 * deltaT)
    {
        return 0;
    }
    else
    {
        ++runTime;
        return 1;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
