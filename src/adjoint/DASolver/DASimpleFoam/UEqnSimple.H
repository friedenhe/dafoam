// Momentum predictor

tmp<fvVectorMatrix> tUEqn(
    fvm::div(phi, U)
    + turbulencePtr_->divDevReff(U));
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
