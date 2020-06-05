// adjoint

daUtilPtr_.reset(new DAUtility);

daOptionPtr_.reset(new DAOption(mesh, pyOptions_));

// need to register thermo and turbulence to mesh 
// before initializing DATurbulenceModel
DARegDbFluidThermo regDbThermo(mesh, thermo);
DARegDbTurbulenceModelCompressible regDbTurbulence(mesh, turbulencePtr_());

daTurbulenceModelPtr_.reset(DATurbulenceModel::New(mesh));

daModelPtr_.reset(new DAModel(mesh));

daRegStatePtr_.reset(DARegState::New(mesh));

daIndexPtr_.reset(new DAIndex(mesh));

daJacConPtr_.reset(DAJacCon::New(mesh));
