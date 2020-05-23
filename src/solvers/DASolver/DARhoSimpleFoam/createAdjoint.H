// adjoint

daUtilPtr_.reset(new DAUtility);

daOptionPtr_.reset(new DAOption(mesh, pyOptions_));

daTurbulenceModelPtr_.reset(DATurbulenceModel::New(mesh));

daModelPtr_.reset(new DAModel(mesh));

daRegStatePtr_.reset(DARegState::New(mesh));

daIndexPtr_.reset(new DAIndex(mesh));
