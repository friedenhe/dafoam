// adjoint
daUtilPtr_.reset(new DAUtility);

daOptionPtr_.reset(new DAOption(mesh,pyOptions_));

autoPtr<DATurbulenceModel> daTurbmodel(
    DATurbulenceModel::New(mesh));

DAModel daModel(mesh);

autoPtr<DARegState> daReg(
    DARegState::New(mesh));

DAIndex daIndex(mesh);
