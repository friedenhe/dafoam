// adjoint

DAUtility daUtil();

DAOption daOption(mesh, pyOptions_);

autoPtr<DATurbulenceModel> daTurbmodel(
    DATurbulenceModel::New(mesh));

DAModel daModel(mesh);

autoPtr<DARegState> daReg(
    DARegState::New(mesh));

DAIndex daIndex(mesh);
