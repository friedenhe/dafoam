// adjoint

DAUtility daUtil();

DAOption daOption(mesh, pyOptions_);

autoPtr<DARegState> daReg(
    DARegState::New(mesh));
