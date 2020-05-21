// adjoint
daUtilPtr_.reset(new DAUtility);

daOptionPtr_.reset(new DAOption(mesh,pyOptions_));

autoPtr<DARegState> daReg(
    DARegState::New(mesh));
