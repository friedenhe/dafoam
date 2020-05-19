// adjoint

DAUtility daUtil;

DAOption daOption(mesh, pyOptions_);

const dictionary& allOptions = daOption.getAllOptions();

autoPtr<DARegState> daReg(
    DARegState::New(mesh));
