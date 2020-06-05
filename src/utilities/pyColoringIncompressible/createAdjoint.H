DAUtility daUtil;

DAOption daOption(mesh, pyOptions_);

// need to register turbulence to mesh
DARegDbSinglePhaseTransportModel regDbTransport(mesh, laminarTransport);
DARegDbTurbulenceModelIncompressible regDbTurbulence(mesh, turbulence());

autoPtr<DATurbulenceModel> daTurbulenceModel(DATurbulenceModel::New(mesh));

DAModel daModel(mesh);

autoPtr<DARegState> daRegState(DARegState::New(mesh));

DAIndex daIndex(mesh);

autoPtr<DAJacCon> daJacCon(DAJacCon::New(mesh));
