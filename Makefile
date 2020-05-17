# make file

default:
	@make opt

opt:
	@echo "******************Compiling Opt Mode******************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allmake_Incompressible_Opt
# compile srs/solvers/incompressible
	cd src/solvers && ./Allmake_Opt
# compile src/utilities
	cd src/utilities && ./Allmake

debug:
	@echo "******************Compiling Debug Mode******************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allmake_Incompressible_Debug
# compile srs/solvers/incompressible
	cd src/solvers && ./Allmake_Debug

clean:
	@echo "********************Cleaning up********************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allclean
# compile srs/solvers/incompressible
	cd src/solvers && ./Allclean
# clean src/utilities
	cd src/utilities && ./Allclean

test:
	@echo "********************Running tests********************"
	cd tests && ./Allrun

