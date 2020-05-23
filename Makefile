# make file

default:
	@make opt

opt:
	@echo "******************Compiling Opt Mode******************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allmake_Incompressible_Opt
# compile src/adjoint with compressible
	cd src/adjoint && ./Allmake_Compressible_Opt
# compile srs/solvers/incompressible
	cd src/pyDASolvers && ./Allmake_Opt
# compile src/utilities
	cd src/utilities && ./Allmake

debug:
	@echo "******************Compiling Debug Mode******************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allmake_Incompressible_Debug
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allmake_Compressible_Debug
# compile srs/solvers/incompressible
	cd src/pyDASolvers && ./Allmake_Debug

clean:
	@echo "********************Cleaning up********************"
# compile src/adjoint with incompressible
	cd src/adjoint && ./Allclean
# compile srs/solvers/incompressible
	cd src/pyDASolvers && ./Allclean
# clean src/utilities
	cd src/utilities && ./Allclean
# clean src/utilities
	cd tests && ./Allclean
# clean python
	rm -rf src/__pycache__/

test:
	@echo "********************Running tests********************"
	cd tests && ./Allclean && ./Allrun

