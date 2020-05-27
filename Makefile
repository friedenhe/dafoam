# make file

default:
	@make opt

opt:
	@echo "******************Compiling Opt Mode******************"
# compile src/adjoint 
	cd src/adjoint && ./Allmake_Opt
# compile src/pySolvers
	cd src/pyDASolvers && ./Allmake_Opt
# compile src/utilities
	cd src/utilities && ./Allmake

debug:
	@echo "******************Compiling Debug Mode******************"
# compile src/adjoint
	cd src/adjoint && ./Allmake_Debug
# compile src/pySolvers
	cd src/pyDASolvers && ./Allmake_Debug

clean:
	@echo "********************Cleaning up********************"
# clean src/adjoint with incompressible
	cd src/adjoint && ./Allclean
# clean src/pySolvers
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

