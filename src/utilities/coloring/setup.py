"""
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Cython setup file for wrapping OpenFOAM libraries and solvers.
        One needs to set include dirs/files and flags according to the
        information in Make/options and Make/files in OpenFOAM libraries
        and solvers. Then, follow the detailed instructions below. The
        python naming convention is to add "py" before the C++ class name
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import petsc4py
import numpy

os.environ["CC"] = "mpicc"
os.environ["CXX"] = "mpicxx"

solverName = "pyColoring"

if os.getenv("DAOF_AD_MODE") == "Passive":
    libSuffix = ""
    definedADMode = "-DDAOF_AD_MODE_Passive"
elif os.getenv("DAOF_AD_MODE") == "T1S":
    libSuffix = "ADF"
    definedADMode = "-DDAOF_AD_MODE_T1S"
elif os.getenv("DAOF_AD_MODE") == "A1S":
    libSuffix = "ADR"
    definedADMode = "-DDAOF_AD_MODE_A1S"
else:
    print("DAOF_AD_MODE not found!")
    exit(1)

# These setup should reproduce calling wmake to compile OpenFOAM libraries and solvers
ext = [
    Extension(
        solverName + libSuffix,
        # All source files, taken from Make/files
        sources=["pyColoring.pyx", "Coloring.C"],
        # All include dirs, refer to Make/options in OpenFOAM
        include_dirs=[
            # These are from Make/options:EXE_INC
            os.getenv("FOAM_SRC") + "/TurbulenceModels/turbulenceModels/lnInclude",
            os.getenv("FOAM_SRC") + "/TurbulenceModels/incompressible/lnInclude",
            os.getenv("FOAM_SRC") + "/TurbulenceModels/compressible/lnInclude",
            os.getenv("FOAM_SRC") + "/transportModels",
            os.getenv("FOAM_SRC") + "/transportModels/incompressible/singlePhaseTransportModel",
            os.getenv("FOAM_SRC") + "/transportModels/compressible/lnInclude",
            os.getenv("FOAM_SRC") + "/thermophysicalModels/basic/lnInclude",
            os.getenv("FOAM_SRC") + "/thermophysicalModels/radiation/lnInclude",
            os.getenv("FOAM_SRC") + "/finiteVolume/lnInclude",
            os.getenv("FOAM_SRC") + "/meshTools/lnInclude",
            os.getenv("FOAM_SRC") + "/sampling/lnInclude",
            os.getenv("FOAM_SRC") + "/fileFormats/lnInclude",
            os.getenv("FOAM_SRC") + "/surfMesh/lnInclude",
            os.getenv("FOAM_SRC") + "/dynamicFvMesh/lnInclude",
            os.getenv("FOAM_SRC") + "/dynamicMesh/lnInclude",
            # These are common for all OpenFOAM executives
            os.getenv("FOAM_SRC") + "/OpenFOAM/lnInclude",
            os.getenv("FOAM_SRC") + "/OSspecific/POSIX/lnInclude",
            os.getenv("FOAM_LIBBIN"),
            # CoDiPack and AdjointMPI
            os.getenv("FOAM_SRC") + "/../CoDiPack/include",
            os.getenv("FOAM_SRC") + "/Pstream/mpi/AdjointMPI/include",
            os.getenv("FOAM_SRC") + "/Pstream/mpi/AdjointMPI/src",
            # DAFoam include
            os.getenv("PETSC_DIR") + "/include",
            petsc4py.get_include(),
            numpy.get_include(),
            os.getenv("PETSC_DIR") + "/" + os.getenv("PETSC_ARCH") + "/include",
            "../../include",
            "../../adjoint/lnInclude",
            "./",
        ],
        # These are from Make/options:EXE_LIBS
        libraries=[
            "turbulenceModels",
            "incompressibleTurbulenceModels",
            "compressibleTurbulenceModels",
            "fluidThermophysicalModels",
            "specie",
            "incompressibleTransportModels",
            "compressibleTransportModels",
            "radiationModels",
            "finiteVolume",
            "sampling",
            "dynamicFvMesh",
            "dynamicMesh",
            "meshTools",
            "fvOptions",
            "DAOption",
            "DAUtility",
            "DACheckMesh",
            "DAStateInfo",
            "DAFvSource",
            "DAModel",
            "DAIndex",
            "DAFunction",
            "DAJacCon",
            "DAColoring",
            "DAResidual",
            "DAField",
            "DAPartDeriv",
            "DALinearEqn",
            "DARegression",
            "DAMisc",
            "DASolver",
            "petsc",
        ],
        # These are pathes of linked libraries
        library_dirs=[
            os.getenv("FOAM_LIBBIN"),
            os.getenv("PETSC_DIR") + "/lib",
            petsc4py.get_include(),
            os.getenv("PETSC_DIR") + "/" + os.getenv("PETSC_ARCH") + "/lib",
        ],
        # All other flags for OpenFOAM, users don't need to touch this
        extra_compile_args=[
            "-std=c++17",
            "-m64",
            "-pthread",
            "-DOPENFOAM=2112",
            "-DDAOF_AD_TOOL_DCO_FOAM",
            #"-Dlinux64",
            #"-DWM_ARCH_OPTION=64",
            "-DWM_DP",
            "-DWM_LABEL_SIZE=32",
            "-Wall",
            "-Wextra",
            "-Wno-deprecated-copy",
            "-Wnon-virtual-dtor",
            "-Wno-unused-parameter",
            "-Wno-invalid-offsetof",
            "-O3",
            "-DNoRepository",
            "-ftemplate-depth-100",
            "-fPIC",
            "-c",
            definedADMode,
        ],
        # Extra link flags for OpenFOAM, users don't need to touch this
        extra_link_args=["-shared", "-Xlinker", "--add-needed", "-Xlinker", "--no-as-needed"],
    )
]


setup(
    name=solverName + libSuffix,
    packages=[solverName + libSuffix],
    description="Cython wrapper for OpenFOAM",
    long_description="Cython wrapper for OpenFOAM",
    ext_modules=cythonize(ext, language_level=3),
)  # languate_level=3 means python3
