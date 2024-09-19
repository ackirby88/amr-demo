#!/bin/bash -e
# Folder structure:
# ================
# project_root
#   install            #install directory
#   builds             #build directory
#   src                #source codes
#     demo
#       amr            #location of source code
#       solvers
#       wake3d
#   scripts            #location of this script


# ======================= #
# project directory paths
# ======================= #
CURRENT_PATH="$(pwd)"
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
PROJECT_ROOT=${MY_PATH}/..

# ====================== #
# folder directory paths
# ====================== #
INSTALL_DIRECTORY=${PROJECT_ROOT}/install
INSTALL_DEMO_DIRECTORY=${INSTALL_DIRECTORY}/demo
BUILD_DIRECTORY=${PROJECT_ROOT}/builds
BUILD_DEMO_DIRECTORY=${BUILD_DIRECTORY}/demo

# ============== #
# demo sources
# ============== #
SOURCES_DIRECTORY=${PROJECT_ROOT}/src
DEMO_DIRECTORY=${SOURCES_DIRECTORY}
DEMO_SOURCES_DIRECTORY=${DEMO_DIRECTORY}

# =============== #
# demo defaults
# =============== #
NFIELDS2D_VAL=4
NFIELDS3D_VAL=5

# ========================= #
# third party library paths
# ========================= #
INSTALL_3RD_PARTY_PATH=${INSTALL_DIRECTORY}/3PL

# p4est path
P4EST_DIRECTORY=${INSTALL_3RD_PARTY_PATH}/p4est

# metis path
METIS_DIRECTORY=${INSTALL_3RD_PARTY_PATH}/metis

# MKL
MKL="OFF"

# ================== #
# compiling defaults
# ================== #
BUILD_DEMO=1
BUILD_TEST=0
BUILD_CLEAN=0
COMPILE_FAIL=0

# ================= #
# compiler defaults
# ================= #
FC=mpif90
CC=mpicc
CXX=mpicxx

C_FLAGS=
CXX_FLAGS=
Fortran_FLAGS=

# ======================== #
# compiler option defaults
# ======================== #
BUILD_SUFFIX="_release"
BUILD_TYPE="Release"

CONN_TYPE_VAL=2

# ======================== #
# make and install command
# ======================== #
MAKE_CMD="make -j4 install"

# ============== #
# print strings
# ============== #
opt_str="[OPTION] "

eC="\x1B[0m"
GC="\x1B[1;32m"
rC="\x1B[0;41m"
gC="\x1B[0;42m"
yC="\x1B[0;33m"
oC="\x1B[3;93m"
aC="\x1B[3;92m"
mC="\x1B[0;43m"

help() {
    echo "Usage: $0 [OPTION]...[COMPILER OPTIONS]...[DEMO OPTIONS]"
    echo " "
    echo "  [OPTION]:"
    echo "    --help     -h      displays this help message"
    echo "    --clean    -c      removes build directory: demo/builds/demo_{}"
    echo "    --release  -opt    compile the project in optimized mode"
    echo "    --debug    -deb    compile the project in debug mode"
    echo "    --testsON  -ton    turn on unit tests (google tests)"
    echo " "
    echo "  [COMPILER OPTIONS]:"
    echo "     CC=<arg>     cc=<arg>    sets the C compiler"
    echo "    CXX=<arg>    cxx=<arg>    sets the C++ compiler"
    echo "     FC=<arg>     fc=<arg>    sets the Fortran compiler"
    echo " "
    echo "      C_FLAGS=<arg>    c_flags=<arg>    sets the C compiler flags"
    echo "    CXX_FLAGS=<arg>  cxx_flags=<arg>    sets the C++ compiler flags"
    echo "     FC_FLAGS=<arg>   fc_flags=<arg>    sets the Fortran compiler flags"
    echo " "
    echo "  [Intel Flags]:"
    echo "     --avx     -avx         sets Intel AVX Instructions"
    echo "     --avx2    -avx2        sets Intel AVX-2 Instructions"
    echo "     --avx512  -avx512      sets Intel AVX-512 Instructions"
    echo "     --mkl     -mkl         sets Intel MKL Library"
    echo " "
    echo -e "  ${aC}Recommended Options:${eC}"
    echo "    Default (-go):"
    echo "      ./build_demo.sh CC=mpicc CXX=mpicxx FC=mpif90"
    echo " "
    echo -e "    Intel MPI (-impi): ${yC}./build_demo.sh -impi${eC}"
    echo "        CC=mpiicc CXX=mpiicpc FC=mpiifort"
}

# ----------------------------- #
# Start the compilation process #
# ----------------------------- #
cd $PROJECT_ROOT

# ============ #
# parse inputs
# ============ #
for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0

  elif [ "$var" == "--clean" -o "$var" == "-clean" -o "$var" == "-c" ]; then
    echo ${opt_str} "Clean and rebuild"
    BUILD_CLEAN=1
    BUILD_DEMO=0

  elif [ "$var" == "--release" -o "$var" == "-release" -o "$var" == "-opt" ]; then
    echo ${opt_str} "Compiling in optimized mode"
    BUILD_SUFFIX="_release"
    BUILD_TYPE="Release"

  elif [ "$var" == "--debug" -o "$var" == "-debug" -o "$var" == "-deb" ]; then
    echo ${opt_str} "Compiling in debug mode"
    BUILD_SUFFIX="_debug"
    BUILD_TYPE="Debug"

  elif [ "$var" == "--testsON" -o "$var" == "-testsON" -o "$var" == "-ton" ]; then
    echo ${opt_str} "Turning on unit testing (google tests)"
    BUILD_TEST=1

  elif [ "${var:0:3}" == "CC=" -o "${var:0:3}" == "cc=" ]; then
    CC=${var:3}
    echo -e "[OPTION]       C Compiler: $yC$CC$eC"

  elif [ "${var:0:4}" == "CXX=" -o "${var:0:4}" == "cxx=" ]; then
    CXX=${var:4}
    echo -e "[OPTION]     CXX Compiler: $yC$CXX$eC"

  elif [ "${var:0:3}" == "FC=" -o "${var:0:3}" == "fc=" ]; then
    FC=${var:3}
    F77=${FC}
    echo -e "[OPTION] Fortran Compiler: $yC$FC$eC"

  elif [ "${var:0:8}" == "C_FLAGS=" -o "${var:0:8}" == "c_flags=" ]; then
    C_FLAGS=${var:8}
    echo -e "[OPTION]       C Compiler Flags: $yC$C_FLAGS$eC"

  elif [ "${var:0:10}" == "CXX_FLAGS=" -o "${var:0:10}" == "cxx_flags=" ]; then
    CXX_FLAGS=${var:10}
    echo -e "[OPTION]     CXX Compiler Flags: $yC$CXX_FLAGS$eC"

  elif [ "${var:0:9}" == "FC_FLAGS=" -o "${var:0:9}" == "fc_flags=" ]; then
    FC=${var:9}
    echo -e "[OPTION] Fortran Compiler Flags: $yC$FC_FLAGS$eC"

  elif [ "$var" == "--AVX" -o "$var" == "-avx" ]; then
    C_FLAGS="${C_FLAGS} -xCore-AVX"
    CXX_FLAGS="${CXX_FLAGS} -xCore-AVX"
    FC_FLAGS="${FC_FLAGS} -xCore-AVX"

  elif [ "$var" == "--AVX2" -o "$var" == "-avx2" ]; then
    C_FLAGS="${C_FLAGS} -xCore-AVX2"
    CXX_FLAGS="${CXX_FLAGS} -xCore-AVX2"
    FC_FLAGS="${FC_FLAGS} -xCore-AVX2"

  elif [ "$var" == "--AVX512" -o "$var" == "-avx512" ]; then
    C_FLAGS="${C_FLAGS} -xCore-AVX512"
    CXX_FLAGS="${CXX_FLAGS} -xCore-AVX512"
    FC_FLAGS="${FC_FLAGS} -xCore-AVX512"

  elif [ "$var" == "-go" ]; then
    CC=mpicc
    CXX=mpicxx
    FC=mpif90

  elif [ "$var" == "-impi" ]; then
    CC=mpiicc
    CXX=mpiicpc
    FC=mpiifort

  fi
done

# ========================= #
# display command line args
# ========================= #
echo " "
echo "$0 $@"

# ----------------------------------------------------- #
# After reading in cmd arg options, set remaining paths #
# ----------------------------------------------------- #
# ====================================== #
# install/build location compiled source
# ====================================== #
COMPILE_INSTALL_DEMO_DIRECTORY="${INSTALL_DEMO_DIRECTORY}${BUILD_SUFFIX}"
COMPILE_BUILD_DEMO_DIRECTORY="${BUILD_DEMO_DIRECTORY}${BUILD_SUFFIX}"

# ============== #
# compiler paths
# ============== #
FC_PATH="`which $FC`"
CC_PATH="`which $CC`"
CXX_PATH="`which $CXX`"
LD_PATH="`which ld`"

# ====================== #
# check source directory
# ====================== #
if [ ! -d "${DEMO_SOURCES_DIRECTORY}" ]; then
  echo " "
  echo "Error:"
  echo "${DEMO_SOURCES_DIRECTORY} does not exist."
  exit 1
fi

# ======================= #
# check install directory
# ======================= #
if [ ! -d "${INSTALL_DIRECTORY}" ]; then
  echo  "${INSTALL_DIRECTORY} does not exist. Making it..."
  mkdir "${INSTALL_DIRECTORY}"
fi

# ====================== #
# check builds directory
# ====================== #
if [ ! -d "${BUILD_DIRECTORY}" ]; then
  echo  "${BUILD_DIRECTORY} does not exist. Making it..."
  mkdir "${BUILD_DIRECTORY}"
fi

UNIT_TEST="OFF"
if [ ${BUILD_TEST} == 1 ]; then
  UNIT_TEST="ON"
  if [ ! -d "${GTEST_DIRECTORY}" ]; then
    echo " "
    echo "==================================="
    echo -e "OPTION: Unit Testing ${rC}OFF${eC}"
    echo "WARNING: Google Test was NOT found."
    echo "         Turning off unit testing."
    echo "==================================="
    UNIT_TEST="OFF"
  else
    echo " "
    echo "==============================="
    echo -e "OPTION: Unit Testing ${gC}ON${eC}"
    echo "        Google Test was found: "
    echo "     ../install/3PL/googletest"
    echo "==============================="
  fi
fi


# ===================== #
# Check 3PL directories #
# ===================== #
if [ ! -d ${METIS_DIRECTORY} ]; then
  echo " "
  echo -e "${rC}ERROR: Required Metis not found!${eC}"
  echo "Metis Location: ${METIS_DIRECTORY}"
  BUILD_DEMO=0
  COMPILE_FAIL=1
fi
if [ ! -d ${P4EST_DIRECTORY} ]; then
  echo " "
  echo -e "${rC}ERROR: Required p4est not found!${eC}"
  echo "p4est Location: ${P4EST_DIRECTORY}"
  BUILD_DEMO=0
  COMPILE_FAIL=1
fi


# =================================================================== #
if [ $BUILD_CLEAN == 1 ]; then
  echo " "
  echo "Clean: removing ${COMPILE_BUILD_DEMO_DIRECTORY}..."
  echo "Clean: removing ${COMPILE_INSTALL_DEMO_DIRECTORY}..."
  echo " "
  rm -rf $COMPILE_BUILD_DEMO_DIRECTORY
  rm -rf $COMPILE_INSTALL_DEMO_DIRECTORY
fi

if [ $BUILD_DEMO == 1 ]; then
  echo " "
  echo -e "${mC} ===================== Building demo ==================== ${eC}"
  echo         "   Compiling Options:"
  echo         "          Build Type: ${BUILD_TYPE}"
  echo         "          Unit Tests: ${UNIT_TEST}"
  echo         " "
  echo         "                  CC: ${CC}"
  echo         "                 CXX: ${CXX}"
  echo         "                  FC: ${FC}"
  echo         " "
  echo         "            CC Flags: ${C_FLAGS}"
  echo         "           CXX Flags: ${CXX_FLAGS}"
  echo         "            FC Flags: ${FC_FLAGS}"
  echo         " "
  echo         "          Build Type: ${BUILD_TYPE}"
  echo         "      Build Location: ${COMPILE_BUILD_DEMO_DIRECTORY}"
  echo         "    Install Location: ${COMPILE_INSTALL_DEMO_DIRECTORY}"
  echo         " Executable Location: ${COMPILE_BUILD_DEMO_DIRECTORY}/bin"
  echo -e "${mC} ========================================================== ${eC}"
  echo " "

  # move to the build directory
  cd $BUILD_DIRECTORY

  if [ ! -d $COMPILE_BUILD_DEMO_DIRECTORY ]; then
    mkdir $COMPILE_BUILD_DEMO_DIRECTORY
  fi
  cd $COMPILE_BUILD_DEMO_DIRECTORY

  cmake -D CMAKE_C_COMPILER=${CC_PATH}                              \
        -D CMAKE_CXX_COMPILER=${CXX_PATH}                           \
        -D CMAKE_C_FLAGS=${C_FLAGS}                                 \
        -D CMAKE_CXX_FLAGS=${CXX_FLAGS}                             \
        -D CMAKE_Fortran_FLAGS=${Fortran_FLAGS}                     \
        -D CMAKE_LINKER=${LD_PATH}                                  \
        -D CMAKE_INSTALL_PREFIX=${COMPILE_INSTALL_DEMO_DIRECTORY}   \
        -D CMAKE_BUILD_TYPE=${BUILD_TYPE}                           \
        -D MKL_LIB=${MKL}                                           \
        -D metis_dir=${METIS_DIRECTORY}                             \
        -D p4est_dir=${P4EST_DIRECTORY}                             \
        -D CONN_TYPE=${CONN_TYPE_VAL}                               \
        -D NFIELDS2D=${NFIELDS2D_VAL}                               \
        -D NFIELDS3D=${NFIELDS3D_VAL}                               \
        -G "Unix Makefiles" ${DEMO_DIRECTORY} | tee cmake_config.out

  ${MAKE_CMD}
  cd ${CURRENT_PATH}

  if [ ! -d "${COMPILE_INSTALL_DEMO_DIRECTORY}" ]; then
    echo "ERROR:"
    echo "${COMPILE_INSTALL_DEMO_DIRECTORY} does not exist."
    COMPILE_FAIL=1
  fi
fi

if [ ${COMPILE_FAIL} == 0 ]; then
  echo " "
  echo -e " ========================================================== "
  echo -e " ${gC}demo build successful! ${eC}"
  echo    "   Compiling Options:"
  echo    "          Build Type: ${BUILD_TYPE}"
  echo    "          Unit Tests: ${UNIT_TEST}"
  echo    " "
  echo    "                  CC: ${CC}"
  echo    "                 CXX: ${CXX}"
  echo    "                  FC: ${FC}"
  echo    " "
  echo    "             C Flags: ${C_FLAGS}"
  echo    "           CXX Flags: ${CXX_FLAGS}"
  echo    "            FC Flags: ${FC_FLAGS}"
  echo    " "
  echo    "          Build Type: ${BUILD_TYPE}"
  echo    "      Build Location: ${COMPILE_BUILD_DEMO_DIRECTORY}"
  echo -e "                    : ${GC}make clean; make -j install${eC} in this directory"
  echo    "    Install Location: ${COMPILE_INSTALL_DEMO_DIRECTORY}"
  echo    " Executable Location: ${COMPILE_BUILD_DEMO_DIRECTORY}/bin"
  echo -e " ========================================================== "
  echo    " "
else
  echo " "
  echo         "======================"
  echo -e "${rC} demo build FAILED! ${eC}"
  echo         "======================"
  echo " "
  exit 1
fi
# =================================================================== #

# Create hyperlink to bin directory
ln -sf ${COMPILE_BUILD_DEMO_DIRECTORY}/bin ../bin

# Create hyperlink to metis lib in bin directory for MACOSX
if [[ "$OSTYPE" == "darwin"* ]]; then
  ln -sf ${INSTALL_DIRECTORY}/3PL/metis/lib/libmetis.dylib ${COMPILE_BUILD_DEMO_DIRECTORY}/bin/libmetis.dylib
fi

# Copy sample input files to the executable directory
cp ${PROJECT_ROOT}/example/input.* ${COMPILE_BUILD_DEMO_DIRECTORY}/bin/.
ln -s ${PROJECT_ROOT}/example/grids ${COMPILE_BUILD_DEMO_DIRECTORY}/bin/.

echo
echo "================================================"
echo -e "${gC} Finished Successfully...${eC}"
echo -e " Executable: ${GC}builds/demo_release/bin/demo?d.mpi${eC}"
echo "================================================"

echo -e "${gC}Demo Build Script Completed Successfully!${eC}"
