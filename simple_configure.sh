#!/bin/bash

# Configure script command line options
function usage() {
    echo "$0 ...options..." 1>&2
    echo "" 1>&2
    echo "All options take the form OPT=<value> unless otherwise specified" 1>&2
    echo "" 1>&2
    echo "FC                   Fortran compiler" 1>&2
    echo "CC                   C compiler" 1>&2
    echo "LAPACK_LIB_DIR       Location of LAPACK libraries" 1>&2
    echo "BIN_DIR              Path for installing executables" 1>&2
    echo "INCLUDE_DIR          Path for installing object files, headers, and modules" 1>&2
    echo "-i|--interactive     Interactive prompts for inputs" 1>&2
    echo "-d|--default         Matt's default values " 1>&2
    echo "                         FC=gfortran" 1>&2
    echo "                         CC=gcc" 1>&2
    echo "                         LAPACK_LIB_DIR=/sw/lib/lapack" 1>&2
    echo "                         BIN_DIR=./bin" 1>&2
    echo "                         INCLUDE_DIR=./include" 1>&2
    echo ""  1>&2
    echo "Note: If trying to install programs in a root directory, may have to use \"sudo $0\""  1>&2
    echo "But remember: with great power must also come great responsibility"  1>&2
    exit 1
}
if [ $# -eq 0 ]
then
    usage
fi

# Parse command line
FC=""
CC=""
LAPACK_LIB_DIR=""
BIN_DIR=""
INCLUDE_DIR=""
INTERACTIVE="N"
while [ "$1" != "" ]
do
    case $1 in
        FC=*)             FC=`echo $1 | sed -e "s/.*=//"`;;
        CC=*)             CC=`echo $1 | sed -e "s/.*=//"`;;
        LAPACK_LIB_DIR=*) LAPACK_LIB_DIR=`echo $1 | sed -e "s/.*=//"`;;
        BIN_DIR=*)        BIN_DIR=`echo $1 | sed -e "s/.*=//"`;;
        INCLUDE_DIR=*)    INCLUDE_DIR=`echo $1 | sed -e "s/.*=//"`;;
        -i|--interactive) INTERACTIVE="Y";;
        -d|--default)     FC="gfortran"
                          CC="gcc"
                          LAPACK_LIB_DIR="/sw/lib/lapack"
                          BIN_DIR="./bin"
                          INCLUDE_DIR="./include";;
        *) echo "!! Error: No option $1" 1>&2; usage;;
    esac
    shift
done

echo "Using variables:"
echo "    FC=$FC"
echo "    CC=$CC"
echo "    LAPACK_LIB_DIR=$LAPACK_LIB_DIR"
echo "    BIN_DIR=$BIN_DIR"
echo "    INCLUDE_DIR=$INCLUDE_DIR"

#####
#	Set up and test Fortran compiler
#####
echo "####################################################################################################"
echo "##########                         WORKING ON FORTRAN STUFF                               ##########"
echo "####################################################################################################"
# Get name of Fortran compiler
if [ -z "$FC" -o $INTERACTIVE == "Y" ]
then
    echo "Enter the name of your Fortran compiler:"
    read FC
fi

echo "Testing $FC installation"

# Look for an executable named $FC
FC_TEST=`which $FC`
if [ -z $FC_TEST ]
then
    echo "!! Error: no executable found named $FC" 1>&2
    echo "!! Check to see if $FC is in your PATH or full path is written correctly" 1>&2
    exit 1
fi

# Compile and run a couple simple Fortran programs
FC_TEST="config_fortran_compiler_test"
if [ -f $FC_TEST ]
then
    rm $FC_TEST
fi

# Fortran 77 test
cat > $FC_TEST.f << EOF
      PROGRAM main
      IMPLICIT none
      write(*,'(A)') '1234567890'
      END
EOF
$FC $FC_TEST.f -o $FC_TEST
chmod +x $FC_TEST
FC_OUTPUT=`./$FC_TEST`
if [ "$FC_OUTPUT" != "1234567890" ]
then
    echo "!! Error: Fortran compiler does not compile Fortran 77 test example!" 1>&2
    exit 1
fi
# Clean up
rm $FC_TEST.f $FC_TEST

# Modern Fortran test
cat > $FC_TEST.f90 << EOF
program main
implicit none
write(*,'(A)') '1234567890'
end
EOF
$FC $FC_TEST.f90 -o $FC_TEST
chmod +x $FC_TEST
FC_OUTPUT=`./$FC_TEST`
if [ "$FC_OUTPUT" != "1234567890" ]
then
    echo "!! Error: Fortran compiler does not compile Fortran 90 test example!" 1>&2
    exit 1
else
    echo "Test complete! Using Fortran compiler: $FC"
fi
# Clean up
rm $FC_TEST.f90 $FC_TEST

echo


#####
#	Set up and test C compiler
#####
echo "####################################################################################################"
echo "##########                            WORKING ON C STUFF                                  ##########"
echo "####################################################################################################"
# Get name of C compiler
if [ -z "$CC" -o $INTERACTIVE == "Y" ]
then
    echo "Enter the name of your C compiler:"
    read CC
fi

echo "Testing $CC installation"

# Look for an executable named $CC
CC_TEST=`which $CC`
if [ -z $CC_TEST ]
then
    echo "!! Error: no executable found named $CC" 1>&2
    echo "!! Check to see if $CC is in your PATH or full path is written correctly" 1>&2
    exit 1
fi

# Compile and run a couple simple C programs
CC_TEST="config_c_compiler_test"
if [ -f $CC_TEST ]
then
    rm $CC_TEST
fi

# C test
cat > $CC_TEST.c << EOF
#include <stdio.h>
#include <stdlib.h>
int main(void)
{
  puts("1234567890");
  return EXIT_SUCCESS;
}

EOF
$CC $CC_TEST.c -o $CC_TEST
chmod +x $CC_TEST
CC_OUTPUT=`./$CC_TEST`
if [ "$CC_OUTPUT" != "1234567890" ]
then
    echo "!! Error: C compiler does not compile C test example!" 1>&2
    exit 1
else
    echo "Test complete! Using C compiler: $CC"
fi
# Clean up
rm $CC_TEST.c $CC_TEST

echo


#####
#	Set up LAPACK path and test LAPACK library
#####
echo "####################################################################################################"
echo "##########                          WORKING ON LAPACK STUFF                               ##########"
echo "####################################################################################################"
if [ $INTERACTIVE == "Y" ]
then
    echo "If you want to install the software that depends on the LAPACK libraries, enter the directory"
    echo "that contains the LAPACK libraries: (Otherwise, leave blank and just hit \"return\")"
    read LAPACK_LIB_DIR
fi

if [ -z "$LAPACK_LIB_DIR" ]
then
    echo "Installing software without LAPACK dependencies"
    CPP_LAPACK=""
else
    # Make sure the directory exists and has the correct libraries in it
    LAPACK_LIB_DIR_CONTENTS=`ls $LAPACK_LIB_DIR 2>&1`
    EXIST=`echo $LAPACK_LIB_DIR_CONTENTS |\
           awk '{if(/No such file or directory/){print "N";exit}else{print "Y";exit}}'`

    if [ "$EXIST" == "N" ]
    then
        echo "!! Error: Directory $LAPACK_LIB_DIR is non-existent" 1>&2
        exit 1
    fi

    if [ -z "$LAPACK_LIB_DIR_CONTENTS" ]
    then
        echo "!! Error: Directory $LAPACK_LIB_DIR is empty" 1>&2
        exit 1
    else
        LAPACK_LIB_LIST=""
        for LIB in lapack refblas tmglib
        do
            EXIST=`ls $LAPACK_LIB_DIR/lib${LIB}.a 2>&1 |\
                   awk '{if(/No such file or directory/){print "N";exit}else{print "Y";exit}}'`
            LIB2=$LIB

            if [ "$EXIST" == "N" -a "$LIB" == "lapack" ]
            then
                echo "LAPACK library \"lapack\" does not exist; looking for library \"reflapack\""
                EXIST=`ls $LAPACK_LIB_DIR/libreflapack.a 2>&1 |\
                       awk '{if(/No such file or directory/){print "N";exit}else{print "Y";exit}}'`
                LIB2=reflapack
            fi

            if [ "$EXIST" == "N" -a "$LIB" == "refblas" ]
            then
                echo "LAPACK library \"refblas\" does not exist; looking for library \"blas\""
                EXIST=`ls $LAPACK_LIB_DIR/libblas.a 2>&1 |\
                       awk '{if(/No such file or directory/){print "N";exit}else{print "Y";exit}}'`
                LIB2=blas
            fi

            if [ "$EXIST" == "N" ]
            then
                echo "!! Error: LAPACK library \"$LIB2\" is not in $LAPACK_LIB_DIR" 1>&2
                exit 1
            fi

            echo "Found LAPACK library: $LIB2!"
            LAPACK_LIB_LIST="$LAPACK_LIB_LIST -l$LIB2"
        done
    fi
    echo "Test complete! Using LAPACK libraries in directory: $LAPACK_LIB_DIR"
    CPP_LAPACK="-DUSE_LAPACK"
fi
echo


#####
#	Directory to install software
#####
echo "####################################################################################################"
echo "##########                    WORKING ON INSTALL DIRECTORY STUFF                          ##########"
echo "####################################################################################################"
if [ -z "$BIN_DIR" -o $INTERACTIVE == "Y" ]
then
    echo "Enter the path to place the compiled executables (Default: ./bin):"
    read ANSWER
    if [ -z "$ANSWER" ]
    then
        ANSWER="./bin"
    fi
    BIN_DIR=$ANSWER
fi
echo "Putting executables in $BIN_DIR"
if [ ! -d $BIN_DIR ]
then
    echo "Directory $BIN_DIR does not exist....creating it"
    mkdir $BIN_DIR
else
    echo "Directory $BIN_DIR already exists"
fi
echo


#####
#	Directory to install object files and Fortran modules (include) files
#####
echo "####################################################################################################"
echo "##########                    WORKING ON INCLUDE DIRECTORY STUFF                          ##########"
echo "####################################################################################################"
if [ -z "$INCLUDE_DIR" -o $INTERACTIVE == "Y" ]
then
    echo "Enter the path to place the compiled object files (Default: ./include):"
    read ANSWER
    if [ -z "$ANSWER" ]
    then
        ANSWER="./include"
    fi
    INCLUDE_DIR=$ANSWER
fi
echo "Putting object files and modules in $INCLUDE_DIR"
if [ ! -d $INCLUDE_DIR ]
then
    echo "Directory $INCLUDE_DIR does not exist....creating it"
    mkdir $INCLUDE_DIR
else
    echo "Directory $INCLUDE_DIR already exists"
fi
echo


echo "####################################################################################################"
echo "##########                           INSTALL THE CODES!                                   ##########"
echo "####################################################################################################"
echo "You used the options:"
echo "    FC=$FC"
echo "    CC=$CC"
echo "    LAPACK_LIB_DIR=$LAPACK_LIB_DIR"
echo "    BIN_DIR=$BIN_DIR"
echo "    INCLUDE_DIR=$INCLUDE_DIR"
echo "Type \"make\" to install codes"


# Set libraries to be read by the compiler
if [ "$LAPACK_LIB_DIR" != "" ]
then
    LAPACK_LIB_DIR="-L${LAPACK_LIB_DIR}"
fi


#####
#	Generating makefile
#####
cat > Makefile << EOF
##### Compiler variables #####
CC = $CC
FC = $FC
FWARN = -Wall -Wextra -Wunused -fbounds-check -fbacktrace
FOPT  = -O1
FFLAG = \$(FWARN) \$(FOPT)


##### Include directory (.o and .mod files) #####
INCLUDE_DIR = $INCLUDE_DIR


##### Executable directory #####
BIN   = $BIN_DIR


##### External libraries #####
# LAPACK (linear algebra)
LAPACK_LIB_DIR = $LAPACK_LIB_DIR
LAPACK_LIB     = $LAPACK_LIB_LIST
LAPACK         = \$(LAPACK_LIB_DIR) \$(LAPACK_LIB)
CPP_LAPACK     = $CPP_LAPACK

# # Super LU (sparse solvers)
# SUPERLU_LIB_DIR = -Lext/SuperLU_5.2.1/lib
# SUPERLU_LIB     = -lsuperlu_5.1
# SUPERLU         = \$(SUPERLU_LIB_DIR) \$(SUPERLU_LIB)
CPP_SUPERLU       = -DUSE_SUPERLU
CPP_SUPERLU       =

# Pre-processing directives
CPP = -cpp \$(CPP_LAPACK) \$(CPP_SUPERLU)


##### ALL COMPILER FLAGS #####
ALL_FLAGS = \$(FFLAGS) -I\$(INCLUDE_DIR) \$(EXTERNAL_LIBS) \$(CPP)


##### GROUPS OF PROGRAMS #####
all: defm geom misc fits seis scripts other

defm: \\
      \$(BIN)/o92util \\
      \$(BIN)/vec2los \\
      \$(BIN)/wraplos \\
      \$(BIN)/triutil

geom: \\
      \$(BIN)/lola2distaz \\
      \$(BIN)/distaz2lola \\
      \$(BIN)/sphfinrot \\
      \$(BIN)/platemotion \\
      \$(BIN)/utm2geo

misc: \\
      \$(BIN)/colortool \\
      \$(BIN)/dateutil \\
      \$(BIN)/eventfrequency \\
      \$(BIN)/ff2gmt \\
      \$(BIN)/grid \\
      \$(BIN)/perturb \\
      \$(BIN)/readkik \\
      \$(BIN)/eqempirical \\
      \$(BIN)/stereo_project

fits: \\
      \$(BIN)/polyfit \\
      \$(BIN)/polyfit_special \\
      \$(BIN)/multifit \\
      \$(BIN)/fltinv

seis: \\
      \$(BIN)/mtutil \\
      \$(BIN)/readGCMT

scripts: \\
      \$(BIN)/coul_hor.sh \\
      \$(BIN)/coul_dip.sh \\
      \$(BIN)/coul_xsec.sh \\
      \$(BIN)/surf_disp.sh \\
      \$(BIN)/simplify_ffm.sh \\
      \$(BIN)/ternary.sh \\
      \$(BIN)/trg_schem.sh \\
      \$(BIN)/gmtcpt.sh

other: \\
      \$(BIN)/numint

####################################
##### MODULE AND OBJECT RULES ######
####################################
INCLUDE_FILES = \$(INCLUDE_DIR)/annealing.o \\
                \$(INCLUDE_DIR)/earth.o \\
                \$(INCLUDE_DIR)/elast.o \\
                \$(INCLUDE_DIR)/error_exit.o \\
                \$(INCLUDE_DIR)/geom.o \\
                \$(INCLUDE_DIR)/io.o \\
                \$(INCLUDE_DIR)/okada92.o \\
                \$(INCLUDE_DIR)/random.o \\
                \$(INCLUDE_DIR)/test.o \\
                \$(INCLUDE_DIR)/tri_disloc.o \\
                \$(INCLUDE_DIR)/trig.o

include: \$(INCLUDE_FILES)

\$(INCLUDE_DIR)/annealing.o: src/annealing_module.f90 \$(INCLUDE_DIR)/random.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/earth.o: src/earth_module.f90 \$(INCLUDE_DIR)/trig.o \$(INCLUDE_DIR)/io.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/elast.o: src/elast_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/error_exit.o: src/error_exit.c
	\$(CC) \$<  -c -o \$@

\$(INCLUDE_DIR)/geom.o: src/geom_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/io.o: src/io_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/okada92.o: src/okada92_module.f90 \$(INCLUDE_DIR)/test.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/random.o: src/random_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/test.o: src/test_module.f90 \$(INCLUDE_DIR)/error_exit.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/trig.o: src/trig_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/tri_disloc.o: src/tri_disloc_module.f90 \$(INCLUDE_DIR)/trig.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE_DIR) \$< -c -o \$@

\$(INCLUDE_DIR)/%.o: \$(INCLUDE_DIR)/%.mod

##################################
##### COMPILED PROGRAM RULES #####
##################################
\$(BIN)/colortool: src/colortool.f90
	\$(FC) \$(FFLAG) -o  \$@ \$^
	rm colormodule.mod

\$(BIN)/dateutil: src/dateutil.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/distaz2lola: src/distaz2lola.f src/geomsubs.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/eqempirical: src/eqempirical.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/eventfrequency: src/eventfrequency.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/ff2gmt: src/ff2gmt.f
	\$(FC) \$(FFLAG) -o \$@ \$^

FLTINV_MODULES = \\
                 src/fltinv_variable_module.f90 \\
                 src/fltinv_lsqr_module.f90 \\
                 src/fltinv_anneal_module.f90
FLTINV_SUBS = src/okada92subs.f src/geomsubs.f src/randsubs.f src/nnls.f90 \
              src/pnpoly.f
\$(BIN)/fltinv: src/fltinv.f90 \$(FLTINV_MODULES) \$(FLTINV_SUBS) \$(INCLUDE_FILES)
	\$(FC) \$(FFLAG) -c \$(FLTINV_MODULES) \$(LAPACK) \$(CPP) -I\$(INCLUDE_DIR)
	\$(FC) \$(FFLAG) -o \$@ \$< \$(FLTINV_MODULES) \$(FLTINV_SUBS) \$(LAPACK) \$(CPP) -I\$(INCLUDE_DIR) -I. \$(INCLUDE_FILES)
	rm *.o *.mod

\$(BIN)/grid: src/grid.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/lola2distaz: src/lola2distaz.f src/geomsubs.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/mtutil: src/mtutil.f src/mtsubs.f
	\$(FC) \$(FFLAG) -o \$(BIN)/mtutil src/mtutil.f src/mtsubs.f \$(LAPACK) \$(CPP)

\$(BIN)/multifit: src/multifit.f
	\$(FC) \$(FFLAG) -o \$(BIN)/multifit src/multifit.f \$(LAPACK) \$(CPP)

\$(BIN)/numint: src/numint.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/o92util: src/o92util.f src/okada92subs.f src/geomsubs.f src/okada92subs_volume.f
	\$(FC) \$(FFLAG) -o \$(BIN)/o92util src/o92util.f src/okada92subs.f src/geomsubs.f src/okada92subs_volume.f

\$(BIN)/perturb: src/perturb.f src/randsubs.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/platemotion: src/platemotion.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/polyfit: src/polyfit.f src/lsqsubs.f
	\$(FC) \$(FFLAG) -o \$(BIN)/polyfit src/polyfit.f src/lsqsubs.f \$(LAPACK) \$(CPP)

\$(BIN)/polyfit_special: src/polyfit_special.f
	\$(FC) \$(FFLAG) -o \$(BIN)/polyfit_special src/polyfit_special.f \$(LAPACK) \$(CPP)

\$(BIN)/readGCMT: src/readGCMT.f90 src/geomsubs.f
	\$(FC) \$(FFLAG) -o \$(BIN)/readGCMT src/readGCMT.f90 src/geomsubs.f
	rm io.mod event_data.mod

\$(BIN)/readkik: src/readkik.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/sphfinrot: src/sphfinrot.f90
	\$(FC) \$(FFLAG) -o \$(BIN)/sphfinrot src/sphfinrot.f90
	rm sphfinrot*.mod

\$(BIN)/stereo_project: src/stereo_project.f90
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/triutil: src/triutil.f90 src/pnpoly.f src/geomsubs.f \$(INCLUDE_FILES)
	\$(FC) \$(FFLAG) -o \$@ \$^ -I\$(INCLUDE_DIR)
	rm triutil.mod

\$(BIN)/utm2geo: src/utm2geo.f src/geomsubs.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/vec2los: src/vec2los.f
	\$(FC) \$(FFLAG) -o \$@ \$^

\$(BIN)/wraplos: src/wraplos.f
	\$(FC) \$(FFLAG) -o \$@ \$^

##############################
##### SHELL SCRIPT RULES #####
##############################
\$(BIN)/coul_dip.sh: scripts/coul_dip.sh
	cp scripts/coul_dip.sh \$(BIN)/coul_dip.sh

\$(BIN)/coul_hor.sh: scripts/coul_hor.sh
	cp scripts/coul_hor.sh \$(BIN)/coul_hor.sh

\$(BIN)/coul_xsec.sh: scripts/coul_xsec.sh
	cp scripts/coul_xsec.sh \$(BIN)/coul_xsec.sh

\$(BIN)/gmtcpt.sh: scripts/gmtcpt.sh
	cp scripts/gmtcpt.sh \$(BIN)/gmtcpt.sh

\$(BIN)/simplify_ffm.sh: scripts/simplify_ffm.sh
	cp scripts/simplify_ffm.sh \$(BIN)/simplify_ffm.sh

\$(BIN)/surf_disp.sh: scripts/surf_disp.sh
	cp scripts/surf_disp.sh \$(BIN)/surf_disp.sh

\$(BIN)/ternary.sh: scripts/ternary.sh
	cp scripts/ternary.sh \$(BIN)/ternary.sh

\$(BIN)/trg_schem.sh: scripts/trg_schem.sh
	cp scripts/trg_schem.sh \$(BIN)/trg_schem.sh

###########################
##### UNIT TEST RULES #####
###########################
test: \\
      test_tri_disloc \
      test_okada92

test_tri_disloc: src/unit_test_tri_disloc.f90 \$(INCLUDE_FILES)
	@\$(FC) \$^ \$(ALL_FLAGS) -o \$@
	\$@ > test_tri_disloc.log
	@echo "test_tri_disloc passed"
	rm \$@ test_tri_disloc.log
	@echo

test_okada92: src/unit_test_okada92.f90 \$(INCLUDE_FILES)
	@\$(FC) \$^ \$(ALL_FLAGS) -o \$@
	\$@ > test_okada92.log
	@echo "test_okada92 passed"
	rm \$@ test_okada92.log

####################
##### CLEAN UP #####
####################
# Clean bin directory
.PHONY: clean
clean:
	-rm \$(BIN)/*
	-rm \$(INCLUDE_DIR)/*


EOF
