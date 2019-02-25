#!/bin/bash

# Configure script command line options
function usage() {
    echo "$0 [-f=FC] [-l=/PATH/TO/LAPACK/LIBRARIES] [-b=/PATH/TO/EXECUTABLES]"
    echo "   [-i] [-d]"
    echo
    echo "-f|--fortran_compiler=FC                   Fortran compiler"
    echo "-l|--lapack_dir=/PATH/TO/LAPACK/LIBRARIES  Location of LAPACK libraries"
    echo "-b|--bin_dir=/PATH/TO/EXECUTABLES          Path for installing executables"
    echo "-i|--interactive                           Interactive prompts for inputs"
    echo "-d|--default                               Use Matt's default values (FC=gfortran, LAPACK_DIR=/sw/lib/lapack, BIN_DIR=bin)"
    echo
    echo "Note: if working in a root directory or trying to install programs in a root directory, may have to use \"sudo $0\""
    exit 1
}
if [ $# -eq 0 ]
then
    usage
fi

# Parse command line
FC=""
LAPACK_LIB_DIR=""
BIN_DIR=""
INTERACTIVE="N"
while [ "$1" != "" ]
do
    case $1 in
        -f=*|--fortran_compiler=*)FC=`echo $1 | sed -e "s/.*=//"`;;
        -l=*|--lapack_dir=*)LAPACK_LIB_DIR=`echo $1 | sed -e "s/.*=//"`;;
        -b=*|--bin_dir=*)BIN_DIR=`echo $1 | sed -e "s/.*=//"`;;
        -i|--interactive)INTERACTIVE="Y";;
        -d|--default)FC="gfortran";LAPACK_LIB_DIR="/sw/lib/lapack";BIN_DIR="bin";;
        *)echo "!! Error: No option $1"; usage;;
    esac
    shift
done

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
    echo "!! Error: no executable found named $FC"
    echo "!! Check to see if $FC is in your PATH or full path is written correctly"
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
    echo "!! Error: Fortran compiler does not compile Fortran 77 test example!"
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
    echo "!! Error: Fortran compiler does not compile Fortran 90 test example!"
    exit 1
else
    echo "Test complete! Using Fortran compiler: $FC"
fi
# Clean up
rm $FC_TEST.f90 $FC_TEST

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
    LAPACK_CPP_OPTION=""
else
    # Make sure the directory exists and has the correct libraries in it
    LAPACK_LIB_DIR_CONTENTS=`ls $LAPACK_LIB_DIR 2>&1`
    EXIST=`echo $LAPACK_LIB_DIR_CONTENTS |\
           awk '{if(/No such file or directory/){print "N";exit}else{print "Y";exit}}'`

    if [ "$EXIST" == "N" ]
    then
        echo "!! Error: Directory $LAPACK_LIB_DIR is non-existent"
        exit 1
    fi

    if [ -z "$LAPACK_LIB_DIR_CONTENTS" ]
    then
        echo "!! Error: Directory $LAPACK_LIB_DIR is empty"
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
                echo "!! Error: LAPACK library \"$LIB2\" is not in $LAPACK_LIB_DIR"
                exit 1
            fi

            echo "Found LAPACK library: $LIB2!"
            LAPACK_LIB_LIST="$LAPACK_LIB_LIST -l$LIB2"
        done
    fi
    echo "Test complete! Using LAPACK libraries in directory: $LAPACK_LIB_DIR"
    LAPACK_CPP_OPTION="-DUSELAPACK"
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
    echo "Enter the path to place the compiled executables (Default: bin):"
    read ANSWER
    if [ -z "$ANSWER" ]
    then
        ANSWER="bin"
    fi
    BIN_DIR=$ANSWER
fi
echo Putting exectuables in $BIN_DIR
if [ ! -d $BIN_DIR ]
then
    echo "Directory $BIN_DIR does not exist....creating it"
    mkdir $BIN_DIR
else
    echo "Directory $BIN_DIR already exists"
fi
echo
echo

echo "####################################################################################################"
echo "##########                           INSTALL THE CODES!                                   ##########"
echo "####################################################################################################"
echo "You used the options: -f=${FC} -l=${LAPACK_LIB_DIR} -b=${BIN_DIR}"
echo "Type \"make\" to install codes"

if [ "$LAPACK_LIB_DIR" != "" ]
then
    LAPACK_LIB_DIR="-L${LAPACK_LIB_DIR}"
fi

#####
#	Generating makefile
#####
cat > Makefile << EOF
##### Compiler variables #####
FC = -$FC
FC = $FC
FWARN = -Wall -Wextra -Wunused -fbounds-check -fbacktrace
FOPT  = -O1
FFLAG = \$(FWARN) \$(FOPT)

##### Include directory (.o and .mod files) #####
INCLUDE = include

##### Executable directory #####
BIN   = $BIN_DIR

##### External libraries #####
LAPACK_LIB_DIR = $LAPACK_LIB_DIR
LAPACK_LIB     = $LAPACK_LIB_LIST
LAPACK         = \$(LAPACK_LIB_DIR) \$(LAPACK_LIB)
CPP_LAPACK     = $LAPACK_CPP_OPTION

CPP = \$(CPP_LAPACK) -cpp

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
                 src/fltinv_gf_module.f90 \\
                 src/fltinv_lsqr_module.f90 \\
                 src/fltinv_anneal_module.f90
FLTINV_SUBS = src/okada92subs.f src/geomsubs.f src/randsubs.f src/nnls.f90 \
              src/pnpoly.f
FLTINV_INCLUDE = \$(INCLUDE)/tri_disloc.o \\
                 \$(INCLUDE)/elast.o \\
                 \$(INCLUDE)/io.o
# SUPERLU = -Lext/SuperLU_5.2.1/lib -lsuperlu_5.1
\$(BIN)/fltinv: src/fltinv.f90 \$(FLTINV_MODULES) \$(FLTINV_SUBS) \$(FLTINV_INCLUDE)
	\$(FC) \$(FFLAG) -c \$(FLTINV_MODULES) \$(LAPACK) \$(CPP) -I\$(INCLUDE)
	\$(FC) \$(FFLAG) -o \$@ \$< \$(FLTINV_MODULES) \$(FLTINV_SUBS) \$(LAPACK) \$(SUPERLU) \$(CPP) -I\$(INCLUDE)  \$(FLTINV_INCLUDE)
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

\$(BIN)/triutil: src/triutil.f90 src/pnpoly.f src/geomsubs.f \$(INCLUDE)/tri_disloc.o
	\$(FC) \$(FFLAG) -o \$@ \$^ -I\$(INCLUDE)
	rm triutil_module.mod

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

test_tri_disloc: src/tri_disloc_unit_tests.f90 src/pnpoly.f src/geomsubs.f \$(INCLUDE)/tri_disloc.o
	\$(FC) \$(FFLAG) -o \$@ \$^ -I\$(INCLUDE)
	\$@
	rm \$@

test_okada92: src/okada92_module.f90 src/okada92_unit_tests.f90 src/okada92subs.f \$(INCLUDE)/trig.o \$(INCLUDE)/test.o
	\$(FC) \$(FFLAG) -c src/okada92_module.f90 -I\$(INCLUDE)
	\$(FC) \$(FFLAG) -o \$@ \$^ -I\$(INCLUDE)
	rm *.o *.mod
	\$@
	@echo "test_okada92 passed"
	rm \$@

################################
##### MODULES AND OBJECTS ######
################################
\$(INCLUDE)/io.o: src/io_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$<

\$(INCLUDE)/trig.o: src/trig_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$<

\$(INCLUDE)/earth.o: src/earth_module.f90 \$(INCLUDE)/trig.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$< -I\$(INCLUDE)

\$(INCLUDE)/tri_disloc.o: src/tri_disloc_module.f90 \$(INCLUDE)/trig.o
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$< -I\$(INCLUDE)

\$(INCLUDE)/elast.o: src/elast_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$< -I\$(INCLUDE)

\$(INCLUDE)/test.o: src/test_module.f90
	\$(FC) \$(FFLAG) -J\$(INCLUDE) -c -o \$@ \$< -I\$(INCLUDE)

\$(INCLUDE)/%.o: \$(INCLUDE)/%.mod

####################
##### CLEAN UP #####
####################
# Clean bin directory
.PHONY: clean
clean:
	-rm \$(BIN)/*
	-rm \$(INCLUDE)/*


EOF
