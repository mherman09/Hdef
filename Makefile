##### Compilers #####
# C compiler executable
CC = gcc

# Fortran compiler executable
FC = gfortran

# Fortran compiler warnings
FWARN = -Wall -Wextra -Wunused -fbounds-check -fbacktrace
# FWARN =

# Fortran compiler optimization
# FOPT  =
FOPT  = -O1
# FOPT  = -O2

# Fortran compiler flags
FFLAGS = $(FWARN) $(FOPT)


##### Include directory (where to put .o and .mod files) #####
INCLUDE_DIR = ./include
INCLUDE_FLAGS = -I$(INCLUDE_DIR) -J$(INCLUDE_DIR)


##### Executable directory (where to put compiled codes) #####
BIN   = ./bin


##### External libraries #####
# LAPACK (linear algebra)
LAPACK_LIB_DIR = /sw/lib/lapack
LAPACK_LIB     = -lreflapack -lrefblas -ltmglib
# LAPACK_LIB     = -llapack -lblas -ltmglib
LAPACK         = -L$(LAPACK_LIB_DIR) $(LAPACK_LIB)
CPP_LAPACK = -DUSE_LAPACK
# LAPACK =
# CPP_LAPACK =

# SuperLU (sparse linear equation solvers)
SUPERLU_LIB_DIR = ./ext/SuperLU_5.2.1/lib
SUPERLU_LIB     = -lsuperlu_5.1
# SUPERLU         = -L$(SUPERLU_LIB_DIR) $(SUPERLU_LIB)
# CPP_SUPERLU = -DUSE_SUPERLU
SUPERLU =
CPP_SUPERLU =


##### Pre-processing #####
CPP_UNIT_TEST = -DUNIT_TEST
# CPP_UNIT_TEST =

CPP_FLAGS = -cpp $(CPP_UNIT_TEST) $(CPP_LAPACK) $(CPP_SUPERLU)



####################################################################################################
####################################################################################################
################### DO NOT CHANGE BELOW HERE UNLESS YOU KNOW WHAT YOU ARE DOING ####################
####################################################################################################
####################################################################################################



##################################
##### COMPILED PROGRAM RULES #####
##################################
all: \
    $(BIN) \
    $(BIN)/colortool \
    $(BIN)/dateutil \
    $(BIN)/distaz2lola \
    $(BIN)/ff2gmt \
    $(BIN)/fltinv \
    $(BIN)/grid \
    $(BIN)/lola2distaz \
    $(BIN)/mtutil \
    $(BIN)/numint \
    $(BIN)/o92util \
    $(BIN)/platemotion \
    $(BIN)/rangen \
    $(BIN)/readGCMT \
    $(BIN)/readkik \
    $(BIN)/sphfinrot \
    $(BIN)/stereo_project \
    $(BIN)/triutil \
    $(BIN)/vec2los \
    $(BIN)/wraplos \
    scripts

$(BIN):
	mkdir -p $(BIN)

$(BIN)/colortool: src/colortool.f90
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/dateutil: src/dateutil.f90 \
                 $(INCLUDE_DIR)/calendar.o \
                 $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/distaz2lola: src/distaz2lola_prog.f90 \
                    $(INCLUDE_DIR)/algebra.o \
                    $(INCLUDE_DIR)/earth.o \
                    $(INCLUDE_DIR)/geom.o \
                    $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/ff2gmt: src/ff2gmt.f90 \
               $(INCLUDE_DIR)/algebra.o \
               $(INCLUDE_DIR)/ffm.o \
               $(INCLUDE_DIR)/geom.o \
               $(INCLUDE_DIR)/io.o \
               $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/fltinv: src/fltinv.f90 \
               $(INCLUDE_DIR)/algebra.o \
               $(INCLUDE_DIR)/annealing.o \
               $(INCLUDE_DIR)/earth.o \
               $(INCLUDE_DIR)/elast.o \
               $(INCLUDE_DIR)/geom.o \
               $(INCLUDE_DIR)/io.o \
               $(INCLUDE_DIR)/okada92.o \
               $(INCLUDE_DIR)/random.o \
               $(INCLUDE_DIR)/solver.o \
               $(INCLUDE_DIR)/trig.o \
               $(INCLUDE_DIR)/tri_disloc.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $(SUPERLU) $^ -o $@

$(BIN)/grid: src/grid.f90 \
             $(INCLUDE_DIR)/algebra.o \
             $(INCLUDE_DIR)/earth.o \
             $(INCLUDE_DIR)/geom.o \
             $(INCLUDE_DIR)/io.o \
             $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/lola2distaz: src/lola2distaz_prog.f90 \
                    $(INCLUDE_DIR)/algebra.o \
                    $(INCLUDE_DIR)/earth.o \
                    $(INCLUDE_DIR)/geom.o \
                    $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/mtutil: src/mtutil.f90 \
               $(INCLUDE_DIR)/eq.o \
               $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $^ -o $@

$(BIN)/numint: src/numint.f
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/o92util: src/o92util.f90 \
                $(INCLUDE_DIR)/algebra.o \
                $(INCLUDE_DIR)/earth.o \
                $(INCLUDE_DIR)/elast.o \
                $(INCLUDE_DIR)/eq.o \
                $(INCLUDE_DIR)/ffm.o \
                $(INCLUDE_DIR)/geom.o \
                $(INCLUDE_DIR)/io.o \
                $(INCLUDE_DIR)/okada92.o \
                $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $^ -o $@

$(BIN)/platemotion: src/platemotion.f90 \
                    $(INCLUDE_DIR)/algebra.o \
                    $(INCLUDE_DIR)/earth.o \
                    $(INCLUDE_DIR)/io.o \
                    $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/rangen: src/rangen.f90 \
               $(INCLUDE_DIR)/io.o \
               $(INCLUDE_DIR)/random.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/readGCMT: src/readGCMT.f90 \
                 $(INCLUDE_DIR)/algebra.o \
                 $(INCLUDE_DIR)/calendar.o \
                 $(INCLUDE_DIR)/earth.o \
                 $(INCLUDE_DIR)/eq.o \
                 $(INCLUDE_DIR)/geom.o \
                 $(INCLUDE_DIR)/io.o \
                 $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $^ -o $@

$(BIN)/readkik: src/readkik.f
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $^ -o $@

$(BIN)/sphfinrot: src/sphfinrot.f90 \
                  $(INCLUDE_DIR)/algebra.o \
                  $(INCLUDE_DIR)/earth.o \
                  $(INCLUDE_DIR)/io.o \
                  $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/stereo_project: src/stereo_project.f90 \
                       $(INCLUDE_DIR)/algebra.o \
                       $(INCLUDE_DIR)/earth.o \
                       $(INCLUDE_DIR)/io.o \
                       $(INCLUDE_DIR)/map.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/triutil: src/triutil.f90 \
                $(INCLUDE_DIR)/algebra.o \
                $(INCLUDE_DIR)/geom.o \
                $(INCLUDE_DIR)/earth.o \
                $(INCLUDE_DIR)/elast.o \
                $(INCLUDE_DIR)/eq.o \
                $(INCLUDE_DIR)/geom.o \
                $(INCLUDE_DIR)/io.o \
                $(INCLUDE_DIR)/trig.o \
                $(INCLUDE_DIR)/tri_disloc.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $^ -o $@

$(BIN)/vec2los: src/vec2los.f90 \
                $(INCLUDE_DIR)/io.o \
                $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@

$(BIN)/wraplos: src/wraplos.f90 \
                $(INCLUDE_DIR)/io.o \
                $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $^ -o $@


scripts: \
    $(BIN)/coul_dip.sh \
    $(BIN)/coul_hor.sh \
    $(BIN)/coul_xsec.sh \
    $(BIN)/gmtcpt.sh \
    $(BIN)/simplify_ffm.sh \
    $(BIN)/surf_disp.sh \
    $(BIN)/ternary.sh \
    $(BIN)/trg_schem.sh

$(BIN)/coul_dip.sh: ./scripts/coul_dip.sh
	cp $^ $@
$(BIN)/coul_hor.sh: ./scripts/coul_hor.sh
	cp $^ $@
$(BIN)/coul_xsec.sh: ./scripts/coul_xsec.sh
	cp $^ $@
$(BIN)/gmtcpt.sh: ./scripts/gmtcpt.sh
	cp $^ $@
$(BIN)/simplify_ffm.sh: ./scripts/simplify_ffm.sh
	cp $^ $@
$(BIN)/surf_disp.sh: ./scripts/surf_disp.sh
	cp $^ $@
$(BIN)/ternary.sh: ./scripts/ternary.sh
	cp $^ $@
$(BIN)/trg_schem.sh: ./scripts/trg_schem.sh
	cp $^ $@

####################################
##### MODULE AND OBJECT RULES ######
####################################
INCLUDE_FILES = \
    $(INCLUDE_DIR)/algebra.o \
    $(INCLUDE_DIR)/annealing.o \
    $(INCLUDE_DIR)/calendar.o \
    $(INCLUDE_DIR)/earth.o \
    $(INCLUDE_DIR)/elast.o \
    $(INCLUDE_DIR)/eq.o \
    $(INCLUDE_DIR)/error_exit.o \
    $(INCLUDE_DIR)/ffm.o \
    $(INCLUDE_DIR)/geom.o \
    $(INCLUDE_DIR)/io.o \
    $(INCLUDE_DIR)/map.o \
    $(INCLUDE_DIR)/okada92.o \
    $(INCLUDE_DIR)/random.o \
    $(INCLUDE_DIR)/solver.o \
    $(INCLUDE_DIR)/test.o \
    $(INCLUDE_DIR)/trig.o \
    $(INCLUDE_DIR)/tri_disloc.o

include: $(INCLUDE_FILES)

$(INCLUDE_DIR)/algebra.o: src/algebra_module.f90 \
                          $(INCLUDE_DIR)/io.o \
                          $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/annealing.o: src/annealing_module.f90 \
                            $(INCLUDE_DIR)/io.o \
                            $(INCLUDE_DIR)/random.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/calendar.o: src/calendar_module.f90 \
                           $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/earth.o: src/earth_module.f90 \
                        $(INCLUDE_DIR)/algebra.o \
                        $(INCLUDE_DIR)/io.o \
                        $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/elast.o: src/elast_module.f90 \
                        $(INCLUDE_DIR)/algebra.o \
                        $(INCLUDE_DIR)/geom.o \
                        $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $(LAPACK) $< -c -o $@

$(INCLUDE_DIR)/eq.o: src/eq_module.f90 \
                     $(INCLUDE_DIR)/io.o \
                     $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $(LAPACK) $< -c -o $@

$(INCLUDE_DIR)/error_exit.o: src/error_exit.c
	$(CC) $<  -c -o $@

$(INCLUDE_DIR)/ffm.o: src/ffm_module.f90 \
                      $(INCLUDE_DIR)/geom.o \
                      $(INCLUDE_DIR)/io.o \
                      $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/geom.o: src/geom_module.f90 \
                       $(INCLUDE_DIR)/algebra.o \
                       $(INCLUDE_DIR)/io.o \
                       $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/io.o: src/io_module.f90
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/okada92.o: src/okada92_module.f90 \
                          $(INCLUDE_DIR)/io.o \
                          $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/map.o: src/map_module.f90 \
                      $(INCLUDE_DIR)/earth.o \
                      $(INCLUDE_DIR)/io.o \
                      $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/random.o: src/random_module.f90
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $(LAPACK) $< -c -o $@

$(INCLUDE_DIR)/solver.o: src/solver_module.f90 \
                         $(INCLUDE_DIR)/io.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $(LAPACK) $< -c -o $@

$(INCLUDE_DIR)/test.o: src/test_module.f90 \
                       $(INCLUDE_DIR)/error_exit.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/trig.o: src/trig_module.f90
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/tri_disloc.o: src/tri_disloc_module.f90 \
                             $(INCLUDE_DIR)/earth.o \
                             $(INCLUDE_DIR)/geom.o \
                             $(INCLUDE_DIR)/io.o \
                             $(INCLUDE_DIR)/trig.o
	$(FC) $(FFLAGS) -J$(INCLUDE_DIR) $(CPP_FLAGS) $< -c -o $@

$(INCLUDE_DIR)/%.o: $(INCLUDE_DIR)/%.mod


###########################
##### UNIT TEST RULES #####
###########################
test: \
      test_annealing \
      test_calendar \
      test_earth \
      test_eq \
      test_ffm \
      test_geom \
      test_map \
      test_okada92 \
      test_random \
      test_solver \
      test_tri_disloc \
      test_executables

test_annealing: src/unit_test_annealing.f90 \
                $(INCLUDE_DIR)/annealing.o \
                $(INCLUDE_DIR)/random.o \
                $(INCLUDE_DIR)/test.o \
                $(INCLUDE_DIR)/error_exit.o \
                $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_annealing.log
	@echo "test_annealing passed"
	rm $@ test_annealing.log include/driver*.mod
	@echo

test_calendar: src/unit_test_calendar.f90 \
               $(INCLUDE_DIR)/calendar.o \
               $(INCLUDE_DIR)/test.o \
               $(INCLUDE_DIR)/error_exit.o \
               $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_calendar.log
	@echo "test_calendar passed"
	rm $@ test_calendar.log
	@echo

test_earth: src/unit_test_earth.f90 \
            $(INCLUDE_DIR)/algebra.o \
            $(INCLUDE_DIR)/earth.o \
            $(INCLUDE_DIR)/test.o \
            $(INCLUDE_DIR)/error_exit.o \
            $(INCLUDE_DIR)/io.o \
            $(INCLUDE_DIR)/trig.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) -o $@
	$@ > test_earth.log
	@echo "test_earth passed"
	rm $@ test_earth.log
	@echo

test_eq: src/unit_test_eq.f90 \
         $(INCLUDE_DIR)/eq.o \
         $(INCLUDE_DIR)/test.o \
         $(INCLUDE_DIR)/error_exit.o \
         $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) -o $@
	$@ > test_eq.log
	@echo "test_eq passed"
	rm $@ test_eq.log
	@echo

test_executables: all
	cd ./test/ && make test

test_ffm: src/unit_test_ffm.f90 \
          $(INCLUDE_DIR)/algebra.o \
          $(INCLUDE_DIR)/ffm.o \
          $(INCLUDE_DIR)/geom.o \
          $(INCLUDE_DIR)/test.o \
          $(INCLUDE_DIR)/error_exit.o \
          $(INCLUDE_DIR)/io.o
	@$(FC) -ffree-line-length-0 $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_ffm.log
	@echo "test_ffm passed"
	rm $@ test_ffm.log
	@echo

test_geom: src/unit_test_geom.f90 \
           $(INCLUDE_DIR)/algebra.o \
           $(INCLUDE_DIR)/geom.o \
           $(INCLUDE_DIR)/test.o \
           $(INCLUDE_DIR)/error_exit.o \
           $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_geom.log
	@echo "test_geom passed"
	rm $@ test_geom.log
	@echo

test_map: src/unit_test_map.f90 \
           $(INCLUDE_DIR)/test.o \
           $(INCLUDE_DIR)/error_exit.o \
           $(INCLUDE_DIR)/io.o \
           $(INCLUDE_DIR)/map.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_map.log
	@echo "test_map passed"
	rm $@ test_map.log
	@echo

test_okada92: src/unit_test_okada92.f90 \
              $(INCLUDE_DIR)/okada92.o \
              $(INCLUDE_DIR)/test.o \
              $(INCLUDE_DIR)/trig.o \
              $(INCLUDE_DIR)/error_exit.o \
              $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_okada92.log
	@echo "test_okada92 passed"
	rm $@ test_okada92.log
	@echo

test_random: src/unit_test_random.f90 \
             $(INCLUDE_DIR)/random.o \
             $(INCLUDE_DIR)/test.o \
             $(INCLUDE_DIR)/error_exit.o \
             $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_random.log
	@echo "test_random passed"
	rm $@ test_random.log
	@echo

test_solver: src/unit_test_solver.f90 \
             $(INCLUDE_DIR)/solver.o \
             $(INCLUDE_DIR)/test.o \
             $(INCLUDE_DIR)/error_exit.o \
             $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) $(LAPACK) $(SUPERLU) -o $@
	$@ > test_solver.log
	@echo "test_solver passed"
	rm $@ test_solver.log
	@echo

test_tri_disloc: src/unit_test_tri_disloc.f90 \
                 $(INCLUDE_DIR)/algebra.o \
                 $(INCLUDE_DIR)/tri_disloc.o \
                 $(INCLUDE_DIR)/geom.o \
                 $(INCLUDE_DIR)/test.o \
                 $(INCLUDE_DIR)/trig.o \
                 $(INCLUDE_DIR)/error_exit.o \
                 $(INCLUDE_DIR)/io.o
	@$(FC) $^ $(FFLAGS) $(INCLUDE_FLAGS) $(CPP_FLAGS) -o $@
	$@ > test_tri_disloc.log
	@echo "test_tri_disloc passed"
	rm $@ test_tri_disloc.log
	@echo

####################
##### CLEAN UP #####
####################
# Clean bin directory
.PHONY: clean
clean:
	-rm $(BIN)/*
	-rm $(INCLUDE_DIR)/*
