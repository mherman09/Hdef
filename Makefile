##### Compiler variables #####
FC = gfortran
#FC = ifort
#FC = pgfortran
FWARN = -Wall -Wextra -Wunused -fbounds-check -fbacktrace
#FOPT  = -O1
FFLAG = $(FWARN) $(FOPT)

##### Executable directory #####
BIN   = ../bin

##### External libraries #####
LAPACK_LIB_DIR = -L/Users/mherman/Research/lapack-3.5.0
LAPACK_LIB     = -llapack -ltmglib -lrefblas
LAPACK         = $(LAPACK_LIB_DIR) $(LAPACK_LIB)

##### Rules #####
all: \
     $(BIN)/o92util \
     $(BIN)/dateutil \
     $(BIN)/mtutil \
     $(BIN)/lola2distaz \
     $(BIN)/distaz2lola \
     $(BIN)/polyfit \
     $(BIN)/polyfit_special \
     $(BIN)/eventfrequency \
     $(BIN)/vec2los \
     $(BIN)/wraplos \
     $(BIN)/sphfinrot \
     $(BIN)/ff2gmt \
     $(BIN)/grid \
     $(BIN)/pt2fin  \
     $(BIN)/readkik \
     $(BIN)/numint \
     $(BIN)/fltinv \
     $(BIN)/platemotion  \
     $(BIN)/perturb \
     $(BIN)/pole2vel \
     $(BIN)/multifit
#     $(BIN)/trigger $(BIN)/twocol2asc \
#     $(BIN)/smooth

$(BIN)/o92util: o92util.f okada92subs.f geomsubs.f okada92subs_volume.f
	$(FC) $(FFLAG) -o $(BIN)/o92util o92util.f okada92subs.f geomsubs.f okada92subs_volume.f

$(BIN)/dateutil: dateutil.f
	$(FC) $(FFLAG) -o $(BIN)/dateutil dateutil.f

$(BIN)/mtutil: mtutil.f mtsubs.f
	$(FC) $(FFLAG) -o $(BIN)/mtutil mtutil.f mtsubs.f $(LAPACK)

$(BIN)/lola2distaz: lola2distaz.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/lola2distaz lola2distaz.f geomsubs.f

$(BIN)/distaz2lola: distaz2lola.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/distaz2lola distaz2lola.f geomsubs.f

$(BIN)/polyfit: polyfit.f lsqsubs.f
	$(FC) $(FFLAG) -o $(BIN)/polyfit polyfit.f lsqsubs.f $(LAPACK)

$(BIN)/polyfit_special: polyfit_special.f
	$(FC) $(FFLAG) -o $(BIN)/polyfit_special polyfit_special.f $(LAPACK)

$(BIN)/eventfrequency: eventfrequency.f
	$(FC) $(FFLAG) -o $(BIN)/eventfrequency eventfrequency.f

$(BIN)/vec2los: vec2los.f
	$(FC) $(FFLAG) -o $(BIN)/vec2los vec2los.f

$(BIN)/wraplos: wraplos.f
	$(FC) $(FFLAG) -o $(BIN)/wraplos wraplos.f

$(BIN)/sphfinrot: sphfinrot.f
	$(FC) $(FFLAG) -o $(BIN)/sphfinrot sphfinrot.f

$(BIN)/ff2gmt: ff2gmt.f
	$(FC) $(FFLAG) -o $(BIN)/ff2gmt ff2gmt.f

$(BIN)/grid: grid.f
	$(FC) $(FFLAG) -o $(BIN)/grid grid.f

$(BIN)/pt2fin: pt2fin.f
	$(FC) $(FFLAG) -o $(BIN)/pt2fin pt2fin.f

$(BIN)/readkik: readkik.f
	$(FC) $(FFLAG) -o $(BIN)/readkik readkik.f

$(BIN)/numint: numint.f
	$(FC) $(FFLAG) -o $(BIN)/numint numint.f

$(BIN)/fltinv: fltinv.f okada92subs.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/fltinv fltinv.f okada92subs.f geomsubs.f $(LAPACK) randsubs.f 

$(BIN)/platemotion: platemotion.f
	$(FC) $(FFLAG) -o $(BIN)/platemotion platemotion.f

$(BIN)/perturb: perturb.f
	$(FC) $(FFLAG) -o $(BIN)/perturb perturb.f randsubs.f

$(BIN)/pole2vel: pole2vel.f
	$(FC) $(FFLAG) -o $(BIN)/pole2vel pole2vel.f

$(BIN)/multifit: multifit.f
	$(FC) $(FFLAG) -o $(BIN)/multifit multifit.f $(LAPACK)

clean:
	rm $(BIN)/o92util
	rm $(BIN)/dateutil
	rm $(BIN)/mtutil
	rm $(BIN)/lola2distaz
	rm $(BIN)/distaz2lola
	rm $(BIN)/polyfit
	rm $(BIN)/eventfrequency
	rm $(BIN)/vec2los
	rm $(BIN)/wraplos
	rm $(BIN)/sphfinrot
	rm $(BIN)/ff2gmt
	rm $(BIN)/grid
	rm $(BIN)/pt2fin
	rm $(BIN)/readkik
	rm $(BIN)/numint
	rm $(BIN)/fltinv
	rm $(BIN)/platemotion 
	rm $(BIN)/perturb
	rm $(BIN)/multifit


