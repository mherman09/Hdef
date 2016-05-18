##### Compiler variables #####
FC    = gfortran
#FC    = ifort
#FC    = pgfortran
FWARN = -Wall -Wextra -Wunused -fbounds-check -fbacktrace
#FOPT  = -O1
FFLAG = $(FWARN) $(FOPT)

##### Directories #####
BIN   = ../bin

##### External libraries #####
LAPACK_LIB_DIR = -L/Users/mherman/Research/lapack-3.5.0
LAPACK_LIB     = -llapack -ltmglib -lrefblas
LAPACK         = $(LAPACK_LIB_DIR) $(LAPACK_LIB)

##### Rules #####
all: \
     $(BIN)/o92util \
     $(BIN)/dateutil \
     $(BIN)/lola2distaz $(BIN)/distaz2lola \
     $(BIN)/polyfit \
     $(BIN)/eventfrequency \
     $(BIN)/vec2los $(BIN)/wraplos \
     $(BIN)/sphfinrot \
     $(BIN)/ff2gmt $(BIN)/grid $(BIN)/pt2fin  \
     $(BIN)/readkik $(BIN)/numint \
     $(BIN)/fltinv \
     $(BIN)/nodplane2 \
     $(BIN)/mag2mom $(BIN)/mt2dc $(BIN)/mom2mag \
     $(BIN)/platemotion \
     $(BIN)/mtutil
#     $(BIN)/trigger $(BIN)/twocol2asc \
#     $(BIN)/smooth

clean:
	rm -f $(OBJ)/*.o

$(BIN)/o92util: o92util.f okada92subs.f geomsubs.f okada92subs_volume.f
	$(FC) $(FFLAG) -o $(BIN)/o92util o92util.f okada92subs.f geomsubs.f okada92subs_volume.f

$(BIN)/mtutil: mtutil.f mtsubs.f
	$(FC) $(FFLAG) -o $(BIN)/mtutil $(LAPACK) mtutil.f mtsubs.f

$(BIN)/mt2dc: mt2dc.f
	$(FC) $(FFLAG) -o $(BIN)/mt2dc $(LAPACK) mt2dc.f

$(BIN)/mom2mag: mom2mag.f
	$(FC) $(FFLAG) -o $(BIN)/mom2mag mom2mag.f

$(BIN)/mag2mom: mag2mom.f
	$(FC) $(FFLAG) -o $(BIN)/mag2mom mag2mom.f

$(BIN)/dateutil: dateutil.f
	$(FC) $(FFLAG) -o $(BIN)/dateutil dateutil.f

$(BIN)/polyfit: polyfit.f lsqsubs.f
	$(FC) $(FFLAG) -o $(BIN)/polyfit $(LAPACK) polyfit.f lsqsubs.f

$(BIN)/ff2gmt: ff2gmt.f
	$(FC) $(FFLAG) -o $(BIN)/ff2gmt ff2gmt.f

$(BIN)/grid: grid.f
	$(FC) $(FFLAG) -o $(BIN)/grid grid.f

$(BIN)/pt2fin: pt2fin.f
	$(FC) $(FFLAG) -o $(BIN)/pt2fin pt2fin.f

$(BIN)/platemotion: platemotion.f
	$(FC) $(FFLAG) -o $(BIN)/platemotion platemotion.f

$(BIN)/lola2distaz: lola2distaz.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/lola2distaz lola2distaz.f geomsubs.f

$(BIN)/distaz2lola: distaz2lola.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/distaz2lola distaz2lola.f geomsubs.f

$(BIN)/sphfinrot: sphfinrot.f
	$(FC) $(FFLAG) -o $(BIN)/sphfinrot sphfinrot.f

$(BIN)/eventfrequency: eventfrequency.f
	$(FC) $(FFLAG) -o $(BIN)/eventfrequency eventfrequency.f

$(BIN)/vec2los: vec2los.f
	$(FC) $(FFLAG) -o $(BIN)/vec2los vec2los.f

$(BIN)/wraplos: wraplos.f
	$(FC) $(FFLAG) -o $(BIN)/wraplos wraplos.f

$(BIN)/readkik: readkik.f
	$(FC) $(FFLAG) -o $(BIN)/readkik readkik.f

$(BIN)/numint: numint.f
	$(FC) $(FFLAG) -o $(BIN)/numint numint.f

$(BIN)/fltinv: fltinv.f okada92subs.f geomsubs.f
	$(FC) $(FFLAG) -o $(BIN)/fltinv fltinv.f okada92subs.f geomsubs.f $(LAPACK) randsubs.f 

$(BIN)/trigger: trigger.f
	$(FC) $(FFLAG) -o $(BIN)/trigger trigger.f

$(BIN)/twocol2asc: twocol2asc.f
	$(FC) $(FFLAG) -o $(BIN)/twocol2asc twocol2asc.f

$(BIN)/smooth: smooth.f
	$(FC) $(FFLAG) -o $(BIN)/smooth smooth.f

$(BIN)/nodplane2: nodplane2.f
	$(FC) $(FFLAG) -o $(BIN)/nodplane2 nodplane2.f
