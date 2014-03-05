FC   = gfortran
FFLG = -Wall -Wextra

BIN  = ../bin
OBJ  = ./obj

OKADA = $(OBJ)/okada92subs.o
GEOM  = $(OBJ)/geomsubs.o
ELAST = $(OBJ)/elastsubs.o
SVD   = $(OBJ)/svdsubs.o


all: $(BIN)/ff2displacement $(BIN)/ff2coulomb $(BIN)/flt2displacement $(BIN)/flt2coulomb \
     $(BIN)/fltxy2displacement \
     $(BIN)/makefaultstxt $(BIN)/ff2gmt $(BIN)/grid $(BIN)/nez2los $(BIN)/wraplos \
     $(BIN)/sdr2sv $(BIN)/daycount $(BIN)/lola2distaz $(BIN)/distaz2lola \
     $(BIN)/simpleinversion \
     $(BIN)/sphfinrot \
     $(BIN)/readkik \
     $(BIN)/pole2vel \
     $(BIN)/eventfrequency

$(BIN)/ff2displacement:    $(OBJ)/ff2displacement.o $(OKADA) $(GEOM)
	$(FC) -o $(BIN)/ff2displacement $(OBJ)/ff2displacement.o $(OKADA) $(GEOM)

$(BIN)/ff2coulomb:         $(OBJ)/ff2coulomb.o $(OKADA) $(GEOM) $(ELAST)
	$(FC) -o $(BIN)/ff2coulomb $(OBJ)/ff2coulomb.o $(OKADA) $(GEOM) $(ELAST)

$(BIN)/ff2gmt:             $(OBJ)/ff2gmt.o
	$(FC) -o $(BIN)/ff2gmt $(OBJ)/ff2gmt.o

$(BIN)/flt2displacement:   $(OBJ)/flt2displacement.o $(OKADA) $(GEOM)
	$(FC) -o $(BIN)/flt2displacement $(OBJ)/flt2displacement.o $(OKADA) $(GEOM)

$(BIN)/fltxy2displacement: $(OBJ)/fltxy2displacement.o $(OKADA) $(GEOM)
	$(FC) -o $(BIN)/fltxy2displacement $(OBJ)/fltxy2displacement.o $(OKADA) $(GEOM)

$(BIN)/flt2coulomb:        $(OBJ)/flt2coulomb.o $(OKADA) $(GEOM) $(ELAST)
	$(FC) -o $(BIN)/flt2coulomb $(OBJ)/flt2coulomb.o $(OKADA) $(GEOM) $(ELAST)

$(BIN)/testokada:          $(OBJ)/testokada.o $(OKADA)
	$(FC) -o $(BIN)/testokada $(OBJ)/testokada.o $(OKADA)

$(BIN)/grid:               $(OBJ)/grid.o
	$(FC) -o $(BIN)/grid $(OBJ)/grid.o

$(BIN)/nez2los:            $(OBJ)/nez2los.o
	$(FC) -o $(BIN)/nez2los $(OBJ)/nez2los.o

$(BIN)/wraplos:            $(OBJ)/wraplos.o
	$(FC) -o $(BIN)/wraplos $(OBJ)/wraplos.o

$(BIN)/sdr2sv:             $(OBJ)/sdr2sv.o
	$(FC) -o $(BIN)/sdr2sv $(OBJ)/sdr2sv.o

$(BIN)/daycount:           $(OBJ)/daycount.o
	$(FC) -o $(BIN)/daycount $(OBJ)/daycount.o

$(BIN)/lola2distaz:        $(OBJ)/lola2distaz.o $(GEOM)
	$(FC) -o $(BIN)/lola2distaz $(OBJ)/lola2distaz.o $(OBJ)/geomsubs.o

$(BIN)/distaz2lola:        $(OBJ)/distaz2lola.o $(GEOM)
	$(FC) -o $(BIN)/distaz2lola $(OBJ)/distaz2lola.o $(OBJ)/geomsubs.o

$(BIN)/pole2vel:        $(OBJ)/pole2vel.o
	$(FC) -o $(BIN)/pole2vel $(OBJ)/pole2vel.o

$(BIN)/makefaultstxt:      $(OBJ)/makefaultstxt.o
	$(FC) -o $(BIN)/makefaultstxt $(OBJ)/makefaultstxt.o

$(BIN)/simpleinversion:    $(OBJ)/simpleinversion.o $(SVD) $(OKADA) $(GEOM)
	$(FC) -o $(BIN)/simpleinversion $(OBJ)/simpleinversion.o $(SVD) $(OKADA) $(GEOM)

$(BIN)/readkik:            $(OBJ)/readkik.o
	$(FC) -o $(BIN)/readkik $(OBJ)/readkik.o

$(BIN)/sphfinrot:          $(OBJ)/sphfinrot.o
	$(FC) -o $(BIN)/sphfinrot $(OBJ)/sphfinrot.o

$(BIN)/eventfrequency:       $(OBJ)/eventfrequency.o
	$(FC) -o $(BIN)/eventfrequency $(OBJ)/eventfrequency.o

$(OBJ)/%.o: %.f
	$(FC) -c $<
	mv *.o $(OBJ)

clean:
	rm -f *.o
