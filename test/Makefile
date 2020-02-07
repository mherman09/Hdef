BIN_DIR := $(shell ./define_bin_dir.sh)

test: \
    test_clip \
    test_dateutil \
    test_distaz2lola \
    test_ff2gmt \
    test_fitutil \
    test_grid \
    test_lola2distaz \
    test_mtutil \
    test_o92util \
    test_platemotion \
    test_rangen \
    test_readGCMT \
    test_sphfinrot \
    test_stereo_project \
    test_surf_disp \
    test_ternary \
    test_trg_schem \
    test_triutil \
    test_vec2los \
    test_wraplos \
    test_fltinv

test_clip: test_clip.sh $(BIN_DIR)/clip test_values.sh
	./test_clip.sh

test_dateutil: test_dateutil.sh $(BIN_DIR)/dateutil test_values.sh
	./test_dateutil.sh

test_distaz2lola: test_distaz2lola.sh $(BIN_DIR)/distaz2lola test_values.sh
	./test_distaz2lola.sh

test_ff2gmt: test_ff2gmt.sh $(BIN_DIR)/ff2gmt test_values.sh
	./test_ff2gmt.sh

test_fitutil: test_fitutil.sh $(BIN_DIR)/fitutil test_values.sh
	./test_fitutil.sh

test_grid: test_grid.sh $(BIN_DIR)/grid test_values.sh
	./test_grid.sh

test_lola2distaz: test_lola2distaz.sh $(BIN_DIR)/lola2distaz test_values.sh
	./test_lola2distaz.sh

test_mtutil: test_mtutil.sh $(BIN_DIR)/mtutil test_values.sh
	./test_mtutil.sh

test_o92util: test_o92util.sh $(BIN_DIR)/o92util test_values.sh
	./test_o92util.sh

test_platemotion: test_platemotion.sh $(BIN_DIR)/platemotion test_values.sh
	./test_platemotion.sh

test_rangen: test_rangen.sh $(BIN_DIR)/rangen test_values.sh
	./test_rangen.sh

test_readGCMT: test_readGCMT.sh $(BIN_DIR)/readGCMT test_values.sh
	./test_readGCMT.sh

test_sphfinrot: test_sphfinrot.sh $(BIN_DIR)/sphfinrot test_values.sh
	./test_sphfinrot.sh

test_stereo_project: test_stereo_project.sh $(BIN_DIR)/stereo_project test_values.sh
	./test_stereo_project.sh

test_surf_disp: test_surf_disp.sh $(BIN_DIR)/surf_disp.sh test_values.sh
	./test_surf_disp.sh

test_ternary: test_ternary.sh $(BIN_DIR)/ternary test_values.sh
	./test_ternary.sh

test_trg_schem: test_trg_schem.sh $(BIN_DIR)/trg_schem.sh test_values.sh
	./test_trg_schem.sh

test_triutil: test_triutil.sh $(BIN_DIR)/triutil test_values.sh
	./test_triutil.sh

test_vec2los: test_vec2los.sh $(BIN_DIR)/vec2los test_values.sh
	./test_vec2los.sh

test_wraplos: test_wraplos.sh $(BIN_DIR)/wraplos test_values.sh
	./test_wraplos.sh

test_fltinv: test_fltinv.sh $(BIN_DIR)/fltinv test_values.sh
	./test_fltinv.sh
