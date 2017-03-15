SVNTAG = MODEL/UFO
include ../../svn.make

pkglib_LTLIBRARIES = libModelUFO.la

libModelUFO_la_SOURCES = \
	SVN_Info.C \
	UFO_Param_Reader.C \
	UFO_Model.C \
	UFO_Color_Functions.C

localinc_HEADERS = \
	UFO_Param_Reader.H \
	UFO_Model.H \
	UFO_Color_Functions.H

if ENABLE_UFO

dist_bin_SCRIPTS = Sherpa-generate-model

ufodir = $(pythondir)/ufo_interface

dist_ufo_PYTHON = \
        __init__.py \
        c_variable.py \
        color_interface.py \
        lorentz_calc_template.C \
        lorentz_structures.py \
        model_template.C \
	run_card_template \
        py_to_cpp.py \
        s_coupling.py \
        s_lorentz.py \
        s_parameter.py \
        s_particle.py \
        s_vertex.py \
        sconstruct_template \
        templates.py \
        tensor.py \
        test.py \
        ufo_exception.py \
        write_model.py \
	write_run_card.py \
	message.py

endif
