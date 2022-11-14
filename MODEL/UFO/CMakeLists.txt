##### Automatically or semiautomaticaly generated, do not edit.#########
#UFO/Makefile.am
 
 set(GITTAG  MODEL/UFO)
set(UFO_SOURCES  UFO_Param_Reader.C UFO_Model.C UFO_Color_Functions.C)
set(UFO_HEADERS UFO_Param_Reader.H UFO_Model.H UFO_Color_Functions.H)
 if (SHERPA_ENABLE_UFO)
 #dist_bin_SCRIPTS = Sherpa-generate-model)
 #ufodir = $(pythondir)/ufo_interface)
#dist_ufo_PYTHON = __init__.py sym_var.py lorentz_calc_template.C color_calc_template.C lorentz_structures.py color_structures.py code_snippets.py model_template.C run_card_template py_to_cpp.py s_coupling.py s_lorentz.py s_color.py s_parameter.py s_particle.py s_vertex.py sconstruct_template templates.py tensor.py test.py ufo_exception.py write_model.py write_run_card.py message.py
 endif()
sherpa_create_git_info(UFO ../UFO/)
list(TRANSFORM UFO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../UFO/")
list(APPEND UFO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ModelUFO SHARED ${UFO_SOURCES} )
sherpa_handle_shared_library(ModelUFO ModelUFO)
