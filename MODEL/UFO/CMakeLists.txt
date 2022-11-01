
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#UFO/Makefile.am
 
 set(GITTAG  MODEL/UFO)
 #include ../../git.make)
 #pkglib_LTLIBRARIES = libModelUFO.la)
set(UFO_SOURCES  UFO_Param_Reader.C UFO_Model.C UFO_Color_Functions.C)
set(UFO_HEADERS UFO_Param_Reader.H UFO_Model.H UFO_Color_Functions.H)
 if (SHERAP_ENABLE_UFO)
 #dist_bin_SCRIPTS = Sherpa-generate-model)
 #ufodir = $(pythondir)/ufo_interface)
#dist_ufo_PYTHON = __init__.py sym_var.py lorentz_calc_template.C color_calc_template.C lorentz_structures.py color_structures.py code_snippets.py model_template.C run_card_template py_to_cpp.py s_coupling.py s_lorentz.py s_color.py s_parameter.py s_particle.py s_vertex.py sconstruct_template templates.py tensor.py test.py ufo_exception.py write_model.py write_run_card.py message.py
 endif()
creategitinfo(UFO ../UFO/)
list(TRANSFORM UFO_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../UFO/")
list(APPEND UFO_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(UFO_esources )
if (SHERPA_BUILD_SHARED)
add_library(ModelUFO SHARED ${ModelUFO_esources}
                             ${UFO_SOURCES} 
)
target_include_directories(ModelUFO PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ModelUFO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ModelUFO DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ModelUFO        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ModelUFO SOVERSION 0.0.0)
set_target_properties(ModelUFO        PROPERTIES DEFINE_SYMBOL "")
endif()
