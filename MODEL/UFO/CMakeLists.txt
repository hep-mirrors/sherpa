set(UFO_SOURCES UFO_Param_Reader.C
                UFO_Model.C
                UFO_Color_Functions.C)

set(UFO_HEADERS UFO_Param_Reader.H
                UFO_Model.H
                UFO_Color_Functions.H)

sherpa_create_git_info(UFO)
 
add_library(ModelUFO SHARED ${UFO_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ModelUFO ModelUFO)

if (SHERPA_ENABLE_UFO)
set(UFO_SOURCES
    src/ufo_interface/__init__.py
    src/ufo_interface/calculators.py
    src/ufo_interface/color_structures.py
    src/ufo_interface/lorentz_algebra.py
    src/ufo_interface/lorentz_printer.py
    src/ufo_interface/lorentz_simplify.py
    src/ufo_interface/lorentz_wavefunction.py
    src/ufo_interface/lorentz_writer.py
    src/ufo_interface/message.py
    src/ufo_interface/model_writer.py
    src/ufo_interface/parameter.py
    src/ufo_interface/parser.py
    src/ufo_interface/particle.py
    src/ufo_interface/py_to_cpp.py
    src/ufo_interface/s_color.py
    src/ufo_interface/tensor.py
    src/ufo_interface/ufo_exception.py
    src/ufo_interface/vertex.py
)

set(_pip_args "--prefix=${CMAKE_BINARY_DIR}")
add_custom_command(OUTPUT ${SHERPA_Python_SITEARCH}/ufo_interface/__init__.py
    COMMAND ${Python_EXECUTABLE} -m pip -q install ${_pip_args} ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${UFO_SOURCES}
    COMMENT "Installing UFO interface")
add_custom_target(ufo_interface ALL DEPENDS ${SHERPA_Python_SITEARCH}/ufo_interface/__init__.py)
add_dependencies(ModelUFO ufo_interface)

endif()
