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
set(PY_UFO_SOURCES
    src/ufo_interface/__init__.py
    src/ufo_interface/calculators.py
    src/ufo_interface/color_structures.py
    src/ufo_interface/contract.py
    src/ufo_interface/lorentz_algebra.py
    src/ufo_interface/lorentz_printer.py
    src/ufo_interface/lorentz_simplify.py
    src/ufo_interface/lorentz_wavefunction.py
    src/ufo_interface/lorentz_writer.py
    src/ufo_interface/message.py
    src/ufo_interface/model_writer.py
    src/ufo_interface/parameter.py
    src/ufo_interface/particle.py
    src/ufo_interface/py_to_cpp.py
    src/ufo_interface/s_color.py
    src/ufo_interface/tensor.py
    src/ufo_interface/ufo_exception.py
    src/ufo_interface/vertex.py
)

set(PY_UFO_TEMPLATES
    src/ufo_interface/Templates/CMakeLists.txt.in
    src/ufo_interface/Templates/color_calc_template.C
    src/ufo_interface/Templates/lorentz_calc_template.C
    src/ufo_interface/Templates/model_template.C
    src/ufo_interface/Templates/run_card_template
)

file(COPY ${PY_UFO_SOURCES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src/ufo_interface)
file(COPY ${PY_UFO_TEMPLATES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src/ufo_interface/Templates)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/ufo_interface/
        DESTINATION ${SHERPA_Python_SITEARCH}/ufo_interface)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/Sherpa-generate-model
        DESTINATION ${CMAKE_INSTALL_BINDIR})

endif()
