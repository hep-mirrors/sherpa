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

file(COPY src/ufo_interface DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src/)
file(COPY src/opt_einsum/opt_einsum DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/src/)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/
        DESTINATION ${SHERPA_Python_SITEARCH}/)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/Sherpa-generate-model
        DESTINATION ${CMAKE_INSTALL_BINDIR})

if(SHERPA_ENABLE_TESTING)
add_custom_target(ufo_python_tests
        COMMAND ${Python_EXECUTABLE} -m unittest discover -v
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/Tests/UFO)
enable_testing()
add_test(NAME ufo_python_tests COMMAND ${CMAKE_COMMAND} --build . --target ufo_python_tests)
endif()

endif()
