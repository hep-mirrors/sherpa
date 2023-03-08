# install all under Examples/, except CI and ".in" files that will be
# configured and installed below
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        DESTINATION ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} 
        COMPONENT doc 
        PATTERN CI EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE
        REGEX ".*\.in$" EXCLUDE
        REGEX "\..*\.sw.$" EXCLUDE)

# configures and installs an ".in" file; the argument PATH must specify the
# file's path, without the ".in" extension
function(CONFIGURE_AND_INSTALL PATH)
  configure_file(${PATH}.in ${CMAKE_CURRENT_BINARY_DIR}/${PATH} @ONLY)
  cmake_path(GET PATH PARENT_PATH target_dir)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PATH}
          DESTINATION ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING}/Examples/${target_dir})
endfunction()

# configure and install all ".in" files (insofar it makes sense to do so)
set (IN_FILE_PATHS
     Userhook/CMakeLists.txt
     V_plus_Bs/LHC_Wbb/CMakeLists.txt
     V_plus_Jets/LHC_WJets/CMakeLists.txt)
if (SHERPA_ENABLE_PYTHON)
  set (IN_FILE_PATHS ${IN_FILE_PATHS}
       API/ME2-Python/test.py
       API/MPIEvents/test.py
       VBF_Interference/BH_Interface/BlackHat.py
       VBF_Interference/BH_Interface/Comix.py
       VBF_Interference/Tests/Comix.py)
endif()
foreach (PATH ${IN_FILE_PATHS})
  configure_and_install(${PATH})
endforeach()
