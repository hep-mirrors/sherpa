# install all under Examples/, except CI and ".in" files that will be
# configured and installed below
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        DESTINATION ${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING} 
        COMPONENT doc 
        PATTERN CI EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE
        REGEX ".*\.in$" EXCLUDE)
