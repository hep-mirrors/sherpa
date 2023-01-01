find_program(SPHINXBUILD NAMES sphinx-build  REQUIRED)
message(STATUS "SHERPA Manual: SPHINXBUILD=${SPHINXBUILD}")
find_package(LATEX COMPONENTS PDFLATEX MAKEINDEX REQUIRED)
message(STATUS "SHERPA Manual: LATEX_FOUND=${LATEX_FOUND} PDFLATEX_COMPILER=${PDFLATEX_COMPILER} MAKEINDEX_COMPILER=${MAKEINDEX_COMPILER}")
find_program(LATEXMK NAMES latexmk  REQUIRED)
message(STATUS "SHERPA Manual: LATEXMK=${LATEXMK}")

set (sphinx_options  "-D release=${SHERPA_VERSION}")
add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/completion.index  ${CMAKE_CURRENT_BINARY_DIR}/options.index
       COMMAND ${SPHINXBUILD} -M completion  ${CMAKE_CURRENT_SOURCE_DIR}/source . ${sphinx_options}
       COMMAND mv completion/completion.index  .
       COMMAND mv completion/options.index .
       COMMAND rm -r completion
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/Sherpa.1 
       COMMAND ${SPHINXBUILD} -M man  ${CMAKE_CURRENT_SOURCE_DIR}/source . ${sphinx_options}
       COMMAND mv man/Sherpa.1 .
       COMMAND rm -r man
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)
       
add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/sherpamanual_html
       COMMAND ${SPHINXBUILD} -M html  ${CMAKE_CURRENT_SOURCE_DIR}/source . ${sphinx_options}
       COMMAND rm -rf sherpamanual_html
       COMMAND mv html sherpamanual_html
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)

add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.pdf
       COMMAND ${SPHINXBUILD} -M latexpdf  ${CMAKE_CURRENT_SOURCE_DIR}/source  . ${sphinx_options}
       COMMAND mv latex/sherpamanual.pdf .
       COMMAND mv latex/sherpamanual.fls .
       COMMAND rm -r latex
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)

add_custom_command(
       OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.info
       COMMAND ${SPHINXBUILD} -M info  ${CMAKE_CURRENT_SOURCE_DIR}/source  . ${sphinx_options}
       COMMAND mv texinfo/sherpamanual.info .
       COMMAND rm -r texinfo
       WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
       VERBATIM)

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/Sherpa.1  DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.info DESTINATION ${CMAKE_INSTALL_INFODIR})
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR} )

ADD_CUSTOM_TARGET(all-files 
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/completion.index  
${CMAKE_CURRENT_BINARY_DIR}/options.index
${CMAKE_CURRENT_BINARY_DIR}/Sherpa.1 
${CMAKE_CURRENT_BINARY_DIR}/sherpamanual_html
${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.pdf
${CMAKE_CURRENT_BINARY_DIR}/sherpamanual.info 
)
ADD_CUSTOM_TARGET(all-formats ALL)
ADD_DEPENDENCIES(all-formats all-files)
set (DATADIR ${CMAKE_INSTALL_DATADIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sherpa-completion.in ${CMAKE_CURRENT_BINARY_DIR}/sherpa-completion @ONLY)

install( FILES
${CMAKE_CURRENT_BINARY_DIR}/completion.index  
${CMAKE_CURRENT_BINARY_DIR}/options.index
${CMAKE_CURRENT_BINARY_DIR}/sherpa-completion
DESTINATION
${CMAKE_INSTALL_DATADIR}/SHERPA-MC
)
