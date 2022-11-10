
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Math/Makefile.am
 
 set(GITTAG  ATOOLS/Math)
 #include ../../git.make)
 #pkglib_LTLIBRARIES = libToolsMath.la)
set(Math_SOURCES  MathTools.C Term.C Algebra_Interpreter.C Marsaglia.C Sobol.C Random.C Function_Base.C Poincare.C ZAlign.C Vector.C Matrix.C Gauss_Integrator.C Scaling.C Axis.C Histogram.C Histogram_2D.C
 Variable.C Permutation.C PermStore.C Kabbala.C Tensor.C)
set(Math_HEADERS MathTools.H Term.H Algebra_Interpreter.H Marsaglia.H Sobol.H Random.H Function_Base.H Poincare.H ZAlign.H Vec3.H Vec4.H Vector.H MyComplex.H Matrix.H Gauss_Integrator.H Variable.H Scaling.H
 Axis.H Histogram.H Histogram_2D.H Permutation.H PermStore.H Cluster_Algorithm.H Cluster_Algorithm.C Kabbala.H Lorentz_Ten2.H Lorentz_Ten3.H Lorentz_Ten4.H Tensor.H Tensor_Build.H Tensor_Contractions.H
 Tensor_Contractions_Epsilon.H Uncertain.H)
 if (SHERPA_ENABLE_GZIP)
 #Sobol/%.gz: Sobol/%)
 #mkdir -p Sobol)
 #gzip < $< > $@)
 set(SOBOL_EXT   .gz)
 else()
 set(SOBOL_EXT )
 endif()
#nobase_#dist_pkgdata_DATA = Sobol/21201$(SOBOL_EXT)
creategitinfo(Math ../Math/)
list(TRANSFORM Math_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Math/")
list(APPEND Math_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ToolsMath SHARED 
                             ${Math_SOURCES} 
)
target_include_directories(ToolsMath PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ToolsMath)
install(TARGETS ToolsMath DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ToolsMath        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ToolsMath SOVERSION 0.0.0)
set_target_properties(ToolsMath        PROPERTIES DEFINE_SYMBOL "")
