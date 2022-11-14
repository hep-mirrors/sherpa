##### Automatically or semiautomaticaly generated, do not edit.#########
#Zfunctions/Makefile.am
 
 #pkglib_LTLIBRARIES = libZfunctions.la)
 set(GITTAG  AMEGIC++/Amplitude/Zfunctions)
 #include ../../../git.make)
set(Zfunctions_SOURCES  Basic_Func.C Basic_Sfuncs.C Basic_Pfunc.C Basic_Mfunc.C Basic_MassTermfunc.C Basic_Vfunc.C Basic_Xfunc.C Basic_Yfunc.C Basic_Zfunc.C Basic_Epsilonfunc.C Unitarityfunc.C G4_Calc.C
 V4_Calc.C VVV_Calc.C Loops_Calc.C Zfunc_Calc.C Calculator.C Mom.C MHVCalculator.C)
set(Zfunctions_HEADERS Basic_Func.H Basic_Sfuncs.H Basic_Func.icc Zfunc_Calc.H Calculator.H Mom.H MHVCalculator.H)
sherpa_create_git_info(Zfunctions ../Zfunctions/)
list(TRANSFORM Zfunctions_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Zfunctions/")
list(APPEND Zfunctions_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(Zfunctions SHARED ${Zfunctions_SOURCES} )
sherpa_handle_shared_library(Zfunctions Zfunctions)
