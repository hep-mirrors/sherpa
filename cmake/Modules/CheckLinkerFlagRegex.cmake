# CheckLinkerFlagRegex
# ------------------
#
# Checks whether a compiler supports a given linker flag.
#

include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckFortranSourceCompiles)

# FLAG: the linker flag to check
# LANGUAGE: the language to test it on.  C, CXX, or Fortran
# REGEX: the regular expression for linker warnings/errors to match against
# RESULT: return variable, 1 for success

macro(check_linker_flag_regex FLAG LANGUAGE REGEX RESULT)
	if(NOT DEFINED ${RESULT})
		set(OLD_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
		set(CMAKE_REQUIRED_FLAGS "${FLAG}")
		# Normalize locale during test compilation.
		set(_CheckLinkerFlag_LOCALE_VARS LC_ALL LC_MESSAGES LANG)
		foreach(v ${_CheckLinkerFlag_LOCALE_VARS})
			set(_CheckLinkerFlag_SAVED_${v} "$ENV{${v}}")
			set(ENV{${v}} C)
		endforeach()
	   
		if("${LANGUAGE}" STREQUAL C)
			check_c_source_compiles("int main(void) { return 0; }" ${RESULT}
			FAIL_REGEX "${REGEX}"
			)
		elseif("${LANGUAGE}" STREQUAL CXX)
			check_cxx_source_compiles("int main(void) { return 0; }" ${RESULT}
			FAIL_REGEX "${REGEX}"
			)
		elseif("${LANGUAGE}" STREQUAL Fortran)
			check_fortran_source_compiles(
			"program test
			print *, \'Hello, world\'
			end program test" 
			${RESULT}
			FAIL_REGEX "${REGEX}"
			)
		else()
			message(FATAL_ERROR "Invalid LANGUAGE argument ${LANGUAGE}")
		endif()
	    
		foreach(v ${_CheckCCompilerFlag_LOCALE_VARS})
	   		set(ENV{${v}} ${_CheckCCompilerFlag_SAVED_${v}})
	     		unset(_CheckCCompilerFlag_SAVED_${v})
	   	endforeach()
		unset(_CheckCCompilerFlag_LOCALE_VARS)
	 	unset(_CheckCCompilerFlag_COMMON_PATTERNS)
	
		set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
	endif()
endmacro (check_linker_flag_regex)
