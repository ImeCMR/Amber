  # Do tests to determine which Fortran types correspond to NCBYTE, NCSHORT, ...
  # The configure file got these by testing an F77 program, invoking
  # UD_FORTRAN_TYPES defined in acinclude.m4. 
  
  #NOTE: this code is included by both pnetcdf AND netcdf-fortran.  Make sure to test both if you change it.
  
  # check sizes of C types
  check_type_size(int SIZEOF_INT)
  check_type_size(long SIZEOF_LONG)
  check_type_size(float SIZEOF_FLOAT)
  check_type_size(double SIZEOF_DOUBLE)
  check_type_size("signed char" SIZEOF_SIGNED_CHAR)  
  check_type_size(short SIZEOF_SHORT)
    
  # The automake script got all paranoid and checked many different Fortran types to see if they exist
  # I feel like we can assume that the standard types exist.
  SET(NCBYTE_T "byte") # 1 byte
  SET(NCSHORT_T "integer*2") # 2 bytes
  SET(NF_INT1_T "integer*1") # 1 byte
  SET(NF_INT2_T "integer*2") # 2 bytes
  
  # Checks the provided C types to see which ones have the number of bytes passed
  # Roughly equivalent to the Automake function UD_CHECK_CTYPE_FORTRAN
  # creates the result variable: NF_<RESULT_PREFIX>_IS_C_<INT,SHORT,LONG,etc.>.  This is the format accepted by config.h.in.
  macro(find_c_type_with_size SIZE RESULT_PREFIX) # 3rd arg: C types to check
  	set(FOUNDANY FALSE)
  	foreach(CTYPE ${ARGN})
	  # figure out type part of variable name
	  string(TOUPPER ${CTYPE} CTYPE_VAR_NAME)
	  string(REPLACE " " "_" CTYPE_VAR_NAME ${CTYPE_VAR_NAME})
	  
	  set(RESULT_VAR_NAME NF_${RESULT_PREFIX}_IS_C_${CTYPE_VAR_NAME})
	  
	  test(${RESULT_VAR_NAME} (NOT FOUNDANY) AND ${SIZEOF_${CTYPE_VAR_NAME}} EQUAL ${SIZE})
	  
	  if(${RESULT_VAR_NAME})
	  	set(FOUNDANY TRUE)
  	  endif()
	  
	  #message("${RESULT_VAR_NAME}: ${${RESULT_VAR_NAME}}")
	endforeach()
	
	if(NOT FOUNDANY)
		message(FATAL_ERROR "Unable to find a C ${RESULT_PREFIX} equivalent type with ${SIZE} bytes")
	endif()
	
	unset(FOUNDANY)
	
  endmacro(find_c_type_with_size)
  
  find_c_type_with_size(1 INT1 "signed char" short int long)
  find_c_type_with_size(2 INT2 short int long)
  find_c_type_with_size(4 INT int long)
  find_c_type_with_size(4 REAL float double)
  find_c_type_with_size(8 DOUBLEPRECISION float double)