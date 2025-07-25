#if 0
    $Id: nfconfig.in,v 1.1 2009/08/24 23:58:47 rcw Exp $ 
#endif


#ifndef UD_NETCDF_CPP_INC
#define UD_NETCDF_CPP_INC


#if 0
    Do not have C-style comments in here because this file is processed
    by both the FORTRAN compiler (for the nf_test/ stuff) and the C
    compiler (for the FORTRAN-callable interface routines) and some 
    FORTRAN preprocessors do not understand the /*...*/ syntax.
#endif


#if 0
    The following macros define the supplementary FORTRAN arithmetic
    datatypes beyond the standard INTEGER, REAL, and DOUBLEPRECISION --
    ostensibly corresponding to 8-bit and 16-bit integers, respectively.
    For example:

	#define NF_INT1_T        byte
	#define NF_INT2_T        integer*2

    These are the types of the relevant arguments in the NF_*_INT1() and
    NF_*_INT2() netCDF FORTRAN function calls.  The word "ostensibly"
    is used advisedly: on some systems an "integer*2" datatype,
    nevertheless, occupies 64 bits (we are not making this up).

    If your FORTRAN system does not have the respective supplementary
    datatype, then do not define the corresponding macro.
#endif
#define NF_INT1_T ${NF_INT1_T}
#define NF_INT2_T ${NF_INT1_T}


#if 0
   Define the following NF_*_IS_C_* macros appropriatly for your system.
   The "INT1", "INT2" and "INT" after the "NF_" refer to the NF_INT1_T
   FORTRAN datatype, the NF_INT2_T FORTRAN datatype, and the INTEGER
   FORTRAN datatype, respectively.  If the respective FORTRAN datatype
   does not exist, then do not define the corresponding macro.
#endif
#cmakedefine NF_INT1_IS_C_SIGNED_CHAR 1
#cmakedefine NF_INT1_IS_C_SHORT 1
#cmakedefine NF_INT1_IS_C_INT 1
#cmakedefine NF_INT1_IS_C_LONG 1
#cmakedefine NF_INT2_IS_C_SHORT 1
#cmakedefine NF_INT2_IS_C_INT 1
#cmakedefine NF_INT2_IS_C_LONG 1
#cmakedefine NF_INT_IS_C_INT 1
#cmakedefine NF_INT_IS_C_LONG 1
#cmakedefine NF_REAL_IS_C_FLOAT 1
#cmakedefine NF_REAL_IS_C_DOUBLE 1
#cmakedefine NF_DOUBLEPRECISION_IS_C_DOUBLE 1
#cmakedefine NF_DOUBLEPRECISION_IS_C_FLOAT 1


#if 0
   Whether the system uses something besides the IEEE floating-point
   format to represent floating-point values.
#endif
#cmakedefine NO_IEEE_FLOAT 1


#if 0
   END OF CUSTOMIZATION
#endif


#if 0
    FORTRAN data types corresponding to netCDF version 2 "byte" and "short"
    data types (e.g. INTEGER*1, INTEGER*2).  See file "ftest.F" for usage.
#endif
#if !defined(NO_NETCDF_2)
#   define NCBYTE_T ${NCBYTE_T}
#   define NCSHORT_T ${NCBYTE_T}
#endif


#endif
