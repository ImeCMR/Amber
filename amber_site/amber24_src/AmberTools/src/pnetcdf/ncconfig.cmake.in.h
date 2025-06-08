/* src/lib/ncconfig.in.  Generated from configure.in by autoheader.  */

/* Define to one of `_getb67', `GETB67', `getb67' for Cray-2 and Cray-YMP
   systems. This function is required for `alloca.c' support on those systems.
   */
#cmakedefine CRAY_STACKSEG_END

/* Define to 1 if using `alloca.c'. */
#cmakedefine C_ALLOCA 1

/* Define if able to support nonblocking routines */
#cmakedefine ENABLE_NONBLOCKING 1

/* Define if Fortran names are lower case */
#cmakedefine F77_NAME_LOWER 1

/* Define if Fortran names are lower case with two trailing underscore2 */
#cmakedefine F77_NAME_LOWER_2USCORE 1

/* Define if Fortran names are lower case with one trailing underscore */
#cmakedefine F77_NAME_LOWER_USCORE 1

/* Define if Fortran names are uppercase */
#cmakedefine F77_NAME_UPPER 1

/* Define to 1 if you have `alloca', as a function or macro. */
#cmakedefine HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#cmakedefine HAVE_ALLOCA_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H 1

/* Define to 1 if you have the <mpio.h> header file. */
#cmakedefine HAVE_MPIO_H 1

/* available */
#cmakedefine HAVE_MPI_CHARACTER 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_DARRAY 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_DUP 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_F90_COMPLEX 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_F90_INTEGER 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_F90_REAL 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_HINDEXED_INTEGER 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_HVECTOR_INTEGER 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_INDEXED_BLOCK 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_RESIZED 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_STRUCT_INTEGER 1

/* available */
#cmakedefine HAVE_MPI_COMBINER_SUBARRAY 1

/* available */
#cmakedefine HAVE_MPI_COMPLEX16 1

/* available */
#cmakedefine HAVE_MPI_COMPLEX32 1

/* available */
#cmakedefine HAVE_MPI_COMPLEX8 1

/* available */
#cmakedefine HAVE_MPI_DOUBLE_PRECISION 1

/* Define to 1 if you have the `MPI_Info_dup' function. */
#cmakedefine HAVE_MPI_INFO_DUP 1

/* available */
#cmakedefine HAVE_MPI_INTEGER 1

/* available */
#cmakedefine HAVE_MPI_INTEGER1 1

/* available */
#cmakedefine HAVE_MPI_INTEGER16 1

/* available */
#cmakedefine HAVE_MPI_INTEGER2 1

/* available */
#cmakedefine HAVE_MPI_INTEGER4 1

/* available */
#cmakedefine HAVE_MPI_INTEGER8 1

/* available */
#cmakedefine HAVE_MPI_LB 1

/* available */
#cmakedefine HAVE_MPI_REAL 1

/* available */
#cmakedefine HAVE_MPI_REAL16 1

/* available */
#cmakedefine HAVE_MPI_REAL4 1

/* available */
#cmakedefine HAVE_MPI_REAL8 1

/* Define to 1 if you have the `MPI_Request_get_status' function. */
#cmakedefine HAVE_MPI_REQUEST_GET_STATUS 1

/* Define to 1 if you have the `MPI_Type_dup' function. */
#cmakedefine HAVE_MPI_TYPE_DUP 1

/* available */
#cmakedefine HAVE_MPI_UB 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#cmakedefine HAVE_PTRDIFF_T 1

/* Define to 1 if the system has the type `ssize_t'. */
#cmakedefine HAVE_SSIZE_T 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H 1

/* Define to 1 if `st_blksize' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_BLKSIZE 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uchar'. */
#cmakedefine HAVE_UCHAR 1

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H 1

/* Type of NC_BYTE */
#define NCBYTE_T ${NCBYTE_T}

/* Type of NC_SHORT */
#define NCSHORT_T ${NCBYTE_T}

/* Type for Fortran INT1 */
#define NF_INT1_T ${NF_INT1_T}

/* Type for Fortran INT2 */
#define NF_INT2_T ${NF_INT2_T}

/* all possibilities for netcdf-fortran types */
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

/* Does sytem have IEEE FLOAT */
#cmakedefine NO_IEEE_FLOAT 1

/* Define if system lacks strerror */
#cmakedefine NO_STRERROR 1

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ${PACKAGE_BUGREPORT}

/* Define to the full name of this package. */
#define PACKAGE_NAME ${PACKAGE_NAME}

/* Define to the full name and version of this package. */
#define PACKAGE_STRING ${PACKAGE_STRING}

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ${PACKAGE_NAME}

/* Define to the version of this package. */
#define PACKAGE_VERSION ${PACKAGE_VERSION}

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE ${SIZEOF_DOUBLE}

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT ${SIZEOF_FLOAT}

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT ${SIZEOF_INT}

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG ${SIZEOF_LONG}

/* The number of bytes in an MPI_Offset */
#define SIZEOF_MPI_OFFSET ${SIZEOF_MPI_OFFSET}

/* The size of `off_t', as computed by sizeof. */
#define SIZEOF_OFF_T ${SIZEOF_OFF_T}

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT ${SIZEOF_SHORT}

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T ${SIZEOF_SIZE_T}

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
#define STACK_DIRECTION 0


/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
#cmakedefine WORDS_BIGENDIAN 1

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
/* #undef YYTEXT_POINTER 1 */

/* Number of bits in a file offset, on hosts where this is settable. */
#cmakedefine _FILE_OFFSET_BITS 1

/* Define for large files, on AIX-style hosts. */
#cmakedefine _LARGE_FILES 1

/* Define to 1 if type `char' is unsigned and you are not using gcc.  */
#ifndef __CHAR_UNSIGNED__
# undef __CHAR_UNSIGNED__
#endif

/* Define to `long int' if <sys/types.h> does not define. */
#cmakedefine off_t 1

/* Define to `unsigned int' if <sys/types.h> does not define. */
#cmakedefine size_t 1
