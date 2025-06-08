#ifndef NFE_MACROS_H
#define NFE_MACROS_H

#ifdef __GNUC__
#  define NFE_GNUC_NORETURN    __attribute__ ((__noreturn__))
#  define NFE_GNUC_UNUSED      __attribute__ ((__unused__))
#  define NFE_GNUC_CONST       __attribute__ ((__const__))
#  define NFE_GNUC_PRINTF(a,b) __attribute__ ((format(printf,a,b)))
#  define NFE_RESTRICT         __restrict__
#else
#  define NFE_GNUC_NORETURN
#  define NFE_GNUC_UNUSED
#  define NFE_GNUC_CONST
#  define NFE_GNUC_PRINTF(a,b)
#  define NFE_RESTRICT
#endif // __GNUC__

#endif // NFE_MACROS_H
