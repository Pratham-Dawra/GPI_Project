
#ifndef __ENUMS_H__
#define __ENUMS_H__

/*! Enum defining the type of wave equation. 
 *  @note The order defined here must match the order of "char* weq_descr[NWEQ]" 
 *        and "char* weq_verbose[NWEQ]" in file read_par_json.c.
 */
typedef enum WEQTYPE {
    AC_ISO = 0,
    VAC_ISO,
    EL_ISO,
    VEL_ISO,
    EL_VTI,
    VEL_VTI,
    EL_TTI,
    VEL_TTI,
    NWEQ
} WEQTYPE;

#endif
