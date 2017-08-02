
// Declarations for registration of native routines

// Following code is, except where otherwise indicated, the output from
//   tools::package_native_routine_registration_skeleton("E:/Rstuff/lmom","E:/R-reg-lmom.txt")

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cdqagie(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdqagse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(cdfwak)(void *, void *, void *, void *, void *);
extern void F77_NAME(lmrexp)(void *, void *, void *, void *);
extern void F77_NAME(lmrgam)(void *, void *, void *, void *);
extern void F77_NAME(lmrgev)(void *, void *, void *, void *);
extern void F77_NAME(lmrglo)(void *, void *, void *, void *);
extern void F77_NAME(lmrgno)(void *, void *, void *, void *);
extern void F77_NAME(lmrgpa)(void *, void *, void *, void *);
extern void F77_NAME(lmrgum)(void *, void *, void *, void *);
extern void F77_NAME(lmrkap)(void *, void *, void *, void *);
extern void F77_NAME(lmrnor)(void *, void *, void *, void *);
extern void F77_NAME(lmrpe3)(void *, void *, void *, void *);
extern void F77_NAME(lmrwak)(void *, void *, void *, void *);
extern void F77_NAME(pelexp)(void *, void *, void *);
extern void F77_NAME(pelgam)(void *, void *, void *);
extern void F77_NAME(pelgev)(void *, void *, void *);
extern void F77_NAME(pelglo)(void *, void *, void *);
extern void F77_NAME(pelgno)(void *, void *, void *);
extern void F77_NAME(pelgpa)(void *, void *, void *);
extern void F77_NAME(pelgum)(void *, void *, void *);
extern void F77_NAME(pelkap)(void *, void *, void *);
extern void F77_NAME(pelnor)(void *, void *, void *);
extern void F77_NAME(pelpe3)(void *, void *, void *);
extern void F77_NAME(pelwa0)(void *, void *, void *);
extern void F77_NAME(pelwak)(void *, void *, void *);
extern void F77_NAME(samlm)(void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"cdqagie", (DL_FUNC) &cdqagie, 11},
    {"cdqagse", (DL_FUNC) &cdqagse, 11},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"cdfwak", (DL_FUNC) &F77_NAME(cdfwak), 5},
    {"lmrexp", (DL_FUNC) &F77_NAME(lmrexp), 4},
    {"lmrgam", (DL_FUNC) &F77_NAME(lmrgam), 4},
    {"lmrgev", (DL_FUNC) &F77_NAME(lmrgev), 4},
    {"lmrglo", (DL_FUNC) &F77_NAME(lmrglo), 4},
    {"lmrgno", (DL_FUNC) &F77_NAME(lmrgno), 4},
    {"lmrgpa", (DL_FUNC) &F77_NAME(lmrgpa), 4},
    {"lmrgum", (DL_FUNC) &F77_NAME(lmrgum), 4},
    {"lmrkap", (DL_FUNC) &F77_NAME(lmrkap), 4},
    {"lmrnor", (DL_FUNC) &F77_NAME(lmrnor), 4},
    {"lmrpe3", (DL_FUNC) &F77_NAME(lmrpe3), 4},
    {"lmrwak", (DL_FUNC) &F77_NAME(lmrwak), 4},
    {"pelexp", (DL_FUNC) &F77_NAME(pelexp), 3},
    {"pelgam", (DL_FUNC) &F77_NAME(pelgam), 3},
    {"pelgev", (DL_FUNC) &F77_NAME(pelgev), 3},
    {"pelglo", (DL_FUNC) &F77_NAME(pelglo), 3},
    {"pelgno", (DL_FUNC) &F77_NAME(pelgno), 3},
    {"pelgpa", (DL_FUNC) &F77_NAME(pelgpa), 3},
    {"pelgum", (DL_FUNC) &F77_NAME(pelgum), 3},
    {"pelkap", (DL_FUNC) &F77_NAME(pelkap), 3},
    {"pelnor", (DL_FUNC) &F77_NAME(pelnor), 3},
    {"pelpe3", (DL_FUNC) &F77_NAME(pelpe3), 3},
    {"pelwa0", (DL_FUNC) &F77_NAME(pelwa0), 3},
    {"pelwak", (DL_FUNC) &F77_NAME(pelwak), 3},
    {"samlm",  (DL_FUNC) &F77_NAME(samlm),  6},
    {NULL, NULL, 0}
};

void R_init_lmom(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);           // added per R-exts manual sec 5.4.2
}

