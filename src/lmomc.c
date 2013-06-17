
#include <R.h>
#include <Rinternals.h>

//-------------------------------------------------------------------------------
//  C declaration of Fortran routines to be called from C

//         subroutine dqagie(bound,   inf,  epsabs,  epsrel,  limit,result,  abserr,
//                           neval,   ier,  last,    IENV)
extern void F77_NAME(dqagie)(double*, int*, double*, double*, int*, double*, double*,
                             int*,    int*, int*,    SEXP);

//         subroutine dqagse(a,       b,       epsabs,  epsrel,  limit,result,  abserr,
//                           neval,   ier,     last,    IENV)
extern void F77_NAME(dqagse)(double*, double*, double*, double*, int*, double*, double*,
                             int*,    int*,    int*,    SEXP);

//--------------------------------------------------------------------------------
//  C functions to be called from R

SEXP cdqagie(SEXP bound, SEXP inf, SEXP epsabs, SEXP epsrel, SEXP limit,
            SEXP result, SEXP abserr, SEXP neval, SEXP ier, SEXP last, SEXP env) {

  F77_CALL(dqagie)(REAL(bound), INTEGER(inf), REAL(epsabs), REAL(epsrel), INTEGER(limit),
                   REAL(result), REAL(abserr), INTEGER(neval), INTEGER(ier), INTEGER(last),
                   env);

  return R_NilValue;
}

SEXP cdqagse(SEXP a, SEXP b, SEXP epsabs, SEXP epsrel, SEXP limit,
            SEXP result, SEXP abserr, SEXP neval, SEXP ier, SEXP last, SEXP env) {

  F77_CALL(dqagse)(REAL(a), REAL(b), REAL(epsabs), REAL(epsrel), INTEGER(limit),
                   REAL(result), REAL(abserr), INTEGER(neval), INTEGER(ier), INTEGER(last),
                   env);

  return R_NilValue;
}
//-------------------------------------------------------------------------------
//  C function to be called from Fortran

void F77_SUB(f)(double *fout, double *fin, int *n, SEXP env) {
  int i;
  SEXP rin, rout;

  PROTECT(rin=allocVector(REALSXP,*n));
  for (i=0; i<*n; i++) REAL(rin)[i]=fin[i];

  defineVar(install("x"),rin,env);

  PROTECT(rout = eval(findVarInFrame(env, install("expr")), env)) ;

  if (length(rout)!=*n) error("evaluation of integrand gave result of wrong length");
  rout=coerceVector(rout,REALSXP);
  for (i=0; i<*n; i++) {
    fout[i]=REAL(rout)[i];
    if (!R_FINITE(fout[i])) error("non-finite integrand at argument %f",fin[i]) ;
  }

  UNPROTECT(2);
  return;
}

