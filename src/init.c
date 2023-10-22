#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _rolleigen_roll_eigen(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _rolleigen_roll_pcr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
  {"_rolleigen_roll_eigen", (DL_FUNC) &_rolleigen_roll_eigen,  9},
  {"_rolleigen_roll_pcr",   (DL_FUNC) &_rolleigen_roll_pcr,   12},
  {NULL, NULL, 0}
};

void R_init_rolleigen(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}