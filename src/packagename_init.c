#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP compute_p_value_change(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_p_values(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP motif_score(SEXP, SEXP);
extern SEXP test_find_theta(SEXP, SEXP, SEXP, SEXP);
extern SEXP transition_matrix(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"compute_p_value_change", (DL_FUNC) &compute_p_value_change, 10},
    {"compute_p_values",       (DL_FUNC) &compute_p_values,        7},
    {"motif_score",            (DL_FUNC) &motif_score,             2},
    {"test_find_theta",        (DL_FUNC) &test_find_theta,         4},
    {"transition_matrix",      (DL_FUNC) &transition_matrix,       2},
    {NULL, NULL, 0}
};

void R_init_atSNP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}