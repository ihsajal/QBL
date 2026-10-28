#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void mcmc_BiCC(
    int *haplotype_map, int *tot_hap,
    double *yb, double *yc,
    int *N, int *num_haplo_id, int *x_length,
    double *freq, double *D,
    double *mu_b, double *beta_b,
    double *mu_c, double *beta_c,
    double *a, double *b, double *a0, double *b0,
    int *unique_map, int *tot_uniq_mat,
    double *lambda_b, double *lambda_c,
    double *u, double *sigma_u,
    double *sigma_sq_b, double *sigma_sq_c,
    int *NUM_IT, int *BURN_IN,
    double *beta_b_out, double *beta_c_out,
    double *lambda_b_out, double *lambda_c_out,
    double *freq_out, double *D_out,
    double *sigma_u_out,
    double *sigma_sq_b_out, double *sigma_sq_c_out
);

static const R_CMethodDef CEntries[] = {
  {"mcmc_BiCC", (DL_FUNC) &mcmc_BiCC, 36, NULL},
  {NULL, NULL, 0, NULL}
};

void R_init_QBL(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
