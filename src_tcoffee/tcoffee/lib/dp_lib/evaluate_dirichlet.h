typedef struct 
{
  int N_COMPONENT;
  double *double_logB_alpha;
  double *exponant_list;
  double **ALPHA;
  double *DM_Q;
  double *alpha_tot;
  int n_aa;
  int tot_n;
}
Mixture;

double int_logB (int *i, int n);
double float_logB (float *i, int n);
double double_logB (double *x, int n);
double *** make_lup_table ( Mixture *D);
double  double_logB2(int j, double *n,Mixture *D);
double compute_exponant ( double *n, int j, Mixture *D);

double *compute_matrix_p ( double *n);
double *compute_dirichlet_p ( double *n);
void precompute_log_B ( double *table,Mixture *D);
double compute_X (double *n,int i,Mixture *D);
Mixture *read_dirichlet ( char *name);
int dirichlet_code( char aa);
int *aa2dirichlet_code_lu ();
int *dirichlet_code2aa_lu ();

double lgamma2 ( double x);
double lgamma_r(double x, int *signgamp);

double **aln2prf (Alignment *A, int ns, int *ls, int len, double **prf);
double **prf2dmx (double **in, double **prf, int len);
