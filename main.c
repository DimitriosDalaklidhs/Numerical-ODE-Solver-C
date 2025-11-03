#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "f_functionh.c"

#ifdef _WIN32
#include <windows.h>
#endif

/* -------- Utility: round((b - a)/h + 1e-12) -------- */
static long num_steps(double a, double b, double h) {
    double x = (b - a) / h + 1e-12;
    return (long)floor(x + 0.5);
}

/* Fill t[i] = a + i*h, i=0..N */
static void fill_t(double *t, long N, double a, double h) {
    long i;
    for (i = 0; i <= N; ++i)
        t[i] = a + i * h;
}

/* ===============================
   Forward Euler (FE)
   =============================== */
static void forward_euler(double (*ff)(double, double),
                          double a, double b, double y0, double h,
                          double **t_out, double **y_out, long *len_out)
{
    long N = num_steps(a, b, h);
    double *t = (double*)malloc((N + 1) * sizeof(double));
    double *y = (double*)calloc((N + 1), sizeof(double));
    long n;

    if (!t || !y) { perror("malloc"); exit(EXIT_FAILURE); }

    fill_t(t, N, a, h);
    y[0] = y0;
    for (n = 0; n < N; ++n)
        y[n + 1] = y[n] + h * ff(t[n], y[n]);

    *t_out = t; *y_out = y; *len_out = N + 1;
}

/* ===============================
   Central Difference (2-step)
   =============================== */
static void central_difference(double (*ff)(double, double),
                               double a, double b, double y0, double h,
                               double **t_out, double **y_out, long *len_out)
{
    long N = num_steps(a, b, h);
    double *t = (double*)malloc((N + 1) * sizeof(double));
    double *y = (double*)calloc((N + 1), sizeof(double));
    long n;

    if (!t || !y) { perror("malloc"); exit(EXIT_FAILURE); }

    fill_t(t, N, a, h);
    y[0] = y0;
    if (N >= 1)
        y[1] = y[0] + h * ff(t[0], y[0]); /* bootstrap */
    for (n = 1; n < N; ++n)
        y[n + 1] = y[n - 1] + 2.0 * h * ff(t[n], y[n]);

    *t_out = t; *y_out = y; *len_out = N + 1;
}

/* ===============================
   Improved Euler (Heun)
   =============================== */
static void improved_euler(double (*ff)(double, double),
                           double a, double b, double y0, double h,
                           double **t_out, double **y_out, long *len_out)
{
    long N = num_steps(a, b, h);
    double *t = (double*)malloc((N + 1) * sizeof(double));
    double *y = (double*)calloc((N + 1), sizeof(double));
    long n;

    if (!t || !y) { perror("malloc"); exit(EXIT_FAILURE); }

    fill_t(t, N, a, h);
    y[0] = y0;
    for (n = 0; n < N; ++n) {
        double k1 = ff(t[n], y[n]);
        double k2 = ff(t[n] + h, y[n] + h * k1);
        y[n + 1] = y[n] + 0.5 * h * (k1 + k2);
    }

    *t_out = t; *y_out = y; *len_out = N + 1;
}

/* ---------------------- Result struct ---------------------- */
typedef struct {
    double h, err_fe, err_cd, err_ie;
} Result;

/* ===============================
   Run once for a given h
   =============================== */
static Result run_once(double a, double b, double y0, double h, const char *title_prefix)
{
    double *t_fe = NULL, *y_fe = NULL; long L_fe = 0;
    double *t_cd = NULL, *y_cd = NULL; long L_cd = 0;
    double *t_ie = NULL, *y_ie = NULL; long L_ie = 0;
    long i;
    double *t, *y_ex;
    double err_fe = 0.0, err_cd = 0.0, err_ie = 0.0;

    forward_euler(f, a, b, y0, h, &t_fe, &y_fe, &L_fe);
    central_difference(f, a, b, y0, h, &t_cd, &y_cd, &L_cd);
    improved_euler(f, a, b, y0, h, &t_ie, &y_ie, &L_ie);

    if (!(L_fe == L_cd && L_fe == L_ie)) {
        fprintf(stderr, "Error: grid mismatch.\n");
        exit(EXIT_FAILURE);
    }

    t = t_fe;
    y_ex = (double*)malloc(L_fe * sizeof(double));
    if (!y_ex) { perror("malloc"); exit(EXIT_FAILURE); }

    for (i = 0; i < L_fe; ++i)
        y_ex[i] = y_exact(t[i]);

    for (i = 0; i < L_fe; ++i) {
        double dfe = fabs(y_ex[i] - y_fe[i]);
        double dcd = fabs(y_ex[i] - y_cd[i]);
        double die = fabs(y_ex[i] - y_ie[i]);
        if (dfe > err_fe) err_fe = dfe;
        if (dcd > err_cd) err_cd = dcd;
        if (die > err_ie) err_ie = die;
    }

    printf("\n=== %s h = %.3f ===\n", title_prefix ? title_prefix : "", h);
    printf("i   t_i        y_FD         y_CD         y_IE         y_exact     |err_FD|    |err_CD|    |err_IE|\n");
    for (i = 0; i < L_fe; ++i) {
        printf("%-3ld%8.4f  %12.8f  %12.8f  %12.8f  %12.8f  %10.3e  %10.3e  %10.3e\n",
               i, t[i], y_fe[i], y_cd[i], y_ie[i], y_ex[i],
               fabs(y_ex[i]-y_fe[i]), fabs(y_ex[i]-y_cd[i]), fabs(y_ex[i]-y_ie[i]));
    }

    printf("\nMax error FE: %.6e\n", err_fe);
    printf("Max error CD: %.6e\n", err_cd);
    printf("Max error IE: %.6e\n", err_ie);

    { Result r;
      r.h = h; r.err_fe = err_fe; r.err_cd = err_cd; r.err_ie = err_ie;
      free(t_fe); free(y_fe); free(t_cd); free(y_cd); free(t_ie); free(y_ie); free(y_ex);
      return r;
    }
}


int main(void)
{
double a, b, y0, h;
char mode[8];

printf("Enter initial value problem data y'(t)=f(t,y), y(a)=y0\n");
printf("a: "); if (scanf("%lf", &a) != 1) return 1;
printf("b: "); if (scanf("%lf", &b) != 1) return 1;
printf("y0: "); if (scanf("%lf", &y0) != 1) return 1;

printf("Do you want to run for h=0.2 and h=0.1? [y/n]: ");
if (scanf("%7s", mode) != 1) return 1;

if (mode[0] == 'y' || mode[0] == 'Y') {
    Result r1 = run_once(a, b, y0, 0.2, "Result");
    Result r2 = run_once(a, b, y0, 0.1, "Result");

    printf("\nERROR SUMMARY (max |y(t_i)-y_i|):\n");
    printf("h       err_FE        err_CD        err_IE\n");
    printf("%-6.3f  %-12.6e %-12.6e %-12.6e\n", r1.h, r1.err_fe, r1.err_cd, r1.err_ie);
    printf("%-6.3f  %-12.6e %-12.6e %-12.6e\n", r2.h, r2.err_fe, r2.err_cd, r2.err_ie);
} else {
    printf("Enter step size h: ");
    if (scanf("%lf", &h) != 1) return 1;
    (void)run_once(a, b, y0, h, "Result");
}





      

    return 0;
}

