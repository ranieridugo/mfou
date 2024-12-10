#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <omp.h>
#include <string.h>
// GNU Scientific Library components
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

#ifndef PI
#define PI	3.14159265358979323846264338327950288
#endif

double * wjk(int[], double, double, double, int);
double * gammauv(int[], double, double, double, int);
double * cuvj(int, double[2], double[2][2], double[2][2], int, int);
int * seq_by1(int, int);
double fVarFOU(double, double, double);
double fCov2FOU(double, double, double, double, double, double, double, double);
char* name = "2FOU";
char* extension = ".csv";

// Simulation parameters (path)
int iNPaths = 2; // Number of paths for each core
int iSteps = 5; // Number of substeps in each interval
char* dir = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/Desktop/tmp/"; // Output directory
int iNCores = 4; // Number of cores (Linux only)


int main()
{
  int begin, end;
  begin = time(NULL);

 #pragma omp parallel num_threads(iNCores)
  {
    int m, n, p, samp;
    int iM, f, i, j, k, l, start;
    double dVar1, dVar2, dCov;
    double dSteps = (double)iSteps;
    gsl_vector * vMu = gsl_vector_alloc(2);
    int tid = omp_get_thread_num(); // parallel


    // Simulation parameters (process)
        double alpha1 = 0.1, alpha2 = 0.1;
        double nu1 = 1, nu2 = 1;
        double H1 = 0.1, H2 = 0.1;
        double mu1 = 0, mu2 = 0;
        double dRho = 0.5, dEta = 0.1;
        double dDelta = 1/252;
        int iN = 7040; // Path length


    p = 2;
    n =  iN * iSteps;
    start = 0;
    double vH[2] = {H1, H2};
    double mRho[2][2] = {{1, dRho}, {dRho, 1}};
    double mEta[2][2] = {{0, dEta},{- dEta, 0}};
    int iInvDelta = (int)round(1.0/dDelta);
    gsl_vector_set(vMu, 0, mu1);
    gsl_vector_set(vMu, 1, mu2);
    iM = pow(2, floor(log(n)/log(2)) + 2);
    dVar1 = fVarFOU(nu1, alpha1, H1);
    dVar2 = fVarFOU(nu2, alpha2, H2);
    dCov = fCov2FOU(nu1, nu2, alpha1, alpha2, H1, H2, mRho[0][1], mEta[0][1]);
    m = pow(2, floor(log(n)/log(2)) + 2);
    gsl_complex temp;
    double * v;
    v = malloc(2 * m * sizeof(double));

    gsl_matrix_complex **tab_list = malloc(m * sizeof(gsl_matrix_complex*));

    for (i = 0; i < m; i++)
    {
        tab_list[i] = gsl_matrix_complex_alloc(p, p);
    }


    for(i = 0; i < p; i++)
    {
        for(j = 0; j < p; j++)
        {
            double * vpCuv = cuvj(iM, vH, mRho, mEta, i, j);
            for(k = 0; k < m; k++)
            {
                v[2 * k] = vpCuv[k];
                v[2 * k + 1] = 0;

            }
            free(vpCuv);
            gsl_fft_complex_radix2_forward(v, 1, m);

            for(l = 0; l < m; l++)
            {
                GSL_SET_COMPLEX(&temp, v[2 * l], v[2 * l + 1]);
                gsl_matrix_complex_set(tab_list[l], i, j, temp);
            }
        }
    }

    free(v);
    gsl_rng * rng_type = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng_type, time(NULL) + tid * 10);
    gsl_matrix *  mVar2FOU = gsl_matrix_alloc(2, 2);
    gsl_matrix_set(mVar2FOU, 0, 0, dVar1);
    gsl_matrix_set(mVar2FOU, 0, 1, dCov);
    gsl_matrix_set(mVar2FOU, 1, 0, dCov);
    gsl_matrix_set(mVar2FOU, 1, 1, dVar2);
    gsl_linalg_cholesky_decomp(mVar2FOU);
    for(samp = start; samp < iNPaths; samp++)
    {
        // printf("The thread %d  executes i = %d\n", tid, samp);
        gsl_matrix_complex * W_mat = gsl_matrix_complex_alloc(m, p);
        gsl_matrix_complex * Z_mat = gsl_matrix_complex_alloc(m, p);
        gsl_vector_complex * Z_vec = gsl_vector_complex_alloc(p);
        gsl_complex bar, bat, foo1, foo2;
        gsl_vector * vStart2FOU = gsl_vector_alloc(2);
        // Populating matrix Z with "Gaussian complex numbers"

        double U_rand, V_rand;
        for(i = 0; i < 2; i++)
        {
            U_rand = gsl_ran_gaussian(rng_type, 1);
            V_rand = gsl_ran_gaussian(rng_type, 1);
            GSL_SET_COMPLEX(&bar, U_rand/sqrt(m), 0);
            GSL_SET_COMPLEX(&bat, V_rand/sqrt(m), 0);
            gsl_matrix_complex_set(Z_mat, 0, i, bar);
            gsl_matrix_complex_set(Z_mat, m/2, i, bat);
        }


        for(f = 0; f < m; f++)
        {

            if((f>0) & (f < m/2))
            {
                for(k = 0; k < 2; k++)
                {
                    U_rand = gsl_ran_gaussian(rng_type, 1);
                    V_rand = gsl_ran_gaussian(rng_type, 1);

                    GSL_SET_COMPLEX(&foo1, U_rand/sqrt(2 * m),
                                    V_rand/sqrt(2 * m));
                    GSL_SET_COMPLEX(&foo2, U_rand/sqrt(2 * m),
                                    - V_rand/sqrt(2 * m));
                    gsl_matrix_complex_set(Z_mat, f, k, foo1);
                    gsl_matrix_complex_set(Z_mat, m - f, k, foo2);
                }
            }


            gsl_vector_complex * W_vec = gsl_vector_complex_alloc(p);
            gsl_matrix_complex * vecpAstore = gsl_matrix_complex_alloc(p, p);
            gsl_matrix_complex * matA = gsl_matrix_complex_alloc(p, p);
            gsl_permutation * perm_mat = gsl_permutation_alloc(p);
            gsl_complex droppodiagonale;
            gsl_matrix_complex * C_mat = gsl_matrix_complex_alloc(p, p);
            gsl_matrix_complex * D_mat = gsl_matrix_complex_alloc(p, p);
            gsl_vector * vpA = gsl_vector_alloc(p);
            gsl_matrix_complex * diag_sqrt_vpA = gsl_matrix_complex_alloc(p, p);
            gsl_matrix_complex * vecpA = gsl_matrix_complex_alloc(p, p);
            gsl_matrix_complex * vecpAinv = gsl_matrix_complex_alloc(p, p);
            int int_signum;

            // Initializing matrices (following Valgrid + ChatGPT suggestion)
            gsl_matrix_complex_set_zero(vecpAstore);
            gsl_matrix_complex_set_zero(matA);
            gsl_matrix_complex_set_zero(C_mat);
            gsl_matrix_complex_set_zero(D_mat);
            gsl_matrix_complex_set_zero(diag_sqrt_vpA);
            gsl_matrix_complex_set_zero(vecpA);
            gsl_matrix_complex_set_zero(vecpAinv);

            gsl_matrix_complex_memcpy(matA, tab_list[f]);
            gsl_eigen_hermv_workspace * wspace = gsl_eigen_hermv_alloc(p);


            gsl_eigen_hermv(matA, vpA, vecpA, wspace);
            gsl_matrix_complex_memcpy(vecpAstore, vecpA);


            if(gsl_vector_get(vpA, 1) > gsl_vector_get(vpA, 0))
            {
                gsl_matrix_complex_swap_columns(vecpA, 0, 1);
                gsl_matrix_complex_swap_columns(vecpAstore, 0, 1);
                gsl_vector_swap_elements(vpA, 0, 1);
            }


//             Inverting eigenvector matrix (vecpA to vecpAinv)

            gsl_linalg_complex_LU_decomp(vecpA, perm_mat, &int_signum);

            // EMMY DEPRECATES vecpA to nan
            gsl_linalg_complex_LU_invert(vecpA, perm_mat, vecpAinv);
            gsl_permutation_free(perm_mat);

            // Creating diagonal matrix
            for(i = 0; i < p; i++)
            {
                if(gsl_vector_get(vpA, i) < 0)
                {
                    gsl_vector_set(vpA, i, 0);
                }
                GSL_SET_COMPLEX(&droppodiagonale, sqrt(gsl_vector_get(vpA, i)), 0);
                gsl_matrix_complex_set(diag_sqrt_vpA, i, i, droppodiagonale);
            }


            // Performing matrix operations to obtain B
            gsl_complex alpha1, alpha2;
            GSL_REAL(alpha1) = 1.0;
            GSL_IMAG(alpha1) = 0.0;
            GSL_REAL(alpha2) = 0.0;
            GSL_IMAG(alpha2) = 0.0;


            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                           alpha1, vecpAstore, diag_sqrt_vpA,
                           alpha2,  C_mat);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                           alpha1, C_mat, vecpAinv,
                           alpha2,  D_mat);




            for(i = 0; i < p; i++)
            {
                gsl_vector_complex_set(Z_vec, i, gsl_matrix_complex_get(Z_mat, f, i));

            }

            gsl_blas_zgemv(CblasNoTrans, alpha1, D_mat, Z_vec, alpha2, W_vec);
            for(i = 0; i < p; i++)
            {

                gsl_matrix_complex_set(W_mat, f, i, gsl_vector_complex_get(W_vec, i));

            }


            gsl_matrix_complex_free(vecpA);
            gsl_matrix_complex_free(vecpAinv);
            gsl_matrix_complex_free(diag_sqrt_vpA);
            gsl_matrix_complex_free(vecpAstore);
            gsl_matrix_complex_free(matA);
            gsl_matrix_complex_free(C_mat);
            gsl_matrix_complex_free(D_mat);
            gsl_vector_complex_free(W_vec);
            gsl_vector_free(vpA);
            gsl_eigen_hermv_free(wspace);

        }


        double * X_GN, * Y_GN, * vOut;
        X_GN = malloc(2 * m * sizeof(double));
        Y_GN = malloc(2 * m * sizeof(double));
        vOut = malloc(2 * n * sizeof(double));
        for(i = 0; i < m; i++)
        {
            X_GN[2 * i] = GSL_REAL(gsl_matrix_complex_get(W_mat, i, 0));
            X_GN[2 * i + 1] = GSL_IMAG(gsl_matrix_complex_get(W_mat, i, 0));

            Y_GN[2 * i] = GSL_REAL(gsl_matrix_complex_get(W_mat, i, 1));
            Y_GN[2 * i + 1] = GSL_IMAG(gsl_matrix_complex_get(W_mat, i, 1));
        }

        gsl_fft_complex_radix2_forward(X_GN, 1, m);
        gsl_fft_complex_radix2_forward(Y_GN, 1, m);


        gsl_ran_multivariate_gaussian(rng_type, vMu, mVar2FOU, vStart2FOU);
        vOut[0] = gsl_vector_get(vStart2FOU, 0);
        vOut[n] = gsl_vector_get(vStart2FOU , 1);


        for(j = 1; j < n; j++)
        {
            vOut[j] = vOut[j - 1] + alpha1 *
                          (mu1 - vOut[j - 1]) * (1/dSteps) * dDelta + nu1 * pow((1/dSteps) * dDelta, H1) * X_GN[2 * j];
            vOut[n + j] = vOut[n + j - 1] + alpha2 *
                          (mu2 - vOut[n + j - 1]) * (1/dSteps) * dDelta + nu2 * pow((1/dSteps) * dDelta, H2) * Y_GN[2 * j];
        }

        // Writing the output
        char fileSpec[strlen(dir) + strlen(name) + strlen(extension) + 150];
        snprintf(fileSpec, sizeof(fileSpec), "%s%s_H1%.2f_H2%.2f_a1%.2f_a2%.2f_v1%.2f_v2%.2f_rho%.2f_eta%.2f_n%i_d%i_s%i_%i_%i%s",
                 dir, name, vH[0], vH[1], alpha1, alpha2, nu1, nu2, mRho[0][1], mEta[0][1], iN, iInvDelta, iSteps, samp, tid, extension);

        FILE * fpt;
        fpt = fopen(fileSpec, "w+");
         fprintf(fpt, "%f", vOut[0]);  // Print the first element without comma

        for (i = 1; i < 2 * iN; i++) {
            fprintf(fpt, ", %f", vOut[iSteps * i]);  // Print the rest with comma
        }
//        for(i = 0; i < 2 * iN; i++)
//        {
//            fprintf(fpt, "%f, ", vOut[iSteps * i]);
//        }
        fclose(fpt);

        free(X_GN);
        free(Y_GN);
        free(vOut);
        gsl_matrix_complex_free(W_mat);
        gsl_matrix_complex_free(Z_mat);
        gsl_vector_complex_free(Z_vec);
        gsl_vector_free(vStart2FOU);
    }
  gsl_matrix_free(mVar2FOU);
  // Assume that each element in tab_list was allocated dynamically
    for (size_t i = 0; i < m; ++i) {
    if (tab_list[i] != NULL) {
        gsl_matrix_complex_free(tab_list[i]);  // Free each gsl_matrix_complex
    }
}

// Finally, free the array of pointers itself
free(tab_list);
gsl_rng_free(rng_type);
gsl_vector_free(vMu);
  }
    end = time(NULL);
//    printf("\nTotal time employed: %f minutes\n", ((double) (end - begin)/60));

    return 0;
}



double * cuvj(int m, double H[2], double rho_mat[2][2], double eta_mat[2][2], int u, int v)
{
    double * z, Hjk, rhojk, etajk, etakj;
    int * in1, in2, * in3, j;

    Hjk = H[u] + H[v];
    rhojk = rho_mat[u][v];
    etajk = eta_mat[u][v];
    etakj = eta_mat[v][u];

    in1 = seq_by1(1, m/2 - 1);
    in2 = m/2;
    in3 = seq_by1(m - (m/2 + 1), m - (m - 1));

    double * gamma1 = gammauv(in1, rhojk, etajk, Hjk, m/2 - 1);
    double * gamma21 = gammauv(&in2, rhojk, etajk, Hjk, 1);
    double * gamma22 = gammauv(&in2, rhojk, etakj, Hjk, 1);
    double * gamma3 = gammauv(in3, rhojk, etakj, Hjk, m/2 - 1);

    z = (double*)calloc(m, sizeof(double));
    z[0] = rhojk;
    z[m/2] = (gamma21[0] + gamma22[0])/2;
    for(j = 1; j < m/2; j++)
    {
        z[j] = gamma1[j - 1];
        z[m/2 + j] = gamma3[j - 1];
    }

    free(gamma1);
    free(gamma21);
    free(gamma22);
    free(gamma3);
    free(in1);
    free(in3);

    return z;
}

double * gammauv(int * h, double rhojk, double etajk, double Hjk, int length_h)   // Translate to GSL
{
    double * out; //* out_minus, * out, * out_at, * out_plus;
    int * h_minus, * h_plus, i;

//    out_minus = (double*)calloc(length_h, sizeof(double));
    out = (double*)calloc(length_h, sizeof(double));
//    out_plus = (double*)calloc(length_h, sizeof(double));
    h_minus = (int*)calloc(length_h, sizeof(int));
    h_plus = (int*)calloc(length_h, sizeof(int));

    for(i = 0; i < length_h; i++)
    {
        h_minus[i] = h[i] - 1;
        h_plus[i] = h[i] + 1;
    }
    double * out_minus = wjk(h_minus, rhojk, etajk, Hjk, length_h);
    double * out_plus = wjk(h_plus, rhojk, etajk, Hjk, length_h);
    double * out_at = wjk(h, rhojk, etajk, Hjk, length_h);

    for(i = 0; i < length_h; i++)
    {
        out[i] = 0.5*(out_minus[i] - 2 * out_at[i] + out_plus[i]);
    }

    free(h_minus);
    free(h_plus);
    free(out_minus);
    free(out_plus);
    free(out_at);

    return out;
}



double * wjk(int h[], double rhojk, double etajk, double Hjk, int length_h)
{
    double *tmp, xlogx_tmp;
    int sign_e, i;

    tmp = (double*)calloc(length_h, sizeof(double));
    for(i = 0; i < length_h; i++)
    {
        if(Hjk != 1)
        {
            if(h[i] < 0)
            {
                sign_e = -1;
            }
            else if(h[i] > 0)
            {
                sign_e = +1;
            }
            else
            {
                sign_e = 0;
            }
            tmp[i] = (rhojk - etajk * sign_e) * pow(fabs(h[i]), Hjk); // fix sign
        }
        else
        {
            if(h[i] == 0)
            {
                xlogx_tmp = 0;
            }
            else
            {
                xlogx_tmp = h[i] * log(fabs(h[i]));
            }
            tmp[i] = rhojk * fabs(h[i]) - etajk * xlogx_tmp;
        }
    }
    return tmp;
}


int * seq_by1(int from, int to)
{
    int * seq_out, seq_length, i;

    if(to > from)
    {
        seq_length = to - from;
        seq_out = (int*)calloc(seq_length + 1, sizeof(int));
        for(i = 0; i <= seq_length; i++)
        {
            seq_out[i] = from + i;
        }
    }
    else if (to < from)
    {
        seq_length = from - to;
        seq_out = (int*)calloc(seq_length + 1, sizeof(int));
        for(i = 0; i <= seq_length; i++)
        {
            seq_out[i] = from - i;
        }
    }
    return seq_out;
}

double fVarFOU(double nu, double alpha, double H)
{
    double z;
    z = pow(nu, 2)/ (2 * pow(alpha, 2 * H)) * tgamma(1 + 2 * H);
    return z;
}

double fCov2FOU(double nu1, double nu2, double alpha1, double alpha2, double H1, double H2, double rho, double eta)
{
    double val;
    if(H1 + H2 != 1) {
        val = tgamma(H1 + H2 + 1) * nu1 * nu2 * 1/(2 * (alpha1 + alpha2)) *
        ((pow(alpha1, (1 - H1 - H2)) + pow(alpha2, (1 - H1 - H2))) * rho +
         (pow(alpha2, (1 - H1 - H2)) - pow(alpha1, (1 - H1 - H2))) * eta);
    }
    else if(H1 + H2 == 1){
        val = nu1 * nu2 * 1/(alpha1 + alpha2) * (rho + eta/2 * (log(alpha2) - log(alpha1)));
    }
        return val;
}


