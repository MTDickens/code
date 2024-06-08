#include <stdio.h>
#include <math.h>

#define MAX_m 200
#define MAX_n 6

double identity(double x)
{
    return pow(x, 4);
}

double f1(double x)
{
    return sin(x);
}

double f2(double x)
{
    return exp(x);
}

int OPA(double (*f)(double t), int m, double x[], double w[], double c[], double *eps);

void print_results(int n, double c[], double eps)
{
    int i;

    printf("%d\n", n);
    for (i = 0; i <= n; i++)
        printf("%12.4e ", c[i]);
    printf("\n");
    printf("error = %9.2e\n", eps);
    printf("\n");
}

int main()
{
    int m, i, n;
    double x[MAX_m], w[MAX_m], c[MAX_n + 1], eps;

    // // begin debugging
    // // m = 90;
    // m = 5;
    // for (i = 0; i < m; i++)
    // {
    //     x[i] = 1 + i;
    //     w[i] = 1;
    // }
    // eps = 0.001;
    // n = OPA(identity, m, x, w, c, &eps);
    // print_results(n, c, eps);
    // // end debugging

    m = 90;
    for (i = 0; i < m; i++)
    {
        x[i] = 3.1415926535897932 * (double)(i + 1) / 180.0;
        w[i] = i;
    }
    eps = 1e-10;
    n = OPA(f1, m, x, w, c, &eps);
    print_results(n, c, eps);

    m = 200;
    for (i = 0; i < m; i++)
    {
        x[i] = 0.01 * (double)i;
        w[i] = i;
    }
    eps = 1e-10;
    n = OPA(f2, m, x, w, c, &eps);
    print_results(n, c, eps);

    return 0;
}

/* Your function will be put here */
#include <string.h>

double inner_prod(double matrix[][MAX_m], double w[], int l1, int l2, int m)
{
    double ret = 0.0;
    for (int j = 0; j != m; ++j)
    {
        ret += w[j] * matrix[l1][j] * matrix[l2][j];
    }
    return ret;
}

double modulo(double vec[], double w[], int m)
{
    double ret = 0.0;
    for (int j = 0; j != m; ++j)
    {
        ret += w[j] * vec[j] * vec[j];
    }
    return sqrt(ret);
}

int OPA(double (*f)(double t), int m, double x[], double w[], double c[], double *eps)
{
    // y[i] = f(x[i])
    double y[MAX_m];
    // ortho_c[i]
    double ortho_c[MAX_n + 1];
    // ortho_base[i][j] = phi_i(x_{j+1})
    double ortho_base[MAX_n + 1][MAX_m];
    double ortho_polynomial[MAX_n + 1][MAX_n + 1];
    for (int i = 0; i != m; ++i)
    {
        y[i] = f(x[i]);
    }
    // Let it be natural base first
    for (int i = 0; i <= MAX_n; ++i)
        for (int j = 0; j != m; ++j)
        {
            ortho_base[i][j] = pow(x[j], (double)i);
        }

    for (int i = 0; i <= MAX_n; ++i)
        for (int j = 0; j <= MAX_n; ++j)
        {
            ortho_polynomial[i][j] = (i == j ? 1.0 : 0.0);
        }
    // Make it orthogonal
    for (int n = 0; n <= MAX_n; ++n)
    {
        // First, make them orthogonal
        double new_ortho_vec[MAX_m];
        memcpy(new_ortho_vec, ortho_base[n], sizeof(double) * m);
        // vec -= sum_j(<original_vec, phi_j> * phi_j)
        for (int j = 0; j <= n - 1; ++j)
        {
            // vec -= <original_vec, phi_j> * phi_j
            double inner = inner_prod(ortho_base, w, n, j, m);
            for (int k = 0; k != m; ++k)
            {
                new_ortho_vec[k] -= inner * ortho_base[j][k];
            }
            for (int k = 0; k <= n; ++k)
            {
                ortho_polynomial[n][k] -= inner * ortho_polynomial[j][k];
            }
        }

        // vec /= sqrt(<vec, vec>)
        double inner_i = modulo(new_ortho_vec, w, m);
        for (int k = 0; k != m; ++k)
        {
            ortho_base[n][k] = new_ortho_vec[k] / inner_i;
        }
        for (int k = 0; k <= n; ++k)
        {
            ortho_polynomial[n][k] = ortho_polynomial[n][k] / inner_i;
        }

        // calculate the current eps
        // add n-th element for ortho_c
        // printf("@@@@@@%p,%p@@@@@@ ", &y, &ortho_c);
        for (int k = 0; k != m; ++k)
        {
            // printf("@@@@@@%12.4e, %12.4e@@@@@@ ", y[0], ortho_c[0]);
            ortho_c[n] += ortho_base[n][k] * w[k] * y[k];
        }
        double current_eps = 0.0;

        for (int j = 0; j != m; ++j)
        {
            double delta = 0.0;
            for (int k = 0; k <= n; ++k)
            {
                delta += ortho_c[k] * ortho_base[k][j];
            }
            current_eps += (delta - y[j]) * (delta - y[j]) * w[j];
        }
        printf("Round[%d]: %12.4e\n", n, current_eps);
        // if (current_eps < *eps || n == MAX_n)
        // {
        //     *eps = current_eps;
        //     // return c
        //     for (int i = 0; i <= n; ++i)
        //     {
        //         c[i] = 0.0;
        //         for (int j = 0; j <= n; ++j)
        //         {
        //             c[i] += ortho_c[j] * ortho_polynomial[i][j];
        //         }
        //     }
        //     return n;
        // }

        if (current_eps < *eps || n == MAX_n)
        {
            *eps = current_eps;
            // return c
            for (int i = 0; i <= n; ++i)
            {
                c[i] = 0.0;
                for (int j = 0; j <= n; ++j)
                {
                    c[i] += ortho_c[j] * ortho_polynomial[j][i];
                }
            }
            // printf("ortho_base:\n");
            // for (int i = 0; i <= n; ++i)
            // {
            //     for (int j = 0; j < m; ++j)
            //     {
            //         printf("%12.4e ", ortho_base[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("ortho_polynomial:\n");
            // for (int i = 0; i <= n; ++i)
            // {
            //     for (int j = 0; j <= n; ++j)
            //     {
            //         printf("%12.4e ", ortho_polynomial[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("ortho_c:\n");
            // for (int i = 0; i <= n; ++i)
            // {
            //     printf("%12.4e\n", ortho_c[i]);
            // }
            // printf("w:\n");
            // for (int i = 0; i < m; ++i)
            // {
            //     printf("%12.4e ", w[i]);
            // }
            // printf("\n");

            // printf("y:\n");
            // for (int i = 0; i < m; ++i)
            // {
            //     printf("%12.4e ", y[0]);
            // }
            // printf("\n");
            return n;
        }
    }
}