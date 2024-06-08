#include <stdio.h>

#define MAX_N 10

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]);

double S(double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[]);

int main()
{
    int n, Type, m, i;
    double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
    double s0, sn, Fmax, t0, tm, h, t;

    scanf("%d", &n);
    for (i = 0; i <= n; i++)
        scanf("%lf", &x[i]);
    for (i = 0; i <= n; i++)
        scanf("%lf", &f[i]);
    scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);

    Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
    for (i = 1; i <= n; i++)
        printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

    scanf("%lf %lf %d", &t0, &tm, &m);
    h = (tm - t0) / (double)m;
    for (i = 0; i <= m; i++)
    {
        t = t0 + h * (double)i;
        printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
    }

    return 0;
}

/* Your functions will be put here */
#include <math.h>

void step1(int n, double h[], double x[])
{
    for (int i = 0; i < n; i++)
        h[i] = x[i + 1] - x[i];
}

void step2(int n, double a[], double f[])
{
    for (int i = 0; i <= n; i++)
        a[i + 1] = f[i];
}

void step3(int n, double al[], double f[], double h[], double s0, double sn)
{
    al[0] = 3 * (f[1] - f[0]) / h[0] - 3 * s0;
    al[n] = 3 * sn - 3 * (f[n] - f[n - 1]) / h[n - 1];
}
void step4(int n, double al[], double f[], double h[], double l[], double u[], double z[])
{
    for (int i = 1; i < n; i++)
        al[i] = 3 / h[i] * (f[i + 1] - f[i]) - 3 / h[i - 1] * (f[i] - f[i - 1]);
}
void step5(int n, double l[], double x[], double h[], double u[], double z[], double al[])
{
    for (int i = 1; i < n; i++)
    {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (al[i] - h[i - 1] * z[i - 1]) / l[i];
    }
}
void step6(int n, double h[], double c[], double a[], double b[], double d[], double z[], double u[])
{
    for (int j = n - 1; j >= 0; j--)
    {
        c[j + 1] = z[j] - u[j] * c[j + 2];
        b[j + 1] = (a[j + 2] - a[j + 1]) / h[j] - h[j] * (c[j + 2] + 2 * c[j + 1]) / 3;
        d[j + 1] = (c[j + 2] - c[j + 1]) / (3 * h[j]);
    }
}

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[])
{
    double x1[MAX_N];
    double x2[MAX_N];
    double x3[MAX_N];
    double x4[MAX_N];
    double x5[MAX_N];

    if (Type == 1)
    {
        step1(n, x1, x);

        step2(n, a, f);

        step3(n, x2, f, x1, s0, sn);

        step4(n, x2, f, x1, x3, x5, x4);

        x3[0] = 2 * x1[0];
        x4[0] = x2[0] / x3[0];
        x5[0] = 0.5;

        step5(n, x3, x, x1, x5, x4, x2);

        x3[n] = x1[n - 1] * (2 - x5[n - 1]);
        x4[n] = (x2[n] - x1[n - 1] * x4[n - 1]) / x3[n];
        c[n + 1] = x4[n];
        step6(n, x1, c, a, b, d, x4, x5);
    }
    else
    {
        step1(n, x1, x);

        step2(n, a, f);

        step4(n, x2, f, x1, x3, x5, x4);

        x4[0] = s0 / 2;
        x5[0] = 0;
        x3[0] = 1;

        step5(n, x3, x, x1, x5, x4, x2);

        x4[n] = sn / 2;
        x3[n] = 1;
        c[n + 1] = sn / 2;

        step6(n, x1, c, a, b, d, x4, x5);
    }
}

double S(double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[])
{
    if (t < x[0] || t > x[n])
    {
        return Fmax;
    }
    else if (t == x[0])
    {
        return a[1];
    }
    else
    {
        int j = 1;
        while (!(x[j - 1] < t && t <= x[j]))
            ++j;
        return a[j] + b[j] * (t - x[j - 1]) + c[j] * pow(t - x[j - 1], 2) + d[j] * pow(t - x[j - 1], 3);
    }
}
