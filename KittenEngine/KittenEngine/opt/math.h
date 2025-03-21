#pragma once
#include <string>

using namespace std;

double r8_epsilon();
double r8_hypot(double x, double y);
double r8_max(double x, double y);
double r8_min(double x, double y);
double r8_uniform_01(int& seed);
void r8mat_print(int m, int n, double a[], string title);
void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi,
	int jhi, string title);
void r8mat_transpose_in_place(int n, double a[]);
void r8vec_copy(int n, double a1[], double a2[]);
double r8vec_max(int n, double r8vec[]);
double r8vec_min(int n, double r8vec[]);
double r8vec_norm(int n, double a[]);
void r8vec_print(int n, double a[], string title);
void svsort(int n, double d[], double v[]);

double r8_abs(double x);
