#include <functional>
#include "math.h"

double flin ( int n, int j, double l, std::function<double(double[], int)> f,
  double x[], int &nf, double v[], double q0[], double q1[], double &qd0, 
  double &qd1, double &qa, double &qb, double &qc );
void minfit ( int n, double tol, double a[], double q[] );
void minny ( int n, int j, int nits, double &d2, double &x1, double &f1, 
  bool fk, std::function<double(double[], int)> f, double x[], double t, double h,
  double v[], double q0[], double q1[], int &nl, int &nf, double dmin, 
  double ldt, double &fx, double &qa, double &qb, double &qc, double &qd0, 
  double &qd1 );
double praxis ( double t0, double h0, int n, int prin, double x[], 
	std::function<double(double[], int)> f);
void print2 ( int n, double x[], int prin, double fx, int nf, int nl );
void quad ( int n, std::function<double(double[], int)> f, double x[], double t,
  double h, double v[], double q0[], double q1[], int &nl, int &nf, double dmin, 
  double ldt, double &fx, double &qf1, double &qa, double &qb, double &qc, 
  double &qd0, double &qd1 );
