#include <functional>
#include "math.h"

double *compass_search ( std::function<double(int, double[])>, int m, 
  double x0[], double delta_tol, double delta_init, int k_max, double &fx, 
  int &k );
