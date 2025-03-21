#include <functional>
#include "math.h"

double best_nearby ( double delta[], double point[], double prevbest, 
  int nvars, std::function<double(double[], int)>, int *funevals );
int hooke ( int nvars, double startpt[], double endpt[], double rho, double eps, 
  int itermax, std::function<double(double[], int)> );
