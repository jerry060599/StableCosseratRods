#include <functional>

void nelmin ( std::function<double(double[])> fn, int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault );
