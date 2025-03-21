#include "math.h"
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

double r8_epsilon()

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
    const double value = 2.220446049250313E-016;

    return value;
}
//****************************************************************************80

double r8_hypot(double x, double y)

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the arguments.
//
//    Output, double R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
//
{
    double a;
    double b;
    double value;

    if (fabs(x) < fabs(y)) {
        a = fabs(y);
        b = fabs(x);
    }
    else {
        a = fabs(x);
        b = fabs(y);
    }
    //
    //  A contains the larger value.
    //
    if (a == 0.0) {
        value = 0.0;
    }
    else {
        value = a * sqrt(1.0 + (b / a) * (b / a));
    }

    return value;
}
//****************************************************************************80

double r8_max(double x, double y)

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Discussion:
//
//    The C++ math library provides the function fmax() which is preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
    double value;

    if (y < x) {
        value = x;
    }
    else {
        value = y;
    }
    return value;
}
//****************************************************************************80

double r8_min(double x, double y)

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Discussion:
//
//    The C++ math library provides the function fmin() which is preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
    double value;

    if (y < x) {
        value = y;
    }
    else {
        value = x;
    }
    return value;
}
//****************************************************************************80

double r8_uniform_01(int& seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
    const int i4_huge = 2147483647;
    int k;
    double r;

    if (seed == 0) {
        cerr << "\n";
        cerr << "R8_UNIFORM_01 - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }

    k = seed / 127773;

    seed = 16807 * (seed - k * 127773) - k * 2836;

    if (seed < 0) {
        seed = seed + i4_huge;
    }
    r = (double)(seed) * 4.656612875E-10;

    return r;
}
//****************************************************************************80

void r8mat_print(int m, int n, double a[], string title)

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
    r8mat_print_some(m, n, a, 1, 1, m, n, title);

    return;
}
//****************************************************************************80

void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi,
    int jhi, string title)

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    R8MAT_PRINT_SOME prints some of an R8MAT.
    //
    //  Discussion:
    //
    //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    //    in column-major order.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    26 June 2013
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, int M, the number of rows of the matrix.
    //    M must be positive.
    //
    //    Input, int N, the number of columns of the matrix.
    //    N must be positive.
    //
    //    Input, double A[M*N], the matrix.
    //
    //    Input, int ILO, JLO, IHI, JHI, designate the first row and
    //    column, and the last row and column to be printed.
    //
    //    Input, string TITLE, a title.
    //
{
# define INCX 5

    int i;
    int i2hi;
    int i2lo;
    int j;
    int j2hi;
    int j2lo;

    cout << "\n";
    cout << title << "\n";

    if (m <= 0 || n <= 0) {
        cout << "\n";
        cout << "  (None)\n";
        return;
    }
    //
    //  Print the columns of the matrix, in strips of 5.
    //
    for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX) {
        j2hi = j2lo + INCX - 1;
        if (n < j2hi) {
            j2hi = n;
        }
        if (jhi < j2hi) {
            j2hi = jhi;
        }
        cout << "\n";
        //
        //  For each column J in the current range...
        //
        //  Write the header.
        //
        cout << "  Col:    ";
        for (j = j2lo; j <= j2hi; j++) {
            cout << setw(7) << j - 1 << "       ";
        }
        cout << "\n";
        cout << "  Row\n";
        cout << "\n";
        //
        //  Determine the range of the rows in this strip.
        //
        if (1 < ilo) {
            i2lo = ilo;
        }
        else {
            i2lo = 1;
        }
        if (ihi < m) {
            i2hi = ihi;
        }
        else {
            i2hi = m;
        }

        for (i = i2lo; i <= i2hi; i++) {
            //
            //  Print out (up to) 5 entries in row I, that lie in the current strip.
            //
            cout << setw(5) << i - 1 << ": ";
            for (j = j2lo; j <= j2hi; j++) {
                cout << setw(12) << a[i - 1 + (j - 1) * m] << "  ";
            }
            cout << "\n";
        }
    }

    return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_in_place(int n, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_IN_PLACE transposes a square R8MAT in place.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input/output, double A[N*N], the matrix to be transposed.
//
{
    int i;
    int j;
    double t;

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            t = a[i + j * n];
            a[i + j * n] = a[j + i * n];
            a[j + i * n] = t;
        }
    }
    return;
}
//****************************************************************************80

void r8vec_copy(int n, double a1[], double a2[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
    int i;

    for (i = 0; i < n; i++) {
        a2[i] = a1[i];
    }
    return;
}
//****************************************************************************80

double r8vec_max(int n, double r8vec[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
    int i;
    double value;

    value = r8vec[0];

    for (i = 1; i < n; i++) {
        if (value < r8vec[i]) {
            value = r8vec[i];
        }
    }
    return value;
}
//****************************************************************************80

double r8vec_min(int n, double r8vec[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
    int i;
    double value;

    value = r8vec[0];

    for (i = 1; i < n; i++) {
        if (r8vec[i] < value) {
            value = r8vec[i];
        }
    }
    return value;
}
//****************************************************************************80

double r8vec_norm(int n, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM, the L2 norm of A.
//
{
    int i;
    double v;

    v = 0.0;

    for (i = 0; i < n; i++) {
        v = v + a[i] * a[i];
    }
    v = sqrt(v);

    return v;
}
//****************************************************************************80

void r8vec_print(int n, double a[], string title)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
    int i;

    cout << "\n";
    cout << title << "\n";
    cout << "\n";
    for (i = 0; i < n; i++) {
        cout << "  " << setw(8) << i
            << ": " << setw(14) << a[i] << "\n";
    }

    return;
}
//****************************************************************************80

void svsort(int n, double d[], double v[])

//****************************************************************************80
//
//  Purpose:
//
//    SVSORT descending sorts D and adjusts the corresponding columns of V.
//
//  Discussion:
//
//    A simple bubble sort is used on D.
//
//    In our application, D contains singular values, and the columns of V are
//    the corresponding right singular vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 August 2016
//
//  Author:
//
//    Original FORTRAN77 version by Richard Brent.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization with Derivatives,
//    Prentice Hall, 1973,
//    Reprinted by Dover, 2002.
//
//  Parameters:
//
//    Input, int N, the length of D, and the order of V.
//
//    Input/output, double D[N], the vector to be sorted.  
//    On output, the entries of D are in descending order.
//
//    Input/output, double V[N,N], an N by N array to be adjusted 
//    as D is sorted.  In particular, if the value that was in D(I) on input is
//    moved to D(J) on output, then the input column V(*,I) is moved to
//    the output column V(*,J).
//
{
    int i;
    int j1;
    int j2;
    int j3;
    double t;

    for (j1 = 0; j1 < n - 1; j1++) {
        //
        //  Find J3, the index of the largest entry in D[J1:N-1].
        //  MAXLOC apparently requires its output to be an array.
        //
        j3 = j1;
        for (j2 = j1 + 1; j2 < n; j2++) {
            if (d[j3] < d[j2]) {
                j3 = j2;
            }
        }
        //
        //  If J1 != J3, swap D[J1] and D[J3], and columns J1 and J3 of V.
        //
        if (j1 != j3) {
            t = d[j1];
            d[j1] = d[j3];
            d[j3] = t;
            for (i = 0; i < n; i++) {
                t = v[i + j1 * n];
                v[i + j1 * n] = v[i + j3 * n];
                v[i + j3 * n] = t;
            }
        }
    }

    return;
}


double r8_abs(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
    double value;

    if (0.0 <= x) {
        value = x;
    }
    else {
        value = -x;
    }
    return value;
}