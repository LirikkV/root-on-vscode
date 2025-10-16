#include <iostream>
#include <stdio.h>
#include <float.h>
#include <math.h>
using namespace ::std;

double F_hyper(double a, double b, double c, double z)
{
  const double accuracy = 1.0e-15;//DBL_EPSILON
  const double N_break_max = 5000.0;
  
  double N = 1.0;

  double X_0 = 1.0;
  double eps = 1.0;
  double Sum = X_0;

  while (eps > accuracy && N < N_break_max)
  {
    double X_multiplier = (a + N - 1.0) * (b + N - 1.0) * z / ((c + N - 1.0) * N);
    X_0 *= X_multiplier;
    Sum += X_0;
    N += 1.0;
    eps = fabs(X_0);
  }

  cout << "Number of steps: " << N << endl;
  return (Sum);

}

void Task_1_c()
{
  double z = 0.99;

  double x = F_hyper(0.5, 1.0, 1.5, -z*z);
  printf("Hypergeometriz = %20.17f\n", x);
  double y = (fabs(z)<DBL_EPSILON) ? 1.0 : atan(z)/z;
  printf("Exact = %20.17f\n", y);
  double eps = fabs((x - y)/y);
  printf("Eps = %9.1e\n", eps);
  //tgamma(10); rint(); pow(2,3)
}