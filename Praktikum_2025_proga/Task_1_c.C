#include <iostream>
#include <stdio.h>
#include <float.h>
#include <math.h>
using namespace ::std;

double F_hyper(double a, double b, double c, double z)
{
  const double accuracy = 0.01;

  double N_break_max = 5000.0;
  double N = 1.0;

  double X_0 = 1.0;
  double X_multiplier = a * b * z / c;
  double eps = fabs(X_0 * X_multiplier);
  double Sum = X_0;
  while (eps > accuracy && eps > DBL_EPSILON && N < N_break_max)
  {
    X_multiplier = (a + N - 1.0) * (b + N - 1.0) * z / ((c + N - 1.0) * N);
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

  double x = F_hyper(0.5, 1.0, 1.5, -0.8 * 0.8);
  printf("%e\n", x);
}