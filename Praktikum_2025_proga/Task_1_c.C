#include <iostream>
#include <stdio.h>
#include <float.h>
#include <math.h>
using namespace ::std;

double F_hyper(double a, double b, double c, double z)
{
  const double accuracy = 1.0e-15; // DBL_EPSILON
  const double N_break_max = 5000.0;
  const double integer_accuracy = 0.001;

  bool is_abc_int = fabs(rint(c - a - b) - (c - a - b)) < integer_accuracy;

  if (fabs(z) > 1)
  {
    cout << "Lirikk's WARNING! Series does not converge, cause |z| > 1" << endl;
    return (NAN);
  }
  //second check for special case: 
  else if(0.5<z && z<=1.0 && (c) && !is_abc_int)
  {
    //cout << "c - a - b is not integer"<<endl;
    double y = tgamma(c)*tgamma(c-a-b)/(tgamma(c-a)*tgamma(c-b))
                *F_hyper(a,b,a+b-c+1.,1.-z) 
                + pow((1.-z),c-a-b) * tgamma(c)*tgamma(a+b-c)/(tgamma(a)*tgamma(b)) 
                * F_hyper(c-a,c-b,c-a-b+1.,1.-z);
    return (y);
  }
  else if(z<-0.5)
  {
    //cout<<"z<0.5"<<endl;
    double y = pow((1-z), -a) * F_hyper(a,c-b,c,z/(z-1));
    return(y);
  }
  else
  {
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
    //cout << "Number of steps: " << N << endl;
    return (Sum);
  }
}

void Task_1_c()
{

  //for "sqrt" function we can use formula (5) if z<-1/2
  // double Hyper_sqrt = F_hyper(-1./4.,1./4., 1./2., -z*z);
  // double Theor_sqrt = 0.5*(sqrt(sqrt(1+z*z)+z)+sqrt(sqrt(1+z*z)-z));

  //for "sqrt" function we also can use formula (5) if z<-1/2
  const double z_min_ln = -0.95;
  const double z_max_ln = 0.95;
  const double z_step_ln = 0.1;
  int N_steps_ln = (int)((z_max_ln-z_min_ln)/z_step_ln) + 1;
  
  printf("\n%-6s\t%-20s\t%-20s\t%-12s\n", "z", "F21(1,1,2,z)", "-ln(1-z)/z", "( eps )");
  for(int i=0;i<=N_steps_ln;i++)
  {
  double z = z_min_ln+i*z_step_ln;
  double Hyper_ln = F_hyper(1.,1.,2.,z);
  double Theor_ln = (fabs(z) < DBL_EPSILON) ? 1.0 : -log(1-z)/z;
  double eps_ln = fabs((Hyper_ln - Theor_ln)/Theor_ln);
  
  printf("%-6.2f\t%20.17f\t%20.17f\t%-10.3e\n", z, Hyper_ln, Theor_ln, eps_ln);
  }

  //for next "arcsin" function we can use recurrent formula (4) if c-a-b not int

  const double z_min_asin = 0.05;
  const double z_max_asin = 0.95;
  const double z_step_asin = 0.1;
  int N_steps_asin = (int)((z_max_asin-z_min_asin)/z_step_asin) + 1;
  
  printf("\n%-6s\t%-20s\t%-20s\t%-12s\n", "z", "F21(1,1,2,z)", "-ln(1-z)/z", "( eps )");
  for(int i=0;i<=N_steps_asin;i++)
  {
  double z = z_min_asin+i*z_step_asin;
  double Hyper_asin = F_hyper(0.5,0.5,1.5,z*z);
  double Theor_asin = (fabs(z) < DBL_EPSILON) ? 1.0 : asin(z)/z;
  double eps_asin = fabs((Hyper_asin - Theor_asin)/Theor_asin);
  
  printf("%-6.2f\t%20.17f\t%20.17f\t%-10.3e\n", z, Hyper_asin, Theor_asin, eps_asin);
  }
  

  // printf("%-16s = %20.17f\n", "Hyper sqrt", Hyper_sqrt);
  // printf("%-16s = %20.17f\n", "Exact sqrt", Theor_sqrt);
  // printf("%-16s = %20.17f\n", "Hyper arcsin", Hyper_arcsin);
  // printf("%-16s = %20.17f\n", "Exact arcsin", Theor_arcsin);
  
  // printf("Eps = %9.1e\n", eps);
  // tgamma(10); rint(); pow(2,3)
}