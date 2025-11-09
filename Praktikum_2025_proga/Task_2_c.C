#include <iostream>
#include <fstream>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "TH1D.h"
#include "TFile.h"

using namespace ::std;

double my_random_double_0_1() {
    return (double)rand() / RAND_MAX;
}

double my_random_Gauss()
{
  double U=my_random_double_0_1();
  double V=my_random_double_0_1();
  double X = sqrt(-2*log(U))*cos(2*M_PI*V);
  return X;
}

void fill_random_0_1(double array[],int size)
{
  for(int i=0;i<size;i++)
  {
    array[i]=my_random_double_0_1();
  }
}

void fill_random_Gauss(double array[],int size)
{
  for(int i=0;i<size;i++)
  {
    array[i]=my_random_Gauss();
  }
}

double stat_x_mean(double array[], int N)
{
  double X=0.;
  for(int i=0;i<N;i++)
  {
    X+=array[i];
  }
  X=X/((double)N);
  return X;
}
double stat_s(double array[], int N)
{
  double Sum =0;
  double x_mean= stat_x_mean(array,N);
  for(int i=0;i<N;i++)
  {
    Sum+=(array[i]-x_mean)*(array[i]-x_mean);
  }
  return sqrt(Sum/(double)N);
}
double stat_gamma_1(double array[], int N)
{
  double x_mean =stat_x_mean(array,N);
  double Sum=0.;
  for(int i=0;i<N;i++)
  {
    Sum+=(array[i]-x_mean)*(array[i]-x_mean)*(array[i]-x_mean);
  }
  double mu_3=Sum/((double)N);
  double s = stat_s(array,N);

  return mu_3/(s*s*s);
}
double stat_gamma_2(double array[], int N)
{
  double x_mean =stat_x_mean(array,N);
  double Sum=0.;
  for(int i=0;i<N;i++)
  {
    Sum+=(array[i]-x_mean)*(array[i]-x_mean)*(array[i]-x_mean)*(array[i]-x_mean);
  }
  double mu_4=Sum/((double)N);
  double s = stat_s(array,N);

  return mu_4/(s*s*s*s) - 3.;
}

void Task_2_c()
{
  srand((unsigned int)time(NULL));

  /*
  TH1D* hLandau = new TH1D("hLandau","Landau.txt histogramm", 100, -10.,40.);
  TH1D* hMyGauss = new  TH1D("hMyGauss", "My Gauss generation", 100, -5.,5.);

  ifstream file_Landau("Praktikum_2025_proga/Landau.txt");
  if(!file_Landau.is_open())
  {
    cerr<< "File Landau.txt hasn't opened"<<endl;
    return;
  }

  
  const int STOP_MAX = 100000;
  int STOP_check =0;
  double value;
  while(file_Landau>>value && STOP_check<STOP_MAX)
  { 
    hLandau->Fill(value);
    hMyGauss->Fill(my_random_Gauss());
    STOP_check++;
  }
  file_Landau.close();

  TFile* file_output = new TFile("Praktikum_2025_proga/Output_task_2.root", "RECREATE");
  hLandau->Write();
  hMyGauss->Write();
  file_output->Close();
  */

  double TEST_ARRAY[4]={0.,1.,2.,0.};
  double a = stat_x_mean(TEST_ARRAY,4);
  double b = stat_s(TEST_ARRAY,4);
  double c =stat_gamma_1(TEST_ARRAY,4);
  double d = stat_gamma_2(TEST_ARRAY,4);

  int N_arr_size = 0;
  cout<<"Enter array's size:"<<endl;
  cin>>N_arr_size;

  double* Arr_Gauss = new double[N_arr_size]();
  double* Arr_0_1 = new double[N_arr_size]();

  fill_random_0_1(Arr_0_1,N_arr_size);
  fill_random_Gauss(Arr_Gauss,N_arr_size);

  cout<<Arr_Gauss[2]<<endl;

  delete[] Arr_Gauss;
  delete[] Arr_0_1;
}