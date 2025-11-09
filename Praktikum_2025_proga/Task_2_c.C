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

  //First part:
  int N_arr_size = 0;
  cout<<"Enter array's size:"<<endl;
  cin>>N_arr_size;

  double* Arr_Gauss = new double[N_arr_size]();
  double* Arr_0_1 = new double[N_arr_size]();

  fill_random_0_1(Arr_0_1,N_arr_size);
  fill_random_Gauss(Arr_Gauss,N_arr_size);

  double x_mean_0_1 = stat_x_mean(Arr_0_1,N_arr_size);
  double s_0_1 = stat_s(Arr_0_1,N_arr_size);
  double gamma_1_0_1 =stat_gamma_1(Arr_0_1,N_arr_size);
  double gamma_2_0_1 = stat_gamma_2(Arr_0_1,N_arr_size);

  double x_mean_Gauss = stat_x_mean(Arr_Gauss,N_arr_size);
  double s_Gauss = stat_s(Arr_Gauss,N_arr_size);
  double gamma_1_Gauss =stat_gamma_1(Arr_Gauss,N_arr_size);
  double gamma_2_Gauss = stat_gamma_2(Arr_Gauss,N_arr_size);

  //output for [0,1]:
  printf("\n%-30s","Uniformly distribution on [0,1]:");
  printf("\n%-12s\t%-10s\t%-10s\t%-10s\t%-10s\n","", "x_mean", "s", "gamma 1", "gamma 2");
  printf("%-12s\t%-14.8f\t%-14.8f\t%-14.8f\t%-14.8f\n","experimental:", x_mean_0_1, s_0_1, gamma_1_0_1, gamma_2_0_1);
  printf("%-12s\t%-14.8f\t%-14.8f\t%-14.8f\t%-14.8f\n","theoretical:", 1/2., sqrt(1/12.), 0., -1.2);
  //output for Gauss:
  printf("\n%-30s","Gauss distribution:");
  printf("\n%-12s\t%-10s\t%-10s\t%-10s\t%-10s\n","", "x_mean", "s", "gamma 1", "gamma 2");
  printf("%-12s\t%-14.8f\t%-14.8f\t%-14.8f\t%-14.8f\n","experimental:", x_mean_Gauss, s_Gauss, gamma_1_Gauss, gamma_2_Gauss);
  printf("%-12s\t%-14.8f\t%-14.8f\t%-14.8f\t%-14.8f\n","theoretical:", 0., 1., 0., 0.);

  delete[] Arr_Gauss;
  delete[] Arr_0_1;

  //Second part:

  // TH1D* hLandau = new TH1D("hLandau","Landau.txt histogramm", 100, -10.,40.);
  // TH1D* hMyGauss = new  TH1D("hMyGauss", "My Gauss generation", 100, -5.,5.);

  // ifstream file_Landau("Praktikum_2025_proga/Landau.txt");
  // if(!file_Landau.is_open())
  // {
  //   cerr<< "File Landau.txt hasn't opened"<<endl;
  //   return;
  // }

  
  // const int STOP_MAX = 100000;
  // int STOP_check =0;
  // double value;
  // while(file_Landau>>value && STOP_check<STOP_MAX)
  // { 
  //   hLandau->Fill(value);
  //   hMyGauss->Fill(my_random_Gauss());
  //   STOP_check++;
  // }
  // file_Landau.close();

  // TFile* file_output = new TFile("Praktikum_2025_proga/Output_task_2.root", "RECREATE");
  // hLandau->Write();
  // hMyGauss->Write();

  // file_output->Close();
  
}