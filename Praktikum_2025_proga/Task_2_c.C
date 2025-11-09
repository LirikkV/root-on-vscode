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

void Task_2_c()
{
  
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
  


  srand((unsigned int)time(NULL));
  cout << my_random_double_0_1() <<endl;
  cout<< my_random_Gauss()<<endl;

}