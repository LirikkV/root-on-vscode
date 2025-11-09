#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TFile.h"

using namespace ::std;

void Task_2_c()
{
  
  TH1D* hLandau = new TH1D("hLandau","Landau.txt histogramm", 100, -10.,40.);

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
    STOP_check++;
  }
  file_Landau.close();

  TFile* file_output = new TFile("Praktikum_2025_proga/Output_task_2.root", "RECREATE");
  hLandau->Write();
  file_output->Close();

}