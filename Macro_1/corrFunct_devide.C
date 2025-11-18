#include "TFile.h"
#include "TH1D.h"

void corrFunct_devide()
{
TFile *f = TFile::Open("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1.root", "READ");

//for Pi+Pi+:
TH1D *hA_Plus = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL");
TH1D *hB_Plus = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL");
TH1D *hA_Plus_norm = (TH1D*)hA_Plus->Clone("hA_Plus_norm");
TH1D *hB_Plus_norm = (TH1D*)hB_Plus->Clone("hB_Plus_norm");
hA_Plus_norm->Scale(1.0/hA_Plus_norm->Integral());
hB_Plus_norm->Scale(1.0/hB_Plus_norm->Integral());
TH1D *CF_Pi_Plus_norm = (TH1D*)hA_Plus_norm->Clone("CF_Pi_Plus_nrom");
TH1D *CF_Pi_Plus = (TH1D*)hA_Plus->Clone("CF_Pi_Plus");
CF_Pi_Plus_norm->Divide(hB_Plus_norm);
CF_Pi_Plus_norm->SetTitle("Corr.Funct Pi+ Pi+ after normalization");
CF_Pi_Plus->Divide(hB_Plus);
CF_Pi_Plus->SetTitle("Corr.Funct Pi+ Pi+");


//for Pi-Pi-:
TH1D *hA_Minus = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL");
TH1D *hB_Minus = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL");
TH1D *hA_Minus_norm = (TH1D*)hA_Minus->Clone("hA_Minus_norm");
TH1D *hB_Minus_norm = (TH1D*)hB_Minus->Clone("hB_Minus_norm");
hA_Minus_norm->Scale(1.0/hA_Minus_norm->Integral());
hB_Minus_norm->Scale(1.0/hB_Minus_norm->Integral());
TH1D *CF_Pi_Minus_norm = (TH1D*)hA_Minus_norm->Clone("CF_Pi_Minus_norm");
TH1D *CF_Pi_Minus = (TH1D*)hA_Minus->Clone("CF_Pi_Minus");
CF_Pi_Minus_norm->Divide(hB_Minus_norm);
CF_Pi_Minus_norm->SetTitle("Corr.Funct Pi- Pi- after normalization");
CF_Pi_Minus->Divide(hB_Minus);
CF_Pi_Minus->SetTitle("Corr.Funct Pi- Pi-");


TFile *f_out = new TFile("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1_pr_1.root", "RECREATE");
CF_Pi_Plus_norm->Write();
CF_Pi_Minus_norm->Write();
CF_Pi_Plus->Write();
CF_Pi_Minus->Write();


f_out->Close();
f->Close();
}