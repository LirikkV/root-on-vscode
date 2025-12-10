#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

void corrFunct_devide()
{
TFile *f = TFile::Open("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1.root", "READ");

//for Pi+Pi+:
TH1D *hA_Plus = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL");
TH1D *hB_Plus = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL");
TH1D *CF_Pi_Plus = (TH1D*)hA_Plus->Clone("CF_Pi_Plus");
CF_Pi_Plus->Divide(hB_Plus);
CF_Pi_Plus->SetTitle("Corr.Funct Pi+ Pi+");
CF_Pi_Plus->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus->GetYaxis()->SetTitle("CF");

//for Pi-Pi-:
TH1D *hA_Minus = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL");
TH1D *hB_Minus = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL");
TH1D *CF_Pi_Minus = (TH1D*)hA_Minus->Clone("CF_Pi_Minus");
CF_Pi_Minus->Divide(hB_Minus);
CF_Pi_Minus->SetTitle("Corr.Funct Pi- Pi-");
CF_Pi_Minus->GetXaxis()->SetTitle("q_inv");
CF_Pi_Minus->GetYaxis()->SetTitle("CF");

TFile *f_out = new TFile("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1_pr.root", "RECREATE");
CF_Pi_Plus->Write();
CF_Pi_Minus->Write();

TCanvas *c1 = new TCanvas("c1", "Canvas",800,600);
hA_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hA_Plus.pdf");
hB_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hB_Plus.pdf");
CF_Pi_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hCF_Pi_Plus.pdf");



// CF_Pi_Plus->GetXaxis()->SetRange
// CF_Pi_Plus->Draw();
// c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hCF_Pi_Plus.pdf");

f_out->Close();
f->Close();
}