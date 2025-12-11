#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "iostream"

void corrFunct_devide()
{
TFile *f = TFile::Open("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1.root", "READ");

//for Pi+Pi+:
TH1D *hA_Plus = (TH1D*)f->Get("hA_Pi_Plus_q_inv_ALL");
TH1D *hB_Plus = (TH1D*)f->Get("hB_Pi_Plus_q_inv_ALL");

//for Pi-Pi-:
TH1D *hA_Minus = (TH1D*)f->Get("hA_Pi_Minus_q_inv_ALL");
TH1D *hB_Minus = (TH1D*)f->Get("hB_Pi_Minus_q_inv_ALL");


//let's normalize CF:
Double_t q_inv_min = 0.2;
Double_t q_inv_max = 0.9;

Double_t A_Pi_Plus_integral = hA_Plus->Integral(hA_Plus->FindBin(q_inv_min), hA_Plus->FindBin(q_inv_max));
Double_t B_Pi_Plus_integral = hB_Plus->Integral(hB_Plus->FindBin(q_inv_min),hB_Plus->FindBin(q_inv_max));
Double_t Scale_factor_Pi_Plus = B_Pi_Plus_integral/A_Pi_Plus_integral;

Double_t A_Pi_Minus_integral = hA_Minus->Integral(hA_Minus->FindBin(q_inv_min), hA_Minus->FindBin(q_inv_max));
Double_t B_Pi_Minus_integral = hB_Minus->Integral(hB_Minus->FindBin(q_inv_min),hB_Minus->FindBin(q_inv_max));
Double_t Scale_factor_Pi_Minus = B_Pi_Minus_integral/A_Pi_Minus_integral;

//For hists scales:
Double_t x1_range = 0.0;
Double_t x2_range = 0.2;
Double_t y1_range = 0.4;
Double_t y2_range = 1.6;
/*
Method #1: 
1)Get relation of A/B in range of q_inv where no correlations 
2)Scale A by relation from 1)
3)CF=A_normilized/B
*/
TH1D *hA_Plus_normalized = (TH1D*)hA_Plus->Clone("CF_Pi_Plus");
hA_Plus_normalized->Scale(Scale_factor_Pi_Plus);
TH1D *CF_Pi_Plus_Meth_1 = (TH1D*)hA_Plus_normalized->Clone("CF_Pi_Plus_Meth_1");
CF_Pi_Plus_Meth_1->Divide(hB_Plus);

CF_Pi_Plus_Meth_1->SetTitle("Corr.Funct Pi+ Pi+ scale A then A/B");
CF_Pi_Plus_Meth_1->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus_Meth_1->GetYaxis()->SetTitle("CF");
CF_Pi_Plus_Meth_1->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Plus_Meth_1->GetXaxis()->SetRangeUser(x1_range,x2_range);


TH1D *hA_Minus_normalized = (TH1D*)hA_Minus->Clone("CF_Pi_Minus");
hA_Minus_normalized->Scale(Scale_factor_Pi_Minus);
TH1D *CF_Pi_Minus_Meth_1 = (TH1D*)hA_Minus_normalized->Clone("CF_Pi_Minus_Meth_1");
CF_Pi_Minus_Meth_1->Divide(hB_Minus);

CF_Pi_Minus_Meth_1->SetTitle("Corr.Funct Pi- Pi- scale A then A/B");
CF_Pi_Minus_Meth_1->GetXaxis()->SetTitle("q_inv");
CF_Pi_Minus_Meth_1->GetYaxis()->SetTitle("CF");
CF_Pi_Minus_Meth_1->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Minus_Meth_1->GetXaxis()->SetRangeUser(x1_range,x2_range);

/*
Method #2: 
1)Get relation of A/B in range of q_inv where no correlations (same as in Method #1)
2)CF_non_normalized = A/B
3)Scale CF_non_normalized by relation from #1
*/
TH1D *CF_Pi_Plus_Meth_2 = (TH1D*)hA_Plus->Clone("CF_Pi_Plus_Meth_2");
CF_Pi_Plus_Meth_2->Divide(hB_Plus);
CF_Pi_Plus_Meth_2->Scale(Scale_factor_Pi_Plus);

CF_Pi_Plus_Meth_2->SetTitle("Corr.Funct Pi+ Pi+ A/B then scale CF");
CF_Pi_Plus_Meth_2->GetXaxis()->SetTitle("q_inv");
CF_Pi_Plus_Meth_2->GetYaxis()->SetTitle("CF");
CF_Pi_Plus_Meth_2->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Plus_Meth_2->GetXaxis()->SetRangeUser(x1_range,x2_range);


TH1D *CF_Pi_Minus_Meth_2 = (TH1D*)hA_Minus->Clone("CF_Pi_Minus_Meth_2");
CF_Pi_Minus_Meth_2->Divide(hB_Minus);
CF_Pi_Minus_Meth_2->Scale(Scale_factor_Pi_Minus);

CF_Pi_Minus_Meth_2->SetTitle("Corr.Funct Pi+ Pi+ A/B then scale CF");
CF_Pi_Minus_Meth_2->GetXaxis()->SetTitle("q_inv");
CF_Pi_Minus_Meth_2->GetYaxis()->SetTitle("CF");
CF_Pi_Minus_Meth_2->GetYaxis()->SetRangeUser(y1_range,y2_range);
CF_Pi_Minus_Meth_2->GetXaxis()->SetRangeUser(x1_range,x2_range);



TFile *f_out = new TFile("/home/lirikk/root-on-vscode/Output_Data_1/Output_file_1_pr.root", "RECREATE");
CF_Pi_Plus_Meth_1->Write();
CF_Pi_Minus_Meth_1->Write();
CF_Pi_Plus_Meth_2->Write();
CF_Pi_Minus_Meth_2->Write();

TCanvas *c1 = new TCanvas("c1", "Canvas",800,600);
hA_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hA_Plus.pdf");
hB_Plus->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hB_Plus.pdf");
CF_Pi_Plus_Meth_1->Draw();
c1->SaveAs("/home/lirikk/root-on-vscode/Output_Data_1/hCF_Pi_Plus.pdf");

f_out->Close();
f->Close();
}