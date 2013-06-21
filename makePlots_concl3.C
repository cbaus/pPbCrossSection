#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TVectorD.h>

#include<iostream>

using namespace std;

void SetAttributes(TH1D* theGraph, int colour, int marker)
{
  theGraph->SetMarkerSize(1.5);
  theGraph->SetLineWidth(1.8);
  theGraph->SetFillColor(kWhite);
  theGraph->SetLineColor(colour);
  theGraph->SetMarkerColor(colour);
  theGraph->SetMarkerStyle(marker);
}

void makePlots_concl3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  gStyle->SetPadTopMargin(0.07);

  ///READ IN VALUES
  TFile f("plots/final_values.root");
  TVectorD* vec_sigma_vis    = NULL;
  TVectorD* vec_sigma_vis_e  = NULL;
  TVectorD* vec_sigma_had    = NULL;
  TVectorD* vec_sigma_had_e  = NULL;
  TVectorD* vec_sigma_inel   = NULL;
  TVectorD* vec_sigma_inel_e = NULL;
  vec_sigma_vis    = (TVectorD*)f.Get("vec_sigma_vis");
  vec_sigma_vis_e  = (TVectorD*)f.Get("vec_sigma_vis_e");
  vec_sigma_had    = (TVectorD*)f.Get("vec_sigma_had");
  vec_sigma_had_e  = (TVectorD*)f.Get("vec_sigma_had_e");
  vec_sigma_inel   = (TVectorD*)f.Get("vec_sigma_inel");
  vec_sigma_inel_e = (TVectorD*)f.Get("vec_sigma_inel_e");
  cout << "Inel - Had - Vis:" << endl;
  vec_sigma_inel->Print();
  vec_sigma_had->Print();
  vec_sigma_vis->Print();
  f.Close();
  ///!READ IN VALUES

  TH1D* h_data = new TH1D("h_data","CMS;cross section;#sigma [b]",6,-0.5,5.5);
  h_data->GetXaxis()->SetBinLabel(1,"inelastic");
  h_data->GetXaxis()->SetBinLabel(2,"hadronic");
  h_data->GetXaxis()->SetBinLabel(3,"visible");
  h_data->GetXaxis()->SetBinLabel(4,"inelastic");
  h_data->GetXaxis()->SetBinLabel(5,"hadronic");
  h_data->GetXaxis()->SetBinLabel(6,"visible");
  TH1D * h_epos = new TH1D("h_epos","EPOS-LHC;;#sigma [b]",6,-0.5,5.5);
  h_epos->GetXaxis()->SetBinLabel(1,"inelastic");
  h_epos->GetXaxis()->SetBinLabel(2,"hadronic");
  h_epos->GetXaxis()->SetBinLabel(3,"visible");
  h_epos->GetXaxis()->SetBinLabel(4,"inelastic");
  h_epos->GetXaxis()->SetBinLabel(5,"hadronic");
  h_epos->GetXaxis()->SetBinLabel(6,"visible");
  TH1D * h_qgsjet = new TH1D("h_qgsjet","QGSJETII-04;;#sigma [b]",6,-0.5,5.5);
  h_qgsjet->GetXaxis()->SetBinLabel(1,"inelastic");
  h_qgsjet->GetXaxis()->SetBinLabel(2,"hadronic");
  h_qgsjet->GetXaxis()->SetBinLabel(3,"visible");
  h_qgsjet->GetXaxis()->SetBinLabel(4,"inelastic");
  h_qgsjet->GetXaxis()->SetBinLabel(5,"hadronic");
  h_qgsjet->GetXaxis()->SetBinLabel(6,"visible");

  h_data->SetBinContent(1,(*vec_sigma_inel)[1]);
  h_data->SetBinContent(2,(*vec_sigma_had)[0]);
  h_data->SetBinContent(3,(*vec_sigma_vis)[0]);
  h_data->SetBinContent(4,(*vec_sigma_inel)[2]);
  h_data->SetBinContent(5,(*vec_sigma_had)[1]);
  h_data->SetBinContent(6,(*vec_sigma_vis)[1]);
  h_data->SetBinError(1,(*vec_sigma_inel_e)[1]);
  h_data->SetBinError(2,(*vec_sigma_had_e)[0]);
  h_data->SetBinError(3,(*vec_sigma_vis_e)[0]);
  h_data->SetBinError(4,(*vec_sigma_inel_e)[2]);
  h_data->SetBinError(5,(*vec_sigma_had_e)[1]);
  h_data->SetBinError(6,(*vec_sigma_vis_e)[1]);

  double eff_epos_single = 0.958;
  double eff_epos_double = 0.929;
  double eff_qgsjet_single = 0.959;
  double eff_qgsjet_double = 0.937;

  h_epos->SetBinContent(1,2.085703);
  h_epos->SetBinContent(2,2.085703*eff_epos_single);
  h_epos->SetBinContent(4,2.085703);
  h_epos->SetBinContent(5,2.085703*eff_epos_double);
  h_qgsjet->SetBinContent(1,2.176422);
  h_qgsjet->SetBinContent(2,2.176422*eff_qgsjet_single);
  h_qgsjet->SetBinContent(4,2.176422);
  h_qgsjet->SetBinContent(5,2.176422*eff_qgsjet_double);

  SetAttributes(h_data,kRed,21);
  SetAttributes(h_epos,kGreen-2,22);
  SetAttributes(h_qgsjet,kBlue,34);

  TCanvas* can1 = new TCanvas;
  h_data->Draw("P");
  h_epos->Draw("SAME P");
  h_qgsjet->Draw("SAME P");

  double min = 1.85;
  double max = 2.8;

  h_data->GetYaxis()->SetRangeUser(min,max);

  TLine* line = new TLine(2.5,min,2.5,max);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend* leg1 = new TLegend(0.26,0.7,0.53,0.9);
  SetLegAtt(leg1);
  leg1->AddEntry(h_data,"CMS","p");
  leg1->AddEntry(h_epos,"EPOS-LHC","p");
  leg1->AddEntry(h_qgsjet,"QGSJetII-04","p");
  leg1->SetBorderSize(1);
  leg1->Draw();


  TPaveText* single_arm = new TPaveText(0.34,0.6,0.44,0.65,"NDC b t l");
  single_arm->SetBorderSize(0);
  single_arm->SetFillStyle(0);
  single_arm->AddText("single-arm");
  single_arm->SetTextSize(0.04);
  single_arm->Draw();

  TPaveText* double_arm = new TPaveText(0.72,0.6,0.82,0.65,"NDC b t l");
  double_arm->SetBorderSize(0);
  double_arm->SetFillStyle(0);
  double_arm->AddText("double-arm");
  double_arm->SetTextSize(0.04);
  double_arm->Draw();

  DataText(0,1);

  can1->SaveAs((string("plots/concl_3")+string(".pdf")).c_str());

}



