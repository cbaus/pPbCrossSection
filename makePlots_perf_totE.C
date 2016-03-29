//total energy plot

#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"

// #ifndef __CINT__
// #include "style.h"
// #endif

#include <iostream>
#include <sstream>
using namespace std;
void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type);

void makePlots_perf_totE()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  // TFile* file = TFile::Open("histos.root");
  // TFile* file2 = TFile::Open("histos.root");
  // TH1D* data=(TH1D*)file->Get("data259163/data259163_h_perf_hf_totE_ZBSingleTrack");
  // TH1D* mc1=(TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_ZBSingleTrack");
  // TH1D* mc2=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_hf_totE_ZBSingleTrack");
  // TH1D* mc3=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_hf_totE_ZBSingleTrack");

  //Show(a,b,c,d,"SingleTrack");

  TFile* file = TFile::Open("histos.root");
  TH1D* data=(TH1D*)file->Get("data259163/data259163_h_perf_hf_totE_single_3gev");
  TH1D* mc1=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_hf_totE_single_3gev");
  TH1D* mc2=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_hf_totE_single_3gev");
  TH1D* mc3=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_hf_totE_single_3gev");

  Show(data,mc1,mc2,mc3,"single");

  TFile* file = TFile::Open("histos.root");
  TH1D* e=(TH1D*)file->Get("data259163/data259163_h_perf_hf_totE_double_1dot5gev");
  TH1D* f=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_hf_totE_double_1dot5gev");
  TH1D* g=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_hf_totE_double_1dot5gev");
  TH1D* h=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_hf_totE_double_1dot5gev");

  Show(data,mc1,mc2,mc3,"double");
}

void Show(TH1D* data,TH1D* mc1,TH1D* mc2,TH1D* mc3, string type)
{
  const int normbin = data->FindBin(0.5);
  const int normendbin = data->GetNbinsX();
  data->Scale(1./double(data->Integral()));
  if(mc1) mc1->Scale(data->Integral(normbin,normendbin)/mc1->Integral(normbin,normendbin));
  if(mc2) mc2->Scale(data->Integral(normbin,normendbin)/mc2->Integral(normbin,normendbin));
  if(mc3) mc3->Scale(data->Integral(normbin,normendbin)/mc3->Integral(normbin,normendbin));

  data->SetMarkerSize(1.2);
  data->SetLineWidth(2.5);
  if(mc1) mc1->SetLineWidth(2.5);
  if(mc2) mc2->SetLineWidth(2.5);
  if(mc3) mc3->SetLineWidth(2.5);


  data->SetMarkerColor(kBlack);
  if(mc1) mc1->SetMarkerColor(kRed);
  if(mc2) mc2->SetMarkerColor(kBlue);
  if(mc3) mc3->SetMarkerColor(kGreen+2);

  data->SetLineColor(kBlack);
  if(mc1) mc1->SetLineColor(kRed);
  if(mc2) mc2->SetLineColor(kBlue);
  if(mc3) mc3->SetLineColor(kGreen+2);

  data->SetTitle("unbiased trigger");
  if(mc1) mc1->SetTitle("Pythia6 Z2*");
  if(mc2) mc2->SetTitle("Pythia8 Monash Tune");
  if(mc3) mc3->SetTitle("Pythia8 MBR Tune");

  data->GetXaxis()->SetRangeUser(0,3.5);
  data->GetYaxis()->SetRangeUser(1e-6,1e2);
  data->GetXaxis()->SetNdivisions(505);
  data->GetXaxis()->SetTitle("_{#sumE_{towers} [TeV]}");
  data->GetYaxis()->SetTitle("events (normalised)");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  TCanvas* c1 = new TCanvas;
  data->Draw("");
  if(mc1) mc1->Draw("HIST SAME");
  if(mc2) mc2->Draw("HIST SAME");
  if(mc3) mc3->Draw("HIST SAME");
  TLegend* leg = c1->BuildLegend(0.22,0.72,0.62,0.92);
  SetLegAtt(leg);
  leg->Draw();

  CMSText(3,0,1,type=="single"?"Single-arm selection":"Double-arm selection");

  c1->SetLogy();
  c1->SaveAs((string("plots/hf_perf_totE_")+type+string(".pdf")).c_str());
  c1->SaveAs((string("plots/hf_perf_totE_")+type+string(".eps")).c_str());
}
