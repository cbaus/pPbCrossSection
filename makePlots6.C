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

#ifndef __CINT__
#include "style.h"
#endif

#include <iostream>
#include <sstream>
using namespace std;
void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type);

void makePlots6()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=(TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_single_3gev");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_single_3gev");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_single_3gev");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_single_3gev");

  Show(a,b,c,d,"single");

  TFile* file = TFile::Open("histos_old.root");
  TH1D* e=(TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_double_1dot5gev");
  TH1D* f=(TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_double_1dot5gev");
  TH1D* g=(TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_double_1dot5gev");
  TH1D* h=(TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_double_1dot5gev");

  Show(e,f,g,h,"double");
}

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type)
{
  const int normbin1 = a->FindBin(1000);
  const int normbin2 = a->FindBin(1500);
  a->Scale(1./double(a->GetEntries()));
  b->Scale(a->Integral(normbin1,normbin2)/b->Integral(normbin1,normbin2));
  c->Scale(a->Integral(normbin1,normbin2)/c->Integral(normbin1,normbin2));
  d->Scale(a->Integral(normbin1,normbin2)/d->Integral(normbin1,normbin2));

  a->SetMarkerSize(1.1);
  a->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  d->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);

  a->SetLineColor(kBlack);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);

  a->SetTitle("Data");
  b->SetTitle("HIJING 1.383");
  c->SetTitle("EPOS-LHC");
  d->SetTitle("QGSJetII-4");

  a->GetXaxis()->SetRange(2,250);
  a->GetYaxis()->SetRangeUser(0.5,5000);
  a->GetXaxis()->SetTitle("_{#sumE_{HF} [GeV]}");
  a->GetYaxis()->SetTitle("events (normalised)");

  a->GetXaxis()->SetNdivisions(505);

  TCanvas* c1 = new TCanvas;
  a->Draw("");
  b->Draw("HIST SAME");
  c->Draw("HIST SAME");
  d->Draw("HIST SAME");
  TLegend* leg = new TLegend(0.1,0.1,0.2,0.2);
  leg->SetX1(0.55);
  leg->SetX2(0.85);
  leg->SetY1(0.67);
  leg->SetY2(0.87);
  leg->AddEntry(a,"","p");
  leg->AddEntry(b,"","l");
  leg->AddEntry(c,"","l");
  leg->AddEntry(d,"","l");
  SetLegAtt(leg);
  DataText(1,1);
  leg->Draw();


  TPaveText* text = new TPaveText(0.23,0.23,0.58,0.28,"NDC b t l");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->AddText(type=="single"?"E > 6 GeV (single-arm)":"E > 2 GeV (double-arm)");
  text->SetTextSize(0.033);
  text->Draw();

  c1->SetLogy();
  c1->SaveAs((string("plots/hf_perf_3_")+type+string(".pdf")).c_str());
}
