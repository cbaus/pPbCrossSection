///HF plots E over eta

#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMarker.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>

void Show(TH1D* a,TH1D* b,TH1D* c, TH1D* epossl, TH1D* d, string type);
int BinToIeta(int bin);
int IetaToRing(int ieta);
int RingToIeta(int ring);

TVectorD hf_calibration(24);


void makePlots_perf_ring_corr_factors()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  for(int ring=0; ring<24; ring++)
    {
      hf_calibration[ring]=1;
    }
  cout << "\\hline Ieta & calibration factor\\\\\\hline\\hline" << endl; 

  TFile* file = TFile::Open("histos_nohfcalib.root");
  TFile* file2 = TFile::Open("histos_nohfcalib.root");
  TH1D* e= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_lev_m");
  TH1D* f= (TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_eta_lev_m");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_m");
  TH1D* epossl= (TH1D*)file->Get("Epos_SL/Epos_SL_h_perf_hf_totE_eta_lev_m");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_m");

  Show(e,f,g,epossl,h,"lev_minus");

  TFile* file = TFile::Open("histos_nohfcalib.root");
  TFile* file2 = TFile::Open("histos_nohfcalib.root");
  TH1D* e= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_lev_p");
  TH1D* f= (TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_eta_lev_p");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_p");
  TH1D* epossl= (TH1D*)file->Get("Epos_SL/Epos_SL_h_perf_hf_totE_eta_lev_p");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_p");

  Show(e,f,g,epossl,h,"lev_plus");
  
  hf_calibration.Print();
  TFile calibfile("plots/hf_calibration_data.root","recreate");
  hf_calibration.Write("hf_calibration");
  calibfile.Close();
}

void Show(TH1D* a, TH1D* b, TH1D* c, TH1D* epossl, TH1D* d, string type)
{
  //a->Scale(a->GetBinContent(2)/a->GetBinContent(2);
  int bin = 2;
  if (type=="lev_minus")
    bin = 12;
  b->Scale(a->GetBinContent(bin)/b->GetBinContent(bin));
  c->Scale(a->GetBinContent(bin)/c->GetBinContent(bin));
  epossl->Scale(a->GetBinContent(bin)/epossl->GetBinContent(bin));
  d->Scale(a->GetBinContent(bin)/d->GetBinContent(bin));

  a->SetMarkerSize(1.2);
  a->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  epossl->SetLineWidth(2.5);
  d->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  epossl->SetMarkerColor(kCyan);
  d->SetMarkerColor(kGreen+2);

  a->SetLineColor(kBlack);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  epossl->SetLineColor(kCyan);
  d->SetLineColor(kGreen+2);

  a->SetTitle("zero bias");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  epossl->SetTitle("EPOS (SL)");
  d->SetTitle("QGSJetII");

  double maximumy = 1.2 * TMath::Max(TMath::Max(a->GetMaximum(),b->GetMaximum()),TMath::Max(c->GetMaximum(),d->GetMaximum()));
  a->GetYaxis()->SetRangeUser(0,60);
  a->GetXaxis()->SetTitle("#eta");
  a->GetYaxis()->SetTitle("E [GeV]");

  TCanvas* c1 = new TCanvas;
  a->Draw("HIST P");
  b->Draw("HIST L SAME");
  c->Draw("HIST L SAME");
  epossl->Draw("HIST L SAME");
  d->Draw("HIST L SAME");
  TLegend* leg = new TLegend(0.3,0.7,0.7,0.9);
  leg->AddEntry(a,"","P");
  leg->AddEntry(b,"","L");
  leg->AddEntry(c,"","L");
  leg->AddEntry(epossl,"","L");
  leg->AddEntry(d,"","L");
  SetLegAtt(leg);
  leg->SetFillColor(kWhite);
  leg->Draw();


  for(int bin=2;bin<=12;bin++)
    {
      double avg = b->GetBinContent(bin) + c->GetBinContent(bin) + d->GetBinContent(bin);
      avg/=3.;
      
      int ieta = BinToIeta(type=="lev_minus"?(-bin):bin);
      
      double calib = avg/a->GetBinContent(bin);
      cout << ieta << " & " << setprecision(3) <<  calib << "\\\\" << endl;
      hf_calibration[IetaToRing(ieta)] = calib;
    }
  cout << "\\hline" << endl;
           
}
int BinToIeta(int bin)
///converts internal bin numbebin to original ieta
{
  if(bin < 0)
    return -40 - 2 - bin;
  else
    return bin - 12 + 40;
}

int IetaToRing(int ieta)
///converts ieta to numbers from 0 to 23
{
  //rings range from [-41,29] and [29,41]
  //29 is skipped in trees
  //41 is skipped usually because of HCAL prescription
  if(ieta < 0)
    return ieta + 41;
  else
    return ieta + 12 - 30;
}

int RingToIeta(int ring)
///converts rings from 0 to 23 to ieta
{
  //rings range from [-41,29] and [29,41]
  //29 is skipped in trees
  //41 is skipped usually because of HCAL prescription
  if(ring < 12)
    return ring - 41;
  else
    return ring - 12 + 30;
}
