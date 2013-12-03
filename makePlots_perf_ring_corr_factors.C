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
  //FROM LEV email 22.11

  vector<double> c_lev_m,c_lev_m_e,c_lev_p,c_lev_p_e;
  //TGraphErrors* gr_lev_p = new TGraphErrors(c_lev_m.size(),&c_lev_m.front());
  c_lev_m.push_back(1.07); //-41 //12
  c_lev_m.push_back(0.97);
  c_lev_m.push_back(0.95);
  c_lev_m.push_back(0.99);
  c_lev_m.push_back(0.96);
  c_lev_m.push_back(0.91);
  c_lev_m.push_back(0.92);
  c_lev_m.push_back(0.86);
  c_lev_m.push_back(0.80);
  c_lev_m.push_back(0.72);
  c_lev_m.push_back(0.69);
  c_lev_m.push_back(0.83);
  c_lev_m.push_back(0.73); //-29 //0

  c_lev_p.push_back(1.01);
  c_lev_p.push_back(0.94);
  c_lev_p.push_back(0.91);
  c_lev_p.push_back(0.89);
  c_lev_p.push_back(0.87);
  c_lev_p.push_back(0.92);
  c_lev_p.push_back(0.84);
  c_lev_p.push_back(0.85);
  c_lev_p.push_back(0.83);
  c_lev_p.push_back(0.67);
  c_lev_p.push_back(0.61);
  c_lev_p.push_back(0.75);
  c_lev_p.push_back(0.66);

  c_lev_m_e.push_back(0.23);
  c_lev_m_e.push_back(0.04);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.04);
  c_lev_m_e.push_back(0.06);
  c_lev_m_e.push_back(0.10);
  c_lev_m_e.push_back(0.09);

  c_lev_p_e.push_back(0.21);
  c_lev_p_e.push_back(0.04);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.07);
  c_lev_p_e.push_back(0.09);
  c_lev_p_e.push_back(0.08);

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

  //Correction Factor Histograms
  TH1D* h_c_avg = a->Clone("h_c_avg");
  TH1D* h_c_lev = a->Clone("h_c_lev");
  h_c_avg->GetYaxis()->SetRangeUser(0,1.2);
  h_c_avg->GetYaxis()->SetTitle("C");
  h_c_lev->SetMarkerStyle(21);
  h_c_avg->SetMarkerSize(1.5);
  h_c_lev->SetMarkerSize(1.5);
  h_c_avg->SetMarkerColor(kRed);
  h_c_lev->SetMarkerColor(kBlack);


  for(int bin=2;bin<=12;bin++)
    {
      double avg = b->GetBinContent(bin) + c->GetBinContent(bin) + d->GetBinContent(bin);
      avg/=3.;
      double calib = avg/a->GetBinContent(bin);
      double c_lev = type=="lev_minus"?c_lev_m[12-(bin-1)]:c_lev_p[bin-1];
      double c_lev_e = type=="lev_minus"?c_lev_m_e[12-(bin-1)]:c_lev_p_e[bin-1];
      
      h_c_avg->SetBinContent(bin,1./calib);
      h_c_avg->SetBinError(bin,0);
      h_c_lev->SetBinContent(bin,c_lev);
      h_c_lev->SetBinError(bin,c_lev_e);
      
      int ieta = BinToIeta(type=="lev_minus"?(-bin):bin);
      
      cout << ieta << " & " << fixed << setprecision(2) <<  1./calib << " & " << c_lev << "$\\pm$" << c_lev_e << "\\\\" << endl;
      hf_calibration[IetaToRing(ieta)] = calib;
    }
  cout << "\\hline" << endl;


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

  TCanvas* c2 = new TCanvas;
  h_c_avg->Draw("P");
  h_c_lev->Draw("P SAME");
  

  leg = new TLegend(0.2,0.3,0.6,0.4);
  leg->AddEntry(h_c_lev,"scale factor (pp 2013/2011 energy flow)","PL");
  leg->AddEntry(h_c_avg,"1/scale factor (pPb MC)","P");
  SetLegAtt(leg);
  leg->Draw();

  CMSText(1,0,0);

  c2->SaveAs((string("plots/hf_perf_corr_vs_")+type+string(".pdf")).c_str());
           
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
