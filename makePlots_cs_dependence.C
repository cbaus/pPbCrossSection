//checks dependence on threshold energy

#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMarker.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"

#ifndef __CINT__
#include "style.h"
#include "makePlots_cs_eff.C" //yes not very nice
#include "makePlots_cs.C"
#endif

#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

using namespace std;

#define modfactor 1.0

void makePlots_cs_dependence(string filename = "histos.root")
{
  TH1::SetDefaultSumw2();
  gROOT->ProcessLine(" .L style.cc+");
  gROOT->ProcessLine(" .L makePlots_cs.C+");
#ifdef __CINT__
  style();
#endif

  TFile* file = TFile::Open(filename.c_str());

  vector<double> cut_energies_single;
  vector<double> cut_energies_double;
  for(double x=1; x<=10; x+=0.5)
    {
      cut_energies_single.push_back(x*2);
      cut_energies_double.push_back(x);
    }
  const int array_size = int(cut_energies_single.size());

  TGraphErrors* g_crosscheck_single = new TGraphErrors(array_size);
  TGraphErrors* g_crosscheck_double = new TGraphErrors(array_size);
  TGraphErrors* g_crosscheck_had_single = new TGraphErrors(array_size);
  TGraphErrors* g_crosscheck_had_double = new TGraphErrors(array_size);
  g_crosscheck_single->SetMarkerColor(kBlack);
  g_crosscheck_double->SetMarkerColor(kBlack);
  g_crosscheck_single->SetMarkerStyle(21);
  g_crosscheck_double->SetMarkerStyle(21);
  g_crosscheck_had_single->SetMarkerColor(g_crosscheck_single->GetMarkerColor());
  g_crosscheck_had_double->SetMarkerColor(g_crosscheck_double->GetMarkerColor());
  g_crosscheck_had_single->SetMarkerStyle(25);
  g_crosscheck_had_double->SetMarkerStyle(25);
  g_crosscheck_single->SetTitle("single-arm;selection threshold [GeV];#sigma [b]");
  g_crosscheck_double->SetTitle("double-arm;selection threshold [GeV];#sigma [b]");

  for(int i=0; i<array_size; i++)
    {
      double cut_single = cut_energies_single[i];
      double cut_double = cut_energies_double[i];

      cout << "Cross checking round " << i << "  single:" << cut_single << " double:" << cut_double << endl << endl;
      bool batch_status = gROOT->IsBatch();
      gROOT->SetBatch(kTRUE);
      makePlots_cs(0,cut_single,cut_double,modfactor,filename);
      gROOT->SetBatch(batch_status);

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

      g_crosscheck_single->SetPoint(i,cut_single,(*vec_sigma_inel)[1]);
      g_crosscheck_double->SetPoint(i,cut_double,(*vec_sigma_inel)[2]); //inel has three
      g_crosscheck_had_single->SetPoint(i,cut_single,(*vec_sigma_had)[0]); //vis has two elements
      g_crosscheck_had_double->SetPoint(i,cut_double,(*vec_sigma_had)[1]);

      g_crosscheck_single->SetPointError(i,0,(*vec_sigma_inel_e)[1]);
      g_crosscheck_double->SetPointError(i,0,(*vec_sigma_inel_e)[2]); //inel has three
      g_crosscheck_had_single->SetPointError(i,0,(*vec_sigma_had_e)[0]); //vis has two elements
      g_crosscheck_had_double->SetPointError(i,0,(*vec_sigma_had_e)[1]);
    }

  //Draw the modified monte carlo lines
  TCanvas* can1 = new TCanvas;
  TH1D* eposin=(TH1D*)file->Get("Epos/Epos_h_hf_cut_single");
  TH1D* eposdiff2in=(TH1D*)file->Get("DPMJet/DPMJet_h_hf_cut_single");
  TH1D* qgsjetin=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_cut_single");
  //Normalise to standard cut values. The lines all cross at this point
  eposin->Scale(g_crosscheck_had_single->Eval(8.)/eposin->Interpolate(8./modfactor));
  eposdiff2in->Scale(g_crosscheck_had_single->Eval(8.)/eposdiff2in->Interpolate(8./modfactor));
  qgsjetin->Scale(g_crosscheck_had_single->Eval(8.)/qgsjetin->Interpolate(8./modfactor));
  TGraph* epos = new TGraph(eposin->GetNbinsX());
  TGraph* eposdiff2 = new TGraph(eposin->GetNbinsX());
  TGraph* qgsjet = new TGraph(qgsjetin->GetNbinsX());
  for(int j=0; j<eposin->GetNbinsX(); j++)
    {
      int bin=j+1;
      epos->SetPoint(j,eposin->GetBinCenter(bin)*modfactor,eposin->GetBinContent(bin));
      eposdiff2->SetPoint(j,eposdiff2in->GetBinCenter(bin)*modfactor,eposdiff2in->GetBinContent(bin));
      qgsjet->SetPoint(j,qgsjetin->GetBinCenter(bin)*modfactor,qgsjetin->GetBinContent(bin));
    }
  epos->SetTitle("EPOS-LHC");
  eposdiff2->SetTitle("DPMJET3.06");
  qgsjet->SetTitle("QGSJetII-04");
  qgsjet->SetLineWidth(3);
  eposdiff2->SetLineWidth(3);
  epos->SetLineWidth(3);
  epos->SetLineColor(kBlue);
  eposdiff2->SetLineColor(kBlue);
  qgsjet->SetLineColor(kRed);
  eposdiff2->SetLineStyle(6);
  qgsjet->SetLineStyle(9);

  g_crosscheck_single->GetYaxis()->SetRangeUser(1.5,2.5);
  g_crosscheck_single->GetXaxis()->SetRangeUser(0,20);
  g_crosscheck_single->SetFillColor(kCyan);
  g_crosscheck_single->SetFillStyle(3002);
  g_crosscheck_single->Draw("A3");
  g_crosscheck_single->Draw("PX");
  g_crosscheck_had_single->Draw("P");
  TLegend* leg1 = new TLegend(0.55,0.71,0.82,0.92);
  leg1->AddEntry(g_crosscheck_single,"inelastic","pf");
  leg1->AddEntry(g_crosscheck_had_single,"hadronic","p");
  leg1->AddEntry(epos,epos->GetTitle(),"l");
  leg1->AddEntry(eposdiff2,eposdiff2->GetTitle(),"l");
  leg1->AddEntry(qgsjet,qgsjet->GetTitle(),"l");
  SetLegAtt(leg1,1.2);
  leg1->Draw();
  epos->Draw("L SAME");
  eposdiff2->Draw("L SAME");
  qgsjet->Draw("L SAME");
  CMSText(3,1,1,"single-arm selection");
  can1->SaveAs((string("plots/cross_check_1_single")+string(".pdf")).c_str());

  TCanvas* can2 = new TCanvas;
  eposin=(TH1D*)file->Get("Epos/Epos_h_hf_cut_double");
  eposdiff2in=(TH1D*)file->Get("DPMJet/DPMJet_h_hf_cut_double");
  qgsjetin=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_cut_double");
  eposin->Scale(g_crosscheck_had_double->Eval(4)/eposin->Interpolate(4/modfactor));
  eposdiff2in->Scale(g_crosscheck_had_double->Eval(4)/eposdiff2in->Interpolate(4/modfactor));
  qgsjetin->Scale(g_crosscheck_had_double->Eval(4)/qgsjetin->Interpolate(4/modfactor));
  TGraph* epos = new TGraph(eposin->GetNbinsX());
  TGraph* eposdiff2 = new TGraph(eposdiff2in->GetNbinsX());
  TGraph* qgsjet = new TGraph(qgsjetin->GetNbinsX());
  for(int j=0; j<eposin->GetNbinsX(); j++)
    {
      int bin=j+1;
      epos->SetPoint(j,eposin->GetBinCenter(bin)*modfactor,eposin->GetBinContent(bin));
      eposdiff2->SetPoint(j,eposdiff2in->GetBinCenter(bin)*modfactor,eposdiff2in->GetBinContent(bin));
      qgsjet->SetPoint(j,qgsjetin->GetBinCenter(bin)*modfactor,qgsjetin->GetBinContent(bin));
    }
  epos->SetTitle("EPOS-LHC");
  eposdiff2->SetTitle("DPMJET3.06");
  qgsjet->SetTitle("QGSJetII-04");
  qgsjet->SetLineWidth(3);
  eposdiff2->SetLineWidth(3);
  epos->SetLineWidth(3);
  epos->SetLineColor(kBlue);
  eposdiff2->SetLineColor(kBlue);
  qgsjet->SetLineColor(kRed);
  eposdiff2->SetLineStyle(6);
  qgsjet->SetLineStyle(9);

  g_crosscheck_double->GetYaxis()->SetRangeUser(1.5,2.5);
  g_crosscheck_double->SetFillColor(kCyan);
  g_crosscheck_double->SetFillStyle(3002);
  g_crosscheck_double->Draw("A3");
  g_crosscheck_double->Draw("PX");
  g_crosscheck_had_double->Draw("P");
  TLegend* leg2 = new TLegend(0.55,0.71,0.82,0.92);
  leg2->AddEntry(g_crosscheck_double,"inelastic","pf");
  leg2->AddEntry(g_crosscheck_had_double,"hadronic","p");
  leg2->AddEntry(epos,epos->GetTitle(),"l");
  leg2->AddEntry(eposdiff2,eposdiff2->GetTitle(),"l");
  leg2->AddEntry(qgsjet,qgsjet->GetTitle(),"l");
  SetLegAtt(leg2,1.2);
  leg2->Draw();
  epos->Draw("L SAME");
  eposdiff2->Draw("L SAME");
  qgsjet->Draw("L SAME");
  CMSText(3,1,1,"double-arm selection");
  can2->SaveAs((string("plots/cross_check_1_double")+string(".pdf")).c_str());
}
