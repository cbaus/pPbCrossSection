///makes plots for largest rapidity gap

#include <TColor.h>
#include <iomanip>
#include "TEfficiency.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TMarker.h"
#include "TPad.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"

#ifndef __CINT__
#include "style.h"
#endif


#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

void makePlots_pt_cuts()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<string> list;
  list.push_back(string("Epos"));
  //list.push_back(string("QGSJetII"));
  //  list.push_back(string("Hijing"));
  //list.push_back(string("DPMJet"));
  //  list.push_back(string("DPMJet_pp"));
//   list.push_back(string("EposDiffWeight150"));
//   list.push_back(string("EposDiffWeight200"));
//   list.push_back(string("EposDiffWeight299"));
//   list.push_back(string("QGSJetIIDiffWeight150"));
//   list.push_back(string("QGSJetIIDiffWeight200"));
//   list.push_back(string("QGSJetIIDiffWeight452"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  //name.push_back(string("QGSJetII-04"));
  //  name.push_back(string("HIJING 1.383"));
  // name.push_back(string("DPMJet 3.06"));
  //  name.push_back(string("DPMJet 3.04 (pp@5.02TeV)"));
//   name.push_back(string("Epos #sigma_{diff}x2"));
//   name.push_back(string("Epos #sigma_{diff}x2.4"));
//   name.push_back(string("Epos #sigma_{diff}x2.4"));
//   name.push_back(string("Epos #sigma_{diff}x2.4"));
//   name.push_back(string("Epos #sigma_{diff}x2.4"));
//   name.push_back(string("Epos #sigma_{diff}x2.4"));


  TGraph graphFindDiffWeightEpos(4);
  TGraph graphFindDiffWeightQgsjet(4);

  for(int i=0; i<int(list.size()); i++)
    {
      cout << i+1 << "/" << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");
      file->cd();

      TH1D* hsingle_all=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("h_mc_pt_single")).c_str());
      TH1D* hdouble_all=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("h_mc_pt_double")).c_str());
      TH1D* hsingle_sel=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("h_mc_pt_single_sel")).c_str());
      TH1D* hdouble_sel=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("h_mc_pt_double_sel")).c_str());

      TEfficiency *eff_single = new TEfficiency(*hsingle_sel,*hsingle_all);
      TEfficiency *eff_double = new TEfficiency(*hdouble_sel,*hdouble_all);
      eff_single->SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson
      eff_double->SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson


      ///////////////DRAWING/////////

      hsingle_sel->Scale(1./double(hsingle_all->GetEntries()));
      hdouble_sel->Scale(1./double(hdouble_all->GetEntries()));
      hsingle_all->Scale(1./double(hsingle_all->GetEntries()));
      hdouble_all->Scale(1./double(hdouble_all->GetEntries()));
      
      hsingle_all->SetFillColor(kGray);
      hdouble_all->SetFillColor(kGray);

          hsingle_sel->SetLineWidth(1.5);
          hdouble_sel->SetLineWidth(1.5);

          hdouble_sel->SetMarkerColor(kMagenta);
          
          hdouble_sel->SetLineColor(hdouble_sel->GetMarkerColor());

          hsingle_sel->SetMarkerStyle(23);
          hdouble_sel->SetMarkerStyle(5);

          TCanvas* can1 = new TCanvas;
          hsingle_all->Draw("HIST");
          hsingle_sel->Draw("SAME");

          hsingle_sel->SetTitle(";p_{T}; Events (normalised)");
          hsingle_sel->GetXaxis()->SetTitleOffset(hsingle_sel->GetXaxis()->GetTitleOffset()*1.15);
          hsingle_sel->SetMinimum(5e-5);
          hsingle_sel->SetMaximum(100.); //off by factor 1.3???
          can1->SetLogx();

          TLegend* leg1;
          leg1 = new TLegend(0.25,0.80,0.45,0.93);
          leg1->AddEntry(hsingle_all,"All","F");
          leg1->AddEntry(hsingle_sel,"Single-arm selection","P");
#ifdef __CINT__
          SetLegAtt(leg1);
#endif
          leg1->Draw();
#ifdef __CINT__
          CMSText(0,0,1,name[i]);
#endif


          can1->SaveAs((string("plots/diff_lrg_")+list[i]+string(".pdf")).c_str());

          TCanvas* can2 = new TCanvas;
          eff_single->SetTitle(";p_{T} [GeV];efficiency");
          eff_single->SetMarkerStyle(hsingle_sel->GetMarkerStyle());
          eff_double->SetMarkerStyle(hdouble_sel->GetMarkerStyle());
          eff_single->SetLineColor(hsingle_sel->GetLineColor());
          eff_double->SetLineColor(hdouble_sel->GetLineColor());
          eff_single->SetMarkerColor(hsingle_sel->GetLineColor());
          eff_double->SetMarkerColor(hdouble_sel->GetLineColor());
          eff_single->Draw();
          gPad->Update();
          eff_single->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.5);
          eff_double->Draw("SAME");
          can2->SetLogx();
          TLegend* leg2;
          leg2 = new TLegend(0.25,0.80,0.45,0.93);
          leg2->AddEntry(eff_single,"Single-arm selection","p");
          leg2->AddEntry(eff_double,"Double-arm selection","p");
#ifdef __CINT__
          SetLegAtt(leg2);
#endif
          leg2->Draw();
#ifdef __CINT__
          CMSText(0,0,1,name[i]);
#endif

          can2->SaveAs((string("plots/diff_lrg_eff_")+list[i]+string(".pdf")).c_str());

    } //list loop
    
}
