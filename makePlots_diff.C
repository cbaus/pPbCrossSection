///makes plots for largest rapidity gap

#include <TColor.h>
#include <iomanip>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraph.h"
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

void makePlots_diff()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("QGSJetII"));
  list.push_back(string("Hijing"));
  list.push_back(string("DPMJet"));
  //  list.push_back(string("DPMJet_pp"));
//   list.push_back(string("EposDiffWeight150"));
//   list.push_back(string("EposDiffWeight200"));
//   list.push_back(string("EposDiffWeight299"));
//   list.push_back(string("QGSJetIIDiffWeight150"));
//   list.push_back(string("QGSJetIIDiffWeight200"));
//   list.push_back(string("QGSJetIIDiffWeight452"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  name.push_back(string("QGSJetII-04"));
  name.push_back(string("HIJING 1.383"));
  name.push_back(string("DPMJet 3.06"));
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
      TFile* file = TFile::Open("histos_deleteme.root");
      file->cd();

      TH1D* sd1=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_SD1")).c_str());
      TH1D* sd2=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_SD2")).c_str());
      TH1D* dd=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_DD")).c_str());
      TH1D* cd=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_CD")).c_str());
      TH1D* nd=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_ND")).c_str());
      TH1D* all=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_all")).c_str());
      TH1D* hsingle=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_single")).c_str());
      TH1D* hdouble=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_lrg_double")).c_str());

      TEfficiency *eff_single = new TEfficiency(*hsingle,*all);
      TEfficiency *eff_double = new TEfficiency(*hdouble,*all);
      eff_single->SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson
      eff_double->SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson


      ///////////////DRAWING/////////

          sd1->Scale(1./double(all->GetEntries()));
          sd2->Scale(1./double(all->GetEntries()));
          dd->Scale(1./double(all->GetEntries()));
          cd->Scale(1./double(all->GetEntries()));
          nd->Scale(1./double(all->GetEntries()));
          all->Scale(1./double(all->GetEntries()));
          hsingle->Scale(1./double(all->GetEntries()));
          hdouble->Scale(1./double(all->GetEntries()));

          sd1->SetLineWidth(2.5);
          sd2->SetLineWidth(2.5);
          dd->SetLineWidth(2.5);
          cd->SetLineWidth(2.5);
          nd->SetLineWidth(2.5);
          all->SetLineWidth(2.5);
          hsingle->SetLineWidth(1.5);
          hdouble->SetLineWidth(1.5);

          cd->SetMarkerColor(TColor::GetColor(120,255,120));
          dd->SetMarkerColor(TColor::GetColor(255,90,90));
          sd2->SetMarkerColor(TColor::GetColor(60,60,255));
          sd1->SetMarkerColor(TColor::GetColor(00,00,150));
          nd->SetMarkerColor(kGray);
          hdouble->SetMarkerColor(kMagenta);
          
          sd1->SetLineColor(sd1->GetMarkerColor());
          sd2->SetLineColor(sd2->GetMarkerColor());
          dd->SetLineColor(dd->GetMarkerColor());
          cd->SetLineColor(cd->GetMarkerColor());
          nd->SetLineColor(nd->GetMarkerColor());
          hdouble->SetLineColor(hdouble->GetMarkerColor());

          sd1->SetFillColor(sd1->GetMarkerColor());
          sd2->SetFillColor(sd2->GetMarkerColor());
          dd->SetFillColor(dd->GetMarkerColor());
          cd->SetFillColor(cd->GetMarkerColor());
          nd->SetFillColor(nd->GetMarkerColor());

          hsingle->SetMarkerStyle(23);
          hdouble->SetMarkerStyle(25);

          TCanvas* can1 = new TCanvas;
          THStack* hs = new THStack("hs","EPOS");
          hs->Add(sd1);
          hs->Add(sd2);
          hs->Add(dd);
          hs->Add(cd);
          hs->Add(nd);
          hs->Draw("HIST");
          hsingle->Draw("SAME");
          hdouble->Draw("SAME");

          hs->SetTitle(";#Deltay^{F}; Events (normalised)");
          hs->GetXaxis()->SetTitleOffset(hs->GetXaxis()->GetTitleOffset()*1.15);
          hs->SetMinimum(5e-5);
          hs->SetMaximum(100.); //off by factor 1.3???
          can1->SetLogy();
          //can1->SetLogx();

          TLegend* leg1;
          if(list[i].find("Hijing") == string::npos && list[i].find("DPMJet") == string::npos)
            {
              leg1 = new TLegend(0.35,0.63,0.55,0.93);
              leg1->AddEntry(sd1,"SD1","F");
              leg1->AddEntry(sd2,"SD2","F");
              leg1->AddEntry(dd,"DD","F");
              leg1->AddEntry(cd,"CD","F");
              leg1->AddEntry(nd,"ND","F");
              leg1->AddEntry(hsingle,"Single-arm selection","P");
              leg1->AddEntry(hdouble,"Double-arm selection","P");
            }
          else
            {
              leg1 = new TLegend(0.25,0.80,0.45,0.93);
              leg1->AddEntry(nd,"All","F");
              leg1->AddEntry(hsingle,"Single-arm selection","P");
              leg1->AddEntry(hdouble,"Double-arm selection","P");
            }
#ifdef __CINT__
          SetLegAtt(leg1);
#endif
          leg1->Draw();
#ifdef __CINT__
          CMSText(0,0,1,name[i]);
#endif


          can1->SaveAs((string("plots/diff_lrg_")+list[i]+string(".pdf")).c_str());

          TCanvas* can2 = new TCanvas;
          eff_single->SetTitle(";#Deltay^{F};efficiency");
          eff_single->SetMarkerStyle(hsingle->GetMarkerStyle());
          eff_double->SetMarkerStyle(hdouble->GetMarkerStyle());
          eff_single->SetLineColor(hsingle->GetLineColor());
          eff_double->SetLineColor(hdouble->GetLineColor());
          eff_single->SetMarkerColor(hsingle->GetLineColor());
          eff_double->SetMarkerColor(hdouble->GetLineColor());
          eff_single->Draw();
          gPad->Update();
          eff_single->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.5);
          eff_double->Draw("SAME");
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
