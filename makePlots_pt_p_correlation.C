///makes plots for correlation between p or p_t on gen leven with E_HF on det level

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

void makePlots_pt_p_correlation()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("QGSJetII"));
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
  name.push_back(string("QGSJetII-04"));
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
      cout << i+1 << "/" << int(list.size()) << endl; TFile* file = TFile::Open("histos_deleteme.root");
      file->cd();

      TH2D* cor_p_single=(TH2D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_p_Ehf_correlation_single")).c_str());
      TH2D* cor_p_double=(TH2D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_p_Ehf_correlation_double")).c_str());
      TH2D* cor_pt_single=(TH2D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_pt_Ehf_correlation_single")).c_str());
      TH2D* cor_pt_double=(TH2D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_pt_Ehf_correlation_double")).c_str());

      ///////////////DRAWING/////////
      gStyle->SetPadLeftMargin(0.17); //0.400
      gStyle->SetPadRightMargin(0.14); //0.400
      gROOT->ForceStyle();
      ostringstream desc;
      TCanvas* can1 = 0;

      can1 = new TCanvas;
      cor_p_single->Draw("COLZ");
      can1->SetLogz();
      cor_p_single->SetTitle(";p (Gen) [GeV]; E_{HF} (Det) [GeV]");
      desc.str(""); desc << "Cor. Factor: " << setprecision(2) << fixed << cor_p_single->GetCorrelationFactor();
#ifdef __CINT__
      CMSText(0,0,1,"single-arm selection",desc.str(),name[i]);
#endif
      can1->SaveAs((string("plots/corr_p_ehf_single_")+list[i]+string(".png")).c_str());

      can1 = new TCanvas;
      cor_p_double->Draw("COLZ");
      can1->SetLogz();
      cor_p_double->SetTitle(";p (Gen) [GeV]; E_{HF} (Det) [GeV]");
      desc.str(""); desc << "Cor. Factor: " << setprecision(2) << fixed << cor_p_double->GetCorrelationFactor();
#ifdef __CINT__
      CMSText(0,0,1,"double-arm selection",desc.str(),name[i]);
#endif
      can1->SaveAs((string("plots/corr_p_ehf_double_")+list[i]+string(".png")).c_str());

      can1 = new TCanvas;
      cor_pt_single->Draw("COLZ");
      can1->SetLogz();
      cor_pt_single->SetTitle(";p_{T} (Gen) [GeV]; E_{HF} (Det) [GeV]");
      desc.str(""); desc << "Cor. Factor: " << setprecision(2) << fixed << cor_pt_single->GetCorrelationFactor();
#ifdef __CINT__
      CMSText(0,0,1,"single-arm selection",desc.str(),name[i]);
#endif
      can1->SaveAs((string("plots/corr_pt_ehf_single_")+list[i]+string(".png")).c_str());

      can1 = new TCanvas;
      cor_pt_double->Draw("COLZ");
      can1->SetLogz();
      cor_pt_double->SetTitle(";p_{T} (Gen) [GeV]; E_{HF} (Det) [GeV]");
      desc.str(""); desc << "Cor. Factor: " << setprecision(2) << fixed << cor_pt_double->GetCorrelationFactor();
#ifdef __CINT__
      CMSText(0,0,1,"double-arm selection",desc.str(),name[i]);
#endif
      can1->SaveAs((string("plots/corr_pt_ehf_double_")+list[i]+string(".png")).c_str());

    } //list loop
    
}
