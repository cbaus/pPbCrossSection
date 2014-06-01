///PLEASE MAKE SURE SIGMA_HAD FOR DATA RATIO IS UP-TO-DATE!!!!
///makes table for diffractive ratios of MC

#include <TColor.h>
#include <iomanip>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TMarker.h"
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

void makePlots_pt_cuts_diff()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

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
  vec_sigma_had->Print();
  f.Close();
  double single_double_weight_data = (*vec_sigma_had)[1]/(*vec_sigma_had)[0]*100;

  vector<string> type;
  type.push_back(string("single"));
  type.push_back(string("double"));

  vector<string> list;
  //list.push_back(string("DPMJet"));
  list.push_back(string("Epos"));
  //list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  name.push_back(string("HIJING 1.383"));
  name.push_back(string("QGSJetII-04"));
  name.push_back(string("DPMJet 3.06"));

  cout << " & SD [$\\%$] & DD [$\\%$] & CD [$\\%$] & ND [$\\%$] & Sum [$\\%$] & "
       << "Ratio_{MC} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "Ratio_{Data} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "\\\\\\hline" << endl;

  for(int i=0; i<int(list.size()); i++)
    {
      TFile* file = TFile::Open("histos_deleteme.root");
      file->cd();

      TH1D* sd1_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_SD1")).c_str());
      TH1D* sd2_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_SD2")).c_str());
      TH1D* dd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_DD")).c_str());
      TH1D* cd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_CD")).c_str());
      TH1D* nd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_ND")).c_str());
      TH1D* all_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_single_all")).c_str());

      TH1D* sd1_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_SD1")).c_str());
      TH1D* sd2_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_SD2")).c_str());
      TH1D* dd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_DD")).c_str());
      TH1D* cd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_CD")).c_str());
      TH1D* nd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_ND")).c_str());
      TH1D* all_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_pt_double_all")).c_str());

      int nbins = all_single->GetNbinsX();

      int nbinsel_single = all_single->FindBin(0.61);
      double n_sd1_single = sd1_single->Integral(0,nbins+1); //includes overflow with nbins+1
      double n_sd2_single = sd2_single->Integral(0,nbins+1);
      double n_dd_single = dd_single->Integral(0,nbins+1);
      double n_cd_single = cd_single->Integral(0,nbins+1);
      double n_nd_single = nd_single->Integral(0,nbins+1);
      double n_all_single = all_single->Integral(0,nbins+1);

      double deltalowedge = (0.61 - all_single->GetBinLowEdge(nbinsel_single))/all_single->GetBinWidth(nbinsel_single);
      double n_sel_sd1_single = sd1_single->Integral(nbinsel_single,nbins+1) - deltalowedge*sd1_single->GetBinContent(nbinsel_single);
      double n_sel_sd2_single = sd2_single->Integral(nbinsel_single,nbins+1) - deltalowedge*sd2_single->GetBinContent(nbinsel_single);;
      double n_sel_dd_single = dd_single->Integral(nbinsel_single,nbins+1) - deltalowedge*dd_single->GetBinContent(nbinsel_single);;
      double n_sel_cd_single = cd_single->Integral(nbinsel_single,nbins+1) - deltalowedge*cd_single->GetBinContent(nbinsel_single);;
      double n_sel_nd_single = nd_single->Integral(nbinsel_single,nbins+1) - deltalowedge*nd_single->GetBinContent(nbinsel_single);;
      double n_sel_all_single = all_single->Integral(nbinsel_single,nbins+1) - deltalowedge*all_single->GetBinContent(nbinsel_single);;


      int nbinsel_double = all_double->FindBin(0.41);
      double n_sd1_double = sd1_double->Integral(0,nbins+1);
      double n_sd2_double = sd2_double->Integral(0,nbins+1);
      double n_dd_double = dd_double->Integral(0,nbins+1);
      double n_cd_double = cd_double->Integral(0,nbins+1);
      double n_nd_double = nd_double->Integral(0,nbins+1);
      double n_all_double = all_double->Integral(0,nbins+1);

      double deltalowedge = (0.41 - all_double->GetBinLowEdge(nbinsel_double))/all_double->GetBinWidth(nbinsel_double);
      double n_sel_sd1_double = sd1_double->Integral(nbinsel_double,nbins+1) - deltalowedge*sd1_double->GetBinContent(nbinsel_double);
      double n_sel_sd2_double = sd2_double->Integral(nbinsel_double,nbins+1) - deltalowedge*sd2_double->GetBinContent(nbinsel_double);;
      double n_sel_dd_double = dd_double->Integral(nbinsel_double,nbins+1) - deltalowedge*dd_double->GetBinContent(nbinsel_double);;
      double n_sel_cd_double = cd_double->Integral(nbinsel_double,nbins+1) - deltalowedge*cd_double->GetBinContent(nbinsel_double);;
      double n_sel_nd_double = nd_double->Integral(nbinsel_double,nbins+1) - deltalowedge*nd_double->GetBinContent(nbinsel_double);;
      double n_sel_all_double = all_double->Integral(nbinsel_double,nbins+1) - deltalowedge*all_double->GetBinContent(nbinsel_double);;

      cout << fixed << setprecision(1) //SAME FOR BOTH: EXISTS ONCE
           << "No Selection & "
           << (n_sd1_single +n_sd2_single)/n_all_single*100. << " & "
           << n_dd_single/n_all_single*100. << " & "
           << n_cd_single/n_all_single*100. << " & "
           << n_nd_single/n_all_single*100. << " & "
           << "100" << "\\\\" << endl;
      
      double single_double_weight_mc = (n_sel_all_double/n_all_double)/(n_sel_all_single/n_all_single)*100;

      cout << fixed << setprecision(1)
           << "Single-arm & "
           << (n_sel_sd1_single + n_sel_sd2_single)/n_all_single*100. << " & "
           << n_sel_dd_single/n_all_single*100. << " & "
           << n_sel_cd_single/n_all_single*100. << " & "
           << n_sel_nd_single/n_all_single*100. << " & "
           << n_sel_all_single/n_all_single*100.<< " & "
           << "\\multirow{2}{*}{" << single_double_weight_mc << "}"
           << "\\multirow{2}{*}{" << single_double_weight_data << "}" //from data. sigma_had
           << "\\\\" << endl;
      cout << fixed << setprecision(1)
           << "Double-arm & "
           << (n_sel_sd1_double + n_sel_sd2_double)/n_all_double*100. << " & "
           << n_sel_dd_double/n_all_double*100. << " & "
           << n_sel_cd_double/n_all_double*100. << " & "
           << n_sel_nd_double/n_all_double*100. << " & "
           << n_sel_all_double/n_all_double*100. << " & "
           << " " //multicolumn
           << "\\\\\\hline" << endl;

      /////////////////DO FRACTIONS IN HISTOGRAMS////////////
      for (int ibin=1; ibin<=nbins; ++ibin)
        {
          sd1_single->SetBinContent(ibin,sd1_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1)); //includes overflow with nbins+1
          sd2_single->SetBinContent(ibin,sd2_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1));
          dd_single->SetBinContent(ibin,dd_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1));
          cd_single->SetBinContent(ibin,cd_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1));
          nd_single->SetBinContent(ibin,nd_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1));
          all_single->SetBinContent(ibin,all_single->Integral(ibin,nbins+1)/all_single->Integral(ibin,nbins+1));


          sd1_double->SetBinContent(ibin,sd1_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));
          sd2_double->SetBinContent(ibin,sd2_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));
          dd_double->SetBinContent(ibin,dd_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));
          cd_double->SetBinContent(ibin,cd_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));
          nd_double->SetBinContent(ibin,nd_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));
          all_double->SetBinContent(ibin,all_double->Integral(ibin,nbins+1)/all_double->Integral(ibin,nbins+1));

        }

      ///////////////DRAWING/////////
      for(int cur=0; cur<int(type.size()); cur++)
        {
          TH1D* sd1 = type[cur]=="single"?sd1_single:sd1_double;
          TH1D* sd2 = type[cur]=="single"?sd2_single:sd2_double;
          TH1D* dd = type[cur]=="single"?dd_single:dd_double;
          TH1D* cd = type[cur]=="single"?cd_single:cd_double;
          TH1D* nd = type[cur]=="single"?nd_single:nd_double;
          TH1D* all = type[cur]=="single"?all_single:all_double;

          sd1->SetLineWidth(2.5);
          sd2->SetLineWidth(2.5);
          dd->SetLineWidth(2.5);
          cd->SetLineWidth(2.5);
          nd->SetLineWidth(2.5);
          all->SetLineWidth(2.5);

          cd->SetMarkerColor(TColor::GetColor(120,255,120));
          dd->SetMarkerColor(TColor::GetColor(255,90,90));
          sd2->SetMarkerColor(TColor::GetColor(60,60,255));
          sd1->SetMarkerColor(TColor::GetColor(00,00,150));
          nd->SetMarkerColor(kGray);

          sd1->SetLineColor(sd1->GetMarkerColor());
          sd2->SetLineColor(sd2->GetMarkerColor());
          dd->SetLineColor(dd->GetMarkerColor());
          cd->SetLineColor(cd->GetMarkerColor());
          nd->SetLineColor(nd->GetMarkerColor());

          sd1->SetFillColor(sd1->GetMarkerColor());
          sd2->SetFillColor(sd2->GetMarkerColor());
          dd->SetFillColor(dd->GetMarkerColor());
          cd->SetFillColor(cd->GetMarkerColor());
          nd->SetFillColor(nd->GetMarkerColor());

          ostringstream txt;
          txt.str(""); txt << "SD1";// (Underflow: " << fixed;
          //if (sd1->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << sd1->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          sd1->SetTitle(txt.str().c_str());

          txt.str(""); txt << "SD2";// (Underflow: " << fixed;
          //if (sd2->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << sd2->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          sd2->SetTitle(txt.str().c_str());

          txt.str(""); txt << "DD";// (Underflow: " << fixed;
          //if (dd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << dd->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          dd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "CD";// (Underflow: " << fixed;
          //if (cd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << cd->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          cd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "ND";// (Underflow: " << fixed;
          //if (nd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << nd->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          nd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "all";// (Underflow: " << fixed;
          //if (all->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << all->GetBinContent(0)*1e3 << "x10^{-3})";
          //else txt << setprecision(0) << 0 << ")";
          all->SetTitle(txt.str().c_str());

          // txt.str(""); txt << "SD1"; sd1->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "SD2"; sd2->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "DD"; dd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "CD"; cd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "ND"; nd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "all)"; all->SetTitle(txt.str().c_str());

          TCanvas* can1 = new TCanvas;
          THStack* hs = new THStack("hs","EPOS");
          hs->Add(sd1);
          hs->Add(sd2);
          hs->Add(dd);
          hs->Add(cd);
          //hs->Add(nd);
          hs->Draw("HIST");

          hs->SetTitle(";p_{T,HF}^{GEN} [GeV];fraction");
          hs->GetXaxis()->SetTitleOffset(hs->GetXaxis()->GetTitleOffset()*1.15);
          hs->GetXaxis()->SetRangeUser(0,3);
          hs->SetMinimum(0);
          hs->SetMaximum(0.2); //off by factor 1.3???
          //can1->SetLogy();
          //can1->SetLogx();

          TLegend* leg1;
          if(list[i]!="Hijing")
            {
              leg1 = new TLegend(0.25,0.63,0.45,0.93);
              leg1->AddEntry(sd1,"","F");
              leg1->AddEntry(sd2,"","F");
              leg1->AddEntry(dd,"","F");
              leg1->AddEntry(cd,"","F");
              //leg1->AddEntry(nd,"","F");
            }
          else
            {
              leg1 = new TLegend(0.25,0.87,0.45,0.93);
              leg1->AddEntry(nd,all->GetTitle(),"F");
            }
#ifdef __CINT__
          SetLegAtt(leg1);
#endif
          leg1->Draw();
#ifdef __CINT__
          CMSText(0,0,1,name[i],type[cur]=="single"?"single-arm selection":"double-arm selection");
#endif

          TLine* line = new TLine(type[cur]=="single"?8:4,2e-5,type[cur]=="single"?8:4,4e-2);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          //can1->SaveAs((string("plots/diff_energy_")+type[cur]+string("_")+list[i]+string(".pdf")).c_str());
        } //type loop
    } //list loop
}
