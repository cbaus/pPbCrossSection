///PLEASE MAKE SURE SIGMA_HAD FOR DATA RATIO IS UP-TO-DATE!!!!
///makes table for diffractive ratios of MC
///makes two plots, how ratio_S/D depends on scaling of sigma_diff

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

void makePlots_diff2()
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
  double single_double_weight_data = (*vec_sigma_had)[1]/(*vec_sigma_had)[0];
  //uncertainty on ratio from EM, EvenSel, Noise = sqrt(0.2**2+0.6**2+0.2**2+1.2**2+0.2**2) = 1.39%
  double ratio_unc = 0.0139;
  double ratio_unc2 = ratio_unc*2.;


  vector<string> type;
  type.push_back(string("single"));
  type.push_back(string("double"));

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));
  list.push_back(string("DPMJet"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  name.push_back(string("HIJING 1.383"));
  name.push_back(string("QGSJetII-04"));
  name.push_back(string("DPMJet 3.06"));

  cout << " & SD [$\\%$] & DD [$\\%$] & CD [$\\%$] & ND [$\\%$] & Sum [$\\%$] & "
       << "Ratio_{MC} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "Ratio_{Data} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "\\\\\\hline" << endl;

  const int nDiffWeight = 300;
  TGraph graphFindDiffWeightHijing(1);
  TGraph graphFindDiffWeightDpmjet(1);
  TGraph graphEffDiffWeightHijing(1);
  TGraph graphEffDiffWeightDpmjet(1);
  TGraph graphFindDiffWeightEpos(nDiffWeight);
  TGraph graphEffDiffWeightEpos(nDiffWeight);
  TGraph graphFindDiffWeightQgsjet(nDiffWeight);
  TGraph graphEffDiffWeightQgsjet(nDiffWeight);

  for(int i=0; i<int(list.size()); i++)
    {
      cout << i+1 << "/" << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");
      file->cd();

      TH1D* sd1_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_SD1")).c_str());
      TH1D* sd2_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_SD2")).c_str());
      TH1D* dd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_DD")).c_str());
      TH1D* cd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_CD")).c_str());
      TH1D* nd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_ND")).c_str());
      TH1D* all_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_all")).c_str());

      TH1D* sd1_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_SD1")).c_str());
      TH1D* sd2_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_SD2")).c_str());
      TH1D* dd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_DD")).c_str());
      TH1D* cd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_CD")).c_str());
      TH1D* nd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_ND")).c_str());
      TH1D* all_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_all")).c_str());

      TH1D* eff_sd1_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_SD1_single")).c_str());
      TH1D* eff_sd2_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_SD2_single")).c_str());
      TH1D* eff_dd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_DD_single")).c_str());
      TH1D* eff_cd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_CD_single")).c_str());
      TH1D* eff_nd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_ND_single")).c_str());
      TH1D* eff_all_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_single")).c_str());
      eff_nd_single->Scale(1./double(eff_all_single->GetBinContent(1)));
      eff_sd1_single->Scale(1./double(eff_all_single->GetBinContent(1)));
      eff_sd2_single->Scale(1./double(eff_all_single->GetBinContent(1)));
      eff_cd_single->Scale(1./double(eff_all_single->GetBinContent(1)));
      eff_dd_single->Scale(1./double(eff_all_single->GetBinContent(1)));
      eff_all_single->Scale(1./double(eff_all_single->GetBinContent(1)));

      TH1D* eff_sd1_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_SD1_double")).c_str());
      TH1D* eff_sd2_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_SD2_double")).c_str());
      TH1D* eff_dd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_DD_double")).c_str());
      TH1D* eff_cd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_CD_double")).c_str());
      TH1D* eff_nd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_ND_double")).c_str());
      TH1D* eff_all_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_hf_cut_double")).c_str());
      eff_nd_double->Scale(1./double(eff_all_double->GetBinContent(1)));
      eff_sd1_double->Scale(1./double(eff_all_double->GetBinContent(1)));
      eff_sd2_double->Scale(1./double(eff_all_double->GetBinContent(1)));
      eff_cd_double->Scale(1./double(eff_all_double->GetBinContent(1)));
      eff_dd_double->Scale(1./double(eff_all_double->GetBinContent(1)));
      eff_all_double->Scale(1./double(eff_all_double->GetBinContent(1)));

      int cut_single = 8.;
      int cut_double = 4.;

      double n_sd1_single = eff_sd1_single->GetBinContent(1);
      double n_sd2_single = eff_sd2_single->GetBinContent(1);
      double n_dd_single = eff_dd_single->GetBinContent(1);
      double n_cd_single = eff_cd_single->GetBinContent(1);
      double n_nd_single = eff_nd_single->GetBinContent(1);
      double n_all_single = eff_all_single->GetBinContent(1);

      double n_sel_sd1_single = eff_sd1_single->Interpolate(cut_single);
      double n_sel_sd2_single = eff_sd2_single->Interpolate(cut_single);
      double n_sel_dd_single = eff_dd_single->Interpolate(cut_single);
      double n_sel_cd_single = eff_cd_single->Interpolate(cut_single);
      double n_sel_nd_single = eff_nd_single->Interpolate(cut_single);
      double n_sel_all_single = eff_all_single->Interpolate(cut_single);

      double n_sd1_double = eff_sd1_double->GetBinContent(1);
      double n_sd2_double = eff_sd2_double->GetBinContent(1);
      double n_dd_double = eff_dd_double->GetBinContent(1);
      double n_cd_double = eff_cd_double->GetBinContent(1);
      double n_nd_double = eff_nd_double->GetBinContent(1);
      double n_all_double = eff_all_double->GetBinContent(1);

      double n_sel_sd1_double = eff_sd1_double->Interpolate(cut_double);
      double n_sel_sd2_double = eff_sd2_double->Interpolate(cut_double);
      double n_sel_dd_double = eff_dd_double->Interpolate(cut_double);
      double n_sel_cd_double = eff_cd_double->Interpolate(cut_double);
      double n_sel_nd_double = eff_nd_double->Interpolate(cut_double);
      double n_sel_all_double = eff_all_double->Interpolate(cut_double);

      double single_double_weight_mc = n_sel_all_double / n_sel_all_single;

      cout << fixed << setprecision(1) //SAME FOR BOTH: EXISTS ONCE
           << "No Selection & "
           << (n_sd1_single +n_sd2_single)*100. << " & "
           << n_dd_single*100. << " & "
           << n_cd_single*100. << " & "
           << n_nd_single*100. << " & "
           << "100" << "\\\\" << endl;

      cout << fixed << setprecision(1)
           << "Single-arm & "
           << (n_sel_sd1_single + n_sel_sd2_single)*100. << " & "
           << n_sel_dd_single*100. << " & "
           << n_sel_cd_single*100. << " & "
           << n_sel_nd_single*100. << " & "
           << n_sel_all_single*100.<< " & " << setprecision(3)
           << "\\multirow{2}{*}{" << single_double_weight_mc << "}"
           << "\\multirow{2}{*}{" << single_double_weight_data << "\\pm" << single_double_weight_data*ratio_unc <<"}" //from data. sigma_had
           << "\\\\" << endl;
      cout << fixed << setprecision(1)
           << "Double-arm & "
           << (n_sel_sd1_double + n_sel_sd2_double)*100. << " & "
           << n_sel_dd_double*100. << " & "
           << n_sel_cd_double*100. << " & "
           << n_sel_nd_double*100. << " & "
           << n_sel_all_double*100. << " & "
           << " " //multicolumn
           << "\\\\\\hline" << endl;

      ///////////////Find optimal diffractive cs weight//////////////////////////
      double diff_frac = (n_sd1_single + n_sd2_single + n_dd_single + n_cd_single); //same for double
      double diff_frac_sel_single = (n_sel_sd1_single + n_sel_sd2_single + n_sel_dd_single + n_sel_cd_single);
      double diff_frac_sel_double = (n_sel_sd1_double + n_sel_sd2_double + n_sel_dd_double + n_sel_cd_double);
      double maxDiffWeight = 3;
      TGraph* graphFindDiffWeight = 0;
      TGraph* graphEffDiffWeight = 0;
      if(list[i]=="Epos")
        {
          graphFindDiffWeight = &graphFindDiffWeightEpos;
          graphEffDiffWeight = &graphEffDiffWeightEpos;
        }
      else if(list[i]=="QGSJetII")
        {
          graphFindDiffWeight = &graphFindDiffWeightQgsjet;
          graphEffDiffWeight = &graphEffDiffWeightQgsjet;
        }
      else if(list[i]=="Hijing")
        {
          graphEffDiffWeightHijing.SetPoint(0,1, n_sel_all_double);
          graphFindDiffWeightHijing.SetPoint(0,1, n_sel_all_double / n_sel_all_single);
        }
      else if(list[i]=="DPMJet")
        {
          graphEffDiffWeightDpmjet.SetPoint(0,1, n_sel_all_double);
          graphFindDiffWeightDpmjet.SetPoint(0,1, n_sel_all_double / n_sel_all_single);
        }

      if(graphFindDiffWeight)
        for (int j = 0; j < nDiffWeight; ++j)
          {
            double diffWeight = maxDiffWeight/double(nDiffWeight)*double(j);
            double single_double_weight =
              (diffWeight * diff_frac_sel_double + n_sel_nd_double) / (1.+diff_frac*(diffWeight-1.)) /
              (diffWeight * diff_frac_sel_single + n_sel_nd_single) / (1.+diff_frac*(diffWeight-1.));
            graphFindDiffWeight->SetPoint(j,diffWeight,single_double_weight);
            graphEffDiffWeight->SetPoint(j,diffWeight,(diffWeight * diff_frac_sel_double + n_sel_nd_double) / (1.+diff_frac*(diffWeight-1.)));
          }


      ///////////////SCALING DIFF TO 25%///////////////diff*x/(1+diff(x-1)=0.25 -> x=(0.25/diff+0.25)/0.75
      double to25 = diff_frac!=0?(0.25/diff_frac-0.25)/0.75:0;
      cout << endl << "fraction of diffraction: " << diff_frac*100. << " weight to 25%: " << setprecision(3) << to25*100. << endl << endl;

      ///////////////DRAWING/////////
      for(int cur=0; cur<int(type.size()); cur++)
        {
          TH1D* sd1 = type[cur]=="single"?sd1_single:sd1_double;
          TH1D* sd2 = type[cur]=="single"?sd2_single:sd2_double;
          TH1D* dd = type[cur]=="single"?dd_single:dd_double;
          TH1D* cd = type[cur]=="single"?cd_single:cd_double;
          TH1D* nd = type[cur]=="single"?nd_single:nd_double;
          TH1D* all = type[cur]=="single"?all_single:all_double;

          sd1->Scale(1./double(all->GetEntries()));
          sd2->Scale(1./double(all->GetEntries()));
          dd->Scale(1./double(all->GetEntries()));
          cd->Scale(1./double(all->GetEntries()));
          nd->Scale(1./double(all->GetEntries()));
          all->Scale(1./double(all->GetEntries()));

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
          txt.str(""); txt << "SD1 (Underflow: " << fixed;
          if (sd1->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << sd1->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
          sd1->SetTitle(txt.str().c_str());

          txt.str(""); txt << "SD2 (Underflow: " << fixed;
          if (sd2->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << sd2->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
          sd2->SetTitle(txt.str().c_str());

          txt.str(""); txt << "DD (Underflow: " << fixed;
          if (dd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << dd->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
          dd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "CD (Underflow: " << fixed;
          if (cd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << cd->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
          cd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "ND (Underflow: " << fixed;
          if (nd->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << nd->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
          nd->SetTitle(txt.str().c_str());

          txt.str(""); txt << "all (Underflow: " << fixed;
          if (all->GetBinContent(0)*1e3 >= 0.05) txt << setprecision(1) << all->GetBinContent(0)*1e3 << "x10^{-3})";
          else txt << setprecision(0) << 0 << ")";
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
          hs->Add(nd);
          hs->Draw("HIST");

          hs->SetTitle(";E_{HF} [GeV];events (normalised)");
          hs->GetXaxis()->SetTitleOffset(hs->GetXaxis()->GetTitleOffset()*1.15);
          hs->SetMinimum(5e-5);
          hs->SetMaximum(10.); //off by factor 1.3???
          can1->SetLogy();
          can1->SetLogx();

          TLegend* leg1;
          if(list[i]!="Hijing")
            {
              leg1 = new TLegend(0.25,0.63,0.45,0.93);
              leg1->AddEntry(sd1,"","F");
              leg1->AddEntry(sd2,"","F");
              leg1->AddEntry(dd,"","F");
              leg1->AddEntry(cd,"","F");
              leg1->AddEntry(nd,"","F");
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

          can1->SaveAs((string("plots/diff_energy_")+type[cur]+string("_")+list[i]+string(".pdf")).c_str());
        } //type loop
    } //list loop

  //////////Draw opt diff plot///////////
  TCanvas* can1 = new TCanvas;
  graphFindDiffWeightQgsjet.SetMarkerStyle(25);
  graphFindDiffWeightQgsjet.SetMarkerStyle(25);
  graphFindDiffWeightQgsjet.SetMarkerSize(1.3);
  graphFindDiffWeightQgsjet.SetLineWidth(2);
  graphFindDiffWeightQgsjet.SetLineColor(kGreen+2);
  graphFindDiffWeightQgsjet.SetLineStyle(9);
  graphFindDiffWeightQgsjet.SetTitle(";#sigma_{diff} scale factor;R_{d/s}");
  graphFindDiffWeightQgsjet.Draw("AC");
  graphFindDiffWeightQgsjet.GetXaxis()->SetRangeUser(0.5,1.5);
  graphFindDiffWeightQgsjet.GetYaxis()->SetRangeUser(0.85,1.10);
  graphFindDiffWeightEpos.SetMarkerSize(1.3);
  graphFindDiffWeightEpos.SetLineWidth(2);
  graphFindDiffWeightEpos.SetLineColor(kBlue);
  graphFindDiffWeightEpos.Draw("C");

  graphFindDiffWeightDpmjet.SetMarkerColor(kMagenta);
  graphFindDiffWeightHijing.SetMarkerColor(kRed);
  graphFindDiffWeightDpmjet.SetMarkerStyle(26);
  graphFindDiffWeightHijing.SetMarkerStyle(4);
  //  graphFindDiffWeightDpmjet.Draw("P");
  //  graphFindDiffWeightHijing.Draw("P");

  TLine* line0 = new TLine(0.5,single_double_weight_data,1.5,single_double_weight_data);
  line0->SetLineStyle(1);
  line0->Draw("SAME");
  TLine* line1a = new TLine(0.5,single_double_weight_data*(1.+ratio_unc),1.5,single_double_weight_data*(1.+ratio_unc));
  line1a->SetLineStyle(3);
  line1a->Draw("SAME");
  TLine* line2a = new TLine(0.5,single_double_weight_data*(1.+ratio_unc2),1.5,single_double_weight_data*(1.+ratio_unc2));
  line2a->SetLineStyle(2);
  line2a->Draw("SAME");
  TLine* line1b = new TLine(0.5,single_double_weight_data/(1.+ratio_unc),1.5,single_double_weight_data/(1.+ratio_unc));
  line1b->SetLineStyle(3);
  line1b->Draw("SAME");
  TLine* line2b = new TLine(0.5,single_double_weight_data/(1.+ratio_unc2),1.5,single_double_weight_data/(1.+ratio_unc2));
  line2b->SetLineStyle(2);
  line2b->Draw("SAME");
  TLegend* leg_opt = new TLegend(0.25,0.2,0.56,0.4);
#ifdef __CINT__
  SetLegAtt(leg_opt);
#endif
  leg_opt->AddEntry(line0,"Data","l");
  leg_opt->AddEntry(line1a,"Data #pm 1 std. dev.","l");
  leg_opt->AddEntry(line2a,"Data #pm 2 std. dev.","l");
  leg_opt->AddEntry(&graphFindDiffWeightEpos,"EPOS-LHC","l");
  leg_opt->AddEntry(&graphFindDiffWeightQgsjet,"QGSJetII","l");
  //  leg->AddEntry(&graphFindDiffWeightHijing,"Hijing","P");
  //  leg->AddEntry(&graphFindDiffWeightDpmjet,"Dpmjet","P");
  leg_opt->Draw("SAME");

  CMSText(3,0,1);
  can1->SaveAs((string("plots/diff_optimal_weight.pdf")).c_str());

  double x_opt=1;
  cout << endl << "-------Optimal value for sigma_diff scale factor--------" << endl;

  x_opt=0.5;
  while(graphFindDiffWeightQgsjet.Eval(x_opt)>single_double_weight_data * (1.+ratio_unc2))
    x_opt += 0.001;
  cout << "QGSJetII-04" << " " << x_opt;
  x_opt=0.5;
  while(graphFindDiffWeightQgsjet.Eval(x_opt)>single_double_weight_data)
    x_opt += 0.001;
  cout << " " << x_opt;
  x_opt=0.5;
  while(graphFindDiffWeightQgsjet.Eval(x_opt)>single_double_weight_data / (1.+ratio_unc2))
    x_opt += 0.001;
  cout << " " << x_opt << endl;

  x_opt=0.5;
  while(graphFindDiffWeightEpos.Eval(x_opt)>single_double_weight_data * (1.+ratio_unc2))
    x_opt += 0.001;
  cout << "EPOS-LHC" << " " << x_opt;
  x_opt=0.5;
  while(graphFindDiffWeightEpos.Eval(x_opt)>single_double_weight_data)
    x_opt += 0.001;
  cout << " " << x_opt;
  x_opt=0.5;
  while(graphFindDiffWeightEpos.Eval(x_opt)>single_double_weight_data / (1.+ratio_unc2))
    x_opt += 0.001;
  cout << " " << x_opt << endl;

  //////////Draw diff scaling - eff plot///////////
  TCanvas* can2 = new TCanvas;
  graphEffDiffWeightQgsjet.SetMarkerStyle(25);
  graphEffDiffWeightQgsjet.SetMarkerSize(1.3);
  graphEffDiffWeightQgsjet.SetLineWidth(2);
  graphEffDiffWeightQgsjet.SetLineStyle(9);
  graphEffDiffWeightQgsjet.SetLineColor(kGreen+2);
  graphEffDiffWeightQgsjet.SetTitle(";#sigma_{diff} scale factor;double-arm #epsilon_{acc}");
  graphEffDiffWeightQgsjet.Draw("AC");
  graphEffDiffWeightQgsjet.GetXaxis()->SetRangeUser(0.5,1.5);
  graphEffDiffWeightQgsjet.GetYaxis()->SetRangeUser(0.8,1);
  graphEffDiffWeightEpos.SetMarkerSize(1.3);
  graphEffDiffWeightEpos.SetLineWidth(2);
  graphEffDiffWeightEpos.SetLineColor(kBlue);
  graphEffDiffWeightEpos.Draw("C");
  graphEffDiffWeightDpmjet.SetMarkerColor(kMagenta);
  graphEffDiffWeightHijing.SetMarkerColor(kRed);
  graphEffDiffWeightDpmjet.SetMarkerStyle(26);
  graphEffDiffWeightHijing.SetMarkerStyle(4);
  //  graphEffDiffWeightDpmjet.Draw("P");
  //  graphEffDiffWeightHijing.Draw("P");

  CMSText(3,0,1);
  TLegend* leg = new TLegend(0.25,0.2,0.6,0.3);
#ifdef __CINT__
  SetLegAtt(leg);
#endif
  leg->AddEntry(&graphFindDiffWeightEpos,"EPOS-LHC","l");
  leg->AddEntry(&graphFindDiffWeightQgsjet,"QGSJetII","l");
  leg->Draw("SAME");
  can2->SaveAs((string("plots/diff_eff_diffweight.pdf")).c_str());

}
