#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVectorD.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iomanip>
#include <iostream>
#include <numeric>      // std::accumulate
#include <sstream>
#include <set>
#include <map>
#include <string>
#include <utility>

#define _LumiCorrpPb 1.142 //only use if trees don't contain vdm calibration factor
#define _LumiCorrPbp 1.138

//#define _diff_qgs_reweight 1.202
//#define _diff_epos_reweight 1.132
#define _diff_qgs_reweight 0.841
#define _diff_epos_reweight 0.877

using namespace std;

TVectorD corr_fac_em(2);
TVectorD corr_fac_eme(2);
TVectorD corr_fac_mc(4);
TVectorD corr_fac_mce(4);
TVectorD corr_fac_epos(2);
TVectorD corr_fac_qgsjet(2);
TVectorD corr_fac_hijing(2);
TVectorD corr_fac_dpmjet(2);

typedef map<string,TH1D*> histomap;
map< string,histomap > histoman;

void makePlots_cs_eff(bool draw=1, double cut_value_single = 8., double cut_value_double = 4.,string filename = "histos.root");
double getEffDiffWeight(string model, double cut, double diffWeight);

//available cuts single 5, 6.4, 7.2, 8, 8.8, 9.6, 10, 15, 20
//available cuts double 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5

void makePlots_cs_eff(bool draw, double cut_value_single, double cut_value_double,string filename)
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<string> type; type.push_back("single"); type.push_back("double");
  vector<string> models;
  models.push_back("DPMJet");
  models.push_back("Hijing");
  models.push_back("Epos");
  models.push_back("QGSJetII");
  
  for(int n=0; n<int(type.size()); n++)
    {
      TFile* file = TFile::Open(filename.c_str());
      TFile* file2 = TFile::Open("histos_deleteme.root");

      TH1D* zb=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n]).c_str());
      TH1D* noise=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n] + string("_noise")).c_str());
      TH1D* h_events=(TH1D*)file->Get(string("data210885/data210885_h_lumi").c_str()); //zb
      TH1D* h_lumi=(TH1D*)file->Get(string("data210885/data210885_h_run_events_lumi").c_str());
      TH1D* sl1=(TH1D*)file->Get((string("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_") + type[n]).c_str());
      TH1D* sl2=(TH1D*)file->Get((string("Starlight_Pythia/Starlight_Pythia_h_hf_cut_") + type[n]).c_str());

      TH1D* eposrew=(TH1D*)file->Get((string("EposDiffWeightOpt/EposDiffWeightOpt_h_hf_cut_") + type[n]).c_str());
      TH1D* qgsrew=(TH1D*)file->Get((string("QGSJetIIDiffWeightOpt/QGSJetIIDiffWeightOpt_h_hf_cut_") + type[n]).c_str());

      histoman.clear();
      for (int i=0; i<int(models.size()); ++i)
        {
          string model=models[i];
          histoman[model]["ALL"]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_") + type[n]).c_str());
          histoman[model]["SD1"]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_SD1_") + type[n]).c_str());
          histoman[model]["SD2"]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_SD2_") + type[n]).c_str());
          histoman[model]["CD" ]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_CD_") + type[n]).c_str());
          histoman[model]["DD" ]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_DD_") + type[n]).c_str());
          histoman[model]["ND" ]=(TH1D*)file->Get((model+string("/")+model+string("_h_hf_cut_ND_") + type[n]).c_str());
        }
      
      TH1D* dpm=histoman["DPMJet"]["ALL"];
      TH1D* hijing=histoman["Hijing"]["ALL"];
      TH1D* epos=histoman["Epos"]["ALL"];
      TH1D* qgs=histoman["QGSJetII"]["ALL"];

      /////
      //LUMI
      double lumi_integral = 0;
      double events_integral = 0;
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
          const double lumicorr = _LumiCorrpPb; //might be different for other run than 210885
          const double lumiPerLS=h_lumi->GetBinContent(i) * lumicorr;
          const double lumiPerLS_error=h_lumi->GetBinError(i) * lumicorr; //not 100% correct since from profile but has no contribution
          if (lumiPerLS<0.) {cerr << "lumi neg: " << i << endl; return;}
          else if (lumiPerLS==0.) {continue;}
          lumi_integral += lumiPerLS;
          events_integral += h_events->GetBinContent(i);
        }

      //zb->Scale(1./double(zb->GetBinContent(1)));
      //noise->Scale(1./double(noise->GetBinContent(1)));
      hijing->Scale(1./double(hijing->GetBinContent(1)));
      epos->Scale(1./double(epos->GetBinContent(1)));
      qgs->Scale(1./double(qgs->GetBinContent(1)));
      dpm->Scale(1./double(dpm->GetBinContent(1)));
      sl1->Scale(1./double(sl1->GetBinContent(1)));
      sl2->Scale(1./double(sl2->GetBinContent(1))/ 195. * 122.); //cross section starlight samples
      if(eposrew) eposrew->Scale(1./double(eposrew->GetBinContent(1)));
      if(qgsrew) qgsrew->Scale(1./double(qgsrew->GetBinContent(1)));

      zb->SetLineWidth(3);
      noise->SetLineWidth(3);
      hijing->SetLineWidth(3);
      epos->SetLineWidth(3);
      qgs->SetLineWidth(3);
      dpm->SetLineWidth(3);
      sl1->SetLineWidth(3);
      sl2->SetLineWidth(3);
      if(eposrew) eposrew->SetLineWidth(3);
      if(qgsrew) qgsrew->SetLineWidth(3);
      noise->SetLineColor(kBlack);
      hijing->SetLineColor(kGreen+2);
      epos->SetLineColor(kBlue);
      qgs->SetLineColor(kRed);
      dpm->SetLineColor(kOrange);
      sl1->SetLineColor(kRed);
      sl2->SetLineColor(kBlue);
      if(eposrew) eposrew->SetLineColor(kMagenta);
      if(qgsrew) qgsrew->SetLineColor(kMagenta);
      noise->SetMarkerColor(noise->GetLineColor());
      hijing->SetMarkerColor(hijing->GetLineColor());
      epos->SetMarkerColor(epos->GetLineColor());
      qgs->SetMarkerColor(qgs->GetLineColor());
      dpm->SetMarkerColor(dpm->GetLineColor());
      sl1->SetMarkerColor(sl1->GetLineColor());
      sl2->SetMarkerColor(sl2->GetLineColor());
      if(eposrew) eposrew->SetMarkerColor(sl2->GetLineColor());
      if(qgsrew) qgsrew->SetMarkerColor(sl2->GetLineColor());
      hijing->SetLineStyle(7);
      qgs->SetLineStyle(9);
      sl2->SetLineStyle(9);
      if(eposrew) eposrew->SetLineStyle(10);
      if(qgsrew) qgsrew->SetLineStyle(10);

      zb->SetTitle("Data");
      noise->SetTitle("Noise");
      if(eposrew) eposrew->SetTitle("EPOS-LHC (#sigma_{diff}x1.12)");
      if(qgsrew) qgsrew->SetTitle("QGSJETII-04 (#sigma_{diff}x1.50)");
      hijing->SetTitle("Hijing 1.383");
      epos->SetTitle("EPOS-LHC");
      qgs->SetTitle("QGSJetII-04");
      dpm->SetTitle("DPMJet 3.0-6");
      sl1->SetTitle("#gamma-p (STARLIGHT+DPMJet)");
      sl2->SetTitle("#gamma-p (STARLIGHT+Pythia)");


      // zb->GetXaxis()->SetLimits(zb->GetBinLowEdge(3),zb->GetBinLowEdge(zb->GetNbinsX())); //cut away first bin
      // noise->GetXaxis()->SetLimits(2,noise->GetBinLowEdge(noise->GetNbinsX())); //cut away first bin
      // hijing->GetXaxis()->SetLimits(hijing->GetBinLowEdge(3),hijing->GetBinLowEdge(hijing->GetNbinsX())); //cut away first bin
      // epos->GetXaxis()->SetLimits(epos->GetBinLowEdge(3),epos->GetBinLowEdge(epos->GetNbinsX())); //cut away first bin
      // qgs->GetXaxis()->SetLimits(qgs->GetBinLowEdge(3),qgs->GetBinLowEdge(qgs->GetNbinsX())); //cut away first bin
      // sl1->GetXaxis()->SetLimits(3,sl1->GetBinLowEdge(sl1->GetNbinsX())); //cut away first bin
      // sl2->GetXaxis()->SetLimits(4,sl2->GetBinLowEdge(sl2->GetNbinsX())); //cut away first bin

      hijing->GetXaxis()->SetRange(2,hijing->GetNbinsX()*(type[n]=="double"?0.5:0.75));
      noise->GetXaxis()->SetRange(2,noise->GetNbinsX()*(type[n]=="double"?0.5:0.75));

      hijing->GetYaxis()->SetRangeUser(0.8,1.001);
      noise->GetYaxis()->SetRangeUser(2e-5,100);//type[n]=="double"?1e-5:1e-5,1.01);

      hijing->GetXaxis()->SetTitle("E_{HF} [GeV]");
      hijing->GetYaxis()->SetTitle("efficiency");
      hijing->GetXaxis()->SetLabelSize(hijing->GetXaxis()->GetLabelSize()*1.2);
      hijing->GetYaxis()->SetLabelSize(hijing->GetYaxis()->GetLabelSize()*1.2);
      hijing->GetXaxis()->SetTitleSize(hijing->GetXaxis()->GetTitleSize()*1.1);
      hijing->GetYaxis()->SetTitleSize(hijing->GetYaxis()->GetTitleSize()*1.1);
      hijing->GetXaxis()->SetTitleOffset(hijing->GetXaxis()->GetTitleOffset()*1.1);
      hijing->GetYaxis()->SetTitleOffset(hijing->GetYaxis()->GetTitleOffset()*1.1);

      const double cut_value = type[n]=="single"?cut_value_single:cut_value_double;

      if(draw)
        {
          TCanvas* can1 = new TCanvas;
          hijing->Draw("HIST L");
          epos->Draw("HIST L SAME");
          qgs->Draw("HIST L SAME");
          dpm->Draw("HIST L SAME");

          TLine* line = new TLine(cut_value,type[n]=="single"?0.8:0.88,cut_value,1.001);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          TLegend* leg = new TLegend(0.1,0.1,0.2,0.2);
          if(type[n]=="single")
            {
              leg->SetX1(0.22);
              leg->SetX2(0.56);
              leg->SetY1(0.35);
              leg->SetY2(0.58);
#ifdef __CINT__
              CMSText(1,1,0,"single-arm selection");
#endif
            }
          if(type[n]=="double")
            {
              leg->SetX1(0.42);
              leg->SetX2(0.75);
              leg->SetY1(0.18);
              leg->SetY2(0.41);
#ifdef __CINT__
              CMSText(1,0,1,"double-arm selection");
#endif
            }

          leg->AddEntry(hijing,"","l");
          leg->AddEntry(epos,"","l");
          leg->AddEntry(qgs,"","l");
          leg->AddEntry(dpm,"","l");
#ifdef __CINT__
          SetLegAtt(leg,1.3);
#endif
          leg->Draw();

          zb->GetYaxis()->SetRangeUser(0,1.01);
          zb->GetXaxis()->SetTitle("E_{HF} [GeV]");
          zb->GetYaxis()->SetTitle("efficiency");

          can1->SaveAs((string("plots/full_p_space_eff_PAS_")+type[n]+string(".pdf")).c_str());
        }

      //////////////////////////////////////////////////////////////////////////////////////////////
      //LOOPING AND OUTPUT//

      double f_eme = 0;
      double f_mce = 0;
      double f_mcesys = 0;
      int count = 0;
      for(int i=zb->FindBin(2); i<=zb->FindBin(10); i++)
        {
          std::set<double> mc_values;
          mc_values.insert(epos->GetBinContent(i));
          mc_values.insert(qgs->GetBinContent(i));
          //mc_values.insert(dpm->GetBinContent(i));

          ++count;
          f_eme += fabs(sl1->GetBinContent(i) - sl2->GetBinContent(i));
          f_mce += fabs(*mc_values.rbegin() - *mc_values.begin());
          f_mcesys += fabs(getEffDiffWeight("Epos", cut_value, _diff_epos_reweight) - getEffDiffWeight("QGSJetII", cut_value, _diff_qgs_reweight));
        }
      f_eme /= double(count);
      f_mce /= double(count);
      f_mcesys /= double(count);


      for(int i=1; i<=zb->GetNbinsX(); i++)
        {
          std::set<double> mc_values;
          mc_values.insert(epos->GetBinContent(i));
          mc_values.insert(qgs->GetBinContent(i));
          //mc_values.insert(dpm->GetBinContent(i));

          const double f_em     = 0.5 * (sl1->GetBinContent(i) + sl2->GetBinContent(i));
          const double f_mc     = std::accumulate(mc_values.begin(),mc_values.end(),0.0) / double(mc_values.size());
          //const double f_mc     = 1./2. * (epos->GetBinContent(i) + qgs->GetBinContent(i));
          const double f_mcsys  = 0.5 * (getEffDiffWeight("Epos", cut_value, _diff_epos_reweight) + getEffDiffWeight("QGSJetII", cut_value, _diff_qgs_reweight));
          const double f_noise  = noise->GetBinContent(i)/noise->GetBinContent(1);
          const double n_sel_zb = zb->GetBinContent(i);
          //const double n_zb     = 1;

          ///////////////////////////////////////////////////////////////////////
          //Number of events
          const double n_zb = (events_integral/lumi_integral);
          const double n_noise = f_noise * n_zb;
          const double n_em = f_em * 0.195; //these are not n but already n/lumi
          if(i!=1)
            {
              noise->SetBinContent(i,n_noise);
              sl1->SetBinContent(i,sl1->GetBinContent(i)*0.195);
              sl2->SetBinContent(i,sl2->GetBinContent(i)*0.195);
            }
          if(i==zb->FindBin(cut_value_single) && type[n]==string("single"))
            {
              // if((i-1)%10 == 0)
              cout
                << setprecision(3)
                << endl << i << "(" << zb->GetBinCenter(i) << ")"
                << endl << "f_mc= " << f_mc << " ± " << f_mce << " ( " << f_mce/f_mc*100. << "%)"
                << endl << "f_mc_scaled= " << f_mcsys << " ± " << f_mcesys << " ( " << f_mcesys/f_mcsys*100. << "%)"
                << endl << "f_em= " << f_em << " ± " << f_eme << " ( " << f_eme/f_em*100. << "%)"
                << endl << "n_em= " << n_em << " ± " << f_eme/f_em*n_em << " ( " << f_eme/f_em*100. << "%)"
                << endl << "for f_noise and n_noise consult makePlots_cs.C" << endl << endl;

              corr_fac_em[0] = f_em;
              corr_fac_mc[0] = f_mc;
              corr_fac_mc[2] = f_mcsys;
              corr_fac_eme[0] = f_eme;
              corr_fac_mce[0] = f_mce;
              corr_fac_mce[2] = f_mcesys;

              corr_fac_epos[0] = epos->GetBinContent(i);
              corr_fac_qgsjet[0] = qgs->GetBinContent(i);
              corr_fac_hijing[0] = hijing->GetBinContent(i);
              corr_fac_dpmjet[0] = dpm->GetBinContent(i);
            }
          if(i==zb->FindBin(cut_value_double) && type[n]==string("double"))
            {
              cout
                << setprecision(3)
                << endl << i << "(" << zb->GetBinCenter(i) << ")"
                << endl << "f_mc= " << f_mc << " ± " << f_mce << " ( " << f_mce/f_mc*100. << "%)"
                << endl << "f_mc_scaled= " << f_mcsys << " ± " << f_mcesys << " ( " << f_mcesys/f_mcsys*100. << "%)"
                << endl << "f_em= " << f_em << " ± " << f_eme << " ( " << f_eme/f_em*100. << "%)"
                << endl << "n_em= " << n_em << " ± " << f_eme/f_em*n_em << " ( " << f_eme/f_em*100. << "%)"
                << endl << "for f_noise and n_noise consult makePlots_cs.C" << endl << endl;

              corr_fac_em[1] = f_em;
              corr_fac_mc[1] = f_mc;
              corr_fac_eme[1] = f_eme;
              corr_fac_mce[1] = f_mce;
              corr_fac_mc[3] = f_mcsys;
              corr_fac_mce[3] = f_mcesys;

              corr_fac_epos[1] = epos->GetBinContent(i);
              corr_fac_qgsjet[1] = qgs->GetBinContent(i);
              corr_fac_hijing[1] = hijing->GetBinContent(i);
              corr_fac_dpmjet[1] = dpm->GetBinContent(i);
            }
        }

      //////////////////////////////////////////////////////////////////////////////////////////////

      if(draw)
        {
          TCanvas* can2 = new TCanvas;
          noise->GetXaxis()->SetTitle("E_{HF} [GeV]");
          noise->GetYaxis()->SetTitle("N/#it{L} [b]");
          noise->GetXaxis()->SetLabelSize(noise->GetXaxis()->GetLabelSize()*1.2);
          noise->GetYaxis()->SetLabelSize(noise->GetYaxis()->GetLabelSize()*1.2);
          noise->GetXaxis()->SetTitleSize(noise->GetXaxis()->GetTitleSize()*1.1);
          noise->GetYaxis()->SetTitleSize(noise->GetYaxis()->GetTitleSize()*1.1);
          noise->GetXaxis()->SetTitleOffset(noise->GetXaxis()->GetTitleOffset()*1.1);
          noise->GetYaxis()->SetTitleOffset(noise->GetYaxis()->GetTitleOffset()*1.0);
          noise->Draw("HIST l");
          sl1->Draw("HIST l SAME");
          sl2->Draw("HIST l SAME");
          TLegend* leg2 = new TLegend(0.1,0.1,0.2,0.2);

          if(type[n]=="single")
            {
              leg2->SetX1(0.22);
              leg2->SetY1(0.32);
              leg2->SetX2(0.56);
              leg2->SetY2(0.52);
#ifdef __CINT__
              CMSText(1,1,0,"single-arm selection");
#endif
            }
          if(type[n]=="double")
            {
              leg2->SetX1(0.42);
              leg2->SetX2(0.93);
              leg2->SetY1(0.60);
              leg2->SetY2(0.80);
#ifdef __CINT__
              CMSText(1,0,1,"double-arm selection");
#endif
            }
          leg2->Draw();

          leg2->AddEntry(noise,"","l");
          leg2->AddEntry(sl1,"","l");
          leg2->AddEntry(sl2,"","l");
#ifdef __CINT__
          SetLegAtt(leg2,1.1);
#endif

          TLine* line = new TLine(cut_value,type[n]=="single"?1e-2:0,cut_value,type[n]=="single"?100:2e-2);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          can2->SetLogy();
          can2->SaveAs((string("plots/starlight_eff_")+type[n]+string(".pdf")).c_str());

          //file->Close();
          //file2->Close();
        }
    }

  TFile outfile("plots/corr_factors.root","recreate");
  corr_fac_em.Write("corr_fac_em");
  corr_fac_mc.Write("corr_fac_mc");
  corr_fac_eme.Write("corr_fac_eme");
  corr_fac_mce.Write("corr_fac_mce");
  corr_fac_epos.Write("corr_fac_epos");
  corr_fac_qgsjet.Write("corr_fac_qgsjet");
  corr_fac_hijing.Write("corr_fac_hijing");
  corr_fac_dpmjet.Write("corr_fac_dpmjet");
  outfile.Close();
}

double getEffDiffWeight(string model, double cut, double diffWeight)
{
  double all =
    histoman[model]["SD1"]->GetBinContent(1) +
    histoman[model]["SD2"]->GetBinContent(1) +
    histoman[model]["DD"]->GetBinContent(1) +
    histoman[model]["CD"]->GetBinContent(1) +
    histoman[model]["ND"]->GetBinContent(1);

  double diff_frac = (
    histoman[model]["SD1"]->GetBinContent(1) +
    histoman[model]["SD2"]->GetBinContent(1) +
    histoman[model]["DD"]->GetBinContent(1) +
    histoman[model]["CD"]->GetBinContent(1) ) /
    all;

  double eff_diff = (
    histoman[model]["SD1"]->Interpolate(cut) +
    histoman[model]["SD2"]->Interpolate(cut) +
    histoman[model]["DD"]->Interpolate(cut) +
    histoman[model]["CD"]->Interpolate(cut) ) /
    all;

  double eff_nd =
    histoman[model]["ND"]->Interpolate(cut) /
    all;

  double mod_eff =
    (diffWeight * eff_diff + eff_nd) / (1.+diff_frac*(diffWeight-1.));

  return mod_eff;
}
