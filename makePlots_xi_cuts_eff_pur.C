//make plots that show the dependence on gen xi cuts. 1 plot shows eff and pur as function of xi. the other combines all models and shows the same

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
#include <sstream>
#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <iterator>     // std::distance

void makePlots_xi_cuts_eff_pur(bool draw=1, string filename = "histos_new.root");


void makePlots_xi_cuts_eff_pur(bool draw, string filename)
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif
  set<double> cut_value_single; cut_value_single.insert(0.5);
  set<double> cut_value_double; cut_value_double.insert(0.5);
  set<double> cut_value2_single; cut_value2_single.insert(0.4); cut_value2_single.insert(0.6);
  set<double> cut_value2_double; cut_value2_double.insert(0.4); cut_value2_double.insert(0.6);

  vector<string> type; type.push_back("single"); type.push_back("double");
  TCanvas* canvases[2][2];
  canvases[0][0] = new TCanvas ; canvases[1][0] = new TCanvas; canvases[0][1] = new TCanvas ; canvases[1][1] = new TCanvas;
  TLegend* legends[2]; legends[0] = new TLegend(0.4,0.28,0.74,0.53) ; legends[1] = new TLegend(0.4,0.28,0.74,0.53);
  TVectorD corr_fac_had(2);
  TVectorD corr_fac_had_xi_epos(2);
  TVectorD corr_fac_had_xi_phythiamonash(2);
  TVectorD corr_fac_had_xi_dpmjet(2);
  TVectorD corr_fac_had_xi_hijing(2);
  TVectorD corr_fac_had_xi_qgsjet(2);
  TVectorD had_xi_value(2);
  TVectorD corr_fac_had_e(2);

  for(int n=0; n<int(type.size()); n++)
    TFile* file = TFile::Open(filename.c_str()); //cannot move outside of loop
    {
      int nBins = 0;

      vector<string> models;
      vector<int> colors;
      vector<double> styles;
      vector<string> titles;

      models.push_back("Epos");   colors.push_back(kBlue);    styles.push_back(25); titles.push_back("EPOS-LHC");
      models.push_back("PythiaMonash"); colors.push_back(kMagenta); styles.push_back(22);titles.push_back("Pythia8 (Monash)");
      models.push_back("PythiaMBR"); colors.push_back(kGreen+2); styles.push_back(4); titles.push_back("Pythia8 (MBR)");
      models.push_back("PythiaZ2Star"); colors.push_back(kOrange); styles.push_back(6); titles.push_back("Pythia6 (Z2*)");
      models.push_back("QGSJetII"); colors.push_back(kRed);   styles.push_back(28);titles.push_back("QGSJetII-04");
      //models.push_back("Starlight_DPMJet"); colors.push_back(kGray); styles.push_back(22);titles.push_back("Starlight dpm"); //not needed i think

      TH1D* interpolateEff = 0;
      TH1D* interpolatePur = 0;
      const int imax = int(models.size());

      //loop for calculating p and epsilon for each model
      for (int i = 0; i<imax; ++i)
        {
          TH1D* histeff = (TH1D*)file->Get(((models[i] + "/" + models[i] + "_h_mc_xi_") + type[n] + "_cut").c_str());
          TH1D* histpur = (TH1D*)file->Get(((models[i] + "/" + models[i] + "_h_mc_xi_") + type[n] + "_bg").c_str());
          cout << "h" << histeff << " " << histpur << endl;
          histeff = (TH1D*)histeff->Clone((string("hist_eff_")+models[i]).c_str()); //original histos needed for later
          histpur = (TH1D*)histpur->Clone((string("hist_pur_")+models[i]).c_str());

          if (!nBins) nBins = histeff->GetNbinsX();

          if(i == 0)
            {
              ostringstream name;
              name.str(""); name << "interpolateEff_" << n;
              interpolateEff = new TH1D(name.str().c_str(),"",histeff->GetNbinsX(), histeff->GetBinLowEdge(1), histeff->GetBinLowEdge(histeff->GetNbinsX()));
              name.str(""); name << "interpolatePur_" << n;
              interpolatePur = new TH1D(name.str().c_str(),"",histpur->GetNbinsX(), histpur->GetBinLowEdge(1), histpur->GetBinLowEdge(histpur->GetNbinsX()));
            }

          for (int bin = 2; bin < histeff->GetNbinsX(); ++bin)
            {
              //this looks a bit strange but is correct when you see what is in those histograms
              double effi = histeff->GetBinContent(bin) / (histeff->GetBinContent(bin) + histpur->GetBinContent(bin)); //N(VIS | HADLEVEL) / N(HADLEVEL)
              double puri = histeff->GetBinContent(bin) / double(histeff->GetBinContent(1)); // N(VIS | HADLEVEL) / N(VIS)
              histpur->SetBinContent(bin, puri);
              histpur->SetBinError(bin, 0); //WARNING change
              histeff->SetBinContent(bin, effi);
              histeff->SetBinError(bin, 0); //WARNING change
              interpolateEff->SetBinContent(bin,interpolateEff->GetBinContent(bin) + effi/double(imax)); //average
              interpolatePur->SetBinContent(bin,interpolatePur->GetBinContent(bin) + puri/double(imax)); // average
              interpolateEff->SetBinError(bin,0.);
              interpolatePur->SetBinError(bin,0.);
            }
          //histeff->Scale(1./double(histeff->GetBinContent(1)));


        }

      //loop for combining all four models
      for (int i = 0; i<imax; ++i)
        {
          TH1D* histeff = (TH1D*)file->Get((string("hist_eff_")+models[i]).c_str()); //see above the Clone(). but gets written to file
          TH1D* histpur = (TH1D*)file->Get((string("hist_pur_")+models[i]).c_str());

          for (int bin = 1; bin < histeff->GetNbinsX(); ++bin)
            {
              double deltaeff = fabs(interpolateEff->GetBinContent(bin) - histeff->GetBinContent(bin)); //inter.. contains average. histeff is eff of curr model
              double deltapur = fabs(interpolatePur->GetBinContent(bin) - histpur->GetBinContent(bin));
              interpolateEff->SetBinError(bin,interpolateEff->GetBinError(bin) + pow(deltaeff,2)/double(imax-1)); //calculate variance (sum (x-mu)**2)/(n-1)
              interpolatePur->SetBinError(bin,interpolatePur->GetBinError(bin) + pow(deltapur,2)/double(imax-1)); //imax-1 = n-1 (unbiased)
            }

          histeff->SetLineWidth(3);
          histpur->SetLineWidth(3);

          histeff->SetLineColor(colors[i]);
          histpur->SetLineColor(colors[i]);

          histeff->SetMarkerColor(histeff->GetLineColor());
          histpur->SetMarkerColor(histpur->GetLineColor());

          histeff->SetLineStyle(1);
          histpur->SetLineStyle(5);

          histeff->SetTitle(titles[i].c_str());
          histpur->SetTitle(titles[i].c_str());

          histeff->GetXaxis()->SetRange(2,histeff->FindBin(1));
          histpur->GetXaxis()->SetRange(2,histeff->FindBin(1));

          histeff->GetYaxis()->SetRangeUser(0.6,1.001);
          histpur->GetYaxis()->SetRangeUser(0.6,1.001);

          histeff->GetXaxis()->SetTitle("log_{10}(#xi_{SD}^{th})");
          histpur->GetXaxis()->SetTitle("log_{10}(#xi_{SD}^{th})");

          histeff->GetYaxis()->SetTitle("");
          histpur->GetYaxis()->SetTitle("");
          /*
            b->GetXaxis()->SetLabelSize(b->GetXaxis()->GetLabelSize()*1.2);
            b->GetYaxis()->SetLabelSize(b->GetYaxis()->GetLabelSize()*1.2);
            b->GetXaxis()->SetTitleSize(b->GetXaxis()->GetTitleSize()*1.1);
            b->GetYaxis()->SetTitleSize(b->GetYaxis()->GetTitleSize()*1.1);
            b->GetXaxis()->SetTitleOffset(b->GetXaxis()->GetTitleOffset()*1.1);
            b->GetYaxis()->SetTitleOffset(b->GetYaxis()->GetTitleOffset()*1.1);
          */
          canvases[n][0]->cd();

          histeff->Draw(i==0?"HIST L":"SAME HIST L");
          histpur->Draw("SAME L");
          TLegend* leg = legends[n];
          if(i==0)
            {
#ifdef __CINT__
              CMSText(3,1,0,type[n]=="single"?"Single-arm":"Double-arm");
#endif
              helper1 = (TH1D*)histeff->Clone("helper1"); helper1->SetLineColor(kBlack);
              helper2 = (TH1D*)histpur->Clone("helper2"); helper2->SetLineColor(kBlack);
              leg->AddEntry(helper1, "efficiency","l");
              leg->AddEntry(helper2, "purity","l");
            }
          leg->AddEntry(histeff,titles[i].c_str(),"l");

          if(i==imax-1)
            {
#ifdef __CINT__
              SetLegAtt(leg,1.1);
#endif
              leg->Draw("");
            }

        } // model

      //get std dev
      for (int bin = 1; bin < nBins; ++bin)
        {
          interpolateEff->SetBinError(bin,sqrt(interpolateEff->GetBinError(bin)));
          interpolatePur->SetBinError(bin,sqrt(interpolatePur->GetBinError(bin)));
        }

      canvases[n][0]->SaveAs((string("plots/hadron_xicuts_eff_pur_allmodels_") + type[n] + string(".eps")).c_str());

      ////////////////averaged plots
      interpolateEff->SetLineWidth(3);
      interpolatePur->SetLineWidth(3);

      interpolateEff->SetLineColor(kBlack);
      interpolatePur->SetLineColor(kBlue);

      interpolateEff->SetMarkerColor(interpolateEff->GetLineColor());
      interpolatePur->SetMarkerColor(interpolatePur->GetLineColor());

      interpolateEff->SetFillColor(interpolateEff->GetLineColor());
      interpolatePur->SetFillColor(interpolatePur->GetLineColor());

      interpolateEff->SetLineStyle(1);
      interpolatePur->SetLineStyle(5);

      interpolateEff->SetFillStyle(3003);
      interpolatePur->SetFillStyle(3003);

      interpolateEff->SetTitle("efficiency");
      interpolatePur->SetTitle("purity");

      interpolateEff->GetXaxis()->SetRange(2,interpolateEff->FindBin(1));
      interpolatePur->GetXaxis()->SetRange(2,interpolateEff->FindBin(1));

      interpolateEff->GetYaxis()->SetRangeUser(0.6,1.001);
      interpolatePur->GetYaxis()->SetRangeUser(0.6,1.001);

      interpolateEff->GetXaxis()->SetTitle("log_{10}(#xi_{SD}^{th})");
      interpolatePur->GetXaxis()->SetTitle("log_{10}(#xi_{SD}^{th})");

      interpolateEff->GetYaxis()->SetTitle("");
      interpolatePur->GetYaxis()->SetTitle("");

      canvases[n][1]->cd();
      interpolateEff->Draw("E3");
      interpolatePur->Draw("E3 SAME");

      TLegend* leg = new TLegend(0.21,0.33,0.55,0.42);
#ifdef __CINT__
      CMSText(3,1,0,type[n]=="single"?"Single-arm":"Double-arm");
#endif
      leg->AddEntry(interpolateEff,"efficiency","FP");
      leg->AddEntry(interpolatePur,"purity","FP");

#ifdef __CINT__
      SetLegAtt(leg,1.1);
#endif
      leg->Draw();

      ////Find xi cut

      double curreff=1.;
      double currpur=0;
      const double step = 0.01;
      for(double j=-8; j<-1; j+=step)
        {
          curreff = interpolateEff->Interpolate(j);
          currpur = interpolatePur->Interpolate(j);
          double currfactor = (interpolatePur->Interpolate(j)) / interpolateEff->Interpolate(j);
          double currerror = sqrt(pow(interpolateEff->GetBinError(interpolateEff->FindBin(j)),2) + pow(interpolatePur->GetBinError(interpolateEff->FindBin(j)),2)); //correct for rel error
          if(curreff>currpur) //eff crossing purity
            {
              if (fabs(corr_fac_had[n] - 1) > fabs(currfactor-1)) // if new factor closer to unity
                {
                  had_xi_value[n] = j;
                  corr_fac_had[n] = currfactor;
                  corr_fac_had_e[n] = currerror;
                }
              break;
            }
          had_xi_value[n] = j - step;
          corr_fac_had[n] = currfactor;
          corr_fac_had_e[n] = currerror;
        }
      cout << "Found xi_cut for type " << n << ": xi=" << had_xi_value[n] << " 1-Pur / eff = " << (interpolatePur->Interpolate(j)) << "/" << interpolateEff->Interpolate(j) << " = " << corr_fac_had[n] << " pm " << corr_fac_had_e[n] << endl;

      //loop to save the efficiency according to had lvl cut for each model
   //    for (int i = 0; i<imax; ++i)
   //      {
      //     TVectorD *model_fac = 0;
      //     if(models[i]=="Epos")
      //       model_fac = &corr_fac_had_xi_epos;
      //     if(models[i]=="QGSJetII")
      //       model_fac = &corr_fac_had_xi_qgsjet;
      //     if(models[i]=="PythiaMBR")
      //       model_fac = &corr_fac_had_xi_hijing;
      //     if(models[i]=="DPMJet")
      //       model_fac = &corr_fac_had_xi_dpmjet;
      //     if(model_fac==0)
      // {
      //   cerr << "cannot find model_fac for model" << endl;
      //       exit(-1);
      // }
      //     TH1D* h1=(TH1D*)file->Get(string(models[i]+string("/")+models[i]+string("_h_mc_xi_")+type[n]+string("_cut")).c_str());
      //     TH1D* h2=(TH1D*)file->Get(string(models[i]+string("/")+models[i]+string("_h_mc_xi_")+type[n]+string("_bg")).c_str());
      //     (*model_fac)[n] = (h1->Interpolate(had_xi_value[n])+h2->Interpolate(had_xi_value[n])) / (h1->GetBinContent(1)+h2->GetBinContent(1));
      //     cout << models[i] << " " << " " << had_xi_value[n] << " " <<  (*model_fac)[n] << endl;
 //        }

      TLine* line = new TLine(had_xi_value[n],0.6,had_xi_value[n],1);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->Draw("SAME");

      canvases[n][1]->SaveAs((string("plots/hadron_xicuts_eff_pur_combined_") + type[n] + string(".eps")).c_str());

    } // type

  // TFile outfile("plots/corr_factors_hadron.root","recreate");
  // corr_fac_had.Write("corr_fac_had");
  // had_xi_value.Write("had_xi_value");
  // corr_fac_had_e.Write("corr_fac_had_e");
  // corr_fac_had_xi_epos.Write("corr_fac_had_epos");
  // corr_fac_had_xi_dpmjet.Write("corr_fac_had_dpmjet");
  // corr_fac_had_xi_hijing.Write("corr_fac_had_hijing");
  // corr_fac_had_xi_qgsjet.Write("corr_fac_had_qgsjet");
  // outfile.Close();

}
