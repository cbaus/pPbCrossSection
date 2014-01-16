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
#include <iterator>     // std::distance

#define _LumiCorrpPb 1.142 //only use if trees don't contain vdm calibration factor
#define _LumiCorrPbp 1.138

TVectorD corr_fac_em(2);
TVectorD corr_fac_eme(2);
TVectorD corr_fac_mc(4);
TVectorD corr_fac_mce(4);
TVectorD corr_fac_epos(2);
TVectorD corr_fac_qgsjet(2);

void makePlots_cs_purvseff(bool draw=1, string filename = "histos.root");


void makePlots_cs_purvseff(bool draw, string filename)
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  //available cuts single 5, 6.4, 7.2, 8, 8.8, 9.6, 10, 15, 20
  //available cuts double 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5
  set<double> cut_value_single; cut_value_single.insert(8);
  set<double> cut_value_double; cut_value_double.insert(4);
  set<double> cut_value2_single; cut_value2_single.insert(6.4); cut_value2_single.insert(10);
  set<double> cut_value2_double; cut_value2_double.insert(3); cut_value2_double.insert(5);

  vector<string> type; type.push_back("single"); type.push_back("double");
  TGraph* pur_single = new TGraph(0);
  TGraph* pur_double = new TGraph(0);
  TGraph* pur_selec_single = new TGraph(cut_value_single.size());
  TGraph* pur_selec_double = new TGraph(cut_value_single.size());
  TGraph* pur_selec2_single = new TGraph(cut_value2_single.size());
  TGraph* pur_selec2_double = new TGraph(cut_value2_single.size());

  for(int n=0; n<int(type.size()); n++)
    {
      TFile* file = TFile::Open(filename.c_str());
      TFile* file2 = TFile::Open("histos_test3.root");

      TH1D* a=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n]).c_str());
      TH1D* a2=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n] + string("_noise")).c_str());
      TH1D* h_events=(TH1D*)file->Get(string("data210885/data210885_h_lumi").c_str()); //zb
      TH1D* h_lumi=(TH1D*)file->Get(string("data210885/data210885_h_run_events_lumi").c_str());
      TH1D* eposrew=(TH1D*)file2->Get((string("EposDiffWeightOpt/EposDiffWeightOpt_h_hf_cut_") + type[n]).c_str());
      TH1D* qgsrew=(TH1D*)file2->Get((string("QGSJetIIDiffWeightOpt/QGSJetIIDiffWeightOpt_h_hf_cut_") + type[n]).c_str());
      TH1D* b=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_cut_") + type[n]).c_str());
      TH1D* c=(TH1D*)file->Get((string("Epos/Epos_h_hf_cut_")+ type[n]).c_str());
      TH1D* d=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_cut_") + type[n]).c_str());
      TH1D* e=(TH1D*)file->Get((string("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_") + type[n]).c_str());
      TH1D* f=(TH1D*)file->Get((string("Starlight_Pythia/Starlight_Pythia_h_hf_cut_") + type[n]).c_str());

      TGraph* pur = type[n]==string("single")?pur_single:pur_double;
      TGraph* pur_selec = type[n]==string("single")?pur_selec_single:pur_selec_double;
      TGraph* pur_selec2 = type[n]==string("single")?pur_selec2_single:pur_selec2_double;

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

      //a->Scale(1./double(a->GetBinContent(1)));
      //a2->Scale(1./double(a2->GetBinContent(1)));
      b->Scale(1./double(b->GetBinContent(1)));
      c->Scale(1./double(c->GetBinContent(1)));
      d->Scale(1./double(d->GetBinContent(1)));
      e->Scale(1./double(e->GetBinContent(1)));
      f->Scale(1./double(f->GetBinContent(1))/ 195. * 122.); //cross section starlight samples
      eposrew->Scale(1./double(eposrew->GetBinContent(1)));
      qgsrew->Scale(1./double(qgsrew->GetBinContent(1)));

      a->SetLineWidth(3);
      a2->SetLineWidth(3);
      b->SetLineWidth(3);
      c->SetLineWidth(3);
      d->SetLineWidth(3);
      e->SetLineWidth(3);
      f->SetLineWidth(3);
      eposrew->SetLineWidth(3);
      qgsrew->SetLineWidth(3);
      a2->SetLineColor(kBlack);
      b->SetLineColor(kGreen+2);
      c->SetLineColor(kBlue);
      d->SetLineColor(kRed);
      e->SetLineColor(kRed);
      f->SetLineColor(kBlue);
      eposrew->SetLineColor(kMagenta);
      qgsrew->SetLineColor(kMagenta);
      a2->SetMarkerColor(a2->GetLineColor());
      b->SetMarkerColor(b->GetLineColor());
      c->SetMarkerColor(c->GetLineColor());
      d->SetMarkerColor(d->GetLineColor());
      e->SetMarkerColor(e->GetLineColor());
      f->SetMarkerColor(f->GetLineColor());
      eposrew->SetMarkerColor(f->GetLineColor());
      qgsrew->SetMarkerColor(f->GetLineColor());
      b->SetLineStyle(7);
      d->SetLineStyle(9);
      f->SetLineStyle(9);
      eposrew->SetLineStyle(10);
      qgsrew->SetLineStyle(10);

      a->SetTitle("Data");
      a2->SetTitle("Noise");
      eposrew->SetTitle("EPOS-LHC (#sigma_{diff}x1.12)");
      qgsrew->SetTitle("QGSJETII-04 (#sigma_{diff}x1.50)");
      b->SetTitle("Hijing 1.383");
      c->SetTitle("EPOS-LHC");
      d->SetTitle("QGSJetII-04");
      e->SetTitle("#gamma-p (STARLIGHT+DPMJet)");
      f->SetTitle("#gamma-p (STARLIGHT+Pythia)");


      // a->GetXaxis()->SetLimits(a->GetBinLowEdge(3),a->GetBinLowEdge(a->GetNbinsX())); //cut away first bin
      // a2->GetXaxis()->SetLimits(2,a2->GetBinLowEdge(a2->GetNbinsX())); //cut away first bin
      // b->GetXaxis()->SetLimits(b->GetBinLowEdge(3),b->GetBinLowEdge(b->GetNbinsX())); //cut away first bin
      // c->GetXaxis()->SetLimits(c->GetBinLowEdge(3),c->GetBinLowEdge(c->GetNbinsX())); //cut away first bin
      // d->GetXaxis()->SetLimits(d->GetBinLowEdge(3),d->GetBinLowEdge(d->GetNbinsX())); //cut away first bin
      // e->GetXaxis()->SetLimits(3,e->GetBinLowEdge(e->GetNbinsX())); //cut away first bin
      // f->GetXaxis()->SetLimits(4,f->GetBinLowEdge(f->GetNbinsX())); //cut away first bin

      b->GetXaxis()->SetRange(2,b->GetNbinsX()*(type[n]=="double"?0.5:0.75));
      a2->GetXaxis()->SetRange(2,a2->GetNbinsX()*(type[n]=="double"?0.5:0.75));

      b->GetYaxis()->SetRangeUser(0.8,1.001);
      a2->GetYaxis()->SetRangeUser(2e-5,100);//type[n]=="double"?1e-5:1e-5,1.01);

      b->GetXaxis()->SetTitle("E_{HF} [GeV]");
      b->GetYaxis()->SetTitle("efficiency");
      b->GetXaxis()->SetLabelSize(b->GetXaxis()->GetLabelSize()*1.2);
      b->GetYaxis()->SetLabelSize(b->GetYaxis()->GetLabelSize()*1.2);
      b->GetXaxis()->SetTitleSize(b->GetXaxis()->GetTitleSize()*1.1);
      b->GetYaxis()->SetTitleSize(b->GetYaxis()->GetTitleSize()*1.1);
      b->GetXaxis()->SetTitleOffset(b->GetXaxis()->GetTitleOffset()*1.1);
      b->GetYaxis()->SetTitleOffset(b->GetYaxis()->GetTitleOffset()*1.1);


      //////////////////////////////////////////////////////////////////////////////////////////////
      //LOOPING AND OUTPUT//

      double f_eme = 0;
      double f_mce = 0;
      double f_mcesys = 0;
      int count = 0;
      for(int i=a->FindBin(2); i<=a->FindBin(10); i++)
        {
          ++count;
          f_eme += fabs(e->GetBinContent(i) - f->GetBinContent(i));
          f_mce += fabs(c->GetBinContent(i) - d->GetBinContent(i));
          f_mcesys += fabs(eposrew->GetBinContent(i) - qgsrew->GetBinContent(i));
        }
      f_eme /= double(count);
      f_mce /= double(count);
      f_mcesys /= double(count);

      pur_single->Expand(a->GetNbinsX());
      pur_double->Expand(a->GetNbinsX());

      for(int i=1; i<=a->GetNbinsX(); i++)
        {
          const double f_em     = 0.5 * (e->GetBinContent(i) + f->GetBinContent(i));
          const double f_mc     = 0.5 * (c->GetBinContent(i) + d->GetBinContent(i));
          const double f_mcsys  = 0.5 * (eposrew->GetBinContent(i) + qgsrew->GetBinContent(i));
          const double f_noise  = a2->GetBinContent(i)/a2->GetBinContent(1);
          const double n_sel_zb = a->GetBinContent(i);

          ///////////////////////////////////////////////////////////////////////
          //Number of events
          const double n_zb = (events_integral/lumi_integral);
          const double n_noise = f_noise * n_zb;
          const double n_em = f_em * 0.195; //these are not n but already n/lumi


          ///////////////////////////////////////////////////////////////////////
          //Purity
          double purity = 2.06 / (2.06 + n_noise + n_em);
          pur->SetPoint(i,purity,f_mc);

          if(i!=1)
            {
              a2->SetBinContent(i,n_noise);
              e->SetBinContent(i,e->GetBinContent(i)*0.195);
              f->SetBinContent(i,f->GetBinContent(i)*0.195);
            }

          set<double>& cut_value = type[n]==string("single")?cut_value_single:cut_value_double;
          set<double>& cut_value2 = type[n]==string("single")?cut_value2_single:cut_value2_double;
          bool selected = false;
          set<double>::iterator it = cut_value.begin();
          for (int point = 0; it != cut_value.end(); ++it, point++)
            if(i==a->FindBin(*it))
              {
                selected = true;
                pur_selec->SetPoint(point,purity,f_mc);
              }

          it = cut_value2.begin();
          for (int point = 0; it != cut_value2.end(); ++it, point++)
            if(i==a->FindBin(*it))
              {
                selected = true;
                pur_selec2->SetPoint(point,purity,f_mc);
              }

          if(selected)
          {
              // if((i-1)%10 == 0)
              cout
                << setprecision(3)
                << endl << i << "(" << a->GetBinCenter(i) << ")"
                << endl << "f_mc= " << f_mc << " ± " << f_mce << " ( " << f_mce/f_mc*100. << "%)"
                << endl << "f_mc_scaled= " << f_mcsys << " ± " << f_mcesys << " ( " << f_mcesys/f_mcsys*100. << "%)"
                << endl << "f_em= " << f_em << " ± " << f_eme << " ( " << f_eme/f_em*100. << "%)"
                << endl << "n_em= " << n_em << " ± " << f_eme/f_em*n_em << " ( " << f_eme/f_em*100. << "%)"
                << endl << "n_noise= " << n_noise << " ± ? ( " << f_eme/f_em*100. << "%)"
                << endl << "purity= " << purity
                << endl << "for f_noise and n_noise consult makePlots_cs.C" << endl << endl;
            }
        }
    }
      //////////////////////////////////////////////////////////////////////////////////////////////

  if(draw)
        {
          TCanvas* can2 = new TCanvas;
          pur_single->Draw("AC");
          double xmin = 0.81, xmax = 1.08, ymin = 0.81, ymax = 1.08;
          pur_single->GetYaxis()->SetRangeUser(ymin,ymax);
          pur_single->GetXaxis()->SetLimits(xmin,xmax);
          pur_single->GetXaxis()->SetNdivisions(505);
          pur_single->GetYaxis()->SetNdivisions(505);

          pur_double->Draw("C");
          pur_single->SetTitle(";purity;efficiency");
          pur_double->SetTitle(";purity;efficiency");
          pur_single->SetLineColor(kRed);
          pur_double->SetLineColor(kBlue);
          pur_single->SetLineWidth(3);
          pur_double->SetLineWidth(3);

          pur_selec_single->SetMarkerColor(pur_single->GetLineColor());
          pur_selec_double->SetMarkerColor(pur_double->GetLineColor());
          pur_selec_single->SetMarkerSize(1.8);
          pur_selec_double->SetMarkerSize(1.8);
          pur_selec_single->SetMarkerStyle(21);
          pur_selec_double->SetMarkerStyle(21);
          pur_selec_single->Draw("P");
          pur_selec_double->Draw("P");

          pur_selec2_single->SetMarkerColor(pur_single->GetLineColor());
          pur_selec2_double->SetMarkerColor(pur_double->GetLineColor());
          pur_selec2_single->SetMarkerSize(1.8);
          pur_selec2_double->SetMarkerSize(1.8);
          pur_selec2_single->Draw("P");
          pur_selec2_double->Draw("P");

          // a2->GetXaxis()->SetLabelSize(a2->GetXaxis()->GetLabelSize()*1.2);
          // a2->GetYaxis()->SetLabelSize(a2->GetYaxis()->GetLabelSize()*1.2);
          // a2->GetXaxis()->SetTitleSize(a2->GetXaxis()->GetTitleSize()*1.1);
          // a2->GetYaxis()->SetTitleSize(a2->GetYaxis()->GetTitleSize()*1.1);
          // a2->GetXaxis()->SetTitleOffset(a2->GetXaxis()->GetTitleOffset()*1.1);
          // a2->GetYaxis()->SetTitleOffset(a2->GetYaxis()->GetTitleOffset()*1.0);

          TLegend* leg2 = new TLegend(0.21,0.73,0.55,0.83);

              leg2->SetX1(0.21);
              leg2->SetY1(0.73);
              leg2->SetX2(0.55);
              leg2->SetY2(0.83);
#ifdef __CINT__
              CMSText(2,1,1);
#endif


          leg2->Draw();

          leg2->AddEntry(pur_single,"single-arm selection","l");
          leg2->AddEntry(pur_double,"double-arm selection","l");
#ifdef __CINT__
          SetLegAtt(leg2,1.1);
#endif

          TLine* line = new TLine(1,ymin,1,ymax);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");
          line = new TLine(xmin,1,xmax,1);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          for(int n=0; n<int(type.size()); n++)of
            {

              pur_selec = type[n]==string("single")?pur_selec_single:pur_selec_double;
              pur_selec2 = type[n]==string("single")?pur_selec2_single:pur_selec2_double;
              set<double>& cut_value = type[n]==string("single")?cut_value_single:cut_value_double;
              set<double>& cut_value2 = type[n]==string("single")?cut_value2_single:cut_value2_double;
              set<double>::iterator cut_ite;

              for (cut_ite = cut_value.begin(), int point = 0; cut_ite != cut_value.end(); ++cut_ite, point++)
                {
                  const double offsetx = type[n]==string("single")?0.01:0.02;
                  const double offsety = type[n]==string("single")?0.01:-0.005;

                  double x;
                  double y;
                  pur_selec->GetPoint(point,x,y);
                  x+=offsetx;
                  y+=offsety;
                  ostringstream string; string << "#color[" << (type[n]==string("single")?2:4) << "]{" << *cut_ite << " GeV}";
                  TLatex* text = new TLatex(x,y,string.str().c_str());
                  text->Draw();
                }
              for (cut_ite = cut_value2.begin(), int point = 0; cut_ite != cut_value2.end(); ++cut_ite, point++)
                {
                  const double offsetx = type[n]==string("single")?0.01:0.02;
                  const double offsety = type[n]==string("single")?0.01:-0.005;

                  double x;
                  double y;
                  pur_selec2->GetPoint(point,x,y);
                  x+=offsetx;
                  y+=offsety;
                  ostringstream string; string << "#color[" << (type[n]==string("single")?2:4) << "]{" << *cut_ite << " GeV}";
                  TLatex* text = new TLatex(x,y,string.str().c_str());
                  text->Draw();
                }
            }


          can2->SaveAs((string("plots/purity_vs_eff")+string(".pdf")).c_str());
        }

}
