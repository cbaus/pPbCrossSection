//makes plot for paper fig. 2. shows optimal working point for e_hf threshold

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
  set<double> cut_value_single; cut_value_single.insert(5);
  set<double> cut_value_double; cut_value_double.insert(2.5);
  set<double> cut_value2_single; cut_value2_single.insert(3); cut_value2_single.insert(5);
  set<double> cut_value2_double; cut_value2_double.insert(2); cut_value2_double.insert(3);

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
      TFile* file2 = TFile::Open("histos.root");

      TH1D* a=(TH1D*)file->Get((string("data247324/data247324_h_hf_cut_") + type[n]).c_str());
      TH1D* a2=(TH1D*)file->Get((string("data247324/data247324_h_hf_cut_") + type[n] + string("_noise")).c_str());
      TH1D* h_events=(TH1D*)file->Get(string("data247324/data247324_h_lumi").c_str()); //zb
      TH1D* h_lumi=(TH1D*)file->Get(string("data247324/data247324_h_run_events_lumi").c_str());
      //TH1D* eposrew=(TH1D*)file2->Get((string("EposDiffWeightOpt/EposDiffWeightOpt_h_hf_cut_") + type[n]).c_str());
      //TH1D* qgsrew=(TH1D*)file2->Get((string("QGSJetIIDiffWeightOpt/QGSJetIIDiffWeightOpt_h_hf_cut_") + type[n]).c_str());
      TH1D* b=(TH1D*)file->Get((string("PythiaMBR/PythiaMBR_h_hf_cut_") + type[n]).c_str());
      TH1D* c=(TH1D*)file->Get((string("Epos/Epos_h_hf_cut_")+ type[n]).c_str());
      TH1D* d=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_cut_") + type[n]).c_str());
      TH1D* e=(TH1D*)file->Get((string("PythiaZ2Star/PythiaZ2Star_h_hf_cut_") + type[n]).c_str());
      TH1D* f=(TH1D*)file->Get((string("PythiaMonash/PythiaMonash_h_hf_cut_") + type[n]).c_str());

      TGraph* pur = type[n]==string("single")?pur_single:pur_double;
      TGraph* pur_selec = type[n]==string("single")?pur_selec_single:pur_selec_double;
      TGraph* pur_selec2 = type[n]==string("single")?pur_selec2_single:pur_selec2_double;

      /////
      //LUMI
      double lumi_integral = 0;
      double events_integral = 0;
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
          const double lumicorr = _LumiCorrpPb; //might be different for other run than 247324
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
      f->Scale(1./double(f->GetBinContent(1)));
      //eposrew->Scale(1./double(eposrew->GetBinContent(1)));
      //qgsrew->Scale(1./double(qgsrew->GetBinContent(1)));




      // a->GetXaxis()->SetLimits(a->GetBinLowEdge(3),a->GetBinLowEdge(a->GetNbinsX())); //cut away first bin
      // a2->GetXaxis()->SetLimits(2,a2->GetBinLowEdge(a2->GetNbinsX())); //cut away first bin
      // b->GetXaxis()->SetLimits(b->GetBinLowEdge(3),b->GetBinLowEdge(b->GetNbinsX())); //cut away first bin
      // c->GetXaxis()->SetLimits(c->GetBinLowEdge(3),c->GetBinLowEdge(c->GetNbinsX())); //cut away first bin
      // d->GetXaxis()->SetLimits(d->GetBinLowEdge(3),d->GetBinLowEdge(d->GetNbinsX())); //cut away first bin
      // e->GetXaxis()->SetLimits(3,e->GetBinLowEdge(e->GetNbinsX())); //cut away first bin
      // f->GetXaxis()->SetLimits(4,f->GetBinLowEdge(f->GetNbinsX())); //cut away first bin


      //////////////////////////////////////////////////////////////////////////////////////////////
      //LOOPING AND OUTPUT//

      double f_mce = 0;
      //double f_mcesys = 0;
      int count = 0;
      for(int i=a->FindBin(2); i<=a->FindBin(10); i++)
        {
          ++count;
          f_mce += fabs(c->GetBinContent(i) - d->GetBinContent(i));
          //f_mcesys += fabs(eposrew->GetBinContent(i) - qgsrew->GetBinContent(i));
        }
      f_mce /= double(count);
      //f_mcesys /= double(count);

      pur_single->Expand(a->GetNbinsX());
      pur_double->Expand(a->GetNbinsX());

      cerr << "fix lumi" << endl;
      for(int i=1; i<=a->GetNbinsX(); i++)
        {
          const double f_mc     = (b->GetBinContent(i) + c->GetBinContent(i) + d->GetBinContent(i) + e->GetBinContent(i) + f->GetBinContent(i)) / 5.;
          //const double f_mcsys  = 0.5 * (eposrew->GetBinContent(i) + qgsrew->GetBinContent(i));
          const double f_noise  = a2->GetBinContent(i)/a2->GetBinContent(1);
          const double f_sel_zb = a->GetBinContent(i)/a->GetBinContent(1);

          ///////////////////////////////////////////////////////////////////////
          //Number of events
          const double n_sel = f_sel_zb * events_integral;//(events_integral/lumi_integral);
          const double n_noise = f_noise * events_integral;
          cout << events_integral << " " << f_noise  << " " << f_sel_zb << endl;


          ///////////////////////////////////////////////////////////////////////
          //Purity
          double purity = 2.06 / (2.06 + n_noise);
          pur->SetPoint(i,purity,f_mc);

          if(i!=1)
            {
              a2->SetBinContent(i,n_noise);
              //e->SetBinContent(i,e->GetBinContent(i)*0.195);
              //f->SetBinContent(i,f->GetBinContent(i)*0.195);
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
                //<< endl << "f_mc_scaled= " << f_mcsys << " ± " << f_mcesys << " ( " << f_mcesys/f_mcsys*100. << "%)"
                << endl << "n_noise= " << n_noise << " ± ? ( " << "?" << "%)"
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

          TLine* line = new TLine(1,ymin,1,ymax);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");
          line = new TLine(xmin,1,xmax,1);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          pur_single->SetTitle(";Purity;Acceptance #it{#epsilon}_{acc}");
          pur_double->SetTitle(";Purity;Acceptance #it{#epsilon}_{acc}");
          pur_single->SetLineColor(kRed);
          pur_double->SetLineColor(kBlue);
          pur_single->SetLineWidth(3);
          pur_double->SetLineWidth(3);
          pur_single->SetLineStyle(5);

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

          TLegend* leg2 = new TLegend(0.21,0.72,0.55,0.82);

              leg2->SetX1(0.21);
              leg2->SetY1(0.72);
              leg2->SetX2(0.55);
              leg2->SetY2(0.82);
#ifdef __CINT__
              CMSText(2,1,1);
#endif


          leg2->Draw();

          leg2->AddEntry(pur_single,"Single-arm selection","l");
          leg2->AddEntry(pur_double,"Double-arm selection","l");
#ifdef __CINT__
          SetLegAtt(leg2,1.1);
#endif

          pur_single->Draw("C");
          pur_double->Draw("C");

          for(int n=0; n<int(type.size()); n++)of
            {

              pur_selec = type[n]==string("single")?pur_selec_single:pur_selec_double;
              pur_selec2 = type[n]==string("single")?pur_selec2_single:pur_selec2_double;
              set<double>& cut_value = type[n]==string("single")?cut_value_single:cut_value_double;
              set<double>& cut_value2 = type[n]==string("single")?cut_value2_single:cut_value2_double;
              set<double>::iterator cut_ite;

              for (cut_ite = cut_value.begin(), int point = 0; cut_ite != cut_value.end(); ++cut_ite, point++)
                {
                  const double offsetx = type[n]==string("single")?0.001:0.023;
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
                  const double offsetx = type[n]==string("single")?0.001:0.023;
                  const double offsety = type[n]==string("single")?0.01:-0.005;

                  if (*cut_ite==3.) offsetx+=0.005;

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
