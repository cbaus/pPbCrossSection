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
#include "TVectorD.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

TVectorD corr_fac_em(2);
TVectorD corr_fac_eme(2);
TVectorD corr_fac_mc(4);
TVectorD corr_fac_mce(4);

void makePlots(bool draw=1, double cut_value_single = 8., double cut_value_double = 2.5);

//available cuts single 5, 6.4, 7.2, 8, 8.8, 9.6, 10, 15, 20
//available cuts double 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5

void makePlots(bool draw, double cut_value_single, double cut_value_double)
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<string> type; type.push_back("single"); type.push_back("double");

  for(int n=0; n<int(type.size()); n++)
    {

      TFile* file2 = TFile::Open("histos_old.root");
      TFile* file = TFile::Open("histos_mc.root");
      TFile* file3 = TFile::Open("plots/hf_cuts_noise.root");

      TH1D* a=(TH1D*)file2->Get((string("data210885/data210885_h_hf_cut_") + type[n]).c_str());
      TH1D* a2=(TH1D*)file2->Get((string("data210885/data210885_h_hf_cut_") + type[n] + string("_noise")).c_str());
      TH1D* eposrew=(TH1D*)file->Get((string("EposSDWeight2/EposSDWeight2_h_hf_cut_") + type[n]).c_str());
      TH1D* eposrewnew=(TH1D*)file->Get((string("EposSDWeight2/EposSDWeight2_h_hf_new_cut_")+ type[n]).c_str());
      TH1D* b=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_cut_") + type[n]).c_str());
      TH1D* b2=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_new_cut_")+ type[n]).c_str());
      TH1D* c=(TH1D*)file->Get((string("Epos/Epos_h_hf_cut_")+ type[n]).c_str());
      TH1D* c2=(TH1D*)file->Get((string("Epos/Epos_h_hf_new_cut_")+ type[n]).c_str());
      TH1D* d=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_cut_") + type[n]).c_str());
      TH1D* d2=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_new_cut_")+ type[n]).c_str());
      TH1D* e=(TH1D*)file->Get((string("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_") + type[n]).c_str());
      TH1D* f=(TH1D*)file->Get((string("Starlight_Pythia/Starlight_Pythia_h_hf_cut_") + type[n]).c_str());
      TVectorD* hf_cuts_equivalent = (TVectorD*)file3->Get("hf_cuts_equivalent");

      a->Scale(1./double(a->GetBinContent(1)));
      a2->Scale(1./double(a2->GetBinContent(1)));
      b->Scale(1./double(b->GetBinContent(1)));
      b2->Scale(1./double(b2->GetBinContent(1)));
      c->Scale(1./double(c->GetBinContent(1)));
      c2->Scale(1./double(c2->GetBinContent(1)));
      d->Scale(1./double(d->GetBinContent(1)));
      d2->Scale(1./double(d2->GetBinContent(1)));
      e->Scale(1./double(e->GetBinContent(1)));
      f->Scale(1./double(f->GetBinContent(1))/ 195. * 122.);
      eposrew->Scale(1./double(eposrew->GetBinContent(1)));
      eposrewnew->Scale(1./double(eposrewnew->GetBinContent(1)));

      a->SetLineWidth(3);
      a2->SetLineWidth(3);
      b->SetLineWidth(3);
      c->SetLineWidth(3);
      d->SetLineWidth(3);
      e->SetLineWidth(3);
      f->SetLineWidth(3);
      eposrew->SetLineWidth(3);
      a2->SetLineColor(kBlack);
      b->SetLineColor(kGreen+2);
      c->SetLineColor(kBlue);
      d->SetLineColor(kRed);
      e->SetLineColor(kRed);
      f->SetLineColor(kBlue);
      eposrew->SetLineColor(kMagenta);
      a2->SetMarkerColor(a2->GetLineColor());
      b->SetMarkerColor(b->GetLineColor());
      c->SetMarkerColor(c->GetLineColor());
      d->SetMarkerColor(d->GetLineColor());
      e->SetMarkerColor(e->GetLineColor());
      f->SetMarkerColor(f->GetLineColor());
      eposrew->SetMarkerColor(f->GetLineColor());
      b->SetLineStyle(7);
      d->SetLineStyle(9);
      f->SetLineStyle(9);
      eposrew->SetLineStyle(10);

      a->SetTitle("Data");
      a2->SetTitle("Noise");
      eposrew->SetTitle("EPOS-LHC (#sigma_{SD}x2)");
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
      a2->GetYaxis()->SetRangeUser(type[n]=="double"?1e-5:1e-5,1.01);

      b->GetXaxis()->SetTitle("E_{HF} tower threshold [GeV]");
      b->GetYaxis()->SetTitle("event fraction");
      b->GetXaxis()->SetLabelSize(b->GetXaxis()->GetLabelSize()*1.2);
      b->GetYaxis()->SetLabelSize(b->GetYaxis()->GetLabelSize()*1.2);
      b->GetXaxis()->SetTitleSize(b->GetXaxis()->GetTitleSize()*1.1);
      b->GetYaxis()->SetTitleSize(b->GetYaxis()->GetTitleSize()*1.1);
      b->GetXaxis()->SetTitleOffset(b->GetXaxis()->GetTitleOffset()*1.1);
      b->GetYaxis()->SetTitleOffset(b->GetYaxis()->GetTitleOffset()*1.1);


      TGraphErrors* b3 = new TGraphErrors(3);
      TGraphErrors* c3 = new TGraphErrors(3);
      TGraphErrors* d3 = new TGraphErrors(3);
      for (int j=0; j<3; j++)
        {
          b3->SetPoint(j,(*hf_cuts_equivalent)[j],b2->GetBinContent(j+2));
          c3->SetPoint(j,(*hf_cuts_equivalent)[j],c2->GetBinContent(j+2));
          d3->SetPoint(j,(*hf_cuts_equivalent)[j],d2->GetBinContent(j+2));
        }
      b3->SetMarkerColor(b->GetLineColor());
      c3->SetMarkerColor(c->GetLineColor());
      d3->SetMarkerColor(d->GetLineColor());
      b3->SetMarkerSize(1.5);
      c3->SetMarkerSize(1.5);
      d3->SetMarkerSize(1.5);

      const double cut_value = type[n]=="single"?cut_value_single:cut_value_double;

      if(draw)
        {
          TCanvas* can1 = new TCanvas;
          b->Draw("HIST L");
          c->Draw("HIST L SAME");
          d->Draw("HIST L SAME");

          TLine* line = new TLine(cut_value,0.8,cut_value,1.001);
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
              leg->SetX1(0.41);
              leg->SetX2(0.74);
              leg->SetY1(0.21);
              leg->SetY2(0.44);
#ifdef __CINT__
              CMSText(1,0,1,"double-arm selection");
#endif
            }

          leg->AddEntry(b,"","l");
          leg->AddEntry(c,"","l");
          leg->AddEntry(d,"","l");
#ifdef __CINT__
          SetLegAtt(leg,1.3);
#endif
          leg->Draw();

          a->GetYaxis()->SetRangeUser(0,1.01);
          a->GetXaxis()->SetTitle("E_{HF} tower threshold [GeV]");
          a->GetYaxis()->SetTitle("event fraction");

          can1->SaveAs((string("plots/full_p_space_eff_PAS_")+type[n]+string(".pdf")).c_str());

          b3->Draw("P");
          c3->Draw("P");
          d3->Draw("P");
          can1->SaveAs((string("plots/full_p_space_eff_")+type[n]+string(".pdf")).c_str());
        }

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
          f_mcesys += fabs(eposrew->GetBinContent(i) - d->GetBinContent(i));
        }
      f_eme /= double(count);
      f_mce /= double(count);
      f_mcesys /= double(count);


      for(int i=1; i<=a->GetNbinsX(); i++)
        {
          const double f_em     = 0.5 * (e->GetBinContent(i) + f->GetBinContent(i));
          const double f_mc     = 0.5 * (c->GetBinContent(i) + d->GetBinContent(i));
          const double f_noise  = a2->GetBinContent(i)/a2->GetBinContent(1);
          const double n_sel_zb = a->GetBinContent(i);
          //const double n_zb     = 1;


          if(i==a->FindBin(cut_value_single) && type[n]==string("single"))
            {
              // if((i-1)%10 == 0)
              cout
                << setprecision(3)
                << endl << i << "(" << a->GetBinCenter(i) << ")"
                << endl << "f_mc= " << f_mc << " ± " << f_mce << " ( " << f_mce/f_mc*100. << "%)"
                << endl << "f_em= " << f_em << " ± " << f_eme << " ( " << f_eme/f_em*100. << "%)"
                << endl << "f_noise= " << f_noise
                << endl << "n_sel_zb= " << n_sel_zb << endl << endl;

              corr_fac_em[0] = f_em;
              corr_fac_mc[0] = f_mc;
              corr_fac_mc[2] = 0.5 * (eposrew->GetBinContent(i) + d->GetBinContent(i));
              corr_fac_eme[0] = f_eme;
              corr_fac_mce[0] = f_mce;
              corr_fac_mce[2] = f_mcesys;
              cout << endl << endl << "Trigger & \\EPOS & \\HIJING & \\QGSJET \\\\\\hline\\hline" << endl;
              cout <<  " single-arm & "
                   << setprecision(3)
                   << c->GetBinContent(i)/c->GetBinContent(1) << " & "
                   << b->GetBinContent(i)/b->GetBinContent(1) << " & "
                   << d->GetBinContent(i)/d->GetBinContent(1)
                   << " \\\\\\hline" << endl;
            }
          if(i==a->FindBin(cut_value_double) && type[n]==string("double"))
            {
              cout
                << setprecision(3)
                << endl << i << "(" << a->GetBinCenter(i) << ")"
                << endl << "f_mc= " << f_mc << " ± " << f_mce << " ( " << f_mce/f_mc*100. << "%)"
                << endl << "f_em= " << f_em << " ± " << f_eme << " ( " << f_eme/f_em*100. << "%)"
                << endl << "f_noise= " << f_noise
                << endl << "n_sel_zb= " << n_sel_zb << endl << endl;

              corr_fac_em[1] = f_em;
              corr_fac_mc[1] = f_mc;
              corr_fac_eme[1] = f_eme;
              corr_fac_mce[1] = f_mce;
              corr_fac_mc[3] = 0.5 * (eposrew->GetBinContent(i) + d->GetBinContent(i));
              corr_fac_mce[3] = f_mcesys;
              cout <<  " double-arm & "
                   << setprecision(3)
                   << c->GetBinContent(i)/c->GetBinContent(1) << " & "
                   << b->GetBinContent(i)/b->GetBinContent(1) << " & "
                   << d->GetBinContent(i)/d->GetBinContent(1)
                   << " \\\\\\hline" << endl;
            }
        }

      //////////////////////////////////////////////////////////////////////////////////////////////

      if(draw)
        {
          TCanvas* can2 = new TCanvas;
          a2->GetXaxis()->SetTitle("E_{HF} tower threshold [GeV]");
          a2->GetYaxis()->SetTitle("event fraction");
          a2->GetXaxis()->SetLabelSize(a2->GetXaxis()->GetLabelSize()*1.2);
          a2->GetYaxis()->SetLabelSize(a2->GetYaxis()->GetLabelSize()*1.2);
          a2->GetXaxis()->SetTitleSize(a2->GetXaxis()->GetTitleSize()*1.1);
          a2->GetYaxis()->SetTitleSize(a2->GetYaxis()->GetTitleSize()*1.1);
          a2->GetXaxis()->SetTitleOffset(a2->GetXaxis()->GetTitleOffset()*1.1);
          a2->GetYaxis()->SetTitleOffset(a2->GetYaxis()->GetTitleOffset()*1.05);
          a2->Draw("HIST l");
          e->Draw("HIST l SAME");
          f->Draw("HIST l SAME");
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

          leg2->AddEntry(a2,"","l");
          leg2->AddEntry(e,"","l");
          leg2->AddEntry(f,"","l");
#ifdef __CINT__
          SetLegAtt(leg2,1.1);
#endif

          TLine* line = new TLine(cut_value,type[n]=="single"?15e-4:0,cut_value,1.01);
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
  outfile.Close();
}
