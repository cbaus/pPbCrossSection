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
#include "TPaveText.h"
#include "TROOT.h"
#include "TVectorD.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

void FindSections(TH1D* h_rate, vector<double> &sectionBegin, vector<double> &sectionEnd, const double& begin, const double& end, const double offset = 0.);

TCanvas* canFind;

void makePlots_vdm_ls()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  TFile* file = TFile::Open("histos_vdm.root");
  TH1D* h_rate = (TH1D*)file->Get((string("pPb/h_rate_ls")).c_str());
  h_rate->Scale(1./h_rate->GetBinWidth(1));
  h_rate->SetTitle(";orbitNb;dN/dt_{Orbit}");
  h_rate->GetXaxis()->SetNdivisions(505);
  h_rate->GetYaxis()->SetTitleOffset(h_rate->GetYaxis()->GetTitleOffset()*0.7);

  canFind = new TCanvas;
  canFind->Divide(1,2);

  vector<string> type;
  type.push_back("X");
  type.push_back("Y");

  const double skip = 5 * 11246;

  vector<double> sectionXBegin;
  vector<double> sectionXEnd;
  vector<double> sectionXTruth;
  // sectionXBegin.push_back(1.93500e8); sectionXEnd.push_back(1.94577e8); sectionXTruth.push_back(100);
  // sectionXBegin.push_back(1.94817e8); sectionXEnd.push_back(1.95273e8); sectionXTruth.push_back(50);
  // sectionXBegin.push_back(1.95451e8); sectionXEnd.push_back(1.95817e8); sectionXTruth.push_back(0);
  // sectionXBegin.push_back(1.95963e8); sectionXEnd.push_back(1.96540e8); sectionXTruth.push_back(-50);
  // sectionXBegin.push_back(1.96897e8); sectionXEnd.push_back(1.97813e8); sectionXTruth.push_back(-100);
  // sectionXBegin.push_back(1.97954e8); sectionXEnd.push_back(1.98557e8); sectionXTruth.push_back(-50);
  // sectionXBegin.push_back(1.98663e8); sectionXEnd.push_back(1.98956e8); sectionXTruth.push_back(0);
  // sectionXBegin.push_back(1.99136e8); sectionXEnd.push_back(1.99595e8); sectionXTruth.push_back(50);
  // sectionXBegin.push_back(1.99747e8); sectionXEnd.push_back(2.00287e8); sectionXTruth.push_back(100);
  // sectionXBegin.push_back(2.00490e8); sectionXEnd.push_back(2.02000e8); sectionXTruth.push_back(0);

  sectionXTruth.push_back(100); //0
  sectionXTruth.push_back(50);  //1
  sectionXTruth.push_back(0);   //2
  sectionXTruth.push_back(-50); //3
  sectionXTruth.push_back(-100);//4
  sectionXTruth.push_back(-50); //5
  sectionXTruth.push_back(0);   //6
  sectionXTruth.push_back(50);  //7
  sectionXTruth.push_back(100); //8
  sectionXTruth.push_back(0);   //9
  canFind->cd(1);
  FindSections(h_rate,sectionXBegin, sectionXEnd, 1.932e8,2.025e8,skip);
#ifdef __CINT__
  CMSText(1,1,1,string("#DeltaX length scale"));
#endif
  if(sectionXTruth.size() != sectionXBegin.size())
    {
      cerr << "Sections error: not same size. Truth/Found: " << sectionXTruth.size() << " " <<  sectionXBegin.size() << endl;
      //exit(-1);
    }

  vector<double> sectionYBegin;
  vector<double> sectionYEnd;
  vector<double> sectionYTruth;
  // sectionYBegin.push_back(2.04171e8); sectionYEnd.push_back(2.04940e8); sectionYTruth.push_back(100);
  // sectionYBegin.push_back(2.05127e8); sectionYEnd.push_back(2.05583e8); sectionYTruth.push_back(50);
  // sectionYBegin.push_back(2.05769e8); sectionYEnd.push_back(2.06234e8); sectionYTruth.push_back(0);
  // sectionYBegin.push_back(2.06344e8); sectionYEnd.push_back(2.06852e8); sectionYTruth.push_back(-50);
  // sectionYBegin.push_back(2.06970e8); sectionYEnd.push_back(2.07367e8); sectionYTruth.push_back(-100);
  // sectionYBegin.push_back(2.07452e8); sectionYEnd.push_back(2.07900e8); sectionYTruth.push_back(-50);
  // sectionYBegin.push_back(2.08035e8); sectionYEnd.push_back(2.08475e8); sectionYTruth.push_back(0);
  // sectionYBegin.push_back(2.08572e8); sectionYEnd.push_back(2.09193e8); sectionYTruth.push_back(50);
  // sectionYBegin.push_back(2.09362e8); sectionYEnd.push_back(2.09923e8); sectionYTruth.push_back(100);
  // sectionYBegin.push_back(2.10075e8); sectionYEnd.push_back(2.11025e8); sectionYTruth.push_back(0);

  sectionYTruth.push_back(100);
  sectionYTruth.push_back(50); 
  sectionYTruth.push_back(0);  
  sectionYTruth.push_back(-50);
  sectionYTruth.push_back(-100);
  sectionYTruth.push_back(-50);
  sectionYTruth.push_back(0);  
  sectionYTruth.push_back(50); 
  sectionYTruth.push_back(100);
  sectionYTruth.push_back(0);
  canFind->cd(2);
  FindSections(h_rate,sectionYBegin, sectionYEnd, 2.038e8,2.12e8,skip);
#ifdef __CINT__
  CMSText(1,1,1,string("#DeltaY length scale"));
#endif
  if(sectionYTruth.size() != sectionYBegin.size())
    {
      cerr << "Sections error: not same size. Truth/Found: " << sectionYTruth.size() << " " <<  sectionYBegin.size() << endl;
      exit(-1);
    }
  


  TCanvas* can1 = new TCanvas;
  can1->Divide(1,2);

  for(int n=0; n<int(type.size()); n++)
    {
      can1->cd(type[n]=="X"?1:2);
      
      vector<double>& sectionBegin = type[n]=="X"?sectionXBegin:sectionYBegin;
      vector<double>& sectionEnd = type[n]=="X"?sectionXEnd:sectionYEnd;
      vector<double>& sectionTruth = type[n]=="X"?sectionXTruth:sectionYTruth;
        
      const int nSec = int(sectionBegin.size());

      TH1D* h_length_scale = (TH1D*)file->Get((string("pPb/h_length_scale_") + type[n]).c_str());
      //      for(int i= 0; i<h_length_scale->GetNbinsX()+1;i++) cout << h_length_scale->Print
      TGraphErrors* h_truth_fit_up = new TGraphErrors(4);
      TGraphErrors* h_truth_fit_down = new TGraphErrors(4);
      TGraphErrors* h_fit_chi2 = new TGraphErrors(nSec);
      ostringstream profiletitle; profiletitle << "truth_fit" << type[n];
      h_truth_fit_down->SetName(profiletitle.str().c_str());
      profiletitle << "_down";
      h_truth_fit_up->SetName(profiletitle.str().c_str());
      profiletitle << "_chi2";
      h_fit_chi2->SetName(profiletitle.str().c_str());

      vector<TF1*> sectionFunc; sectionFunc.resize(nSec);
      vector<TFitResultPtr> sectionFuncPtr; sectionFuncPtr.resize(nSec);

      ostringstream histotitle; histotitle << ";orbit number; vertex " << type[n] << " position [mm]";
      h_length_scale->SetTitle(histotitle.str().c_str());
      h_length_scale->SetMarkerSize(0.5);
      h_length_scale->GetXaxis()->SetNdivisions(505);
      h_length_scale->GetYaxis()->SetTitleOffset(h_length_scale->GetYaxis()->GetTitleOffset()*0.7);
      if(type[n]=="X")
        h_length_scale->GetXaxis()->SetRangeUser(193000000,202000000);
      else
        h_length_scale->GetXaxis()->SetRangeUser(203000000,211000000);
      h_length_scale->GetYaxis()->SetRangeUser(0.03,0.12);
      h_length_scale->Draw();

#ifdef __CINT__
      CMSText(1,1,1,string("#Delta")+type[n]+string(" length scale"));
#endif
      //FIT EACH SECTION
      for (int i=0; i<nSec; i++)
        {
          ostringstream funcTitle;
          funcTitle << "section" << type[n] << "Func_" << i;
          sectionFunc[i] = new TF1(funcTitle.str().c_str(),"pol0",sectionBegin[i],sectionEnd[i]);
          sectionFunc[i]->SetLineWidth(2);
          sectionFunc[i]->SetLineColor(kRed);
          sectionFunc[i]->SetLineStyle(1);
              
          sectionFuncPtr[i] = h_length_scale->Fit(sectionFunc[i],"QSN","",sectionBegin[i],sectionEnd[i]);
          sectionFunc[i]->Draw("SAME");
          
        }
      
      //DETERMINE THE MEAN WHERE TRUTH IS AT NOMINAL POSITION
      double xw = 0;
      double w = 0;
      for(int i=0; i<nSec; i++)
        {
          if(sectionTruth[i]==0)
            {
              double weight = 1./pow(sectionFuncPtr[i]->ParError(0),2);
              xw += sectionFuncPtr[i]->Parameter(0) * weight;
              w += weight;
              cout << sectionFuncPtr[i]->Parameter(0) << " weight: " << weight << endl;
            }
        }
      double sigma = sqrt(1./w);
      double y0 = xw/w;
      cout << "y0 = " << y0 << "+-" << sigma << endl;

      //DRAW SPECIAL OBJECTS AND FILL TRUTH VS TRACKER PLOT
      for (int i=0; i<nSec; i++)
        {
          ostringstream text;
          double size = 0.01;
          double y = sectionFuncPtr[i]->Parameter(0);
          double yum = (y-y0)*10000.;
          double yumerror = (sectionFuncPtr[i]->ParError(0))*10000.;
          double chi2 = sectionFuncPtr[i]->Chi2()/double(sectionFuncPtr[i]->Ndf());
          text << "#Delta" << type[n] << "_{Fit}=" << fixed << setprecision(1) << yum;
          TPaveText* txt = new TPaveText(sectionBegin[i],y+1.5*size,sectionEnd[i],y+2.5*size,"b t l");
          txt->SetTextFont(62);
          txt->SetFillStyle(0);
          txt->SetTextColor(kRed);
          txt->SetTextSize(0.04);
          txt->SetBorderSize(0);
          txt->AddText(text.str().c_str());
          txt->Draw("SAME");
          if(i>=1 && i<=4)
            {
              cout << "up " << i << " " << sectionTruth[i]*(type[n]=="X"?-1:1) << endl;
              h_truth_fit_up->SetPoint(i-1,sectionTruth[i]*(type[n]=="X"?-1:1), yum); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_up->SetPointError(i-1,0, yumerror);
            }
          if(i>=5 && i<=8) //anti hysteresis
            {
              cout << "down " << i << " " << sectionTruth[i]*(type[n]=="X"?-1:1) << endl;
              h_truth_fit_down->SetPoint(i-5,sectionTruth[i]*(type[n]=="X"?-1:1), yum); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_down->SetPointError(i-5,0, yumerror);
            }
          h_fit_chi2->SetPoint(i,i,chi2);
          
        }
      h_truth_fit_up->Print("ALL");
      h_truth_fit_down->Print("ALL");

      
      TCanvas* can2 = new TCanvas;
      ostringstream titleup; titleup << ";#Delta" << type[n] << " (LHC) [#mum]; #Delta" << type[n] << " (CMS) [#mum]";
      h_truth_fit_up->SetTitle(titleup.str().c_str());
      h_truth_fit_up->SetMarkerSize(1.3);
      h_truth_fit_up->SetMarkerColor(kRed);
      h_truth_fit_up->SetLineColor(kRed);
      h_truth_fit_up->SetLineWidth(2);
      h_truth_fit_up->GetXaxis()->SetLimits(-110,110);
      h_truth_fit_up->GetYaxis()->SetRangeUser(-110,110);
      h_truth_fit_up->Draw("AP");
      h_truth_fit_up->Fit("pol1");
      h_truth_fit_up->GetFunction("pol1")->SetLineWidth(2);
      h_truth_fit_up->GetFunction("pol1")->SetLineColor(kRed);
      h_truth_fit_down->SetMarkerColor(kBlue);
      h_truth_fit_down->SetMarkerStyle(25);
      h_truth_fit_down->SetMarkerSize(1.3);
      h_truth_fit_down->SetLineColor(kBlue);
      h_truth_fit_down->SetLineWidth(2);
      h_truth_fit_down->Draw("P");
      h_truth_fit_down->Fit("pol1");
      h_truth_fit_down->GetFunction("pol1")->SetLineWidth(2);
      h_truth_fit_down->GetFunction("pol1")->SetLineColor(kBlue);
      TLegend* leg = new TLegend(0.22,0.72,0.5,0.82);
      leg->AddEntry(h_truth_fit_up,(string("#Delta")+type[n]+string(" upwards")).c_str(),"lp");
      leg->AddEntry(h_truth_fit_down,(string("#Delta")+type[n]+string(" downwards")).c_str(),"lp");
#ifdef __CINT__
      SetLegAtt(leg);
      CMSText(1,1,1);
#endif
      leg->Draw("SAME");
      can2->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_2.pdf")).c_str());

      TCanvas* can3 = new TCanvas;
      h_fit_chi2->SetTitle(";section; #chi^{2}/NDF");
      h_fit_chi2->SetMarkerSize(2);
      h_fit_chi2->Draw("AP");
#ifdef __CINT__
      CMSText(1,1,1);
#endif
      can3->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_3.pdf")).c_str());
    }

  can1->SaveAs((string("plots/vdm_length_scale")+string("_1.pdf")).c_str());
  canFind->SaveAs((string("plots/vdm_length_scale")+string("_sections.pdf")).c_str());
}

void FindSections(TH1D* h_rate, vector<double> &sectionBegin, vector<double> &sectionEnd, const double& begin, const double& end, const double offset)
{
  if(begin < h_rate->GetBinLowEdge(1) || end > h_rate->GetBinLowEdge(h_rate->GetNbinsX()+1))
    {
      cerr << "Borders outside of histogram " << begin << " " <<  h_rate->GetBinLowEdge(1) << " " << end << " " << h_rate->GetBinLowEdge(h_rate->GetNbinsX()+1) << endl;
      exit(-1);
    }

  const double cut = h_rate->GetMaximum() / 3.;
  cout << "Using cut value = " << cut << endl;

  const double step = h_rate->GetBinWidth(1) / 100.;

  double prevValue = -1e9;

  for(double current = begin; current < end; current += step)
    {
      const double value = h_rate->Interpolate(current);
      if (prevValue < cut && value >= cut)
        sectionBegin.push_back(current+3*11246+offset);
      if (prevValue >= cut && value < cut)
        sectionEnd.push_back(current-3*11246);
      prevValue = value;
    }

  if(sectionBegin.size() == sectionEnd.size()+1)
    sectionEnd.push_back(end);
    
  if(sectionBegin.size() != sectionEnd.size())
    {
      cerr << "sections size error" << endl;
      exit(-1);
    }

  ostringstream name; name << "rate" << begin << end;
  TH1D* draw = (TH1D*)(h_rate->Clone(name.str().c_str()));
  draw->GetXaxis()->SetRangeUser(begin,end);
  draw->Draw("HIST L");
  for (int i=0; i < int(sectionBegin.size()); ++i)
    {
      TBox* box = new TBox(sectionBegin[i],cut/2.,sectionEnd[i],cut);
      box->SetFillColor(kRed);
      box->SetFillStyle(3001);
      box->DrawClone();
    }
  return;
}
