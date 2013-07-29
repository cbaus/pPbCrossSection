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
void makePlots_vdm_ls()
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  TFile* file = TFile::Open("histos_vdm.root");

  vector<string> type;
  type.push_back("X");
  type.push_back("Y");

  vector<double> sectionXBegin;
  vector<double> sectionXEnd;
  vector<double> sectionXTruth;
  sectionXBegin.push_back(1.93500e8); sectionXEnd.push_back(1.94577e8); sectionXTruth.push_back(100);
  sectionXBegin.push_back(1.94817e8); sectionXEnd.push_back(1.95273e8); sectionXTruth.push_back(50);
  sectionXBegin.push_back(1.95451e8); sectionXEnd.push_back(1.95817e8); sectionXTruth.push_back(0);
  sectionXBegin.push_back(1.95963e8); sectionXEnd.push_back(1.96540e8); sectionXTruth.push_back(-50);
  sectionXBegin.push_back(1.96897e8); sectionXEnd.push_back(1.97813e8); sectionXTruth.push_back(-100);
  sectionXBegin.push_back(1.97954e8); sectionXEnd.push_back(1.98557e8); sectionXTruth.push_back(-50);
  sectionXBegin.push_back(1.98663e8); sectionXEnd.push_back(1.98956e8); sectionXTruth.push_back(0);
  sectionXBegin.push_back(1.99136e8); sectionXEnd.push_back(1.99595e8); sectionXTruth.push_back(50);
  sectionXBegin.push_back(1.99747e8); sectionXEnd.push_back(2.00287e8); sectionXTruth.push_back(100);
  sectionXBegin.push_back(2.00490e8); sectionXEnd.push_back(2.02000e8); sectionXTruth.push_back(0);

  vector<double> sectionYBegin;
  vector<double> sectionYEnd;
  vector<double> sectionYTruth;
  sectionYBegin.push_back(2.04171e8); sectionYEnd.push_back(2.04940e8); sectionYTruth.push_back(100);
  sectionYBegin.push_back(2.05127e8); sectionYEnd.push_back(2.05583e8); sectionYTruth.push_back(50);
  sectionYBegin.push_back(2.05769e8); sectionYEnd.push_back(2.06234e8); sectionYTruth.push_back(0);
  sectionYBegin.push_back(2.06344e8); sectionYEnd.push_back(2.06852e8); sectionYTruth.push_back(-50);
  sectionYBegin.push_back(2.06970e8); sectionYEnd.push_back(2.07367e8); sectionYTruth.push_back(-100);
  sectionYBegin.push_back(2.07452e8); sectionYEnd.push_back(2.07900e8); sectionYTruth.push_back(-50);
  sectionYBegin.push_back(2.08035e8); sectionYEnd.push_back(2.08475e8); sectionYTruth.push_back(0);
  sectionYBegin.push_back(2.08572e8); sectionYEnd.push_back(2.09193e8); sectionYTruth.push_back(50);
  sectionYBegin.push_back(2.09362e8); sectionYEnd.push_back(2.09923e8); sectionYTruth.push_back(100);
  sectionYBegin.push_back(2.10075e8); sectionYEnd.push_back(2.11025e8); sectionYTruth.push_back(0);
  
  TCanvas* can1 = new TCanvas;
  can1->Divide(1,2);

  for(int n=0; n<int(type.size()); n++)
    {
      can1->cd(type[n]=="X"?1:2);
      
      vector<double>& sectionBegin = type[n]=="X"?sectionXBegin:sectionYBegin;
      vector<double>& sectionEnd = type[n]=="X"?sectionXEnd:sectionYEnd;
      vector<double>& sectionTruth = type[n]=="X"?sectionXTruth:sectionYTruth;
        
      const int nSec = int(sectionBegin.size());

      TH1D* h_length_scale =(TH1D*)file->Get((string("pPb/h_length_scale_") + type[n]).c_str());
      TGraphErrors* h_truth_fit = new TGraphErrors(nSec);
      ostringstream profiletitle; profiletitle << "truth_fit" << type[n];
      h_truth_fit->SetName(profiletitle.str().c_str());

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
              
          sectionFuncPtr[i] = h_length_scale->Fit(sectionFunc[i],"SQN","",sectionBegin[i],sectionEnd[i]);
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
          text << "#Delta" << type[n] << "_{Fit}=" << fixed << setprecision(1) << yum;
          TPaveText* txt = new TPaveText(sectionBegin[i],y+1.5*size,sectionEnd[i],y+2.5*size,"b t l");
          txt->SetTextFont(62);
          txt->SetFillStyle(0);
          txt->SetTextColor(kRed);
          txt->SetTextSize(0.04);
          txt->SetBorderSize(0);
          txt->AddText(text.str().c_str());
          txt->Draw("SAME");
          h_truth_fit->SetPoint(i,sectionTruth[i], yum);
          h_truth_fit->SetPointError(i,0, yumerror);
        }

      
      TCanvas* can2 = new TCanvas;
      h_truth_fit->SetTitle(";#DeltaX (LHC) [#mum]; #DeltaX (CMS) [#mum]");
      h_truth_fit->Draw("AP");
      h_truth_fit->Fit("pol1");
#ifdef __CINT__
      CMSText(1,0,1);
#endif
      can2->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_2.pdf")).c_str());
    }

  can1->SaveAs((string("plots/vdm_length_scale")+string("_1.pdf")).c_str());
}
