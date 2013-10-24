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
#include "TMinuit.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TVectorD.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>   
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#define BEGINX 80.8e6
#define ENDX 86.4e6
#define BEGINY 88.4e6
#define ENDY 93.8e6

void FindSections(TH1D* h_rate, vector<double> &sectionBegin, vector<double> &sectionEnd, const double& begin, const double& end, const double offset = 3*11246., const double skip = 15.*11246.);

TCanvas* canFind;

void makePlots_vdm_ls(double skip = 10.*11246.)
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  TFile* file = TFile::Open("histos_vdmpp.root");
  TH1D* h_rate = (TH1D*)file->Get((string("pp/h_rate_ls")).c_str());
  h_rate->Scale(1./h_rate->GetBinWidth(1));
  h_rate->SetTitle(";orbitNb;dN/dT_{Orbit}");
  h_rate->GetXaxis()->SetNdivisions(505);
  h_rate->GetYaxis()->SetTitleOffset(h_rate->GetYaxis()->GetTitleOffset()*0.7);

  canFind = new TCanvas;
  canFind->Divide(1,2);

  vector<string> type;
  type.push_back("X");
  type.push_back("Y");

  vector<double> sectionXBegin;
  vector<double> sectionXEnd;
  vector<double> sectionXTruth;
  vector<double> sectionXMean;
  vector<double> sectionXMeanE;
  vector<TGraphErrors*> sectionXChi;

//   sectionXTruth.push_back(100); //0
//   sectionXTruth.push_back(50);  //1
//   sectionXTruth.push_back(0);   //2
//   sectionXTruth.push_back(-50); //3
//   sectionXTruth.push_back(-100);//4
//   sectionXTruth.push_back(-50); //5
//   sectionXTruth.push_back(0);   //6
//   sectionXTruth.push_back(50);  //7
//   sectionXTruth.push_back(100); //8
//   sectionXTruth.push_back(0);   //9

  //sectionXTruth.push_back(-70); //0
  sectionXTruth.push_back(-370);  //1
  sectionXTruth.push_back(-220);   //2
  sectionXTruth.push_back(-70); //3
  sectionXTruth.push_back(80);//4
  sectionXTruth.push_back(230); //5
  sectionXTruth.push_back(80);   //6
  sectionXTruth.push_back(-70);  //7
  sectionXTruth.push_back(-220); //8
  sectionXTruth.push_back(-370);   //9
  canFind->cd(1);
  FindSections(h_rate,sectionXBegin, sectionXEnd, BEGINX,ENDX,3.*11246.,skip);
  //FindSections(h_rate,sectionXBegin, sectionXEnd, 1.932e8,2.025e8,3.*11246.,skip);
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
  vector<double> sectionYMean;
  vector<double> sectionYMeanE;
  vector<TGraphErrors*> sectionYChi;

//   sectionYTruth.push_back(100);
//   sectionYTruth.push_back(50); 
//   sectionYTruth.push_back(0);  
//   sectionYTruth.push_back(-50);
//   sectionYTruth.push_back(-100);
//   sectionYTruth.push_back(-50);
//   sectionYTruth.push_back(0);  
//   sectionYTruth.push_back(50); 
//   sectionYTruth.push_back(100);
//   sectionYTruth.push_back(0);


  //sectionYTruth.push_back(70); //0
  sectionYTruth.push_back(-370);  //1
  sectionYTruth.push_back(-220);   //2
  sectionYTruth.push_back(-70); //3
  sectionYTruth.push_back(80);//4
  sectionYTruth.push_back(230); //5
  sectionYTruth.push_back(80);   //6
  sectionYTruth.push_back(-70);  //7
  sectionYTruth.push_back(-220); //8
  sectionYTruth.push_back(-370);   //9
  canFind->cd(2);
  //FindSections(h_rate,sectionYBegin, sectionYEnd, 2.038e8,2.12e8,3.*11246.,skip);
  FindSections(h_rate,sectionYBegin, sectionYEnd, BEGINY,ENDY,3.*11246.,skip);
#ifdef __CINT__
  CMSText(1,1,1,string("#DeltaY length scale"));
#endif
  if(sectionYTruth.size() != sectionYBegin.size())
    {
      cerr << "Sections error: not same size. Truth/Found: " << sectionYTruth.size() << " " <<  sectionYBegin.size() << endl;
      //exit(-1);
    }
  


  TCanvas* can1 = new TCanvas;
  can1->Divide(1,2);

  //GENERAL LOOP OVER SECTIONS
  for(int n=0; n<int(type.size()); n++)
    {
      can1->cd(type[n]=="X"?1:2);
      
      vector<double>& sectionBegin      = type[n]=="X"?sectionXBegin:sectionYBegin;
      vector<double>& sectionEnd        = type[n]=="X"?sectionXEnd:sectionYEnd;
      vector<double>& sectionTruth      = type[n]=="X"?sectionXTruth:sectionYTruth;
      vector<double>& sectionMean       = type[n]=="X"?sectionXMean:sectionYMean;
      vector<double>& sectionMeanE      = type[n]=="X"?sectionXMeanE:sectionYMeanE;
      vector<TGraphErrors*>& sectionChi = type[n]=="X"?sectionXChi:sectionYChi;
        
      const int nSec = int(sectionBegin.size());

      sectionMean.resize(nSec);
      sectionMeanE.resize(nSec);

      TH2D* h_length_scale = (TH2D*)file->Get((string("pp/h_length_scale_") + type[n]).c_str());
      //      for(int i= 0; i<h_length_scale->GetNbinsX()+1;i++) cout << h_length_scale->Print
      TGraphErrors* h_truth_fit_up = new TGraphErrors(0);
      TGraphErrors* h_truth_fit_down = new TGraphErrors(0);
      TGraphErrors* h_truth_fit_crosscheck_up = new TGraphErrors(0);
      TGraphErrors* h_truth_fit_crosscheck_down = new TGraphErrors(0);
      ostringstream profiletitle; profiletitle << "truth_fit_" << type[n];
      ostringstream crosschecktitle; crosschecktitle << "truth_fit_crosscheck_" << type[n];
      h_truth_fit_down->SetName(profiletitle.str().c_str());
      h_truth_fit_crosscheck_down->SetName(crosschecktitle.str().c_str());
      profiletitle << "_down";
      crosschecktitle << "_down";
      h_truth_fit_up->SetName(profiletitle.str().c_str());
      h_truth_fit_crosscheck_up->SetName(crosschecktitle.str().c_str());

      vector<TF1*> sectionFunc; sectionFunc.resize(nSec);
      vector<TFitResultPtr> sectionFuncPtr; sectionFuncPtr.resize(nSec);

      ostringstream histotitle; histotitle << ";orbit number; vertex " << type[n] << " position [mm]";
      h_length_scale->SetTitle(histotitle.str().c_str());
      h_length_scale->SetMarkerSize(0.5);
      h_length_scale->GetXaxis()->SetNdivisions(505);
      h_length_scale->GetYaxis()->SetTitleOffset(h_length_scale->GetYaxis()->GetTitleOffset()*0.7);
      if(type[n]=="X")
        h_length_scale->GetXaxis()->SetRangeUser(BEGINX,ENDX);
      else
        h_length_scale->GetXaxis()->SetRangeUser(BEGINY,ENDY);
//       if(type[n]=="X")
//         h_length_scale->GetXaxis()->SetRangeUser(193000000,202000000);
//       else
//         h_length_scale->GetXaxis()->SetRangeUser(203000000,211000000);
//      h_length_scale->GetYaxis()->SetRangeUser(0.00,0.15);
      h_length_scale->Draw("COLZ");

#ifdef __CINT__
      CMSText(1,1,1,string("#Delta")+type[n]+string(" length scale"));
#endif
      //FIT EACH SECTION
      TCanvas* canFit = new TCanvas;
      canFit->Divide(TMath::Nint(nSec/2.),2);
      for (int i=0; i<nSec; i++)
        {
          int binStartNotSkipped,binStart, binEnd, helper;
          double x1,x2,xend;

          h_length_scale->GetBinXYZ(h_length_scale->FindBin(sectionBegin[i]),binStartNotSkipped,helper,helper); binStartNotSkipped++;
          h_length_scale->GetBinXYZ(h_length_scale->FindBin(sectionBegin[i]+skip),binStart,helper,helper); binStart++;
          h_length_scale->GetBinXYZ(h_length_scale->FindBin(sectionEnd[i]),binEnd,helper,helper); binEnd--;
          x1=h_length_scale->GetBinLowEdge(binStartNotSkipped);
          x2=h_length_scale->GetBinLowEdge(binStart);
          xend=h_length_scale->GetBinLowEdge(binEnd+1);

          if(binStart >= binEnd)
            {
              cerr << "MAJOR WARNING Chosen skipping value too large for section " << i << "(" << binStart << "," << binEnd << ")" << endl;
            }

          ostringstream funcTitle;
          funcTitle << "section" << type[n] << "Func_" << i;
          sectionFunc[i] = new TF1(funcTitle.str().c_str(),"gaus");
          sectionFunc[i]->SetLineWidth(2);
          sectionFunc[i]->SetLineColor(kRed);
          sectionFunc[i]->SetLineStyle(1);

          funcTitle << "_Prof";
          TH1D* helperProfile = h_length_scale->ProjectionY(funcTitle.str().c_str(),binStart,binEnd);
          canFit->cd(i+1);
          helperProfile->Draw();
          sectionFuncPtr[i] = helperProfile->Fit(sectionFunc[i],"QS");


          if(Int_t(sectionFuncPtr[i]) !=0)
            {
              cerr << " !!! MAJOR WARNING." << endl
                   << "in section " << i << " fit did not converge: " << gMinuit->fCstatu.Data() << endl;
            }

          sectionFunc[i]->SetLineWidth(2);
          sectionFunc[i]->Draw("SAME");

          sectionMean[i] = helperProfile->GetMean();
          sectionMeanE[i] = helperProfile->GetMeanError();

          if(x1<x2)
            {
              //can1->cd();
              TBox* box1 = new TBox(x1,0.1,x2,0.105);
              box1->SetFillColor(kRed);
              box1->SetFillStyle(3001);
              box1->DrawClone();
              TBox* box2 = new TBox(x2,0.1,xend,0.105);
              box2->SetFillColor(kGreen);
              box2->SetFillStyle(3001);
              box2->DrawClone();
            }
        }
      
//       //DETERMINE THE MEAN WHERE TRUTH IS AT NOMINAL POSITION -> now done below
//       double xw = 0;
//       double w = 0;
//       for(int i=0; i<nSec; i++)
//         {
//           if(sectionTruth[i]==70 || sectionTruth[i]==-70)
//             {
//               double weight = 1./pow(sectionFuncPtr[i]->ParError(1),2);
//               xw += sectionFuncPtr[i]->Parameter(1) * weight;
//               w += weight;
//               //cout << sectionFuncPtr[i]->Parameter(1) << " weight: " << weight << endl;
//             }
//         }
//       double sigma = sqrt(1./w);
//       double y0 = xw/w;
//       cout << "y0 = " << y0 << "+-" << sigma << endl;

      //CHECK CHI2 OPTIMISATION
      TCanvas* can3 = new TCanvas;
      sectionChi.resize(nSec);
      string drawoption("AL");
      TLegend* legchi = new TLegend(0.55,0.5,0.85,0.9);
      for(int i=0; i<nSec; i++)
        {
          if(i<1 || i >8)
            continue;
          
          const int nSkipBins=10;
          ostringstream graphTitle;
          graphTitle << "section_" << type[n] << "_chi2_Graph_" << i;
          sectionChi[i] = (new TGraphErrors(0));
          sectionChi[i]->SetName(graphTitle.str().c_str());
          for (int skipBin=0; skipBin<nSkipBins; skipBin++)
            {
              int binStart, binEnd, helper;
              h_length_scale->GetBinXYZ(h_length_scale->FindBin(sectionBegin[i]),binStart,helper,helper); binStart++;
              h_length_scale->GetBinXYZ(h_length_scale->FindBin(sectionEnd[i]),binEnd,helper,helper); binEnd--;

              binStart += skipBin;

              if(binStart >= binEnd)
                continue;

              double x = double(h_length_scale->GetBinLowEdge(binStart) - h_length_scale->GetBinLowEdge(binStart-skipBin)) / 11246.;
              //cout << "skip bin: " << skipBin << "(" << binStart << ")" << endl;

              ostringstream funcTitle;
              funcTitle << "section" << type[n] << "Func_" << i << "_" << skipBin;

              TF1* helperFunc = new TF1(funcTitle.str().c_str(),"gaus");
              funcTitle << "_Prof";
              TH1D* helperProfile = h_length_scale->ProjectionY(funcTitle.str().c_str(),binStart,binEnd);
              TFitResultPtr helperPtr = helperProfile->Fit(helperFunc,"QSN");

              if(Int_t(helperPtr) !=0)
                {
                  cout << "fit failed" << endl;
                  cout << "skip bin: " << skipBin << "(" << binStart << ")" << endl;
                  continue;
                }
              
              double chi2 = helperPtr->Chi2()/double(helperPtr->Ndf());
              sectionChi[i]->SetPoint(sectionChi[i]->GetN(),x,chi2);
              //cout << skipBin << " " << x << " " << chi2 << endl;
            }
          //sectionChi[i]->Print("ALL");
          sectionChi[i]->SetTitle(";skip interval [s]; #chi^{2}/NDF");
          sectionChi[i]->GetYaxis()->SetRangeUser(0,6);
          sectionChi[i]->GetXaxis()->SetLimits(0,30);
          sectionChi[i]->SetLineWidth(2);
          sectionChi[i]->SetLineColor(i);
          sectionChi[i]->Draw(drawoption.c_str());
          drawoption = string("L");
          
          ostringstream legtext ; legtext << "Section " << i << " (" << sectionTruth[i] << " mm)";
          legchi->AddEntry(sectionChi[i],legtext.str().c_str(),"l");

        }
#ifdef __CINT__
      CMSText(1,1,1,string("#Delta")+type[n]+string(" length scale"));
      SetLegAtt(legchi);
#endif
      legchi->Draw("SAME");
      TLine* line = new TLine(skip/11246.,0,skip/11246.,2.5);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->Draw("SAME");
      can3->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_3.pdf")).c_str());
          
      //Determine reference value y0 to subtract from plotting
      TGraph helperGraph(nSec);
      double y0=0;
      for (int i=0; i<nSec; i++)
        {
          double y = sectionFuncPtr[i]->Parameter(1);
          helperGraph.SetPoint(i,sectionTruth[i]*(type[n]=="X"?-1:1), y);
          y0 = helperGraph.Eval(0);
        }
      cout << "y0 = " << y0 << endl;

      //TRUTH VS TRACKER PLOT
      for (int i=0; i<nSec; i++)
        {
          ostringstream text;
          double size = 0.01;
          double y = sectionFuncPtr[i]->Parameter(1);
          double y_cc = sectionMean[i]; //just used for cross checking (cc) the gaussian fit
          double yum = (y-y0)*1e4; //cm -> µm => 1e-2 -> 1e-6 => 1e4
          double yum_cc = (y_cc-y0)*1e4; //cm -> µm => 1e-2 -> 1e-6 => 1e4
          double yumerror = (sectionFuncPtr[i]->ParError(1))*10000.;
          double yumerror_cc = sectionMeanE[i] * 10000;
          text << "#Delta" << type[n] << "_{Fit}=" << fixed << setprecision(1) << yum;
          TPaveText* txt = new TPaveText(sectionBegin[i],y+1.5*size,sectionEnd[i],y+2.5*size,"b t l");
          txt->SetTextFont(62);
          txt->SetFillStyle(0);
          txt->SetTextColor(kRed);
          txt->SetTextSize(0.04);
          txt->SetBorderSize(0);
          txt->AddText(text.str().c_str());
          txt->Draw("SAME");
          if(i>=1 && i<=4) //WARNING also setpoint needs to be changed
            {
              cout << "up " << i << " " << sectionTruth[i]*(type[n]=="X"?-1:1) << " " << yum << "+-" << yumerror << endl;
              h_truth_fit_up->SetPoint(i-1,sectionTruth[i]*(type[n]=="X"?-1:1), yum); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_up->SetPointError(i-1,0, yumerror);
              h_truth_fit_crosscheck_up->SetPoint(i-1,sectionTruth[i]*(type[n]=="X"?-1:1), yum_cc); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_crosscheck_up->SetPointError(i-1,0, yumerror_cc);
            }
          if(i>=5 && i<=8) //anti hysteresis //WARNING also setpoint needs to be changed
            {
              cout << "down " << i << " " << sectionTruth[i]*(type[n]=="X"?-1:1) << " " << yum << "+-" << yumerror << endl;
              h_truth_fit_down->SetPoint(i-5,sectionTruth[i]*(type[n]=="X"?-1:1), yum); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_down->SetPointError(i-5,0, yumerror);
              h_truth_fit_crosscheck_down->SetPoint(i-5,sectionTruth[i]*(type[n]=="X"?-1:1), yum_cc); // Coordinate Sys different z? X has to be inverted.
              h_truth_fit_crosscheck_down->SetPointError(i-5,0, yumerror_cc);
            }
        }
      canFit->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_Fits.pdf")).c_str());
      
      TCanvas* can2 = new TCanvas;
      ostringstream titleup; titleup << ";#Delta" << type[n] << " (LHC) [#mum]; #Delta" << type[n] << " (CMS) [#mum]";
      h_truth_fit_up->SetTitle(titleup.str().c_str());
      h_truth_fit_up->SetMarkerSize(1.3);
      h_truth_fit_up->SetMarkerColor(kRed);
      h_truth_fit_up->SetLineColor(kRed);
      h_truth_fit_up->SetLineWidth(2);
      h_truth_fit_up->GetXaxis()->SetLimits(-400,400);
      h_truth_fit_up->GetYaxis()->SetRangeUser(-400,400);
      h_truth_fit_up->Draw("AP");
      TFitResultPtr fit_up = h_truth_fit_up->Fit("pol1","QS");
      TFitResultPtr fit_crosscheck_up = h_truth_fit_crosscheck_up->Fit("pol1","QNS");
      h_truth_fit_up->GetFunction("pol1")->SetLineWidth(2);
      h_truth_fit_up->GetFunction("pol1")->SetLineColor(kRed);
      h_truth_fit_down->SetMarkerColor(kBlue);
      h_truth_fit_down->SetMarkerStyle(25);
      h_truth_fit_down->SetMarkerSize(1.3);
      h_truth_fit_down->SetLineColor(kBlue);
      h_truth_fit_down->SetLineWidth(2);
      h_truth_fit_down->Draw("P");
      TFitResultPtr fit_down = h_truth_fit_down->Fit("pol1","QS");
      TFitResultPtr fit_crosscheck_down = h_truth_fit_crosscheck_down->Fit("pol1","QNS");
      h_truth_fit_down->GetFunction("pol1")->SetLineWidth(2);
      h_truth_fit_down->GetFunction("pol1")->SetLineColor(kBlue);
      TLegend* leg = new TLegend(0.22,0.72,0.6,0.82);
      ostringstream legup,legdown;
      legup   << "#Delta" << type[n] << " low-to-high (Slope:" << fixed << setprecision(3) << fit_up  ->Parameter(1) << " )";
      legdown << "#Delta" << type[n] << " high-to-low (Slope:" << fixed << setprecision(3) << fit_down->Parameter(1) << " )";
      leg->AddEntry(h_truth_fit_up  ,legup  .str().c_str(),"lp");
      leg->AddEntry(h_truth_fit_down,legdown.str().c_str(),"lp");
#ifdef __CINT__
      SetLegAtt(leg);
      CMSText(1,1,1);
#endif
      leg->Draw("SAME");
      can2->SaveAs((string("plots/vdm_length_scale_")+type[n]+string("_2.pdf")).c_str());

      double average = (fit_up->Parameter(1) + fit_down->Parameter(1)) / 2.;
      double average_cc = (fit_crosscheck_up->Parameter(1) + fit_crosscheck_down->Parameter(1)) / 2.;
      double stat = sqrt( pow(fit_up->Parameter(1)/2.,2) * pow(fit_down->ParError(1),2) + pow(fit_down->Parameter(1)/2.,2) * pow(fit_up->ParError(1),2));
      double sysupdown = fabs(fit_up->Parameter(1) - fit_down->Parameter(1))/ 2.;
      double sysfit = fabs(average_cc - average) / 2.;
      double sysskip = (type[n]=="X"?0.001:0.003);
      double sys = sqrt ( sysfit*sysfit + sysupdown*sysupdown + sysskip*sysskip);

      cout << endl << endl << endl;
      cout << fixed << setprecision(3)
           << "        Î”" << type[n] << endl
           << "--Correction Factor: " << average << endl
           << "up: " << fit_up->Parameter(1) << endl
           << "down: " << fit_down->Parameter(1) << endl
           << "--Stat: " << stat << endl
           << "--Sys Combined: " << sys << endl
           << "Sys. Up Down: " << sysupdown << endl
           << "Sys. Fitting: " << sysfit << endl
           << "Sys. Skipping: " << sysskip << endl;
      cout << endl << endl << endl;
    }

  can1->SaveAs((string("plots/vdm_length_scale")+string("_1.pdf")).c_str());
  canFind->SaveAs((string("plots/vdm_length_scale")+string("_sections.pdf")).c_str());
}

void FindSections(TH1D* h_rate, vector<double> &sectionBegin, vector<double> &sectionEnd, const double& begin, const double& end, const double offset, const double skip)
{
  if(begin < h_rate->GetBinLowEdge(1) || end > h_rate->GetBinLowEdge(h_rate->GetNbinsX()+1))
    {
      cerr << "Borders outside of histogram " << begin << " " <<  h_rate->GetBinLowEdge(1) << " " << end << " " << h_rate->GetBinLowEdge(h_rate->GetNbinsX()+1) << endl;
      exit(-1);
    }
  
  //Get Median
  int n = h_rate->GetXaxis()->GetNbins(); 
  //std::vector<double>  x(n);
  //h_rate->GetXaxis()->GetCenter(&x[0]);
  double * y = h_rate->GetArray();
  // exclude underflow/overflows from bin content array y
  double median = TMath::Median(n, &y[0]);

  //Define cut values (this is a bit arbitrary)
  const double cut_up = median * 1.2;
  const double cut_down = median * 0.5;
  cout << "Using cut values = " << cut_down << " < " << median << " < " << cut_up << endl;

  //Search in histogram
  const double step = h_rate->GetBinWidth(1) / 100.;
  bool in_section = false;

  for(double current = begin; current < end; current += step)
    {
      const double value = h_rate->Interpolate(current);

      //Start a section
      if (!in_section && (cut_down <= value && value <= cut_up))
        {
          sectionBegin.push_back(current+offset);
          in_section = true;
        }

      //End a section
      if (in_section && (value < cut_down || cut_up < value))
        {
          in_section = false;
          sectionEnd.push_back(current-offset);
        }
    }

  if(in_section)
    sectionEnd.push_back(end);
    
  if(sectionBegin.size() != sectionEnd.size())
    {
      cerr << "sections size error" << endl;
      exit(-1);
    }

  ostringstream name; name << "rate_" << begin << "_" << end;
  TH1D* draw = (TH1D*)(h_rate->Clone(name.str().c_str()));
  draw->GetXaxis()->SetRangeUser(begin,end);
  draw->Draw("HIST L");
  for (int i=0; i < int(sectionBegin.size()); ++i)
    {
      TBox* boxskip = new TBox(sectionBegin[i],cut_down/2.,sectionBegin[i]+skip,cut_down);
      boxskip->SetFillColor(kRed);
      boxskip->SetFillStyle(3001);
      boxskip->DrawClone();
//       TBox* box = new TBox(sectionBegin[i]+skip,cut_down/2.,sectionEnd[i],cut_down);
//       box->SetFillColor(kGreen-2);
//       box->SetFillStyle(3001);
//       box->DrawClone();
      TPaveText* txt = new TPaveText(sectionBegin[i]+skip,cut_down/2.,sectionEnd[i],cut_down,"b t l");
      ostringstream text; text << i;
      txt->AddText(text.str().c_str());
      txt->SetFillStyle(3001);
      txt->SetFillColor(kGreen-2);
      txt->Draw("SAME");
    }
  return;
}
