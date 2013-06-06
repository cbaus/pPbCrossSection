#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iostream>
#include <sstream>

using namespace std;

#define _DoNoisePerLS 0
#define _LumiCorr 1.0925

const double fac_error = 1.0;
const double cut_value_single = 6.*fac_error;
const double cut_value_double = 2.*fac_error;

void makePlots_cs1()
{
  TH1::SetDefaultSumw2();
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
#endif

  vector<int> run_num;
  run_num.push_back(210614);
  run_num.push_back(210885);
  run_num.push_back(210998);
  run_num.push_back(211000);
  run_num.push_back(211001);
  run_num.push_back(211032);
  run_num.push_back(211256);
  run_num.push_back(211371);
  run_num.push_back(211390);
  run_num.push_back(211460);
  run_num.push_back(211532);
  run_num.push_back(211538);
  run_num.push_back(211607);

  //Histogram with all runs
  TH1D* h_runs_single = new TH1D("h_runs_single","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_single");
  TH1D* h_runs_double = new TH1D("h_runs_double","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_double");
  TH1D* h_runs_single_pPb = new TH1D("h_runs_single_pPb","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_single");
  TH1D* h_runs_double_pPb = new TH1D("h_runs_double_pPb","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_double");
  TH1D* h_runs_single_Pbp = new TH1D("h_runs_single_Pbp","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_single");
  TH1D* h_runs_double_Pbp = new TH1D("h_runs_double_Pbp","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_double");


  TFile f("plots/corr_factors.root");
  TVectorD* f_em = NULL;
  TVectorD* f_mc = NULL;
  TVectorD* f_eme = NULL;
  TVectorD* f_mce = NULL;
  f_em  = (TVectorD*)f.Get("corr_fac_em");
  f_mc  = (TVectorD*)f.Get("corr_fac_mc");
  f_eme = (TVectorD*)f.Get("corr_fac_eme");
  f_mce = (TVectorD*)f.Get("corr_fac_mce");
  if(!f_mc || !f_em || !f_mce || !f_eme) {cerr << "error" << endl; return;}
  cout << "MC Eff. Correction:" << endl;
  f_mc->Print();
  f_mce->Print();
  cout << "EM Eff. Correction:" << endl;
  f_em->Print();
  f_eme->Print();

  double w = 0;
  double xw = 0;
  double sigma_em_runs_single=0;
  double sigma_mc_runs_single=0;
  double sigma_pu_runs_single=0;
  double sigma_em_runs_double=0;
  double sigma_mc_runs_double=0;
  double sigma_pu_runs_double=0;

  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos_old.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_single_";
      TH1D* h_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_double_";
      TH1D* h_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_lumi_";
      TH1D* h_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single_noise";
      TH1D* h_noise_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double_noise";
      TH1D* h_noise_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single";
      TH1D* h_zb_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double";
      TH1D* h_zb_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_all_lumi";
      TH1D* h_hf_noise_all_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_selected_single_lumi";
      TH1D* h_hf_noise_selected_single_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_selected_double_lumi";
      TH1D* h_hf_noise_selected_double_lumi=(TH1D*)file->Get(filename_ss.str().c_str());




      //Noise
      double f_noise_single = h_noise_single->GetBinContent(h_noise_single->FindBin(cut_value_single))/h_noise_single->GetBinContent(1);
      double f_noise_double = h_noise_double->GetBinContent(h_noise_double->FindBin(cut_value_double))/h_noise_double->GetBinContent(1);

      if(_DoNoisePerLS)
        {
          h_hf_noise_selected_single_lumi->Divide(h_hf_noise_all_lumi);
          h_hf_noise_selected_single_lumi->Rebin(20);
          h_hf_noise_selected_single_lumi->Scale(1./20.);
          h_hf_noise_selected_double_lumi->Divide(h_hf_noise_all_lumi);
          h_hf_noise_selected_double_lumi->Rebin(20);
          h_hf_noise_selected_double_lumi->Scale(1./20.);
        }

      //ZB
      const double f_zb_single = h_zb_single->GetBinContent(1)/h_zb_single->GetBinContent(h_zb_single->FindBin(cut_value_single));
      const double f_zb_double = h_zb_double->GetBinContent(1)/h_zb_double->GetBinContent(h_zb_double->FindBin(cut_value_double));


      double sigma_em_run_single=0;
      double sigma_mc_run_single=0;
      double sigma_pu_run_single=0;
      double sigma_em_run_double=0;
      double sigma_mc_run_double=0;
      double sigma_pu_run_double=0;
      int n = 0;

      //Loop over LumiSections
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
          //Noise per LumiSection
          if(_DoNoisePerLS)
            {
              cout << "f_noise=" << f_noise_single;
              f_noise_single = h_hf_noise_selected_single_lumi->GetBinContent(h_hf_noise_selected_single_lumi->FindBin(h_lumi->GetBinCenter(i)));
              cout << " f_noise=" << f_noise_single << endl;
              f_noise_double = h_hf_noise_selected_double_lumi->GetBinContent(h_hf_noise_selected_double_lumi->FindBin(h_lumi->GetBinCenter(i)));
            }
          //PILEUP!!
          const double lumiPerLS=h_lumi->GetBinContent(i);
          const double lumiPerLS_error=h_lumi->GetBinError(i);
          const double lambda = lumiPerLS*2.1/11246./23.31/296.;
          const double p1_single = (*f_mc)[0];
          const double p1_double = (*f_mc)[1];
          double f_pileup_single = 1;
          double f_pileup_double = 1;
          for (int k=2; k<100; k++)
            {
              double add_single =  p1_single/(1.-pow(1.-p1_single,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
              double add_double =  p1_double/(1.-pow(1.-p1_double,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
              f_pileup_single += add_single;
              f_pileup_double += add_double;
              if(add_single < 0.0001 && add_double < 0.0001)
                break;
            }
          //cout << "lambda=" << lambda << " " << f_pileup_single << endl;

          //Other Corrections
          if(lumiPerLS>0.)
            {
              const double n_cut_single = h_single->GetBinContent(i) / lumiPerLS / _LumiCorr;
              const double n_zb_single = (f_zb_single*n_cut_single);
              const double n_noise_single = f_noise_single * n_zb_single;
              const double n_em_single = (*f_em)[0] * n_cut_single;

              const double n_cut_double = h_double->GetBinContent(i) / lumiPerLS / _LumiCorr;
              const double n_zb_double = (f_zb_double*n_cut_double);
              const double n_noise_double = f_noise_double * n_zb_double;
              const double n_em_double = (*f_em)[1] * n_cut_double;

              if(i<200)
                {
                }
              const double A_single = (1./f_pileup_single) - f_noise_single;
              const double A_double = (1./f_pileup_double) - f_noise_double;
              if( A_double != 0 && A_single != 0)
                {
                  double n_single = ((n_cut_single/(*f_mc)[0]) - n_noise_single - n_em_single) / A_single;
                  double n_double = ((n_cut_double/(*f_mc)[1]) - n_noise_double - n_em_double) / A_double;
                  double error_single = sqrt( pow(h_single->GetBinContent(i),2)/pow(lumiPerLS,4)*pow(lumiPerLS_error,2) + 1./pow(lumiPerLS,2)*pow(h_single->GetBinError(i),2));
                  double error_double = sqrt( pow(h_double->GetBinContent(i),2)/pow(lumiPerLS,4)*pow(lumiPerLS_error,2) + 1./pow(lumiPerLS,2)*pow(h_double->GetBinError(i),2));
                  double sigma_mc_single = (*f_mce)[0] / A_single * n_cut_single / pow((*f_mc)[0],2);
                  double sigma_mc_double = (*f_mce)[1] / A_double * n_cut_double / pow((*f_mc)[1],2);
                  //should use n_single which is n_inel right?
                  double sigma_em_single = sqrt(pow((*f_eme)[0],2)+pow(0.2*(*f_em)[0],2));  sigma_em_single *= n_single / A_single;
                  double sigma_em_double = sqrt(pow((*f_eme)[1],2)+pow(0.2*(*f_em)[1],2));  sigma_em_double *= n_double / A_double;
                  double sigma_pu_single = sqrt(pow(((n_cut_single/(*f_mc)[0]) - n_noise_single - n_em_single)/pow(f_pileup_single,2),2) * pow(0.2*(f_pileup_single-1),2));
                  double sigma_pu_double = sqrt(pow(((n_cut_double/(*f_mc)[1]) - n_noise_double - n_em_double)/pow(f_pileup_double,2),2) * pow(0.2*(f_pileup_double-1),2));

                  h_single->SetBinContent(i, n_single);
                  h_double->SetBinContent(i, n_double);
                  h_single->SetBinError(i, error_single);
                  h_double->SetBinError(i, error_double);

                  if(i<00)
                    {
                      cout << n_cut_single<< " " << n_cut_single/(*f_mc)[0]<< " " << n_noise_single<< " " << n_em_single<< " " << ((1./f_pileup_single) - f_noise_single) << " " << f_pileup_single << endl;
                      cout << n_zb_double << " " << n_cut_double<< " " << n_cut_double/(*f_mc)[1]<< " " << n_noise_double<< " " << n_em_double<< " " << ((1./f_pileup_double) - f_noise_double) << endl;
                      cout << "sigma_em_single=" << sigma_em_single << " sigma_mc_single=" << sigma_mc_single << endl;
                      cout << "sigma_em_double=" << sigma_em_double << " sigma_mc_double=" << sigma_mc_double << endl;
                      cout << "pileuperror (" << f_pileup_double-1 << ")= " << sigma_pu_double << endl;
                    }

                  n++;
                  sigma_mc_run_single += sigma_mc_single;
                  sigma_mc_run_double += sigma_mc_double;
                  sigma_em_run_single += sigma_em_single;
                  sigma_em_run_double += sigma_em_double;
                  sigma_pu_run_single += sigma_pu_single;
                  sigma_pu_run_double += sigma_pu_double;
                }
              else {cerr << "div 0 in bin: " << i << endl; return;}
            }
          else if (lumiPerLS<0.) {cerr << "lumi neg: " << i << endl; return;}
        }

      sigma_em_run_single /= double(n);
      sigma_em_run_double /= double(n);
      sigma_mc_run_single /= double(n);
      sigma_mc_run_double /= double(n);
      sigma_pu_run_single /= double(n);
      sigma_pu_run_double /= double(n);

      sigma_em_runs_single += sigma_em_run_single;
      sigma_em_runs_double += sigma_em_run_double;
      sigma_mc_runs_single += sigma_mc_run_single;
      sigma_mc_runs_double += sigma_mc_run_double;
      sigma_pu_runs_single += sigma_pu_run_single;
      sigma_pu_runs_double += sigma_pu_run_double;

      cout << "sigma_em_single=" << sigma_em_run_single/2.16*100 << "%"
           << " sigma_mc_single=" << sigma_mc_run_single/2.16*100 << "%"
           << " sigma_pu_single=" << sigma_pu_run_single/2.16*100  << "%"<< endl;
      cout << "sigma_em_double=" << sigma_em_run_double/2.16*100 << "%"
           << " sigma_mc_double=" << sigma_mc_run_double/2.16*100 << "%"
           << " sigma_pu_single=" << sigma_pu_run_double/2.16*100  << "%"<< endl;

      h_double->SetMarkerColor(kRed);
      h_double->SetLineColor(kRed);

      h_single->SetLineWidth(2);
      h_double->SetLineWidth(2);

      h_single->SetTitle("single-arm;lumisection;#sigma_{inel} / b");
      h_double->SetTitle("double-arm;lumisection;#sigma_{inel} / b");

      TH1D* projection_single =  new TH1D("projection_single","",50,0,5);
      TH1D* projection_double =  new TH1D("projection_double","",50,0,5);
      projection_double->SetLineColor(kRed);
      projection_double->SetMarkerColor(kRed);
      TCanvas* c1 = NULL;
      if(run_num[run] == 210885)
        c1 = new TCanvas;

      TFitResultPtr fit_single = h_single->Fit("pol0","QS");
      TFitResultPtr fit_double = h_double->Fit("pol0","QS");

      h_single->GetXaxis()->SetRangeUser(0,1200);
      h_double->GetXaxis()->SetRangeUser(0,1200);
      //h_single->GetYaxis()->SetRangeUser(2.0,2.8);
      //h_double->GetYaxis()->SetRangeUser(2.0,2.8);

      if(run_num[run] == 210885)
        {
          h_single->Draw();
          h_double->Draw("SAME");
#ifdef __CINT__
          CMSPreliminary();
          DataText(false,true,"run #210885");
#endif
          c1->SaveAs((string("plots/CS_run")+string(".pdf")).c_str());

          for(int bin=1; bin<=h_single->GetNbinsX(); bin++)
            {
              if(h_single->GetBinError(bin) > 0)
                projection_single->Fill(h_single->GetBinContent(bin),1./pow(h_single->GetBinError(bin),2));
              if(h_double->GetBinError(bin) > 0)
                projection_double->Fill(h_double->GetBinContent(bin),1./pow(h_double->GetBinError(bin),2));
            }
          c1 = new TCanvas;
          projection_single->Draw("HIST");
          projection_double->Draw("HIST SAME");
          projection_single->SetTitle("single-arm;#sigma_{inel} in b;weighted count");
          projection_double->SetTitle("double-arm;#sigma_{inel} in b;weighted count");
          TLegend* leg = c1->BuildLegend(0.55,0.7,0.85,0.9);
          SetLegAtt(leg);
          leg->Draw();
          c1->SaveAs((string("plots/CS_run_proj")+string(".pdf")).c_str());
        }

      double weight = 1./pow(fit_double->ParError(0),2);
      xw  += fit_double->Parameter(0) * weight;
      w += weight;

      h_runs_single->SetBinContent(run+1,run,fit_single->Parameter(0));
      h_runs_single->SetBinError(run+1,0,fit_single->ParError(0));

      h_runs_double->SetBinContent(run+1,run,fit_double->Parameter(0));
      h_runs_double->SetBinError(run+1,0,fit_double->ParError(0));

      if(run_num[run] <= 211256) //pPb or Pbp?
        {
          h_runs_single_pPb->SetBinContent(run+1,run,fit_single->Parameter(0));
          h_runs_single_pPb->SetBinError(run+1,0,fit_single->ParError(0));

          h_runs_double_pPb->SetBinContent(run+1,run,fit_double->Parameter(0));
          h_runs_double_pPb->SetBinError(run+1,0,fit_double->ParError(0));
        }
      else
        {
          h_runs_single_Pbp->SetBinContent(run+1,run,fit_single->Parameter(0));
          h_runs_single_Pbp->SetBinError(run+1,0,fit_single->ParError(0));

          h_runs_double_Pbp->SetBinContent(run+1,run,fit_double->Parameter(0));
          h_runs_double_Pbp->SetBinError(run+1,0,fit_double->ParError(0));
        }
      h_runs_double_Pbp->GetXaxis()->SetBinLabel(run+1,runname.c_str());
      h_runs_double_pPb->GetXaxis()->SetBinLabel(run+1,runname.c_str());
      cout << "Single: " << fit_single->Parameter(0) << " " << fit_single->ParError(0) << " --- "
           << "Double: " << fit_double->Parameter(0) << " " << fit_double->ParError(0) << endl;

      if(run_num[run] == 210885)
        c1 = new TCanvas; //please use a new one to paint your fit crap
    }
  double sigma = sqrt(1./w);
  double average = xw/w;

  int n = run_num.size();
  sigma_em_runs_single /= double(n);
  sigma_em_runs_double /= double(n);
  sigma_mc_runs_single /= double(n);
  sigma_mc_runs_double /= double(n);
  sigma_pu_runs_single /= double(n);
  sigma_pu_runs_double /= double(n);

  h_runs_double_pPb->GetYaxis()->SetRangeUser(1.8,2.6);
  h_runs_double_pPb->SetLineColor(kRed);
  h_runs_double_pPb->SetMarkerColor(kRed);
  h_runs_double_pPb->SetTitle("double-arm (pPb);run number;#sigma_{inel} / b");

  h_runs_single_pPb->SetLineColor(kRed-1);
  h_runs_single_pPb->SetMarkerColor(kRed-1);
  h_runs_single_pPb->SetMarkerStyle(34);
  h_runs_single_pPb->SetTitle("single-arm (pPb)");

  h_runs_double_Pbp->SetLineColor(kBlue);
  h_runs_double_Pbp->SetMarkerColor(kBlue);
  h_runs_double_Pbp->SetTitle("double-arm (Pbp);run number;#sigma_{inel} / b");

  h_runs_single_Pbp->SetLineColor(kBlue-1);
  h_runs_single_Pbp->SetMarkerColor(kBlue-1);
  h_runs_single_Pbp->SetMarkerStyle(34);
  h_runs_single_Pbp->SetTitle("single-arm (Pbp)");



  TCanvas* can3 = new TCanvas;
  h_runs_double_pPb->Draw("P");
  h_runs_single_pPb->Draw("PSAME");
  h_runs_double_Pbp->Draw("PSAME");
  h_runs_single_Pbp->Draw("PSAME");
  TLegend* leg3 = can3->BuildLegend(0.5,0.65,0.8,0.85);
  SetLegAtt(leg3);
  leg3->Draw();
#ifdef __CINT__
  CMSPreliminary();
  DataText(true,true);
#endif
  can3->SaveAs((string("plots/CS_runs")+string(".pdf")).c_str());

  TCanvas* can4 = new TCanvas;
  h_runs_double_pPb->Draw("P");
  h_runs_double_Pbp->Draw("PSAME");
  TLegend* leg4 = can4->BuildLegend(0.5,0.65,0.8,0.85);
  SetLegAtt(leg4);
  leg4->Draw();
#ifdef __CINT__
  CMSPreliminary();
  DataText(true,true);
#endif
  can4->SaveAs((string("plots/CS_runs_pas")+string(".pdf")).c_str());


  TFitResultPtr fit_runs_single = h_runs_single->Fit("pol0","S");
  TFitResultPtr fit_runs_double = h_runs_double->Fit("pol0","S");
  double sigmainel = (fit_runs_single->Parameter(0)+fit_runs_double->Parameter(0))/2.;
  double runrun = (fit_runs_single->ParError(0)+fit_runs_double->ParError(0))/2.;
  cout << "sigma_inel: " << sigmainel << " b"<< endl;
  cout << "Run-by-run variation: " << runrun << " = " << runrun/sigmainel*100 << "%" << endl;

  if(fabs(average-fit_runs_double->Parameter(0))/(average+fit_runs_double->Parameter(0))*2 > 0.001)
    {
      cerr << "The weighted average is more than 1 per mille different from the fit" < endl;
      return;
    }


  cout << "sigma_em_single=" << sigma_em_runs_single/fit_runs_single->Parameter(0)*100 << "%"
       << " sigma_mc_single=" << sigma_mc_runs_single/fit_runs_single->Parameter(0)*100 << "%"
       << " sigma_pu_single=" << sigma_pu_runs_single/fit_runs_single->Parameter(0)*100  << "%"<< endl;
  cout << "sigma_em_double=" << sigma_em_runs_double/fit_runs_double->Parameter(0)*100 << "%"
       << " sigma_mc_double=" << sigma_mc_runs_double/fit_runs_double->Parameter(0)*100 << "%"
       << " sigma_pu_single=" << sigma_pu_runs_double/fit_runs_double->Parameter(0)*100  << "%"<< endl;
  cout << "sigma_em=" << (sigma_em_runs_single+sigma_em_runs_double)/2./sigmainel*100 << "%"
       << " sigma_mc=" << (sigma_mc_runs_single+sigma_mc_runs_double)/2./sigmainel*100 << "%"
       << " sigma_pu=" << (sigma_pu_runs_single+sigma_pu_runs_double)/2./sigmainel*100  << "%"
       << " sigma_combine=" << fabs(fit_runs_single->Parameter(0)-fit_runs_double->Parameter(0))/2./sigmainel*100 << "%"<< endl;

  //PULL Distribution
  TH1D* h_pull_single = new TH1D("h_pull_single","single-arm;#frac{x-#mu}{#sigma};N_{run}",19,-8,8);
  TH1D* h_pull_double = new TH1D("h_pull_double","double-arm;#frac{x-#mu}{#sigma};N_{run}",19,-8,8);
  for (int run=0; run<int(run_num.size()); run++)
    {
      const double pull_single = (h_runs_single->GetBinContent(run+1) - fit_runs_single->Parameter(0) ) / h_runs_single->GetBinError(run+1);
      const double pull_double = (h_runs_double->GetBinContent(run+1) - fit_runs_double->Parameter(0) ) / h_runs_double->GetBinError(run+1);
      h_pull_single->Fill(pull_single);
      h_pull_double->Fill(pull_double);
    }
  h_pull_double->SetLineColor(kRed);
  h_pull_single->SetLineColor(kBlack);
  h_pull_double->SetMarkerColor(kRed);
  h_pull_single->SetMarkerColor(kBlack);
  TCanvas* can5 = new TCanvas;
  TF1* f_pull_single = new TF1("pull_single","gaus",-20,20);
  TF1* f_pull_double = new TF1("pull_double","gaus",-20,20);
  TFitResultPtr fit_pull_single = h_pull_single->Fit(f_pull_single,"SWL");
  TFitResultPtr fit_pull_double = h_pull_double->Fit(f_pull_double,"SWL");
  h_pull_single->Draw("HIST");
  h_pull_single->GetYaxis()->SetRangeUser(0,h_pull_double->GetMaximum()*1.4);
  h_pull_double->Draw("HIST SAME");
  TLegend* leg5 = can5->BuildLegend(0.6,0.7,0.85,0.9);
  SetLegAtt(leg5);
  f_pull_single->SetLineWidth(2);
  f_pull_double->SetLineWidth(2);
  f_pull_single->SetLineColor(kBlack);
  f_pull_double->SetLineColor(kRed);
  f_pull_single->Draw("SAME L");
  f_pull_double->Draw("SAME L");
  leg5->Draw();
#ifdef __CINT__
  CMSPreliminary();
  DataText(true,true);
#endif
  can5->SaveAs((string("plots/CS_pull")+string(".pdf")).c_str());

}
