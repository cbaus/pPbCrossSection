#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMarker.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

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

  TVectorD vec_sigma_vis(2);
  TVectorD vec_sigma_had(2);
  TVectorD vec_sigma_inel(3);
  TVectorD vec_sigma_vis_e(2);
  TVectorD vec_sigma_had_e(2);
  TVectorD vec_sigma_inel_e(3);

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
  TH1D* h_runs_vis_single = new TH1D("h_runs_vis_single","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_single");
  TH1D* h_runs_vis_double = new TH1D("h_runs_vis_double","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_double");
  TH1D* h_runs_had_single = new TH1D("h_runs_had_single","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_single");
  TH1D* h_runs_had_double = new TH1D("h_runs_had_double","",run_num.size(),-0.5,run_num.size()-0.5);// h_runs_single->SetName("h_runs_double");
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
  double n_noise_runs_single=0;
  double n_em_runs_single=0;
  double A_runs_single=0;
  double A_runs_vis_single=0;
  double n_noise_runs_double=0;
  double n_em_runs_double=0;
  double A_runs_double=0;
  double A_runs_vis_double=0;
  double sigma_em_runs_single=0;
  double sigma_mc_runs_single=0;
  double sigma_pu_runs_single=0;
  double sigma_em_runs_double=0;
  double sigma_mc_runs_double=0;
  double sigma_pu_runs_double=0;

  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos_skiprings.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_single_6";
      TH1D* h_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_double_2";
      TH1D* h_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_lumi";
      TH1D* h_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single_noise";
      TH1D* h_noise_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double_noise";
      TH1D* h_noise_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single";
      TH1D* h_zb_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double";
      TH1D* h_zb_double=(TH1D*)file->Get(filename_ss.str().c_str());

      // filename_ss.str(""); filename_ss << "Epos" << runname << "/Epos" << runname << "_h_hf_cut_single";
      // TH1D* h_epos_single=(TH1D*)file->Get(filename_ss.str().c_str());
      // filename_ss.str(""); filename_ss << "Epos" << runname << "/Epos" << runname << "_h_hf_cut_double";
      // TH1D* h_epos_double=(TH1D*)file->Get(filename_ss.str().c_str());
      // filename_ss.str(""); filename_ss << "QGSJetII" << runname << "/QGSJetII" << runname << "_h_hf_cut_single";
      // TH1D* h_epos_single=(TH1D*)file->Get(filename_ss.str().c_str());
      // filename_ss.str(""); filename_ss << "QGSJetII" << runname << "/QGSJetII" << runname << "_h_hf_cut_double";
      // TH1D* h_epos_double=(TH1D*)file->Get(filename_ss.str().c_str());

      TH1D *h_vis_single = (TH1D*)h_single->Clone("h_vis_single");
      TH1D *h_had_single = (TH1D*)h_single->Clone("h_had_single");
      TH1D *h_vis_double = (TH1D*)h_double->Clone("h_vis_double");
      TH1D *h_had_double = (TH1D*)h_double->Clone("h_had_double");


      ///////////////////////////////////////////////////////////////////////
      //Acceptance
      const double eff_acc_single = (*f_mc)[0];
      const double eff_acc_double = (*f_mc)[1];
      const double eff_acc_single_e = (*f_mce)[0];
      const double eff_acc_double_e = (*f_mce)[1];

      ///////////////////////////////////////////////////////////////////////
      //EM
      const double eff_em_single = (*f_em)[0];
      const double eff_em_double = (*f_em)[1];
      const double eff_em_single_e = (*f_eme)[0];
      const double eff_em_double_e = (*f_eme)[1];

      ///////////////////////////////////////////////////////////////////////
      //Noise
      double f_noise_single = h_noise_single->GetBinContent(h_noise_single->FindBin(cut_value_single))/h_noise_single->GetBinContent(1);
      double f_noise_double = h_noise_double->GetBinContent(h_noise_double->FindBin(cut_value_double))/h_noise_double->GetBinContent(1);

      ///////////////////////////////////////////////////////////////////////
      //ZB
      const double f_zb_single = h_zb_single->GetBinContent(1)/h_zb_single->GetBinContent(h_zb_single->FindBin(cut_value_single));
      const double f_zb_double = h_zb_double->GetBinContent(1)/h_zb_double->GetBinContent(h_zb_double->FindBin(cut_value_double));


      double n_noise_run_single=0;
      double n_em_run_single=0;
      double A_run_single=0;
      double A_run_vis_single=0;
      double n_noise_run_double=0;
      double n_em_run_double=0;
      double A_run_double=0;
      double A_run_vis_double=0;
      double sigma_em_run_single=0;
      double sigma_mc_run_single=0;
      double sigma_pu_run_single=0;
      double sigma_em_run_double=0;
      double sigma_mc_run_double=0;
      double sigma_pu_run_double=0;
      int n = 0;


      ///////////////////////////////////////////////////////////////////////
      //Loop over LumiSections
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
          ///////////////////////////////////////////////////////////////////////
          //PILEUP!!
          const double lumiPerLS=h_lumi->GetBinContent(i) * _LumiCorr;
          const double lumiPerLS_error=h_lumi->GetBinError(i) * _LumiCorr; //not 100% correct since from profile but has no contribution
          if (lumiPerLS<0.) {cerr << "lumi neg: " << i << endl; return;}
          else if (lumiPerLS==0.) {continue;}

          double nbunches = run_num[run]<=210638 ? 206 : 296.;

          const double lambda = lumiPerLS*2.13/11246./23.31/296.;
          // double f_pileup_single = 1;
          // double f_pileup_double = 1;
          // for (int k=2; k<100; k++)
          //   {
          //     double add_single =  eff_acc_single/(1.-pow(1.-eff_acc_single,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
          //     double add_double =  eff_acc_double/(1.-pow(1.-eff_acc_double,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
          //     f_pileup_single += add_single;
          //     f_pileup_double += add_double;
          //     if(add_single < 0.0001 && add_double < 0.0001)
          //       break;
          //   }
          double f_pileup_single = 0;
          double f_pileup_double = 0;
          double nomi_single = 0;
          double denomi_single = 0;
          double nomi_double = 0;
          double denomi_double = 0;
          for (int k=1; k<100; k++)
            {
              double nomi_add_single =  eff_acc_single * k * TMath::Poisson(k,lambda);
              double denomi_add_single =  (1.-pow(1.-eff_acc_single,k)) * 1 * TMath::Poisson(k,lambda);
              double nomi_add_double =  eff_acc_double * k * TMath::Poisson(k,lambda);
              double denomi_add_double =  (1.-pow(1.-eff_acc_double,k)) * 1 * TMath::Poisson(k,lambda);
              nomi_single += nomi_add_single;
              denomi_single += denomi_add_single;
              nomi_double += nomi_add_double;
              denomi_double += denomi_add_double;
              if(nomi_add_single < 0.0001 && nomi_add_double < 0.0001 && denomi_add_single < 0.0001 && denomi_add_double < 0.0001)
                break;
            }
          if(!denomi_single || !denomi_double)
            {
              cerr << "Error in PU calculation" << endl; return;
            }
          f_pileup_single = nomi_single / denomi_single;
          f_pileup_double = nomi_double / denomi_double;
          //cout << "lambda=" << lambda << " " << f_pileup_single << endl;


          ///////////////////////////////////////////////////////////////////////
          //Other Corrections

          const double n_cut_single = h_single->GetBinContent(i) / lumiPerLS;
          const double n_zb_single = (f_zb_single*n_cut_single);
          const double n_noise_single = f_noise_single * n_zb_single;
          const double n_em_single = eff_em_single * 0.195; //these are not n but already n/lumi

          const double n_cut_double = h_double->GetBinContent(i) / lumiPerLS;
          const double n_zb_double = (f_zb_double*n_cut_double);
          const double n_noise_double = f_noise_double * n_zb_double;
          const double n_em_double = eff_em_double * 0.195; //relates to sigma_inel

          const double A_single = (eff_acc_single/f_pileup_single) - f_noise_single;
          const double A_double = (eff_acc_double/f_pileup_double) - f_noise_double;
          const double A_vis_single = (  1.  /f_pileup_single) - f_noise_single;
          const double A_vis_double = (  1.  /f_pileup_double) - f_noise_double;

          if( A_double == 0 || A_single == 0) {cerr << "div 0 in bin: " << i << endl; return;}

          double n_single = (n_cut_single - n_noise_single - n_em_single) / A_single;
          double n_double = (n_cut_double - n_noise_double - n_em_double) / A_double;
          double n_vis_single = (n_cut_single - n_noise_single) / A_vis_single;
          double n_vis_double = (n_cut_double - n_noise_double) / A_vis_double;
          double n_had_single = (n_cut_single - n_noise_single - n_em_single) / A_vis_single;
          double n_had_double = (n_cut_double - n_noise_double - n_em_double) / A_vis_double;

          //stat error for histos
          double error_single = sqrt( pow(h_single->GetBinContent(i),2)/pow(lumiPerLS,4)*pow(lumiPerLS_error,2) + 1./pow(lumiPerLS,2)*pow(h_single->GetBinError(i),2));
          double error_double = sqrt( pow(h_double->GetBinContent(i),2)/pow(lumiPerLS,4)*pow(lumiPerLS_error,2) + 1./pow(lumiPerLS,2)*pow(h_double->GetBinError(i),2));

          //systematic errors according to PAS formulas
          double sigma_mc_single = fabs(-1 * n_single / ( eff_acc_single - f_noise_single * f_pileup_single) * eff_acc_single_e);
          double sigma_mc_double = fabs(-1 * n_double / ( eff_acc_double - f_noise_double * f_pileup_double) * eff_acc_double_e);
          double sigma_em_single = fabs(n_em_single / (eff_acc_single / f_pileup_single - f_noise_single) * eff_em_single_e);
          double sigma_em_double = fabs(n_em_double / (eff_acc_double / f_pileup_double - f_noise_double) * eff_em_double_e);
          double sigma_pu_single = fabs(n_single / (f_pileup_single - pow(f_pileup_single,2) * f_noise_single / eff_acc_single) * 0.2*(f_pileup_single-1));
          double sigma_pu_double = fabs(n_double / (f_pileup_double - pow(f_pileup_double,2) * f_noise_double / eff_acc_double) * 0.2*(f_pileup_double-1));

          h_single->SetBinContent(i, n_single);
          h_double->SetBinContent(i, n_double);
          h_vis_single->SetBinContent(i, n_vis_single);
          h_vis_double->SetBinContent(i, n_vis_double);
          h_had_single->SetBinContent(i, n_had_single);
          h_had_double->SetBinContent(i, n_had_double);
          h_single->SetBinError(i, error_single);
          h_double->SetBinError(i, error_double);
          h_vis_single->SetBinError(i, error_single);
          h_vis_double->SetBinError(i, error_double);
          h_had_single->SetBinError(i, error_single);
          h_had_double->SetBinError(i, error_double);

          if(i<00)
            {
              cout << n_zb_single << " " << n_cut_single<< " " << n_cut_single/eff_acc_single<< " " << n_noise_single<< " " << n_em_single<< " " << f_pileup_single  << " " << pow(f_pileup_single,2) * f_noise_single / eff_acc_single << endl;
              cout << n_zb_double << " " << n_cut_double<< " " << n_cut_double/eff_acc_double<< " " << n_noise_double<< " " << n_em_double<< " " << f_pileup_double  << " " << pow(f_pileup_double,2) * f_noise_double / eff_acc_double << endl;
              cout << "sigma_em_single=" << sigma_em_single << " sigma_mc_single=" << sigma_mc_single << endl;
              cout << "sigma_em_double=" << sigma_em_double << " sigma_mc_double=" << sigma_mc_double << endl;
              cout << "pileuperror (" << f_pileup_double-1 << ")= " << sigma_pu_double << endl;
              cout << "eff_em_single " << eff_em_single << endl;
            }

          n++;
          n_noise_run_single += n_noise_single;
          n_em_run_single    += n_em_single;
          A_run_single       += A_single;
          A_run_vis_single   += A_vis_single;
          n_noise_run_double += n_noise_double;
          n_em_run_double    += n_em_double;
          A_run_double       += A_double;
          A_run_vis_double   += A_vis_double;
          sigma_mc_run_single += sigma_mc_single;
          sigma_mc_run_double += sigma_mc_double;
          sigma_em_run_single += sigma_em_single;
          sigma_em_run_double += sigma_em_double;
          sigma_pu_run_single += sigma_pu_single;
          sigma_pu_run_double += sigma_pu_double;
        }//Loop over LumiSections

      n_noise_run_single /= double(n);
      n_em_run_single    /= double(n);
      A_run_single       /= double(n);
      A_run_vis_single   /= double(n);
      n_noise_run_double /= double(n);
      n_em_run_double    /= double(n);
      A_run_double       /= double(n);
      A_run_vis_double   /= double(n);

      n_noise_runs_single += n_noise_run_single;
      n_em_runs_single    += n_em_run_single;
      A_runs_single       += A_run_single;
      A_runs_vis_single   += A_run_vis_single;
      n_noise_runs_double += n_noise_run_double;
      n_em_runs_double    += n_em_run_double;
      A_runs_double       += A_run_double;
      A_runs_vis_double   += A_run_vis_double;

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

      h_single->SetTitle("E>6 GeV (single-arm);lumisection;#sigma_{inel} [b]");
      h_double->SetTitle("E>2 GeV (double-arm);lumisection;#sigma_{inel} [b]");

      TH1D* projection_single =  new TH1D("projection_single","",50,0,5);
      TH1D* projection_double =  new TH1D("projection_double","",50,0,5);
      projection_double->SetLineColor(kRed);
      projection_double->SetMarkerColor(kRed);
      TCanvas* c1 = NULL;
      if(run_num[run] == 210885)
        c1 = new TCanvas;

      TFitResultPtr fit_single = h_single->Fit("pol0","QS");
      TFitResultPtr fit_double = h_double->Fit("pol0","QS");
      TFitResultPtr fit_vis_single = h_vis_single->Fit("pol0","QS");
      TFitResultPtr fit_vis_double = h_vis_double->Fit("pol0","QS");
      TFitResultPtr fit_had_single = h_had_single->Fit("pol0","QS");
      TFitResultPtr fit_had_double = h_had_double->Fit("pol0","QS");

      h_single->GetXaxis()->SetRangeUser(0,1100);
      h_double->GetXaxis()->SetRangeUser(0,1100);
      //h_single->GetYaxis()->SetRangeUser(2.0,2.8);
      //h_double->GetYaxis()->SetRangeUser(2.0,2.8);

      if(run_num[run] == 210885)
        {
          h_single->Draw();
          h_double->Draw("SAME");
#ifdef __CINT__
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
          projection_single->GetYaxis()->SetRangeUser(0,TMath::Max(projection_single->GetMaximum(),projection_double->GetMaximum())*1.4);
          projection_single->GetXaxis()->SetRangeUser(0,fit_double->Parameter(0)*2.);
          projection_single->SetMarkerStyle(21);
          projection_double->SetMarkerStyle(22);
          projection_single->SetTitle("E>6 GeV (single-arm);#sigma_{inel} [b];weighted count");
          projection_double->SetTitle("E>2 GeV (double-arm);#sigma_{inel} [b];weighted count");

          TLegend* leg = new TLegend(0.4,0.75,0.85,0.9);
          leg->AddEntry(projection_single,"","lp");
          leg->AddEntry(projection_double,"","lp");
#ifdef __CINT__
          DataText(true,true);
          SetLegAtt(leg);
#endif
          TMarker *m_single = new TMarker(fit_single->Parameter(0),0.1,projection_single->GetMarkerStyle());
          TMarker *m_double = new TMarker(fit_double->Parameter(0),0,projection_double->GetMarkerStyle());

          m_single->SetMarkerColor(projection_single->GetMarkerColor());
          m_double->SetMarkerColor(projection_double->GetMarkerColor());

          projection_single->Draw("HIST");
          projection_double->Draw("HIST SAME");
          m_single->Draw("SAME");
          m_double->Draw("SAME");
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

      h_runs_vis_single->SetBinContent(run+1,run,fit_vis_single->Parameter(0));
      h_runs_vis_single->SetBinError(run+1,0,fit_vis_single->ParError(0));

      h_runs_vis_double->SetBinContent(run+1,run,fit_vis_double->Parameter(0));
      h_runs_vis_double->SetBinError(run+1,0,fit_vis_double->ParError(0));

      h_runs_had_single->SetBinContent(run+1,run,fit_had_single->Parameter(0));
      h_runs_had_single->SetBinError(run+1,0,fit_had_single->ParError(0));

      h_runs_had_double->SetBinContent(run+1,run,fit_had_double->Parameter(0));
      h_runs_had_double->SetBinError(run+1,0,fit_had_double->ParError(0));

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

  n_noise_runs_single /= double(n);
  n_em_runs_single /= double(n);
  A_runs_single /= double(n);
  A_runs_vis_single /= double(n);
  n_noise_runs_double /= double(n);
  n_em_runs_double /= double(n);
  A_runs_double /= double(n);
  A_runs_vis_double /= double(n);

  sigma_em_runs_single /= double(n);
  sigma_em_runs_double /= double(n);
  sigma_mc_runs_single /= double(n);
  sigma_mc_runs_double /= double(n);
  sigma_pu_runs_single /= double(n);
  sigma_pu_runs_double /= double(n);

  h_runs_double_pPb->GetYaxis()->SetRangeUser(1.9,2.7);
  h_runs_double_pPb->SetLineColor(kRed);
  h_runs_double_pPb->SetMarkerColor(kRed);
  h_runs_double_pPb->SetTitle("E>2 GeV (double-arm) (pPb);run number;#sigma_{inel} [b]");

  h_runs_single_pPb->SetLineColor(kBlue);
  h_runs_single_pPb->SetMarkerColor(kBlue);
  h_runs_single_pPb->SetMarkerStyle(34);
  h_runs_single_pPb->SetTitle("E>6 GeV (single-arm) (pPb)");

  h_runs_double_Pbp->SetLineColor(kRed-1);
  h_runs_double_Pbp->SetMarkerColor(kRed-1);
  h_runs_double_Pbp->SetTitle("E>2 GeV (double-arm) (Pbp);run number;#sigma_{inel} [b]");

  h_runs_single_Pbp->SetLineColor(kBlue-1);
  h_runs_single_Pbp->SetMarkerColor(kBlue-1);
  h_runs_single_Pbp->SetMarkerStyle(34);
  h_runs_single_Pbp->SetTitle("E>6 GeV (single-arm) (Pbp)");



  ///////////////////////////////////////////////////////////////////////
  //FITTING EVERY RUN
  TFitResultPtr fit_runs_vis_single = h_runs_vis_single->Fit("pol0","S");
  TFitResultPtr fit_runs_had_single = h_runs_had_single->Fit("pol0","S");
  TFitResultPtr fit_runs_single = h_runs_single->Fit("pol0","S");

  TFitResultPtr fit_runs_vis_double = h_runs_vis_double->Fit("pol0","S");
  TFitResultPtr fit_runs_had_double = h_runs_had_double->Fit("pol0","S");
  TFitResultPtr fit_runs_double = h_runs_double->Fit("pol0","S");

  if(fabs(average-fit_runs_double->Parameter(0))/(average+fit_runs_double->Parameter(0))*2 > 0.001 || // = a/((a+b)/2)
     fabs(fit_runs_double->ParError(0)-sigma)/(sigma+fit_runs_double->ParError(0))*2 > 0.001
     )
    {
      cerr << "The weighted average is more than 1 per mille different from the fit " << average << " " << sigma << endl;
      return;
    }

  ///////////////////////////////////////////////////////////////////////
  //RESULTS
  cout.precision(4);
  double sigmainel_double = fit_runs_double->Parameter(0);
  double sigmainel_single = fit_runs_single->Parameter(0);
  double sigmainel = (sigmainel_single + sigmainel_double)/2.;
  double run1=fit_runs_single->ParError(0);
  double run2=fit_runs_double->ParError(0);
  double runrun = (run1 + run2)/2.;
  cout << endl << endl
       << " !!! sigma_inel_single: " << sigmainel_single << " b +- (" << run1 << " = " << run1/sigmainel_single*100 << "% (stat))" << endl
       << " !!! sigma_inel_double: " << sigmainel_double << " b +- (" << run2 << " = " << run2/sigmainel_double*100 << "% (stat))" << endl
       << " !!! sigma_inel       : " << sigmainel << " b +- (" << runrun << " = " << runrun/sigmainel*100 << "% (stat))" << endl << endl;

  cout << " !! corr_acc_single=" << sigmainel_single << " b "
       << " !! corr_acc_double=" << n_em_runs_double/A_runs_double << " b " << endl << endl;
  cout << " !! corr_noise_single=" << n_noise_runs_single/A_runs_single << " b "
       << " !! corr_em_single=" << n_em_runs_single/A_runs_single<< " b " << endl;
  cout << " !! corr_noise_double=" << n_noise_runs_double/A_runs_double << " b "
       << " !! corr_em_double=" << n_em_runs_double/A_runs_double << " b " << endl << endl;


  cout << " ! sigma_em_single=" << sigma_em_runs_single/fit_runs_single->Parameter(0)*100 << "%"
       << " sigma_mc_single=" << sigma_mc_runs_single/fit_runs_single->Parameter(0)*100 << "%"
       << " sigma_pu_single=" << sigma_pu_runs_single/fit_runs_single->Parameter(0)*100  << "%"
       << " sigma_lumi=" << sigmainel_single * 0.08 << " b" << endl << endl;
  cout << " ! sigma_em_double=" << sigma_em_runs_double/fit_runs_double->Parameter(0)*100 << "%"
       << " sigma_mc_double=" << sigma_mc_runs_double/fit_runs_double->Parameter(0)*100 << "%"
       << " sigma_pu_double=" << sigma_pu_runs_double/fit_runs_double->Parameter(0)*100  << "%"
       << " sigma_lumi=" << sigmainel_double * 0.08 << " b" << endl << endl;
  cout << " ! sigma_em=" << (sigma_em_runs_single+sigma_em_runs_double)/2./sigmainel*100 << "%"
       << " sigma_mc=" << (sigma_mc_runs_single+sigma_mc_runs_double)/2./sigmainel*100 << "%"
       << " sigma_pu=" << (sigma_pu_runs_single+sigma_pu_runs_double)/2./sigmainel*100  << "%"
       << " sigma_combine=" << fabs(fit_runs_single->Parameter(0)-fit_runs_double->Parameter(0))/2./sigmainel*100 << "%"
       << " sigma_lumi=" << sigmainel * 0.08 << " b" << endl << endl;


  //Write to files
  vec_sigma_inel[0] = sigmainel;
  vec_sigma_inel[1] = sigmainel_single;
  vec_sigma_inel[2] = sigmainel_double;
  vec_sigma_inel_e[0] = runrun;
  vec_sigma_inel_e[1] = run1;
  vec_sigma_inel_e[2] = run2;

  vec_sigma_vis[0] = fit_runs_vis_single->Parameter(0);
  vec_sigma_vis[1] = fit_runs_vis_double->Parameter(0);
  vec_sigma_vis_e[0] = fit_runs_vis_single->ParError(0);
  vec_sigma_vis_e[1] = fit_runs_vis_double->ParError(0);

  vec_sigma_had[0] = fit_runs_had_single->Parameter(0);
  vec_sigma_had[1] = fit_runs_had_double->Parameter(0);
  vec_sigma_had_e[0] = fit_runs_had_single->ParError(0);
  vec_sigma_had_e[1] = fit_runs_had_double->ParError(0);

  TFile outfile("plots/final_values.root","recreate");
  vec_sigma_vis.Write("vec_sigma_vis");
  vec_sigma_vis_e.Write("vec_sigma_vis_e");
  vec_sigma_had.Write("vec_sigma_had");
  vec_sigma_had_e.Write("vec_sigma_had_e");
  vec_sigma_inel.Write("vec_sigma_inel");
  vec_sigma_inel_e.Write("vec_sigma_inel_e");
  outfile.Close();

  ///////////////////////////////////////////////////////////////////////
  //Drawing cs over run plot
  TCanvas* can3 = new TCanvas;
  h_runs_double_pPb->Draw("P");
  h_runs_single->GetFunction("pol0")->SetLineColor(kBlue);
  h_runs_double->GetFunction("pol0")->SetLineColor(kRed);
  h_runs_single->GetFunction("pol0")->Draw("SAME");
  h_runs_double->GetFunction("pol0")->Draw("SAME");
  h_runs_single_pPb->Draw("PSAME");
  h_runs_double_Pbp->Draw("PSAME");
  h_runs_single_Pbp->Draw("PSAME");
  TLegend* leg3 = new TLegend(0.3,0.60,0.8,0.80);
  leg3->AddEntry(h_runs_single_pPb,"","pl");
  leg3->AddEntry(h_runs_single_Pbp,"","pl");
  leg3->AddEntry(h_runs_double_pPb,"","pl");
  leg3->AddEntry(h_runs_double_Pbp,"","pl");
#ifdef __CINT__
  SetLegAtt(leg3);
#endif
  leg3->Draw();
#ifdef __CINT__
  DataText(true,true);
#endif
  can3->SaveAs((string("plots/CS_runs")+string(".pdf")).c_str());
  for(int bin=1; bin <= h_runs_double_pPb->GetNbinsX(); bin++)
    {
      ostringstream ss_bin; ss_bin << bin;
      h_runs_double_pPb->GetXaxis()->SetBinLabel(bin,ss_bin.str().c_str());
    }
  h_runs_double_pPb->GetXaxis()->SetLabelSize(0.07);
  can3->SaveAs((string("plots/CS_runs_pas")+string(".pdf")).c_str());



  ///////////////////////////////////////////////////////////////////////
  //PULL Distribution
  TH1D* h_pull_single = new TH1D("h_pull_single","E>6 GeV (single-arm);#frac{x-#mu}{#sigma};N_{run}",19,-8,8);
  TH1D* h_pull_double = new TH1D("h_pull_double","E>2 GeV (double-arm);#frac{x-#mu}{#sigma};N_{run}",19,-8,8);
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
  TFitResultPtr fit_pull_single = h_pull_single->Fit(f_pull_single,"QSWL");
  TFitResultPtr fit_pull_double = h_pull_double->Fit(f_pull_double,"QSWL");
  h_pull_single->Draw("HIST");
  h_pull_single->GetYaxis()->SetRangeUser(0,h_pull_double->GetMaximum()*1.45);
  h_pull_double->Draw("HIST SAME");
  TLegend* leg5 = can5->BuildLegend(0.47,0.74,0.92,0.88);
#ifdef __CINT__
  SetLegAtt(leg5);
#endif
  f_pull_single->SetLineWidth(2);
  f_pull_double->SetLineWidth(2);
  f_pull_single->SetLineColor(kBlack);
  f_pull_double->SetLineColor(kRed);
  f_pull_single->Draw("SAME L");
  f_pull_double->Draw("SAME L");
  leg5->Draw();
#ifdef __CINT__
  DataText(true,true);
#endif
  cout << "GAUSSIAN FITS (Pull)" << endl
       << "single: mean=" << fit_pull_single->Parameter(1) << " +- " << fit_pull_single->ParError(1) << " "
       << "sigma=" << fit_pull_single->Parameter(2) << " +- " << fit_pull_single->ParError(2) << endl
       << "double: mean=" << fit_pull_double->Parameter(1) << " +- " << fit_pull_double->ParError(1) << " "
       << "sigma=" << fit_pull_double->Parameter(2) << " +- " << fit_pull_double->ParError(2) << endl;
  can5->SaveAs((string("plots/CS_pull")+string(".pdf")).c_str());

}
