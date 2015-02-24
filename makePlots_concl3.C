//macro for vis - had - inel plot
//important please always update uncertainties by hand!
//makePlots_pt_cuts_eff_pur & makePlots_cs

#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TVectorD.h>

#include<iostream>
#include<iomanip>

using namespace std;

#define WITH_ALICE 1
#define EPOS_CS 2.081702 //ran again with more statistics. this introduced slight changes
#define QGS_CS 2.180836
#define DPM_CS 2.165823

template<typename T>
void SetAttributes(T* theGraph, int colour, int marker)
{
  theGraph->SetMarkerSize(1.5);
  theGraph->SetLineWidth(1.8);
  theGraph->SetFillColor(kWhite);
  theGraph->SetLineColor(colour);
  theGraph->SetMarkerColor(colour);
  theGraph->SetMarkerStyle(marker);
}

void makePlots_concl3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  gStyle->SetPadTopMargin(0.07);

  //Get Hadron Level uncertainty
  TFile f0("plots/corr_factors_hadron.root");
  TVectorD* cor_fac_had           = NULL;
  TVectorD* cor_fac_had_e         = NULL;
  cor_fac_had           = (TVectorD*)f0.Get("corr_fac_had");
  cor_fac_had_e         = (TVectorD*)f0.Get("corr_fac_had_e");
  cout << "Data had level correction:" << endl;
  cor_fac_had->Print();
  cor_fac_had_e->Print();

  //Uncertainty values
  double s_lumi = 3.5, d_lumi = 3.5;
  double s_pu = 0.1,   d_pu = 0.1;
  double s_acc = 0.5,  d_acc = 1.6;
  double s_diff = 0.8, d_diff = 1.1;
  double s_em = 0.6,   d_em = 0.1;
  double s_mod = 1.7,  d_mod = 0.8;
  double s_sel = 0.6,  d_sel = 0.2;
  double s_noi = 1.3,  d_noi = 0.2;
  double s_hadlvl = (*cor_fac_had_e)[0]*100., d_hadlvl = (*cor_fac_had_e)[1]*100.;

  double s_withoutl = sqrt(pow(s_pu,2)+pow(s_acc,2)+pow(s_diff,2)+pow(s_em,2)+pow(s_mod,2)+pow(s_sel,2)+pow(s_noi,2));
  double d_withoutl = sqrt(pow(d_pu,2)+pow(d_acc,2)+pow(d_diff,2)+pow(d_em,2)+pow(d_mod,2)+pow(d_sel,2)+pow(d_noi,2));

  double combined = sqrt( pow(sqrt( pow(s_withoutl,2)+pow(d_withoutl,2) )/2.,2) + pow(3.5,2) );
  double combined_withoutl = sqrt( pow(s_withoutl,2)+pow(d_withoutl,2) )/2.;

  double s_vis  = sqrt(pow(s_lumi,2)+pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2));
  double d_vis  = sqrt(pow(d_lumi,2)+pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2));
  double s_had  = sqrt(pow(s_lumi,2)+pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2) + pow(s_em,2) + pow(s_mod,2)+pow(s_hadlvl,2));
  double d_had  = sqrt(pow(d_lumi,2)+pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2) + pow(d_em,2) + pow(d_mod,2)+pow(d_hadlvl,2));
  double s_inel = sqrt(pow(s_lumi,2)+pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2) + pow(s_em,2) + pow(s_mod,2)+pow(s_acc,2)+pow(s_diff,2));
  double d_inel = sqrt(pow(d_lumi,2)+pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2) + pow(d_em,2) + pow(d_mod,2)+pow(d_acc,2)+pow(d_diff,2));

  double s_vis_withoutl  = sqrt(pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2));
  double d_vis_withoutl  = sqrt(pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2));
  double s_had_withoutl  = sqrt(pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2) + pow(s_em,2) + pow(s_mod,2)+pow(s_hadlvl,2));
  double d_had_withoutl  = sqrt(pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2) + pow(d_em,2) + pow(d_mod,2)+pow(d_hadlvl,2));
  double s_inel_withoutl = sqrt(pow(s_pu,2)+pow(s_sel,2)+pow(s_noi,2) + pow(s_em,2) + pow(s_mod,2)+pow(s_acc,2)+pow(s_diff,2));
  double d_inel_withoutl = sqrt(pow(d_pu,2)+pow(d_sel,2)+pow(d_noi,2) + pow(d_em,2) + pow(d_mod,2)+pow(d_acc,2)+pow(d_diff,2));



  ///READ IN VALUES
  TFile f("plots/final_values.root");
  TVectorD* vec_sigma_vis    = NULL;
  TVectorD* vec_sigma_vis_e  = NULL;
  TVectorD* vec_sigma_had    = NULL;
  TVectorD* vec_sigma_had_e  = NULL;
  TVectorD* vec_sigma_inel   = NULL;
  TVectorD* vec_sigma_inel_e = NULL;
  vec_sigma_vis    = (TVectorD*)f.Get("vec_sigma_vis");
  vec_sigma_vis_e  = (TVectorD*)f.Get("vec_sigma_vis_e");
  vec_sigma_had    = (TVectorD*)f.Get("vec_sigma_had");
  vec_sigma_had_e  = (TVectorD*)f.Get("vec_sigma_had_e");
  vec_sigma_inel   = (TVectorD*)f.Get("vec_sigma_inel");
  vec_sigma_inel_e = (TVectorD*)f.Get("vec_sigma_inel_e");
  cout << "Inel - Had - Vis:" << endl;
  vec_sigma_inel->Print();
  vec_sigma_had->Print();
  vec_sigma_vis->Print();
  f.Close();
  const double sigma_had_single = (*vec_sigma_had)[0]*(*cor_fac_had)[0];
  const double sigma_had_double = (*vec_sigma_had)[1]*(*cor_fac_had)[1];

  cout << "Data sigma had single = " << fixed << setprecision(3) << sigma_had_single << endl;
  cout << "Data sigma had double = " << fixed << setprecision(3) << sigma_had_double << endl;
  cout << endl;

  cout << "Uncertainty vis single    = " << fixed << setprecision(1) << s_vis_withoutl << " +lum unc = " << setprecision(3) << s_vis/100.*(*vec_sigma_vis)[0] << "mb" << endl;
  cout << "Uncertainty vis double    = " << fixed << setprecision(1) << d_vis_withoutl << "  +lum unc = " << setprecision(3) << d_vis/100.*(*vec_sigma_vis)[1] << "mb" << endl;
  cout << "Uncertainty had single    = " << fixed << setprecision(1) << s_had_withoutl << "  +lum unc = " << setprecision(3) << s_had/100.*sigma_had_single << "mb" << endl;
  cout << "Uncertainty had double    = " << fixed << setprecision(1) << d_had_withoutl << "  +lum unc = " << setprecision(3) << d_had/100.*sigma_had_double << "mb" << endl;
  cout << "Uncertainty inel single   = " << fixed << setprecision(1) << s_inel_withoutl << "  +lum unc = " << setprecision(3) << s_inel/100.*(*vec_sigma_inel)[1] << "mb" << endl;
  cout << "Uncertainty inel double   = " << fixed << setprecision(1) << d_inel_withoutl << "  +lum unc = " << setprecision(3) << d_inel/100.*(*vec_sigma_inel)[2] << "mb" << endl;
  cout << endl;

  cout << fixed << setprecision(1) << "single (without lumi)=" << s_withoutl << endl;
  cout << fixed << setprecision(1) << "double (without lumi)=" << d_withoutl << endl;
  cout << fixed << setprecision(1) << combined_withoutl << " = " << setprecision(3) << combined_withoutl/100.*(*vec_sigma_inel)[0] << " mb (syst.)" << endl;

  cout << setprecision(3) << s_lumi/100.*(*vec_sigma_inel)[0] << " mb (lumi.)" << endl;
  cout << "Uncertainty inel combined = " << fixed << setprecision(1) << combined << " = " << setprecision(3) << combined/100.*(*vec_sigma_inel)[0] << "mb" << endl;
  cout << endl << endl;


  TVectorD* corr_fac_epos = NULL;
  TVectorD* corr_fac_qgsjet = NULL;
  TVectorD* corr_fac_dpmjet = NULL;
  corr_fac_epos  = (TVectorD*)f0.Get("corr_fac_had_epos"); //changed to had pt eff
  corr_fac_qgsjet  = (TVectorD*)f0.Get("corr_fac_had_qgsjet");
  corr_fac_dpmjet  = (TVectorD*)f0.Get("corr_fac_had_dpmjet");
  if(!corr_fac_epos || !corr_fac_qgsjet || !corr_fac_dpmjet) {cerr << "error" << endl; return;}
  cout << "EPOS Eff. Correction:" << endl;
  corr_fac_epos->Print();
  cout << "QGSJETII-04 Eff. Correction:" << endl;
  corr_fac_qgsjet->Print();
  cout << "DPMJet Eff. Correction:" << endl;
  corr_fac_dpmjet->Print();
  ///!READ IN VALUES

  TH1D* h_data = new TH1D("h_data","CMS;;#sigma [b]",6,-0.5,5.5);
  h_data->GetXaxis()->SetBinLabel(3,"#splitline{hadronic}{inelastic}");
  h_data->GetXaxis()->SetBinLabel(2,"#splitline{#splitline{}{}}{#splitline{hadron level}{p_{T,S}>0.61 GeV}}");
  h_data->GetXaxis()->SetBinLabel(1,"visible");
  h_data->GetXaxis()->SetBinLabel(6,"#splitline{hadronic}{inelastic}");
  h_data->GetXaxis()->SetBinLabel(5,"#splitline{#splitline{}{}}{#splitline{hadron level}{p_{T,D}>0.41 GeV}}");
  h_data->GetXaxis()->SetBinLabel(4,"visible");
  TH1D * h_epos = new TH1D("h_epos","EPOS-LHC;;#sigma [b]",6,-0.5,5.5);
  TH1D * h_qgsjet = new TH1D("h_qgsjet","QGSJETII-04;;#sigma [b]",6,-0.5,5.5);
  TH1D * h_dpmjet = new TH1D("h_dpmjet","DPMJETII-04;;#sigma [b]",6,-0.5,5.5);

#ifdef WITH_ALICE
  TGraphErrors* h_alice_ppb = new TGraphErrors(1);
  TGraphErrors* h_alice_pbp = new TGraphErrors(1);
  TGraphErrors* h_lhcb = new TGraphErrors(1);

  h_alice_ppb->SetPoint(0,3.1,2.09);
  h_alice_ppb->SetPointError(0,0,0.06);

  h_alice_pbp->SetPoint(0,3.2,2.12);
  h_alice_pbp->SetPointError(0,0,0.06);

  h_lhcb->SetPoint(0,0.1,2.09);
  h_lhcb->SetPointError(0,0,0.12);

  SetAttributes<TGraph>(h_lhcb,kCyan+3,22);
  SetAttributes<TGraph>(h_alice_ppb,kCyan+3,29);
  SetAttributes<TGraph>(h_alice_pbp,kCyan+3,30);

#endif

  h_data->SetBinContent(3,(*vec_sigma_inel)[1]);
  h_data->SetBinContent(2,sigma_had_single);
  h_data->SetBinContent(1,(*vec_sigma_vis)[0]);
  h_data->SetBinContent(6,(*vec_sigma_inel)[2]);
  h_data->SetBinContent(5,sigma_had_double);
  h_data->SetBinContent(4,(*vec_sigma_vis)[1]);
  h_data->SetBinError(3,s_inel/100.*(*vec_sigma_inel)[1]);
  h_data->SetBinError(2,s_had/100.*sigma_had_single);
  h_data->SetBinError(1,s_vis/100.*(*vec_sigma_vis)[0]);
  h_data->SetBinError(6,d_inel/100.*(*vec_sigma_inel)[2]);
  h_data->SetBinError(5,d_had/100.*sigma_had_double);
  h_data->SetBinError(4,d_vis/100.*(*vec_sigma_vis)[1]);

  double eff_epos_single = (*corr_fac_epos)[0];
  double eff_epos_double = (*corr_fac_epos)[1];
  double eff_qgsjet_single = (*corr_fac_qgsjet)[0];
  double eff_qgsjet_double = (*corr_fac_qgsjet)[1];
  double eff_dpmjet_single = (*corr_fac_dpmjet)[0];
  double eff_dpmjet_double = (*corr_fac_dpmjet)[1];

  h_epos->SetBinContent(3,EPOS_CS); //this is cross section from model
  h_epos->SetBinContent(2,EPOS_CS*eff_epos_single);
  h_epos->SetBinContent(6,EPOS_CS);
  h_epos->SetBinContent(5,EPOS_CS*eff_epos_double);
  h_qgsjet->SetBinContent(3,QGS_CS);
  h_qgsjet->SetBinContent(2,QGS_CS*eff_qgsjet_single);
  h_qgsjet->SetBinContent(6,QGS_CS);
  h_qgsjet->SetBinContent(5,QGS_CS*eff_qgsjet_double);
  h_dpmjet->SetBinContent(3,DPM_CS);
  h_dpmjet->SetBinContent(2,DPM_CS*eff_dpmjet_single);
  h_dpmjet->SetBinContent(6,DPM_CS);
  h_dpmjet->SetBinContent(5,DPM_CS*eff_dpmjet_double);

  SetAttributes<TH1D>(h_data,kRed,20);
  SetAttributes<TH1D>(h_epos,kGreen-1,22);
  SetAttributes<TH1D>(h_qgsjet,kBlue,34);
  SetAttributes<TH1D>(h_dpmjet,kMagenta,23);

  TCanvas* can1 = new TCanvas;
  h_data->Draw("P");
  h_epos->Draw("SAME P");
  h_qgsjet->Draw("SAME P");
  h_dpmjet->Draw("SAME P");
#ifdef WITH_ALICE
  h_alice_ppb->Draw("P");
  h_alice_pbp->Draw("P");
  h_lhcb->Draw("P");
#endif

  double min = 1.7;
  double max = 2.7;

  h_data->GetYaxis()->SetRangeUser(min,max);

  h_data->GetXaxis()->SetLabelSize(h_data->GetXaxis()->GetLabelSize()*0.9);
  h_data->GetXaxis()->SetLabelOffset(h_data->GetXaxis()->GetLabelOffset()*1.7);
  //h_data->GetXaxis()->LabelsOption("V");
  h_data->GetXaxis()->SetTitleOffset(h_data->GetXaxis()->GetTitleOffset()*1.3);

  TLine* line = new TLine(2.5,min,2.5,max);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("SAME");

  TLegend* leg1 = new TLegend(0.24,0.7,0.55,0.9);
  SetLegAtt(leg1);
  leg1->AddEntry(h_data,"CMS","p");
  leg1->AddEntry(h_epos,"EPOS-LHC","p");
  leg1->AddEntry(h_qgsjet,"QGSJetII-04","p");
  leg1->AddEntry(h_dpmjet,"DPMJet 3.06","p");
#ifdef WITH_ALICE
  leg1->AddEntry(h_alice_ppb,"ALICE V0 (pPb)","p");
  leg1->AddEntry(h_alice_pbp,"ALICE V0 (Pbp","p");
  leg1->AddEntry(h_lhcb,"LHCb","p");
#endif
  leg1->SetBorderSize(1);
  leg1->Draw();


  TPaveText* single_arm = new TPaveText(0.34,0.6,0.44,0.65,"NDC b t l");
  single_arm->SetBorderSize(0);
  single_arm->SetFillStyle(0);
  single_arm->AddText("single-arm");
  single_arm->SetTextSize(0.04);
  single_arm->Draw();

  TPaveText* double_arm = new TPaveText(0.72,0.6,0.82,0.65,"NDC b t l");
  double_arm->SetBorderSize(0);
  double_arm->SetFillStyle(0);
  double_arm->AddText("double-arm");
  double_arm->SetTextSize(0.04);
  double_arm->Draw();

  CMSText(1,0,1);

  can1->SaveAs((string("plots/concl_3")+string(".pdf")).c_str());

  cout << endl;
  cout << "epos  inel=" << EPOS_CS << " hadsingle=" << EPOS_CS*eff_epos_single  << " haddouble=" << EPOS_CS*eff_epos_double << endl;;
  cout << "qgsjet  inel=" << QGS_CS << " hadsingle=" << QGS_CS*eff_qgsjet_single  << " haddouble=" << QGS_CS*eff_qgsjet_double << endl;
  cout << "dpmjet  inel=" << DPM_CS << " hadsingle=" << DPM_CS*eff_dpmjet_single  << " haddouble=" << DPM_CS*eff_dpmjet_double << endl;

  cout << endl << "Differences for cross sections:" << endl;
  cout << endl << "Single:" << endl
       << "HAD -> INEL:" << fixed << setprecision(1) << (*vec_sigma_inel)[1]*1000 - sigma_had_single*1000 << " mb" << endl
       << "VIS -> HAD: " << fixed << setprecision(1) <<  sigma_had_single*1000 - (*vec_sigma_vis)[0]*1000 << " mb" << endl;
  cout << endl << "Double:" << endl
       << "HAD -> INEL:" << fixed << setprecision(1) <<  (*vec_sigma_inel)[2]*1000 - sigma_had_double*1000 << " mb" << endl
       << "VIS -> HAD: " << fixed << setprecision(1) <<  sigma_had_double*1000 - (*vec_sigma_vis)[1]*1000 << " mb" << endl;
  cout << endl << "Single/Double" << endl
       << "VIS: " << (*vec_sigma_vis)[0]*1000 - (*vec_sigma_vis)[1]*1000 << endl
       << "HAD: " << sigma_had_single*1000 - sigma_had_double*1000 << endl;

}



