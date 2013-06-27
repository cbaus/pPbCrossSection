#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TVectorD.h>

#include <iostream>

using namespace std;

#define CSe 0.082

void SetAttributes(TGraphErrors* theGraph, int colour, int marker)
{
  theGraph->SetMarkerSize(1.8);
  theGraph->SetLineWidth(2.1);
  theGraph->SetFillColor(kWhite);
  theGraph->SetLineColor(colour);
  theGraph->SetMarkerColor(colour);
  theGraph->SetMarkerStyle(marker);
  theGraph->SetLineStyle(marker);
}

double GetSqrtS(double elab)
{
  TLorentzVector s = TLorentzVector(0,0,elab,sqrt(pow(elab,2)+pow(0.938,2))) + TLorentzVector(0,0,0,0.938);
  return s.M();
}


void makePlots_concl2()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  gStyle->SetPadTopMargin(0.07);

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
  ///!READ IN VALUES

  double CS = (*vec_sigma_inel)[0];

  TGraphErrors* g_ihep  = new TGraphErrors(5);
  g_ihep->SetName("g_ihep");
  g_ihep->SetTitle("IHEP;#sqrt{s_{NN}} [GeV];#sigma_{inel} [b]");
  SetAttributes(g_ihep,kBlue,21);
  g_ihep->SetPoint(0,GetSqrtS(20),1.739);  g_ihep->SetPointError(0,0,0.030);
  g_ihep->SetPoint(1,GetSqrtS(30),1.870);  g_ihep->SetPointError(1,0,0.023);
  g_ihep->SetPoint(2,GetSqrtS(40),1.780);  g_ihep->SetPointError(2,0,0.020);
  g_ihep->SetPoint(3,GetSqrtS(50),1.785);  g_ihep->SetPointError(3,0,0.029);
  g_ihep->SetPoint(4,GetSqrtS(60),1.930);  g_ihep->SetPointError(4,0,0.050);

  TGraphErrors* g_cosmic  = new TGraphErrors(3);
  g_cosmic->SetName("g_cosmic");
  g_cosmic->SetTitle("COSMIC;#sqrt{s_{NN}} [GeV];#sigma_{inel} [b]");
  SetAttributes(g_cosmic,kOrange+2,23);
  g_cosmic->SetPoint(0,GetSqrtS(700),1.786);  g_cosmic->SetPointError(0,0,0.096);
  g_cosmic->SetPoint(1,GetSqrtS(1500),1.863);  g_cosmic->SetPointError(1,0,0.105);
  g_cosmic->SetPoint(2,GetSqrtS(3500),1.951);  g_cosmic->SetPointError(2,0,0.125);

  TGraphErrors* g_fnal  = new TGraphErrors(3);
  g_fnal->SetName("g_fnal");
  g_fnal->SetTitle("IHEP;#sqrt{s_{NN}} [GeV];#sigma_{inel} [b]");
  SetAttributes(g_fnal,kGreen-2,22);
  g_fnal->SetPoint(0,GetSqrtS(60),1.730);  g_fnal->SetPointError(0,0,0.052);
  g_fnal->SetPoint(1,GetSqrtS(200),1.765);  g_fnal->SetPointError(1,0,0.053);
  g_fnal->SetPoint(2,GetSqrtS(280),1.752);  g_fnal->SetPointError(2,0,0.053);

  TGraphErrors* g_glauber  = new TGraphErrors(9);
  g_glauber->SetName("g_glauber");
  g_glauber->SetTitle("(COMPETE+TOTEM)+Glauber");
  SetAttributes(g_glauber,kBlack,1);
  g_glauber->SetPoint(0,5,1.87851);
  g_glauber->SetPoint(1,10,1.80454);
  g_glauber->SetPoint(2,30,1.78135);
  g_glauber->SetPoint(3,100,1.82648);
  g_glauber->SetPoint(4,300,1.90241);
  g_glauber->SetPoint(5,1000,1.99123);
  g_glauber->SetPoint(6,3000,2.07849);
  g_glauber->SetPoint(7,5020,2.12448);
  g_glauber->SetPoint(8,8800,2.16933);

  TGraphErrors* g_eposlhc  = new TGraphErrors(7);
  g_eposlhc->SetName("g_eposlhc");
  g_eposlhc->SetTitle("EPOS-LHC");
  SetAttributes(g_eposlhc,kGreen-1,7);
  g_eposlhc->SetPoint(0,GetSqrtS(20),1.646732);
  g_eposlhc->SetPoint(1,GetSqrtS(60),1.663151);
  g_eposlhc->SetPoint(2,100,1.784655);
  g_eposlhc->SetPoint(3,200,1.831986);
  g_eposlhc->SetPoint(4,1000,1.927643);
  g_eposlhc->SetPoint(5,2000,1.991543);
  g_eposlhc->SetPoint(6,5020,2.085703);

  TGraphErrors* g_q4  = new TGraphErrors(7);
  g_q4->SetName("g_q4");
  g_q4->SetTitle("QGSJetII-04");
  SetAttributes(g_q4,kBlue,2);
  g_q4->SetPoint(0,GetSqrtS(20),1.651290);
  g_q4->SetPoint(1,GetSqrtS(60),1.704261);
  g_q4->SetPoint(2,100,1.877070);
  g_q4->SetPoint(3,200,1.925534);
  g_q4->SetPoint(4,1000,2.044743);
  g_q4->SetPoint(5,2000,2.100428);
  g_q4->SetPoint(6,5020,2.176422);

  TGraphErrors* g_cms  = new TGraphErrors(1);
  g_cms->SetName("g_cms");
  g_cms->SetTitle("CMS");
  SetAttributes(g_cms,kRed,20);
  g_cms->SetPoint(0,5020,CS);  g_cms->SetPointError(0,0,CS*CSe);


  TLegend* leg1 = new TLegend(0.25,0.65,0.45,0.90);
  SetLegAtt(leg1);
  leg1->AddEntry(g_cms,"CMS","p");
  leg1->AddEntry(g_ihep,"IHEP","p");
  leg1->AddEntry(g_fnal,"FNAL","p");
  leg1->AddEntry(g_cosmic,"Avakian et al.","p");
  leg1->AddEntry(g_glauber,"(COMPETE+TOTEM)+Glauber","l");
  leg1->AddEntry(g_eposlhc,"EPOS-LHC","l");
  leg1->AddEntry(g_q4,"QGSJetII-04","l");

  TCanvas* can1 = new TCanvas;
  g_ihep->GetXaxis()->SetLimits(5,9000);
  g_ihep->GetYaxis()->SetRangeUser(1.500,2.600);
  g_ihep->GetXaxis()->SetTitleOffset(g_ihep->GetXaxis()->GetTitleOffset()*1.2);

  g_ihep->Draw("AP");
  g_fnal->Draw("P");
  g_cosmic->Draw("P");
  g_glauber->Draw("C");
  g_eposlhc->Draw("C");
  g_q4->Draw("C");
  g_cms->Draw("P");
  can1->SetLogx();
  leg1->Draw();
  CMSText(1,0,1);
  can1->SaveAs((string("plots/concl_2")+string(".pdf")).c_str());

}
