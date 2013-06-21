#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TROOT.h"

#define CSe 0.082
#define CSe2 0.017



void DataText2(const bool left=true, const bool top=true, const std::string txt="pPb, #sqrt{s_{_{NN}}}=5.02 TeV",bool sim=false);

void CMSPreliminary() {
  TLatex* txt = new TLatex();
  txt->SetNDC();
  txt->SetTextFont(62);
  txt->SetTextColor(kGray+1);
  txt->SetTextSize(0.07);
  txt->SetTextAlign(32);
  txt->DrawLatex(1-gStyle->GetPadRightMargin(), 0.95, "CMS PRELIMINARY");
}
void DataText2(const bool left, const bool top, const std::string str, bool sim) {
  const double marg = 0.005, margY = 0.01;
  TPaveText* txt
    = new TPaveText(left ? gStyle->GetPadLeftMargin()+marg : 1-gStyle->GetPadRightMargin()-marg-0.30,
                    top ? 1-gStyle->GetPadTopMargin()-margY-0.08 : gStyle->GetPadBottomMargin() + margY,
                    left ? gStyle->GetPadLeftMargin()+marg+0.60 : 1-gStyle->GetPadRightMargin()-marg,
                    top ? 1-gStyle->GetPadTopMargin()-margY : gStyle->GetPadBottomMargin() + margY + 0.04,
                    "NDC b t l");
    txt->SetTextAlign(10*(left ? 1 : 3) + (top ? 1 : 3));
    txt->SetTextFont(42);
    txt->SetFillStyle(0);
    txt->SetTextColor(kBlack);
    txt->SetTextSize(0.066);
    txt->SetBorderSize(0);
    txt->AddText(sim==1?"CMS Simulation":"CMS Preliminary");
    txt->AddText(str.c_str());
    txt->Draw();
}

void SetAttributes(TGraphErrors* theGraph, int colour, int marker)
{
  theGraph->SetMarkerSize(1.5);
  theGraph->SetLineWidth(1.8);
  theGraph->SetFillColor(kWhite);
  theGraph->SetLineColor(colour);
  theGraph->SetMarkerColor(colour);
  theGraph->SetMarkerStyle(marker);
}

void makePlots_concl1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TCanvas* can1 = new TCanvas;
  can1->Divide(2);

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


  TGraphErrors* pp_epos  = new TGraphErrors(1);
  pp_epos->SetName("pp_epos");
  pp_epos->SetTitle("EPOS");
  pp_epos->SetPoint(0,1.15,0.064783);
  SetAttributes(pp_epos,kBlue,22);
  TGraphErrors* pb_epos  = new TGraphErrors(1);
  pb_epos->SetName("pb_epos");
  pb_epos->SetTitle("EPOS");
  pb_epos->SetPoint(0,1.15,2.096168);
  SetAttributes(pb_epos,kBlue,22);

  TGraphErrors* pp_eposlhc  = new TGraphErrors(1);
  pp_eposlhc->SetName("pp_eposlhc");
  pp_eposlhc->SetTitle("EPOSLHC");
  pp_eposlhc->SetPoint(0,1.1,0.067448);
  SetAttributes(pp_eposlhc,kGreen-2,23);
  TGraphErrors* pb_eposlhc  = new TGraphErrors(1);
  pb_eposlhc->SetName("pb_eposlhc");
  pb_eposlhc->SetTitle("EPOSLHC");
  pb_eposlhc->SetPoint(0,1.1,2.085703);
  SetAttributes(pb_eposlhc,kGreen-2,23);


  TGraphErrors* pp_qgsjetii03  = new TGraphErrors(1);
  pp_qgsjetii03->SetName("pp_qgsjetii03");
  pp_qgsjetii03->SetTitle("QGSJETII-03");
  pp_qgsjetii03->SetPoint(0,1,0.073279);
  SetAttributes(pp_qgsjetii03,kMagenta-2,25);
  TGraphErrors* pb_qgsjetii03  = new TGraphErrors(1);
  pb_qgsjetii03->SetName("pb_qgsjetii03");
  pb_qgsjetii03->SetTitle("QGSJETII-03");
  pb_qgsjetii03->SetPoint(0,1,2.094112);
  SetAttributes(pb_qgsjetii03,kMagenta-2,25);

  TGraphErrors* pp_qgsjetii  = new TGraphErrors(1);
  pp_qgsjetii->SetName("pp_qgsjetii");
  pp_qgsjetii->SetTitle("QGSJETII-04");
  pp_qgsjetii->SetPoint(0,1.05,0.069422);
  SetAttributes(pp_qgsjetii,kMagenta,24);
  TGraphErrors* pb_qgsjetii  = new TGraphErrors(1);
  pb_qgsjetii->SetName("pb_qgsjetii");
  pb_qgsjetii->SetTitle("QGSJETII-04");
  pb_qgsjetii->SetPoint(0,1.05,2.176422);
  SetAttributes(pb_qgsjetii,kMagenta,24);


  TGraphErrors* pp_fit  = new TGraphErrors(1);
  pp_fit->SetName("pp_fit");
  pp_fit->SetTitle("COMPETE and Glauber");
  pp_fit->SetPoint(0,0.95,0.070);
  pp_fit->SetPointError(0,0,0.070*0.04);
  SetAttributes(pp_fit,kBlack,20);
  TGraphErrors* pb_fit  = new TGraphErrors(1);
  pb_fit->SetName("pb_fit");
  pb_fit->SetTitle("COMPETE and Glauber");
  pb_fit->SetPoint(0,0.95,2.130);
  pb_fit->SetPointError(0,0,0.04);
  SetAttributes(pb_fit,kBlack,20);

  TGraphErrors* pb_this  = new TGraphErrors(1);
  pb_this->SetName("pb_this");
  pb_this->SetTitle("THIS WORK");
  pb_this->SetPoint(0,0.9,CS);
  pb_this->SetPointError(0,0,CS*CSe);
  SetAttributes(pb_this,kRed,21);
  TGraphErrors* pb_this_up  = new TGraphErrors(1);
  pb_this_up->SetName("pb_this_up");
  pb_this_up->SetTitle("THIS WORK");
  pb_this_up->SetPoint(0,0.9,CS*(1+CSe2));
  SetAttributes(pb_this_up,kRed,23);
  pb_this_up->SetMarkerSize(1.2);
  TGraphErrors* pb_this_down  = new TGraphErrors(1);
  pb_this_down->SetName("pb_this_down");
  pb_this_down->SetTitle("THIS WORK");
  pb_this_down->SetPoint(0,0.9,CS*(1.-CSe2));
  SetAttributes(pb_this_down,kRed,22);
  pb_this_down->SetMarkerSize(1.2);


  ///RANGE
  pp_epos->GetYaxis()->SetTitle("#sigma_{inel} [b]");
  pp_epos->GetYaxis()->SetTitleSize(pp_epos->GetYaxis()->GetTitleSize()*1.5);
  pp_epos->GetYaxis()->SetTitleOffset(pp_epos->GetYaxis()->GetTitleOffset()/1.5);
  pp_epos->GetXaxis()->SetNdivisions(0);
  pp_epos->GetXaxis()->SetLimits(0.6,1.4);
  pp_epos->GetHistogram()->SetMaximum(0.076);
  pp_epos->GetHistogram()->SetMinimum(0.064);

  pb_epos->GetXaxis()->SetLimits(0.6,1.4);
  pb_epos->GetXaxis()->SetNdivisions(0);
  pb_epos->GetHistogram()->SetMaximum(2.430);
  pb_epos->GetHistogram()->SetMinimum(1.830);
  ///RANGE

  can1->cd(1);
  pp_epos->Draw("AP");
  pp_eposlhc->Draw("P");
  pp_qgsjetii03->Draw("P");
  pp_qgsjetii->Draw("P");
  pp_fit->Draw("P");

  can1->cd(2);
  pb_epos->Draw("AP");
  pb_eposlhc->Draw("P");
  pb_qgsjetii03->Draw("P");
  pb_qgsjetii->Draw("P");
  pb_fit->Draw("P");
  pb_this->Draw("P");
  pb_this_up->Draw("P");
  pb_this_down->Draw("P");

  TLegend* leg1 = new TLegend(0.283,0.157,0.884,0.357);
  SetLegAtt(leg1);
  leg1->AddEntry(pb_this,"CMS","p");
  leg1->AddEntry(pb_fit,"Menon et al. + Glauber","p");
  leg1->AddEntry(pb_epos,"EPOS","p");
  leg1->AddEntry(pb_eposlhc,"EPOS-LHC","p");
  leg1->AddEntry(pb_qgsjetii03,"QGSJetII-03","p");
  leg1->AddEntry(pb_qgsjetii,"QGSJetII-04","p");
  leg1->SetBorderSize(1);

  can1->cd(1);
  DataText2(true,true,"pp, #sqrt{s_{_{NN}}}=5.02 TeV",true);
  can1->cd(2);
  DataText2(true,true);
  leg1->Draw();
  can1->SaveAs((string("plots/concl_1")+string(".pdf")).c_str());

  TCanvas* can2 = new TCanvas;
  pb_epos->GetYaxis()->SetTitle("#sigma_{inel} [b]");
  pb_epos->GetXaxis()->SetLimits(0.8,1.25);
  pb_epos->GetXaxis()->SetNdivisions(0);
  pb_epos->GetHistogram()->SetMaximum(2.40);
  pb_epos->GetHistogram()->SetMinimum(1.95);

  pb_epos->Draw("AP");
  pb_eposlhc->Draw("P");
  pb_qgsjetii03->Draw("P");
  pb_qgsjetii->Draw("P");
  pb_fit->Draw("P");
  pb_this->Draw("P");
  pb_this_up->Draw("P");
  pb_this_down->Draw("P");

  leg1->SetX1NDC(0.5);
  leg1->SetY1NDC(0.72);
  leg1->SetX2NDC(0.915);
  leg1->SetY2NDC(0.943);
  leg1->Draw();
  DataText(0,0);
  can2->SaveAs((string("plots/concl_1_PAS")+string(".pdf")).c_str());


}
