///Plot the stacked hit distribution of (noise,em + data/mc) for HF

#include "THStack.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TChain.h"
#include "TFile.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVectorD.h"

#include <string>
#include <iostream>
#include <sstream>

//#include "style.h"

using namespace std;

TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal);


void cas_noise()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

    TFile* file = TFile::Open("histos_noise.root");
    TH1D* data=(TH1D*)file->Get(string("data210885/cas_with_rapgap_dy3").c_str());
    TH1D* noise=(TH1D*)file->Get(string("data210885/cas_with_rapgap_noise").c_str());
    TH1D* bxexcl=(TH1D*)file->Get(string("data210885/cas_with_rapgap_noise_bxexcl").c_str());
    TH1D* unpaired=(TH1D*)file->Get(string("data210885/cas_with_rapgap_noise_unpaired").c_str());
    TH1D* zb=(TH1D*)file->Get(string("data210885/cas_with_rapgap_zb").c_str());
    TH1D* epos=(TH1D*)file->Get(string("Epos/cas_with_rapgap_dy3").c_str());

    cout << data << " " << noise << " " << bxexcl << " " << unpaired << " " << zb << " " << epos << endl;
  //RESCALING-------------------------------
  // TH1D** toBeRescaled = &b;
  // double fac = 1./1.3;
  // double* newbinning = new double[(*toBeRescaled)->GetNbinsX()+1];
  // for (int i=0; i<=(*toBeRescaled)->GetNbinsX();i++)
  //   newbinning[i]=(*toBeRescaled)->GetBinLowEdge(i+1)*fac;
  // newbinning[0]=(*toBeRescaled)->GetBinLowEdge(1);
  // newbinning[(*toBeRescaled)->GetNbinsX()]=(*toBeRescaled)->GetBinLowEdge((*toBeRescaled)->GetNbinsX()+1);
  // TH1D* bs = new TH1D("rescaled","rescaled",(*toBeRescaled)->GetNbinsX(),newbinning);
  // for (int i=1; i<=bs->GetNbinsX();i++)
  //   {
  //     bs->SetBinContent(i,(*toBeRescaled)->GetBinContent(i));
  //     bs->SetBinError(i,(*toBeRescaled)->GetBinError(i));
  //   }
  //(*toBeRescaled)=bs;
  //RESCALING-------------------------------

//   const int normbin = data->FindBin(300);
//   const int normendbin = data->GetNbinsX();
//   double eposscale=double(epos->Integral())/double(sl1->Integral());
//   double eposscale2=double(epos->Integral())/double(sl2->Integral());
//   if(noise) noise->Scale(1./double(zb->Integral()) * 48.216 / (9.9845*(1.-(84.+296.)/3568.))); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet
//   const double zbRatioLeft = double(data->Integral()) / double(zb->Integral());
// cout << zbRatioLeft << endl;
//  if(noise) noise->Scale(zbRatioLeft);
//   data->Scale(1./double(data->Integral()));

//   if(type == "dy0")
//     {
//       eposnorm=data->Integral(normbin,normendbin)/epos->Integral(normbin,normendbin);
//       qgsnorm=data->Integral(normbin,normendbin)/qgs->Integral(normbin,normendbin);
//       dpmnorm=data->Integral(normbin,normendbin)/dpm->Integral(normbin,normendbin);
//       hijingnorm=data->Integral(normbin,normendbin)/hijing->Integral(normbin,normendbin);
//     }

//   hijing->Scale(hijingnorm);
//   dpm->Scale(dpmnorm);
//   epos->Scale(eposnorm);
//   qgs->Scale(qgsnorm);

//   sl1->Scale(eposscale*195./2085.*eposnorm);
//   sl2->Scale(eposscale2*122./2085.*eposnorm);


    data  ->Scale(1./data->Integral(1,8));
    zb  ->Scale(1./zb->Integral(1,8));
    noise->Scale(1./noise->Integral(1,8));
    bxexcl->Scale(1./bxexcl->Integral(1,8));
    unpaired->Scale(1./unpaired->Integral(1,8));
//double(data->Integral()) * 48.216 / (9.9845*(1.-(84.+296.)/3568.))); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet

  //zb->Add(noise,-1); //subtract noise from data. don't do it
  zb->SetMarkerSize(1);

  // data->SetLineWidth(2);
  // noise->SetLineWidth(1.5);
  // hijing->SetLineWidth(2);
  // dpm->SetLineWidth(2);
  // epos->SetLineWidth(2);
  // qgs->SetLineWidth(2);
  // sl->SetLineWidth(2);


  zb->SetMarkerColor(kBlack);
  if(noise) noise->SetMarkerColor(kGreen+2);
  if(bxexcl) bxexcl->SetMarkerColor(kRed);
  if(unpaired) unpaired->SetMarkerColor(kBlue);

  if(noise) noise->SetMarkerStyle(34);
  if(bxexcl) bxexcl->SetMarkerStyle(26);
  if(unpaired) unpaired->SetMarkerStyle(28);

  zb->SetLineColor(zb->GetMarkerColor());
  if(noise) noise->SetLineColor(noise->GetMarkerColor());
  if(bxexcl) bxexcl->SetLineColor(bxexcl->GetMarkerColor());
  if(unpaired) unpaired->SetLineColor(unpaired->GetMarkerColor());

  zb->SetFillColor(zb->GetMarkerColor());
  if(noise) noise->SetFillColor(noise->GetMarkerColor());
  if(bxexcl) bxexcl->SetFillColor(bxexcl->GetMarkerColor());
  if(unpaired) unpaired->SetFillColor(unpaired->GetMarkerColor());


  zb->SetTitle("unbiased collisions");
  if(noise) noise->SetTitle("Random (no beam)");
  if(bxexcl) bxexcl->SetTitle("Random - distant from filled bunches");
  if(unpaired) unpaired->SetTitle("Unpaired bunches (one beam)");

  zb->GetXaxis()->SetLimits(1,200);
  zb->GetYaxis()->SetRangeUser(1e-5,5e1);
  zb->GetXaxis()->SetTitle("E_{CAS,tot} [GeV]");
  zb->GetYaxis()->SetTitle("dN/dE/L [mb]");
  zb->GetXaxis()->SetTitleOffset(zb->GetXaxis()->GetTitleOffset()*1.1);


  TCanvas* c1 = new TCanvas;
  zb->Draw("P");
  noise->Draw("SAME");
  bxexcl->Draw("SAME");
  unpaired->Draw("SAME");
  zb->Draw("SAME AXIS");

  TLegend* leg = new TLegend(0.23,0.72,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(zb,"","p");
  if(unpaired) leg->AddEntry(unpaired,"","p");
  if(noise) leg->AddEntry(noise,"","p");
  if(bxexcl) leg->AddEntry(bxexcl,"","p");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  string text;
  CMSText(3,0,1,text);

  // TLine* line = new TLine(type=="single"?8:4,1e-6,type=="single"?8:4,0.1);
  // line->SetLineWidth(2);
  // line->SetLineStyle(2);
  // line->Draw("SAME");

  c1->SaveAs((string("plots/cas_noise") + string(".pdf")).c_str());

}

TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal)
{
  TList *list = new TList;
  list->Add(bg1);
  list->Add(bg2);
  list->Add(signal);
  ostringstream newname; newname << signal->GetName() << "_merged";
  TH1D* newsignal = (TH1D*)signal->Clone(newname.str().c_str());
  newsignal->Reset();
  newsignal->Merge(list);
  return newsignal;
}
