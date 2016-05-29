#include "TCanvas.h"
#include "TLine.h"
#include "TH1D.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#ifndef __CINT__
#include "style.h"
#endif

void noise_bxexcl()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_noise.root");
  TH1D* noise_delta = (TH1D*)file->Get((string("data210885/cas_with_noise_deltabx").c_str()));
  if(noise_delta) noise_delta->SetMarkerColor(kBlack);
  if(noise_delta) noise_delta->SetMarkerStyle(34);
  if(noise_delta) noise_delta->SetLineColor(noise_delta->GetMarkerColor());
  if(noise_delta) noise_delta->SetFillColor(noise_delta->GetMarkerColor());
  if(noise_delta) noise_delta->SetTitle("Random trigger");

  noise_delta->GetXaxis()->SetRangeUser(-20,20);
  noise_delta->GetYaxis()->SetRangeUser(0,30);
  noise_delta->GetYaxis()->SetTitle("<E_{CAS,tot}> [GeV]");
  noise_delta->GetXaxis()->SetTitle("#DeltaBX [25 ns]");


  TFitResultPtr fit =noise_delta->Fit("expo","S","",-1.2,6);
      //->Parameter(0)


  TCanvas* c1 = new TCanvas;
  noise_delta->Draw("P");


  TLine* line = new TLine(0,0,0,30);
  line->SetLineWidth(2);
  line->SetLineStyle(3);
  line->Draw("SAME");

  ostringstream fit_title;
fit_title << setprecision(2) << "Fit: exp(" << fit->Parameter(1) <<"#times(x+" << fit->Parameter(0) << "))";

  TLegend* leg = new TLegend(0.21,0.78,0.41,0.93);
  SetLegAtt(leg);
  leg->AddEntry(noise_delta,"","p");
  leg->AddEntry(noise_delta->GetFunction("expo"),fit_title.str().c_str(),"l");
  leg->AddEntry(line,"#DeltaBX=0","l");
  leg->Draw();

  CMSText(3,0,1);

  c1->SaveAs((string("plots/cas_noise_deltabx") + string(".pdf")).c_str());
}
