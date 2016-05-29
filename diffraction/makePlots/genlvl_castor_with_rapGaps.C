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

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);
void ShowStack(TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,string,string);
TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal);

double epos_lumi = 1.;
double qgs_lumi = 1.;
double hijing_lumi = 1.;
double dpm_lumi = 1.;
double sl1_lumi = 1.;
double sl2_lumi = 1.;

void genlvl_castor_with_rapGaps()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> histnames;
  //histnames.push_back("scaling"); //has to go first. only used for scaling
  histnames.push_back("dy0");
  histnames.push_back("dy1");
  histnames.push_back("dy2");
  histnames.push_back("dy3");
  vector<string> sides; sides.push_back("plus"); sides.push_back("minus");
  for(int j = 0; j < int(sides.size()); j++)
    {
      string side = sides[j];
      for(int i = 0; i < int(histnames.size()); i++)
        {
          string histname = histnames[i];
          TFile* file = TFile::Open("histos_new.root");
          TH1D* hijing=(TH1D*)file->Get((string("Hijing/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());
          TH1D* dpm=(TH1D*)file->Get((string("DPMJet/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());
          TH1D* epos=(TH1D*)file->Get((string("Epos/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());
          TH1D* qgs=(TH1D*)file->Get((string("QGSJetII/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());
          TH1D* sl1=(TH1D*)file->Get((string("Starlight_DPMJet/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());
          TH1D* sl2=(TH1D*)file->Get((string("Starlight_Pythia/gen_cas_") + side + string("_with_rapgap_") + histname).c_str());


          //ShowStack(0,0,0,epos,qgs,hijing,dpm,sl1,sl2,histname);
          ShowStack(0,0,0,epos,qgs,hijing,dpm,sl1,sl2,side,histname);
        }
    }
}

void ShowStack(TH1D* data,TH1D* noise,TH1D* zb, TH1D* epos, TH1D* qgs, TH1D* hijing,TH1D* dpm,TH1D* sl1,TH1D* sl2, string side, string type)
  {
  if(type == "dy0")
    {
      epos_lumi  = 2085.7/epos->Integral(0,epos->GetNbinsX()+1);
      qgs_lumi   = 2176.4/qgs->Integral(0,epos->GetNbinsX()+1);
      dpm_lumi   = 2165.8/dpm->Integral(0,epos->GetNbinsX()+1);
      hijing_lumi= 2124.5/hijing->Integral(0,epos->GetNbinsX()+1);
      sl1_lumi   = 195.0/sl1->Integral(0,epos->GetNbinsX()+1);
      sl2_lumi   = 122.0/sl2->Integral(0,epos->GetNbinsX()+1);
      //return;
    }

  epos  ->Scale(epos_lumi);
  qgs   ->Scale(qgs_lumi);
  dpm   ->Scale(dpm_lumi);
  hijing->Scale(hijing_lumi);
  sl1   ->Scale(sl1_lumi);
  sl2   ->Scale(sl2_lumi);

  sl1->SetBit(TH1::kIsAverage);
  sl2->SetBit(TH1::kIsAverage);
  sl1->Add(sl2);
  TH1D* sl = sl1;

  // hijing->SetLineWidth(2);
  // dpm->SetLineWidth(2);
  // epos->SetLineWidth(2);
  // qgs->SetLineWidth(2);
  // sl->SetLineWidth(2);

  hijing->SetMarkerColor(kRed);
  dpm->SetMarkerColor(kMagenta);
  epos->SetMarkerColor(kBlue);
  qgs->SetMarkerColor(kGreen+2);
  sl->SetMarkerColor(kOrange-9);

  hijing->SetMarkerStyle(4);
  dpm->SetMarkerStyle(5);
  epos->SetMarkerStyle(25);
  qgs->SetMarkerStyle(28);
  sl->SetMarkerStyle(22);

  hijing->SetLineColor(hijing->GetMarkerColor());
  dpm->SetLineColor(dpm->GetMarkerColor());
  epos->SetLineColor(epos->GetMarkerColor());
  qgs->SetLineColor(qgs->GetMarkerColor());
  sl->SetLineColor(sl->GetMarkerColor());

  hijing->SetFillColor(hijing->GetMarkerColor());
  dpm->SetFillColor(dpm->GetMarkerColor());
  epos->SetFillColor(epos->GetMarkerColor());
  qgs->SetFillColor(qgs->GetMarkerColor());
  sl->SetFillColor(sl->GetMarkerColor());

  hijing->SetTitle("HIJING 1.383");
  dpm->SetTitle("DPMJet 3.06");
  epos->SetTitle("EPOS-LHC");
  qgs->SetTitle("QGSJetII-04");
  sl->SetTitle("#gamma-p (STARLIGHT+DPMJET/PYTHIA)");

  hijing->GetXaxis()->SetLimits(1,200);
  hijing->GetYaxis()->SetRangeUser(1e-3,1e5);
  hijing->GetXaxis()->SetTitle("E_{CAS,tot} [GeV]");
  hijing->GetYaxis()->SetTitle("dN/L [mb]");
  hijing->GetXaxis()->SetTitleOffset(hijing->GetXaxis()->GetTitleOffset()*1.1);

  TCanvas* c1 = new TCanvas;
  hijing->Draw("P");
  sl->Draw("SAME HIST F");
  hijing->Draw("SAME P");
  dpm->Draw("SAME");
  epos->Draw("SAME");
  qgs->Draw("SAME");

  TLegend* leg = new TLegend(0.23,0.72,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(hijing,"","p");
  leg->AddEntry(dpm,"","p");
  leg->AddEntry(epos,"","p");
  leg->AddEntry(qgs,"","p");
  leg->AddEntry(sl,"","f");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  string text;
  if (type == "dy0") text = "#Deltay_{F,Pb}>0";
  if (type == "dy5") text = "#Deltay_{F,Pb}>5";
  if (type == "dy10") text = "#Deltay_{F,Pb}>10";
  CMSText(3,0,1,text);

  // TLine* line = new TLine(type=="single"?8:4,1e-6,type=="single"?8:4,0.1);
  // line->SetLineWidth(2);
  // line->SetLineStyle(2);
  // line->Draw("SAME");

  c1->SaveAs((string("plots/genlvl_cas_") + side + string("_with_rapgap_") + type + string(".pdf")).c_str());

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
