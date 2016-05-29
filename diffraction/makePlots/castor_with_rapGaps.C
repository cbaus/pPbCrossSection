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
#include "TPolyLine.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVectorD.h"

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

//#include "style.h"

class helperClass {
private:
  std::vector<double> xPl;
  std::vector<double> yPl;

  TPolyLine* pline;
public:
  helperClass() {}
  ~helperClass() {}// delete pline; }

  void histFillGaus(TH1* h, const int& nEvents = 1000) { h->FillRandom("gaus",nEvents); }

  void createPolyLine(TH1* h, const double& Xerr) {
    const int nBins = h->GetNbinsX();

    xPl.resize(2*nBins+1);
    yPl.resize(2*nBins+1);

    // loop over bins of histogram
    for(int ibin=1; ibin<=nBins; ibin++) {
      // for TPolyLine
      // move from the left and use lower error
      double valLowX = h->GetBinCenter(ibin) * (1-Xerr);
      double valLowY = h->GetBinContent(ibin) - h->GetBinErrorLow(ibin);

      // always check if values not hit boundary of the frame
      valLowX = valLowX < h->GetBinCenter(1)     ? h->GetBinCenter(1) :
                valLowX > h->GetBinCenter(nBins) ? h->GetBinCenter(nBins) : valLowX;

      valLowY = valLowY < h->GetMinimum() ? h->GetMinimum() :
                valLowY > h->GetMaximum() ? h->GetMaximum() : valLowY;

      xPl[ibin-1] = valLowX;
      yPl[ibin-1] = valLowY;

      // switch move from the right and use higher error
      double valUpX = h->GetBinCenter(ibin) * (1+Xerr);
      double valUpY = h->GetBinContent(ibin) + h->GetBinErrorUp(ibin);

      // always check if values not hit boundary of the frame
      valUpX = valUpX < h->GetBinCenter(1)     ? h->GetBinCenter(1) :
               valUpX > h->GetBinCenter(nBins) ? h->GetBinCenter(nBins) : valUpX;

      valUpY = valUpY < h->GetMinimum() ? h->GetMinimum() :
                valUpY > h->GetMaximum() ? h->GetMaximum() : valUpY;

      xPl[2*nBins-ibin] = valUpX;
      yPl[2*nBins-ibin] = valUpY;

      // set the end point
      if( ibin == 1 ) {
        xPl[2*nBins] = valLowX;
        yPl[2*nBins] = valLowY;
      }
    }

    pline = new TPolyLine(2*nBins+1,&xPl.front(),&yPl.front());
  }

  void setPolyLineStyle() {
    pline->SetLineColor(kRed);
    pline->SetLineWidth(2);
    pline->SetFillColor(16);
  }

  void drawPolyLine(TH1* h,const double& Xerr,const TString& drawCmd = "") {
    createPolyLine(h,Xerr);
    setPolyLineStyle();

    pline->Draw(drawCmd);
  }

};

using namespace std;

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);
void ShowStack(TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,string);
TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal);

double epos_lumi = 1.;
double qgs_lumi = 1.;
double hijing_lumi = 1.;
double dpm_lumi = 1.;
double sl1_lumi = 1.;
double sl2_lumi = 1.;

void castor_with_rapGaps()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> histnames;
  histnames.push_back("zb"); //has to go first. only used for scaling
  histnames.push_back("dy0");
  histnames.push_back("dy1");
  histnames.push_back("dy2");
  histnames.push_back("dy3");
  for(int i = 0; i < int(histnames.size()); i++)
  {
    TFile* file = TFile::Open("histos.root");
    TH1D* data=(TH1D*)file->Get((string("data210885/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* noise=(TH1D*)file->Get((string("data210885/cas_with_rapgap_noise_bxexcl").c_str()));
    TH1D* zb=(TH1D*)file->Get((string("data210885/cas_with_rapgap_zb").c_str()));
    TH1D* hijing=(TH1D*)file->Get((string("Hijing/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* dpm=(TH1D*)file->Get((string("DPMJet/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* epos=(TH1D*)file->Get((string("Epos/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* qgs=(TH1D*)file->Get((string("QGSJetII/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* sl1=(TH1D*)file->Get((string("Starlight_DPMJet/cas_with_rapgap_") + histnames[i]).c_str());
    TH1D* sl2=(TH1D*)file->Get((string("Starlight_Pythia/cas_with_rapgap_") + histnames[i]).c_str());

    //noise = 0;

    ShowStack(data,noise,zb,epos,qgs,hijing,dpm,sl1,sl2,histnames[i]);
  }
}

void ShowStack(TH1D* data,TH1D* noise,TH1D* zb, TH1D* epos, TH1D* qgs, TH1D* hijing,TH1D* dpm,TH1D* sl1,TH1D* sl2, string type)
  {
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

  if(type == "zb")
    {
      epos_lumi  = 2085.7/epos->Integral(0,epos->GetNbinsX()+1);
      qgs_lumi   = 2176.4/qgs->Integral(0,epos->GetNbinsX()+1);
      dpm_lumi   = 2165.8/dpm->Integral(0,epos->GetNbinsX()+1);
      hijing_lumi= 2124.5/hijing->Integral(0,epos->GetNbinsX()+1);
      sl1_lumi   = 195.0/sl1->Integral(0,epos->GetNbinsX()+1);
      sl2_lumi   = 122.0/sl2->Integral(0,epos->GetNbinsX()+1);
      return;
    }

  data  ->Scale(1./1310000.);
  epos  ->Scale(epos_lumi);
  qgs   ->Scale(qgs_lumi);
  dpm   ->Scale(dpm_lumi);
  hijing->Scale(hijing_lumi);
  sl1   ->Scale(sl1_lumi);
  sl2   ->Scale(sl2_lumi);
  if(noise) noise->Scale(1./noise->Integral(1,8)*data->Integral(1,7));
  //noise->Scale(1./double(data->Integral()) * 48.216 / (9.9845*(1.-(84.+296.)/3568.))); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet

  sl1->SetBit(TH1::kIsAverage);
  sl2->SetBit(TH1::kIsAverage);
  sl1->Add(sl2);
  TH1D* sl = sl1;
  //sl->Reset();

  //data->Add(noise,-1); //subtract noise from data. don't do it
  data->SetMarkerSize(1);

  // data->SetLineWidth(2);
  // noise->SetLineWidth(1.5);
  // hijing->SetLineWidth(2);
  // dpm->SetLineWidth(2);
  // epos->SetLineWidth(2);
  // qgs->SetLineWidth(2);
  // sl->SetLineWidth(2);


  data->SetMarkerColor(kBlack);
  if(noise) noise->SetMarkerColor(kGreen-10);
  hijing->SetMarkerColor(kRed);
  dpm->SetMarkerColor(kMagenta);
  epos->SetMarkerColor(kBlue);
  qgs->SetMarkerColor(kGreen+2);
  sl->SetMarkerColor(kOrange-9);

  if(noise) noise->SetMarkerStyle(34);
  hijing->SetMarkerStyle(4);
  dpm->SetMarkerStyle(5);
  epos->SetMarkerStyle(25);
  qgs->SetMarkerStyle(28);
  sl->SetMarkerStyle(22);

  data->SetLineColor(data->GetMarkerColor());
  if(noise) noise->SetLineColor(noise->GetMarkerColor());
  hijing->SetLineColor(hijing->GetMarkerColor());
  dpm->SetLineColor(dpm->GetMarkerColor());
  epos->SetLineColor(epos->GetMarkerColor());
  qgs->SetLineColor(qgs->GetMarkerColor());
  sl->SetLineColor(sl->GetMarkerColor());

  data->SetFillColor(data->GetMarkerColor());
  if(noise) noise->SetFillColor(noise->GetMarkerColor());
  hijing->SetFillColor(hijing->GetMarkerColor());
  dpm->SetFillColor(dpm->GetMarkerColor());
  epos->SetFillColor(epos->GetMarkerColor());
  qgs->SetFillColor(qgs->GetMarkerColor());
  sl->SetFillColor(sl->GetMarkerColor());


  data->SetTitle("Data");
  if(noise) noise->SetTitle("Noise");
  hijing->SetTitle("HIJING 1.383");
  dpm->SetTitle("DPMJet 3.06");
  epos->SetTitle("EPOS-LHC");
  qgs->SetTitle("QGSJetII-04");
  sl->SetTitle("#gamma-p (STARLIGHT+DPMJET/PYTHIA)");

  data->GetXaxis()->SetLimits(1,200);
  data->GetYaxis()->SetRangeUser(1e-3,1e5);
  if (type == "dy0")
    data->GetYaxis()->SetRangeUser(1e-3,1e6);
  if (type == "dy3")
    data->GetYaxis()->SetRangeUser(1e-3,1e8);
  data->GetXaxis()->SetTitle("E_{CAS,tot} [GeV]");
  data->GetYaxis()->SetTitle("dN/dE/L [mb]");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  THStack* h_s_bg = new THStack("h_s_gb","events");
  if(noise) h_s_bg->Add(noise,"HIST F");
  h_s_bg->Add(sl,"HIST F");

  if(1) //switch to add background to MC
    {
      hijing = merge(noise,sl,hijing);
      dpm = merge(noise,sl,dpm);
      epos   = merge(noise,sl,epos);
      qgs    = merge(noise,sl,qgs);
    }

  TCanvas* c1 = new TCanvas;
  data->Draw("P");
  helperClass hlpc;
  h_s_bg->Draw("SAME");
  hlpc.drawPolyLine(data,0.22,"f same");
  hijing->Draw("SAME");
  dpm->Draw("SAME");
  epos->Draw("SAME");
  qgs->Draw("SAME");
  data->Draw("SAME P");
  data->Draw("SAME AXIS");

  TLegend* leg = new TLegend(0.23,0.68,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(data,"","p");
  leg->AddEntry(hijing,"","p");
  leg->AddEntry(dpm,"","p");
  leg->AddEntry(epos,"","p");
  leg->AddEntry(qgs,"","p");
  leg->AddEntry(sl,"","f");
  if(noise) leg->AddEntry(noise,"","f");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  string text;
  if (type == "dy0") text = "#Deltay_{F,Pb}>0";
  if (type == "dy1") text = "#Deltay_{F,Pb}>2";
  if (type == "dy2") text = "#Deltay_{F,Pb}>8";
  if (type == "dy3") text = "#Deltay_{F,Pb}>10";
  CMSText(3,0,1);

  // TLine* line = new TLine(type=="single"?8:4,1e-6,type=="single"?8:4,0.1);
  // line->SetLineWidth(2);
  // line->SetLineStyle(2);
  // line->Draw("SAME");

  c1->SaveAs((string("plots/cas_with_rapgap_") + type + string(".pdf")).c_str());

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
