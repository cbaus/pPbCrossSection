//make plots that show 2D xix and xiy distributions

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
#include "TROOT.h"
#include "TSystem.h"
#include "TVectorD.h"

#ifndef __CINT__
#include "style.h"
#endif

#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <iterator>     // std::distance

#include "modelInfoClass.h"


void makePlots_xix_xiy(string filename = "histos_new.root")
{
  gROOT->ProcessLine(" .L style.cc+");
#ifdef __CINT__
  style();
  gStyle->SetPadLeftMargin(0.15); //0.20
  gStyle->SetPadRightMargin(0.15); //0.05
  gROOT->ForceStyle();
#endif

  TFile* file = TFile::Open(filename.c_str());
  modelInfoClass modelInfo;
  for (int iModel=0; iModel<modelInfo.GetN(); ++iModel)
    {
      TH1D* hist = 0;
      hist = (TH1D*)file->Get(((modelInfo.models[iModel] + "/" + modelInfo.models[iModel] + "_h_mc_xix_xiy")).c_str());
      if (!hist)
        {
          cerr << "Model " << modelInfo.names[iModel] << " failure" << endl;
          continue;
        }
      TCanvas* c = new TCanvas;

      hist->SetTitle(";log_{10}(#xi_{x});log_{10}(#xi_{y})");
      hist->GetYaxis()->SetTitleOffset(hist->GetYaxis()->GetTitleOffset()*0.7);

      hist->Draw("COLZ");

#ifdef __CINT__
      CMSText(3,0,1,modelInfo.names[iModel]);
#endif
      c->SaveAs((string("plots/hadron_xix_xiy_") + modelInfo.models[iModel] + string(".eps")).c_str());
      // break;
    }
    return;
}
