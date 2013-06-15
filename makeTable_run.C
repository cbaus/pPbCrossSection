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

void makeTable_run()
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


  double totallumi = 0;
  cout << "\\hlineRun number & Intgrated luminosity (nb^-1) & Peak Pileup (\\%)\\\\\\hline\\hline" << endl; 
  for (int run=0; run<int(run_num.size()); run++)
    {
      TFile* file = TFile::Open("histos_old.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;

      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_lumi";
      TH1D* h_lumi=(TH1D*)file->Get(filename_ss.str().c_str());


      const double lumiPerLS=h_lumi->GetMaximum() * _LumiCorr;
      const double lumiIntegral=h_lumi->Integral() * _LumiCorr * 1e-9;
      totallumi += lumiIntegral;

      if (lumiPerLS<=0.) {cerr << "lumi neg: " << i << endl; return;}
      
      const double lambda = lumiPerLS*2.13/11246./23.31/296.;

      if(runname == "211256") cout << "\\hline\\multicolumn{3}{|c|}{Beam Reversal} \\\\\\hline\\hline" << endl;
      cout << runname << " & " << setprecision(3) << lumiIntegral << " & " << lambda*100 << "\\\\\\hline" << endl;
    }
  cout << "\\hline" << endl;
  cout << "Total (" << run_num.size() ")" << " & " << setprecision(3) << totallumi << " & - \\\\\\hline" << endl;
}
      
