#define _MAXEVT -50000
#define _SkipHFRings 1
#define _HFEnergyScale 1.0 //0.8

#include "TChain.h"
#include "TFile.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVectorD.h"

#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <utility>

#include <CastorTreeVariables.h>
#include <ParticleInfo.h>

//#include "style.h"

using namespace std;
int main()
{
  TH1::SetDefaultSumw2();

  vector<string> sample_fname;
  vector<string> sample_name;
  vector<e_type> sample_type;

  //*************************************************************INPUT***********************************************************
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/VDM210986/*.root"); sample_name.push_back("pPb"); sample_type.push_back(DATA);

  //**************************************************************OUTPUT*********************************************************
  TFile* out_file = new TFile("histos.root","RECREATE");

  //****************************************************************LOOP*******************************************************************

  for (int sample=0; sample<int(sample_name.size()); sample++)
    {

      TChain* tree = new TChain("cAnalyzer/ikCastorTree");
      const int nFiles = tree->Add(sample_fname[sample].c_str()); // file name(s)

      if (nFiles == 0) {
        cout << "No tree files have been found \"" << sample_fname[sample] << "\"" << endl;
        return 0;
      }

      if (tree->GetEntries() == 0) {
        cout << "No events found in file \"" << sample_fname[sample] << "\"" << endl;
        return 0;
      }

      AnalysisEvent* event = 0;
      tree->SetBranchAddress("event", &event);
      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];

      
      TH2D* h_length_scale_x;
      h_length_scale_x = new TH2D((add + string("_h_length_scale_x")).c_str(),"",100,195000000,210000000,30,0.06,0.07);


      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);

          const double orbitNb = event->orbitNb;
          const double nVertex = event->nVertex;
          const double vertexX = event->vertexX;
          const double vertexY = event->vertexY;

          if(195000000 < orbitNb && orbitNb<201512861 && nVertex==1 && !(0.06688 < vertexX && vertexX<0.0669))
            h_length_scale_x->Fill(vertexX,orbitNb);

          vector<double> sectionXBegin;
          vector<double> sectionXEnd;
          int nSec = int(sectionXBegin.size());


          vector<double> sectionYBegin;
          vector<double> sectionYEnd;

          vector<TF1*> sectionXFunc; sectionXFunc.resize(sectionXBegin.size());
          vector<TFitResultPtr> sectionXFuncPtr; sectionXFuncPtr.resize(sectionXBegin.size());
          vector<TF1*> sectionYFunc; sectionYFunc.resize(sectionYBegin.size());
          vector<TFitResultPtr> sectionYFuncPtr; sectionYFuncPtr.resize(sectionYBegin.size());


          for (int i=0; i<nSec; i++)
            {
              ostringstream funcTitle;
              funcTitle << "sectionXFunc_" << i;
              sectionXFunc = new TF1(funcTitle.str().c_str(),"pol0",sectionXBegin[i],sectionXEnd[i]);
              sectionXFunc->SetLineWidth(2);
              sectionXFunc->SetLineColor(kRed);
              sectionXFunc->SetLineStyle(9);
              
              TFitResultPtr sectionXFuncPtr = h_length_scale_x->Fit(sectionXFunc[i],"QWSL",sectionXBegin[i],sectionXEnd[i]);
              sectionXFunc[i]->Draw("SAME");

              ostringstream text;
              double size = 0.02;
              double y = sectionXFuncPtr->Parameter(0)
              text << "#DeltaX_{Fit}=" << y;
              TPaveText* txt = new TPaveText(sectionXBegin[i],y+size,sectionXEnd[i],y+2*size,"b t l");
              txt->SetTextFont(42);
              txt->SetFillStyle(0);
              txt->SetTextColor(kRed);
              txt->SetTextSize(0.033);
              txt->SetBorderSize(0);
              txt->AddText(text.c_str());
              txt->Draw("SAME");
            }
        }

      //******************************************AFTER EVENT LOOP*******************************************
      cout << endl << "--- Processed " << n_total << " events." << endl << endl;
    }

  out_file->Write();
  out_file->Save();

  return 0;
}
