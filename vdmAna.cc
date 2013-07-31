#define _MAXEVT -500000
#define _SkipHFRings 1
#define _HFEnergyScale 1.0 //0.8

#include "TChain.h"
#include "TFile.h"
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

#include "corr/CastorCorrFactorpPb2013.h"
#include "CastorTreeVariables.h"
#include "ParticleInfo.h"

//#include "style.h"

using namespace std;
int main()
{
  TH1::SetDefaultSumw2();

  vector<string> sample_fname;
  vector<string> sample_name;

  //*************************************************************INPUT***********************************************************
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/VDM210986/*.root"); sample_name.push_back("pPb");
  sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/VDM"); sample_name.push_back("pPb");

  //**************************************************************OUTPUT*********************************************************
  TFile* out_file = new TFile("histos_vdm.root","RECREATE");

  //****************************************************************LOOP*******************************************************************

  for (int sample=0; sample<int(sample_name.size()); sample++)
    {

      TChain* tree = new TChain("cAnalyzer/ikCastorTree");
      int nFiles = tree->Add((sample_fname[sample]+string("4/*.root")).c_str()); // file name(s)
      nFiles += tree->Add((sample_fname[sample]+string("3/*.root")).c_str()); // file name(s)
      nFiles += tree->Add((sample_fname[sample]+string("2/*.root")).c_str()); // file name(s)
      nFiles += tree->Add((sample_fname[sample]+string("/*.root")).c_str()); // file name(s)

      if (nFiles == 0) {
        cout << "No tree files have been found \"" << sample_fname[sample] << "\"" << endl;
        return 0;
      }

      if (tree->GetEntries() == 0) {
        cout << "No events found in file \"" << sample_fname[sample] << "\"" << endl;
        return 0;
      }

      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("event*", 1);
      tree->SetBranchStatus("orbit*", 1);
      tree->SetBranchStatus("CASTOR.Sector", 1);
      tree->SetBranchStatus("CASTOR.Module", 1);
      tree->SetBranchStatus("CASTOR.Energy", 1);
      tree->SetBranchStatus("*ertex*", 1);

      AnalysisEvent* event = 0;
      tree->SetBranchAddress("event", &event);
      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];

      TH1D* h_rate_ls;
      h_rate_ls = new TH1D("h_rate_ls","",250,192000000,212000000);
      TH2D* h_length_scale_x;
      h_length_scale_x = new TH2D("h_length_scale_X","",250,192000000,212000000,200,-0.02,0.16);
      TH2D* h_length_scale_y;
      h_length_scale_y = new TH2D("h_length_scale_Y","",250,192000000,212000000,200,-0.02,0.16);
      TProfile* h_vdm_castor_e;
      h_vdm_castor_e = new TProfile("h_vdm_castor_e","",1000,135e6,192e6);//option s std dev, no option sigma/sqrt(mu)

      double n_total = double(tree->GetEntries());
      if(_MAXEVT<n_total && _MAXEVT>0)
        n_total = double(_MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 50000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);

          const double orbitNb = event->orbitNb;
          const double nVertex = event->nVertex;
          const double vertexX = event->vertexX;
          const double vertexY = event->vertexY;
          const bool vertexIsFake = event->vertexIsFake;
          
          if(192000000 < orbitNb && orbitNb<212000000 && nVertex==1)
            {
              if(vertexIsFake)
                continue;

              h_length_scale_x->Fill(orbitNb,vertexX);
              h_length_scale_y->Fill(orbitNb,vertexY);
              h_rate_ls->Fill(orbitNb);
            }

          

          if(orbitNb)
            {
              double sumCastorE = 0;
              for (vector<RecHitCASTOR>::const_iterator it = event->CASTOR.begin(); it != event->CASTOR.end(); ++it)
                {
                  const int secNb = it->GetSectorId();
                  const int modNb = it->GetModuleId();
                  if(castor::channelQuality[secNb-1][modNb-1] == false)
                    continue;
                  double corrFactor = castor::channelGainQE[secNb-1][modNb-1] * castor::absEscaleFactor;
                  sumCastorE += it->Energy * corrFactor;
                }
              h_vdm_castor_e->Fill(orbitNb,sumCastorE);
            }
        }

      //******************************************AFTER EVENT LOOP*******************************************
      cout << endl << "--- Processed " << n_total << " events." << endl << endl;
    }

  out_file->Write();
  out_file->Save();

  return 0;
}
