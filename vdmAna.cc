#define _MAXEVT -500000
#define _SkipHFRings 1
#define _HFEnergyScale 1.0 //0.8

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVectorD.h"

#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "corr/CastorCorrFactorpPb2013.h"
#include "CastorTreeVariables.h"
#include "ParticleInfo.h"

/////////////SETTINGS////////
#define STARTORBIT 80e6
#define ENDORBIT 95e6
//#define STARTORBIT 192e6
//#define ENDORBIT 212e6
/////////////////////////////

//#include "style.h"
using namespace std;

void convert(string filenameL, string filenameR, TGraph*& h_vpos, TGraph*& h_hpos);

int main()
{
  TH1::SetDefaultSumw2();

  vector<string> sample_fname;
  vector<string> sample_name;
  TGraph *h_b1_vpos,*h_b1_hpos=0,*h_b2_vpos=0,*h_b2_hpos=0;
  
  //*************************************************************INPUT***********************************************************
  //sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/VDM210986/*.root"); sample_name.push_back("pPb");
  //sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/VDM211561/*.root"); sample_name.push_back("Pbp");
sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/VDM211821/*.root"); sample_name.push_back("pp");

  //**************************************************************OUTPUT*********************************************************
  TFile* out_file = new TFile("histos_vdmpp.root","RECREATE");
  convert("vdmfiles/B1Lpos.csv","vdmfiles/B1Rpos.csv",h_b1_vpos,h_b1_hpos);
  convert("vdmfiles/B2Lpos.csv","vdmfiles/B2Lpos.csv",h_b2_vpos,h_b2_hpos); // WARNUNG twice L position because R is corrupt
  h_b1_vpos->Write();
  h_b2_vpos->Write();
  h_b1_hpos->Write();
  h_b2_hpos->Write();
  TGraph a(1200);
  for (int i=0; i<1200; i++)
    {
      double x = h_b1_hpos->GetX()[i];
      a.SetPoint(i,x,h_b1_hpos->Eval(x)-h_b2_hpos->Eval(x));
    }
  a.Write();

  //****************************************************************LOOP*******************************************************************

  for (int sample=0; sample<int(sample_name.size()); sample++)
    {

      TChain* tree = new TChain("cAnalyzer/ikCastorTree");
      int nFiles = tree->Add(sample_fname[sample].c_str());
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
      tree->SetBranchStatus("orbitNb", 1);
      tree->SetBranchStatus("time", 1);
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
      h_rate_ls = new TH1D("h_rate_ls","",500,STARTORBIT,ENDORBIT);
      TH1D* h_rate_s;
      h_rate_s = new TH1D("h_rate_s","",500,0,100000);
      TH2D* h_length_scale_x;
      h_length_scale_x = new TH2D("h_length_scale_X","",500,STARTORBIT,ENDORBIT,200,-0.02,0.16);
      TH2D* h_length_scale_y;
      h_length_scale_y = new TH2D("h_length_scale_Y","",500,STARTORBIT,ENDORBIT,200,-0.02,0.16);
      TProfile* h_vdm_castor_e;
      h_vdm_castor_e = new TProfile("h_vdm_castor_e","",1000,135e6,192e6);//option s std dev, no option sigma/sqrt(mu)

      double n_total = double(tree->GetEntries());
      if(_MAXEVT<n_total && _MAXEVT>0)
        n_total = double(_MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 100000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << fixed << iEvent << " / " << int(n_total) << endl;
          tree->GetEntry(iEvent);
          
          const long long time = event->time;
          const double timeUnix = (time>>32) + (time-((time>>32)<<32))/1e6;
          const double timeS     = (time>>32)-1.35941088e9 + (time-((time>>32)<<32))/1e6;
          const double orbitNb = event->orbitNb;
          const int nVertex = event->nVertex;
          const double vertexX = double(event->vertexX);
          const double vertexY = double(event->vertexY);
          const double vertexXe = double(event->vertexXe);
          const double vertexYe = double(event->vertexYe);
          const bool vertexIsFake = event->vertexIsFake;
          
          if(STARTORBIT < orbitNb && orbitNb<ENDORBIT && nVertex==1)
            {
              if(vertexIsFake)
                continue;

              h_length_scale_x->Fill(orbitNb,vertexX,1./vertexXe/vertexXe);
              h_length_scale_y->Fill(orbitNb,vertexY,1./vertexYe/vertexYe);
              h_rate_ls->Fill(orbitNb);
              h_rate_s->Fill(timeS);
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
  out_file->Close();

  return 0;
}


void convert(string filenameL, string filenameR, TGraph*& h_vpos, TGraph*& h_hpos)
{
  ostringstream treename;

  double timestampL, hposL, hresL, vposL, vresL;
  double timestampR, hposR, hresR, vposR, vresR;
  const int estimate = 1200;
  h_vpos = new TGraph(estimate);
  cout << "here" << h_vpos << endl;
  treename.str("");  treename << "vpos_" << filenameL;
  h_vpos->SetName(treename.str().c_str());
  h_hpos = new TGraph(estimate);
  treename.str("");  treename << "h_pos_" << filenameL;
  h_hpos->SetName(treename.str().c_str());

  fstream theFileL,theFileR;
  theFileL.open(filenameL.c_str(), std::ios::in);
  theFileR.open(filenameR.c_str(), std::ios::in);
  if (!theFileL.is_open() || !theFileR.is_open())
    throw runtime_error ("cannot open file");
  string line,element;
  int i=0;
  while(getline ( theFileL, line) )
    {
      stringstream ssline(line);
      getline(ssline, element,','); stringstream sselement1(element); sselement1 >> timestampL;
      getline(ssline, element,','); stringstream sselement2(element); sselement2 >> hposL;
      getline(ssline, element,','); stringstream sselement3(element); sselement3 >> hresL;
      getline(ssline, element,','); stringstream sselement4(element); sselement4 >> vposL;
      getline(ssline, element,','); stringstream sselement5(element); sselement5 >> vresL;

      if(!getline ( theFileR, line))
        {
          cout << "read " << i << "line" << endl;
          throw runtime_error ("one file longer than the other");
        }

      stringstream ssline2(line);
      getline(ssline2, element,','); stringstream sselement6(element); sselement6 >> timestampR;
      getline(ssline2, element,','); stringstream sselement7(element); sselement7 >> hposR;
      getline(ssline2, element,','); stringstream sselement8(element); sselement8 >> hresR;
      getline(ssline2, element,','); stringstream sselement9(element); sselement9 >> vposR;
      getline(ssline2, element,','); stringstream sselement0(element); sselement0 >> vresR;
      
      if(timestampL != timestampR)
        {
          cout << timestampL << " " << timestampR << endl;
          throw runtime_error ("files not synchronised");
        }

      h_vpos->SetPoint(i,timestampL*1e-3,(vposL+vposR)/2.);
      h_hpos->SetPoint(i,timestampL*1e-3,(hposL+hposR)/2.);
      i++;
    }
//     TTree *T = new TTree("treename","data from ascii file");
//     Long64_t nlines = T->ReadFile(filename.c_str(),"timestamp/D:hpos:hres:vpos:vres");
//     cout << "!Found  " << nlines << " lines in file " << filename << " of type " << type << endl;
  assert(estimate == i);
  cout << "here" << h_vpos << endl;
  return;
}
