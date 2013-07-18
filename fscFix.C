#include <map>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
using namespace std;
void fscFix()
{
  TFile* file1 = TFile::Open("/afs/cern.ch/work/j/jgomez2/public/ForColin/FSCTrees_PAMinBiasUPCRaw_210614merged.root");
  TFile* file2 = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r210614_autoforest_v51.root");
  TFile* fileout = TFile::Open("fsctree.root","RECREATE");

  TTree* tree0 = 0;
  TTree* tree1 = 0;
  TTree* tree2 = 0;
  TTree* treeout = 0;
  file1->GetObject("FSCDigiTree",tree0);
  file1->GetObject("BeamTree",tree1);

  map<int,int> evmap;
  
  int event1;
  int event2;
  float posFSC1cFtsz[10], posFSC2cFtsz[10], posFSC3cFtsz[10], posFSC4cFtsz[10];
  float negFSC1cFtsz[10], negFSC2cFtsz[10], negFSC3cFtsz[10], negFSC4cFtsz[10];

  tree0->SetBranchStatus("*",0);
  tree0->SetBranchStatus("*FSC*fC*",1);
  tree0->SetBranchAddress("posFSC1fC",posFSC1cFtsz);
  tree0->SetBranchAddress("posFSC2fC",posFSC2cFtsz);
  tree0->SetBranchAddress("posFSC3fC",posFSC3cFtsz);
  tree0->SetBranchAddress("posFSC4fC",posFSC4cFtsz);
  tree0->SetBranchAddress("negFSC1fC",negFSC1cFtsz);
  tree0->SetBranchAddress("negFSC2fC",negFSC2cFtsz);
  tree0->SetBranchAddress("negFSC3fC",negFSC3cFtsz);
  tree0->SetBranchAddress("negFSC4fC",negFSC4cFtsz);

  tree1->SetBranchStatus("*",0);
  tree1->SetBranchStatus("Event",1);
  tree1->SetBranchAddress("Event",&event1);

  const int n1=tree1->GetEntries();
  for(int i=0; i<n1; i++)
    {
      if (i%1000000 == 0) cout << "Entry: " << i << endl;
      tree1->GetEntry(i);
      evmap.insert(pair<int,int>(event1,i));
    }
  cout << "Done" << endl;
  
  ((TDirectory*)file2->Get("hiEvtAnalyzer"))->GetObject("HiTree",tree2);
  
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("evt",1);
  tree2->SetBranchAddress("evt",&event2);

  float fsc_sum_plus=0;
  float fsc_sum_minus=0;
  fileout->cd();
  treeout = new TTree;
  treeout->Branch("fsc_sum_plus",&fsc_sum_plus,"fsc_sum_plus/F");
  treeout->Branch("fsc_sum_minus",&fsc_sum_minus,"fsc_sum_minus/F");


  const int n2=tree2->GetEntries();
  for(int i=0; i<n2; i++)
    {
      if (i%100000 == 0) cout << "Entry: " << i << endl;
      tree2->GetEntry(i);
      map<int,int>::iterator it;
      it=evmap.find(event2);
      if (it == evmap.end() )
        {
          //cerr << "cannot find element" << i << endl;
          fsc_sum_minus = -1000;
          fsc_sum_plus = -1000;
        }
      else
        {
          tree0->GetEntry(it->second);
          tree1->GetEntry(it->second);

          fsc_sum_minus = posFSC1cFtsz[5] + posFSC2cFtsz[5] + posFSC3cFtsz[5] + posFSC4cFtsz[5];
          fsc_sum_minus += posFSC1cFtsz[6] + posFSC2cFtsz[6] + posFSC3cFtsz[6] + posFSC4cFtsz[6];
          fsc_sum_plus = negFSC1cFtsz[5] + negFSC2cFtsz[5] + negFSC3cFtsz[5] + negFSC4cFtsz[5];
          fsc_sum_plus += negFSC1cFtsz[6] + negFSC2cFtsz[6] + negFSC3cFtsz[6] + negFSC4cFtsz[6];
        }
      //cout << i << " : " << event1 << " " << event2 << " " << fsc_sum_minus << " " << fsc_sum_plus << " " << it->second << endl;
      treeout->Fill();
    }
  treeout->Write("fsctree");
  fileout->Close();
}
