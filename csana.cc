#define MAXEVT -50000

#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"
#include "iostream"

#include "style.h"

using namespace std;

void csana()
{
  gROOT->ProcessLine(" .L style.cc++");
  //**********************************************INPUT******************************
  TFile* data_file = NULL; //****FILE
  data_file = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r210658_autoforest_v51.root");
  if (!data_file) cerr << "Cannot find data file" << endl;

  TTree* data_tree = NULL; //****TREE0
  ((TDirectory*)data_file->Get("rechitanalyzer"))->GetObject("castor",data_tree);
  if (!data_tree) cerr << "Cannot find data tree" << endl;

  TTree* data_friend_tree0 = NULL; //****TREE1
  ((TDirectory*)data_file->Get("pfTowers"))->GetObject("hf",data_friend_tree0);
  if (!data_friend_tree0) cerr << "Cannot find data friend tree" << endl;
  data_tree->AddFriend(data_friend_tree0);

  TTree* data_friend_tree1 = NULL; //****TREE2
  ((TDirectory*)data_file->Get("hltanalysis"))->GetObject("HltTree",data_friend_tree1);
  if (!data_friend_tree1) cerr << "Cannot find data friend tree" << endl;
  data_tree->AddFriend(data_friend_tree1);

  TTree* data_friend_tree2 = NULL; //****TREE3
  ((TDirectory*)data_file->Get("pixelTrack"))->GetObject("trackTree",data_friend_tree2);
  if (!data_friend_tree2) cerr << "Cannot find data friend tree" << endl;
  data_tree->AddFriend(data_friend_tree2);



  TFile* mc_file = NULL; //****FILE
  mc_file = TFile::Open("root://eoscms//eos/cms//store/caf/user/dgulhan/pPb_Hijing_MB/HiForest_v03_mergedv02/merged_forest_0.root");
  if (!mc_file) cerr << "Cannot find mc file" << endl;

  TTree* mc_tree = NULL; //****TREE0
  ((TDirectory*)mc_file->Get("rechitanalyzer"))->GetObject("castor",mc_tree);
  if (!mc_tree) cerr << "Cannot find mc tree" << endl;

  TTree* mc_friend_tree0 = NULL; //****TREE1
  ((TDirectory*)mc_file->Get("pfTowers"))->GetObject("hf",mc_friend_tree0);
  if (!mc_friend_tree0) cerr << "Cannot find mc friend tree" << endl;
  mc_tree->AddFriend(mc_friend_tree0);

  TTree* mc_friend_tree1 = NULL; //****TREE2
  ((TDirectory*)mc_file->Get("hltanalysis"))->GetObject("HltTree",mc_friend_tree1);
  if (!mc_friend_tree1) cerr << "Cannot find mc friend tree" << endl;
  mc_tree->AddFriend(mc_friend_tree1);

  TTree* mc_friend_tree2 = NULL; //****TREE3
  ((TDirectory*)mc_file->Get("pixelTrack"))->GetObject("trackTree",mc_friend_tree2);
  if (!mc_friend_tree2) cerr << "Cannot find mc friend tree" << endl;
  mc_tree->AddFriend(mc_friend_tree2);

  const int n_cas_rechits = 224;
  const int n_hf_rechits = 1800;
  
  int hf_n;
  float hf_e[n_hf_rechits];
  int hf_depth[n_hf_rechits];

  float cas_e[n_cas_rechits];
  int cas_sat[n_cas_rechits];
  int cas_depth[n_cas_rechits];
  int cas_iphi[n_cas_rechits];

  int nTrk;

  int zero_bias;
  int bptx_p_nm;
  int bptx_np_m;

  //*********************************OUTPUT************************+++

  TFile* out_file = new TFile("histos.root","RECREATE");

  TH1D* h_zero_count_zb_coll;
  TH1D* h_zero_count_zb_no_coll;
  TH1D* h_hf_towers_zb_coll;
  TH1D* h_hf_towers_zb_no_coll;
  TH1D* h_eff;



  //**********************************LOOP*******************

  for (int type=0; type<2; type++)
    {
      TTree* tree = NULL;

      switch (type)
        {
        case 0:tree = mc_tree; break;
        case 1:tree = data_tree; break;
        default: cerr << "switch op wrong" << endl; break;
        }

      tree->SetBranchAddress("hf.n",&hf_n);
      tree->SetBranchAddress("hf.e",hf_e);
      tree->SetBranchAddress("hf.depth",hf_depth);

      tree->SetBranchAddress("e",cas_e);
      tree->SetBranchAddress("saturation",cas_sat);
      tree->SetBranchAddress("depth",cas_depth);
      tree->SetBranchAddress("iphi",cas_iphi);


      tree->SetBranchAddress("nTrk",&nTrk);

      tree->SetBranchAddress("HLT_PAZeroBias_v1",&zero_bias);
      tree->SetBranchAddress("HLT_PABptxPlusNotBptxMinus_v1",&bptx_p_nm);
      tree->SetBranchAddress("HLT_PABptxMinusNotBptxPlus_v1",&bptx_np_m);

      switch(type)
        {
        case 0:
          out_file->mkdir("mc");
          out_file->cd("mc");
          h_zero_count_zb_coll               = new TH1D("mc_h_zero_count_zb_coll","",100,0,1000);
          h_zero_count_zb_no_coll            = new TH1D("mc_h_zero_count_zb_no_coll","",100,0,1000);
          h_hf_towers_zb_coll                = new TH1D("mc_h_hf_towers_zb_coll","",100,0,200);
          h_hf_towers_zb_no_coll             = new TH1D("mc_h_hf_towers_zb_no_coll","",100,0,200);
          h_eff                              = new TH1D("mc_h_eff","",9,-0.5,8.5);
          break;
        case 1:
          out_file->mkdir("data");
          out_file->cd("data");
          h_zero_count_zb_coll               = new TH1D("data_h_zero_count_zb_coll","",100,0,1000);
          h_zero_count_zb_no_coll            = new TH1D("data_h_zero_count_zb_no_coll","",100,0,1000);
          h_hf_towers_zb_coll                = new TH1D("data_h_hf_towers_zb_coll","",100,0,200);
          h_hf_towers_zb_no_coll             = new TH1D("data_h_hf_towers_zb_no_coll","",100,0,200);
          h_eff                              = new TH1D("data_h_eff","",9,-0.5,8.5);
          break;
        default: cerr << "switch op wrong" << endl; break;
        }
      
      h_eff->GetXaxis()->SetBinLabel(1,"ZB");
      h_eff->GetXaxis()->SetBinLabel(2,"No coll (BBTX)");
      h_eff->GetXaxis()->SetBinLabel(3,"HF#pm zero count < 500");
      h_eff->GetXaxis()->SetBinLabel(4,"HF energy cut >3GeV");
      h_eff->GetXaxis()->SetBinLabel(5,"HF energy cut >4GeV");
      h_eff->GetXaxis()->SetBinLabel(6,"HF energy cut >5GeV");
      h_eff->GetXaxis()->SetBinLabel(7,"Pixel tracks >0");
      h_eff->GetXaxis()->SetBinLabel(8,"Pixel tracks >3");
      h_eff->GetXaxis()->SetBinLabel(9,"Pixel tracks >7");


      for(int iEvent=0; iEvent<tree->GetEntries(); iEvent++)
        {
          if(iEvent==MAXEVT) break;
          string s_type = type?string("DATA"):string("MC");
          if(iEvent % 10000 == 0) cout << s_type.c_str() << " Entry: " << iEvent << " / " << tree->GetEntries() << endl;
          tree->GetEntry(iEvent);
          
          if(type == 1 && zero_bias == 0)
            continue;
          
          double sum_cas_e_em = 0;
          double sum_cas_e_had = 0;
          
          bool no_coll = bptx_p_nm==1 || bptx_np_m==1;
          
          int hf_zero_count = n_hf_rechits - hf_n;
          
          for(int ch_cas=0; ch_cas<n_cas_rechits; ch_cas++) // no ZS
            {break;
              if(cas_depth[ch_cas] < 3)
                sum_cas_e_em += cas_e[ch_cas];
              else
                sum_cas_e_had += cas_e[ch_cas];
            }
          double hf_energy_max = -1;
          for(int ch_hf=0; ch_hf<hf_n; ch_hf++)
            {
              if(hf_e[ch_hf] > hf_energy_max)
                hf_energy_max = hf_e[ch_hf];
            }
          //Filling HISTOS
          if(1)                    h_zero_count_zb_coll->Fill(hf_zero_count);
          if(no_coll)              h_zero_count_zb_no_coll->Fill(hf_zero_count);
          if(1)                    h_hf_towers_zb_coll->Fill(hf_energy_max);
          if(no_coll)              h_hf_towers_zb_no_coll->Fill(hf_energy_max);
          if(zero_bias)            h_eff->Fill(0);
          if(no_coll == 1)         h_eff->Fill(1);
          if(hf_zero_count < 500)  h_eff->Fill(2);
          if(hf_energy_max > 3)    h_eff->Fill(3);
          if(hf_energy_max > 4)    h_eff->Fill(4);
          if(hf_energy_max > 5)    h_eff->Fill(5);
          if(nTrk > 0)             h_eff->Fill(6);
          if(nTrk > 3)             h_eff->Fill(7);
          if(nTrk > 7)             h_eff->Fill(8);
          
        }
      double lumi = 1854.344875;// nb^-1
      double n_total = double(tree->GetEntries());
      double n_zb = h_eff->GetBinContent(1) * (MAXEVT>0?n_total/double(MAXEVT):1.) * 40.;
      double n_mb = h_eff->GetBinContent(4) * (MAXEVT>0?n_total/double(MAXEVT):1.) * 40.;

      cout << endl << "Cross section ZB: " << n_zb/lumi << " --- cross section MB: " << n_mb/lumi << endl;
    }
  
  //******************After LOOP*******************
  

  //******************END************************

  out_file->Write();
  out_file->Save();
}
