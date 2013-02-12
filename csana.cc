#define MAXEVT -50000

#include "TTree.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"
#include <iostream>
#include <vector>
#include <string>

#include "style.h"

enum e_type
{
  MC = 0,
  DATA
};

vector<string> sample_fname;
vector<string> sample_name;
vector<e_type> sample_type;
using namespace std;

void csana()
{
  gROOT->ProcessLine(" .L style.cc+");
  //**********************************************INPUT******************************
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r210614_autoforest_v51.root"); sample_name.push_back("data"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root"); sample_name.push_back("Epos"); sample_type.push_back(MC); 
  sample_fname.push_back("root://eoscms//eos/cms//store/caf/user/dgulhan/pPb_Hijing_MB/HiForest_v03_mergedv02/merged_forest_0.root"); sample_name.push_back("HIJING"); sample_type.push_back(MC); 
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/pPb_5020_QGSJetII_5_3_8_HI_patch2/forest.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC); 

  const int n_cas_rechits = 224;
  const int n_hf_rechits = 1728;
  
  int hf_n;
  float hf_e[n_hf_rechits];
  int hf_depth[n_hf_rechits];
  float hf_eta[n_hf_rechits];

  float cas_e[n_cas_rechits];
  int cas_sat[n_cas_rechits];
  int cas_depth[n_cas_rechits];
  int cas_iphi[n_cas_rechits];

  int nTrk;

  int zero_bias;
  int min_bias;
  int random;
  int bptx_p_m;
  int bptx_p_nm;
  int bptx_np_m;

  //*********************************OUTPUT************************+++

  TFile* out_file = new TFile("histos.root","RECREATE");

  TH1D* h_zero_count_zb_coll;
  TH1D* h_zero_count_zb_no_coll;
  TH1D* h_hf_hits_coll;
  TH1D* h_hf_hits_plus_min;
  TH1D* h_hf_hits_minus_min;
  TH1D* h_hf_hits_noise;
  TH1D* h_hf_hits_beam_gas;
  TH1D* h_castor_hf_diff_3;
  TH1D* h_castor_hf_diff_5;
  TH1D* h_castor_hf_diff_10;
  TH1D* h_eff;



  //**********************************LOOP*******************

  for (int sample=0; sample<int(sample_name.size()); sample++)
    {
      TFile* file = NULL; //****FILE
      file = TFile::Open(sample_fname[sample].c_str());
      if (!file || !file->IsOpen() || file->IsZombie()) cerr << "Cannot find file" << endl;
      else cout << "file opened" << endl;
      
      TTree* tree = NULL; //****TREE0
      ((TDirectory*)file->Get("rechitanalyzer"))->GetObject("castor",tree);
      if (!tree) cerr << "Cannot find tree" << endl;
      
      TTree* friend_tree0 = NULL; //****TREE1
      ((TDirectory*)file->Get("pfTowers"))->GetObject("hf",friend_tree0);
      if (!friend_tree0) cerr << "Cannot find friend tree" << endl;
      tree->AddFriend(friend_tree0);
      
      TTree* friend_tree1 = NULL; //****TREE2
      ((TDirectory*)file->Get("hltanalysis"))->GetObject("HltTree",friend_tree1);
      if (!friend_tree1) cerr << "Cannot find friend tree" << endl;
      tree->AddFriend(friend_tree1);
      
      TTree* friend_tree2 = NULL; //****TREE3
      ((TDirectory*)file->Get("pixelTrack"))->GetObject("trackTree",friend_tree2);
      if (!friend_tree2) cerr << "Cannot find friend tree" << endl;
      tree->AddFriend(friend_tree2);
      
      tree->SetBranchStatus("*",0);

      tree->SetBranchStatus("hf.n",1);
      tree->SetBranchStatus("hf.e",1);
      tree->SetBranchStatus("hf.depth",1);
      tree->SetBranchStatus("hf.eta",1);

      tree->SetBranchStatus("e",1);
      tree->SetBranchStatus("saturation",1);
      tree->SetBranchStatus("depth",1);
      tree->SetBranchStatus("iphi",1);


      tree->SetBranchStatus("nTrk",1);

      tree->SetBranchStatus("HLT_PAZeroBias_v1",1);
      tree->SetBranchStatus("HLT_PAL1Tech53_MB_SingleTrack_v1",1);
      tree->SetBranchStatus("HLT_PARandom_v1",1);
      tree->SetBranchStatus("L1Tech_BPTX_plus_AND_minus.v0",1);
      tree->SetBranchStatus("L1Tech_BPTX_plus_AND_NOT_minus.v0",1);
      tree->SetBranchStatus("L1Tech_BPTX_minus_AND_not_plus.v0",1);
      //________________________________

      tree->SetBranchAddress("hf.n",&hf_n);
      tree->SetBranchAddress("hf.e",hf_e);
      tree->SetBranchAddress("hf.depth",hf_depth);
      tree->SetBranchAddress("hf.eta",hf_eta);

      tree->SetBranchAddress("e",cas_e);
      tree->SetBranchAddress("saturation",cas_sat);
      tree->SetBranchAddress("depth",cas_depth);
      tree->SetBranchAddress("iphi",cas_iphi);


      tree->SetBranchAddress("nTrk",&nTrk);

      tree->SetBranchAddress("HLT_PAZeroBias_v1",&zero_bias);
      tree->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1",&min_bias);
      tree->SetBranchAddress("HLT_PARandom_v1",&random);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0",&bptx_p_m);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0",&bptx_p_nm);
      tree->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0",&bptx_np_m);

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      h_zero_count_zb_coll      = new TH1D((sample_name[sample] + string("_h_zero_count_zb_coll")).c_str(),"",100,728,1728);
      h_zero_count_zb_no_coll   = new TH1D((sample_name[sample] + string("_h_zero_count_zb_no_coll")).c_str(),"",100,728,1728);
      h_hf_hits_coll            = new TH1D((sample_name[sample] + string("_h_hf_hits_coll")).c_str(),"",200,0,200);
      h_hf_hits_plus_min        = new TH1D((sample_name[sample] + string("_h_hf+_hits_min")).c_str(),"",200,0,200);
      h_hf_hits_minus_min       = new TH1D((sample_name[sample] + string("_h_hf-_hits_min")).c_str(),"",200,0,200);
      h_hf_hits_noise           = new TH1D((sample_name[sample] + string("_h_hf_hits_noise")).c_str(),"",200,0,200);
      h_hf_hits_beam_gas        = new TH1D((sample_name[sample] + string("_h_hf_hits_beam_gas")).c_str(),"",200,0,200);
      h_castor_hf_diff_3        = new TH1D((sample_name[sample] + string("_h_castor_hf_diff_3")).c_str(),"",100,0,10000);
      h_castor_hf_diff_5        = new TH1D((sample_name[sample] + string("_h_castor_hf_diff_5")).c_str(),"",100,0,10000);
      h_castor_hf_diff_10       = new TH1D((sample_name[sample] + string("_h_castor_hf_diff_10")).c_str(),"",100,0,10000);
      h_eff                     = new TH1D((sample_name[sample] + string("_h_eff")).c_str(),"",9,-0.5,8.5);
      
      h_eff->GetXaxis()->SetBinLabel(1,"MinBias");
      h_eff->GetXaxis()->SetBinLabel(2,"ZeroBias");
      h_eff->GetXaxis()->SetBinLabel(3,"Noise");
      h_eff->GetXaxis()->SetBinLabel(4,"Beam Gas");
      h_eff->GetXaxis()->SetBinLabel(5,"HF energy cut >3GeV");
      h_eff->GetXaxis()->SetBinLabel(6,"HF energy cut >4GeV");
      h_eff->GetXaxis()->SetBinLabel(7,"HF energy cut >5GeV");
      h_eff->GetXaxis()->SetBinLabel(8,"HF energy cut >3GeV (noise or bg)");
      h_eff->GetXaxis()->SetBinLabel(9,"HF energy cut >4GeV (noise or bg)");


      for(int iEvent=0; iEvent<tree->GetEntries(); iEvent++)
        {
          if(iEvent==MAXEVT) break;
          if(iEvent % 10000 == 0) cout << sample << "/" << sample_name.size() << " " << sample_name[sample].c_str() << " Entry: " << iEvent << " / " << (MAXEVT>0?MAXEVT:tree->GetEntries()) << endl;
          tree->GetEntry(iEvent);
          
          bool coll=0, noise=0, beam_gas=0;

          coll          = zero_bias && bptx_p_m; //double beam
          beam_gas      = (bptx_np_m || bptx_p_nm); // only single beam
          noise         = !bptx_p_m && !bptx_np_m && !bptx_p_nm; //not both and not single beam

          if(sample_type[sample] == MC)
            {
              beam_gas = 0;
              noise = 0;
              coll = 1;
              min_bias = 1;
            }

          if(!coll && !noise && !beam_gas && !min_bias) //not intersted
            continue;
          

          double sum_cas_e_em = 0;
          double sum_cas_e_had = 0;          
          
          int hf_zero_count = n_hf_rechits - hf_n;
          
          for(int ch_cas=0; ch_cas<n_cas_rechits; ch_cas++) // no ZS
            {//break;
              if(cas_depth[ch_cas] < 3)
                sum_cas_e_em += cas_e[ch_cas];
              else
                sum_cas_e_had += cas_e[ch_cas];
            }
          double hf_double_energy_max = -1;
          double hf_m_energy_max = -1;
          double hf_p_energy_max = -1;
          for(int ch_hf=0; ch_hf<hf_n; ch_hf++)
            {
              if(hf_eta[ch_hf] > 0. && hf_e[ch_hf] > hf_p_energy_max)
                hf_p_energy_max = hf_e[ch_hf];
              if(hf_eta[ch_hf] <= 0. && hf_e[ch_hf] > hf_m_energy_max)
                hf_m_energy_max = hf_e[ch_hf];
            }
          hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
          //Filling HISTOS
          if(coll)                                    h_zero_count_zb_coll->Fill(hf_zero_count);
          if(noise || beam_gas)                       h_zero_count_zb_no_coll->Fill(hf_zero_count);

          if(coll)                                    h_hf_hits_coll->Fill(hf_double_energy_max);
          if(min_bias)                                h_hf_hits_plus_min->Fill(hf_p_energy_max);
          if(min_bias)                                h_hf_hits_minus_min->Fill(hf_m_energy_max);
          if(noise)                                   h_hf_hits_noise->Fill(hf_double_energy_max);
          if(beam_gas)                                h_hf_hits_beam_gas->Fill(hf_double_energy_max);

          if(coll && hf_double_energy_max < 3)        h_castor_hf_diff_3->Fill(sum_cas_e_em);
          if(coll && hf_double_energy_max < 5)        h_castor_hf_diff_5->Fill(sum_cas_e_em);
          if(coll && hf_double_energy_max < 10)       h_castor_hf_diff_10->Fill(sum_cas_e_em);

          if(min_bias)                                h_eff->Fill(0);
          if(coll)                                    h_eff->Fill(1);
          if(noise)                                   h_eff->Fill(2);
          if(beam_gas)                                h_eff->Fill(3);
          if(coll && hf_double_energy_max > 3)        h_eff->Fill(4);
          if(coll && hf_double_energy_max > 4)        h_eff->Fill(5);
          if(coll && hf_double_energy_max > 5)        h_eff->Fill(6);
          if((noise || beam_gas) && hf_double_energy_max > 3)       h_eff->Fill(7);
          if((noise || beam_gas) && hf_double_energy_max > 4)       h_eff->Fill(8);
          
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
