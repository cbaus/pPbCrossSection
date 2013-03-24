#define MAXEVT 20000

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"

#include <iostream>
#include <vector>
#include <string>

#include <CastorTreeVariables.h>
#include <ParticleInfo.h>

//#include "style.h"

enum e_type
  {
    MC = 0,
    DATA
  };

using namespace std;

vector<string> sample_fname;
vector<string> sample_name;
vector<e_type> sample_type;

int main()
{
  TH1::SetDefaultSumw2();
  //style();

  //*************************************************************INPUT***********************************************************
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data/*_1*_*.root"); sample_name.push_back("data"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("Epos"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Hijing/*.root"); sample_name.push_back("Hijing"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/QGSJetII/reeMC.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlighDPMjet/treeMC.root"); sample_name.push_back("Starlight_DPMJet");  sample_type.push_back(MC);

  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos.root","RECREATE");

  TH1D* h_zero_count_zb_coll;
  TH1D* h_zero_count_zb_no_coll;
  TH1D* h_hf_hits_coll;
  TH1D* h_hf_hits_plus;
  TH1D* h_hf_hits_minus;
  TH1D* h_hf_hits_noise;
  TH1D* h_hf_hits_beam_gas;
  TH1D* h_castor_hf_diff_3;
  TH1D* h_castor_hf_diff_5;
  TH1D* h_castor_hf_diff_10;
  TH1D* h_castor_gap_hf;
  TH1D* h_castor_nogap_hf;
  TH1D* h_eff;

  TH2D* h_hf_fsc_p;
  TH2D* h_hf_fsc_m;

  TH1D* h_hf_cut_single;
  TH1D* h_hf_cut_single_noise;
  TH1D* h_hf_cut_double;
  TH1D* h_hf_cut_double_noise;

  TProfile* h_hf_hits_coll_lumi;
  TProfile* h_hf_hits_plus_lumi;
  TProfile* h_hf_hits_minus_lumi;
  TProfile* h_hf_hits_noise_lumi;
  TProfile* h_hf_totE_coll_lumi;
  TProfile* h_hf_totE_plus_lumi;
  TProfile* h_hf_totE_minus_lumi;
  TProfile* h_hf_totE_noise_lumi;
  TH1D* h_lumi;
  TH1D* h_lumi_3GeV;

  TH1D* h_perf_hf_rechits_single_3gev;
  TH1D* h_perf_hf_rechits_double_1dot5gev;
  TH1D* h_perf_hf_totE_single_3gev;
  TH1D* h_perf_hf_totE_double_1dot5gev;
  TH1D* h_perf_hf_totE_eta_single_3gev;
  TH1D* h_perf_hf_totE_eta_double_1dot5gev;

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

      //________________________________________________________________________________

      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("event", 1);
      tree->SetBranchStatus("genProcessID", 1);
      tree->SetBranchStatus("*HFtowers*", 1);
      tree->SetBranchStatus("*HFtowers.*", 1);
      tree->SetBranchStatus("CASTOR*", 1);
      tree->SetBranchStatus("*lumi*", 1);

      tree->SetBranchStatus("HLT_PAZeroBias_v1",1);
      tree->SetBranchStatus("HLT_PAL1Tech53_MB_SingleTrack_v1",1);
      tree->SetBranchStatus("HLT_PARandom_v1",1);
      tree->SetBranchStatus("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",1);
      tree->SetBranchStatus("L1Tech_BPTX_plus_AND_NOT_minus.v0_DecisionBeforeMask",1);
      tree->SetBranchStatus("L1Tech_BPTX_minus_AND_not_plus.v0_DecisionBeforeMask",1);

      //________________________________________________________________________________


      AnalysisEvent* event = 0;
      tree->SetBranchAddress("event", &event);
      int genProcessID;
      tree->SetBranchAddress("genProcessID", &genProcessID);

      float fsc_sum_minus;
      float fsc_sum_plus;

      int zero_bias;
      int min_bias;
      int random;
      int bptx_p_m;
      int bptx_p_nm;
      int bptx_np_m;
      tree->SetBranchAddress("HLT_PAZeroBias_v1",&zero_bias);
      tree->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1",&min_bias);
      tree->SetBranchAddress("HLT_PARandom_v1",&random);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",&bptx_p_m);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_DecisionBeforeMask",&bptx_p_nm);
      tree->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_DecisionBeforeMask",&bptx_np_m);

      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];
      h_zero_count_zb_coll      = new TH1D((add + string("_h_zero_count_zb_coll")).c_str(),"",100,728,1728);
      h_zero_count_zb_no_coll   = new TH1D((add + string("_h_zero_count_zb_no_coll")).c_str(),"",100,728,1728);
      h_hf_hits_coll            = new TH1D((add + string("_h_hf_hits_coll")).c_str(),"",200,0,200);
      h_hf_hits_plus            = new TH1D((add + string("_h_hf+_hits_min")).c_str(),"",200,0,200);
      h_hf_hits_minus           = new TH1D((add + string("_h_hf-_hits_min")).c_str(),"",200,0,200);
      h_hf_hits_noise           = new TH1D((add + string("_h_hf_hits_noise")).c_str(),"",200,0,200);
      h_hf_hits_beam_gas        = new TH1D((add + string("_h_hf_hits_beam_gas")).c_str(),"",200,0,200);
      h_castor_hf_diff_3        = new TH1D((add + string("_h_castor_hf_diff_3")).c_str(),"",100,0,10000);
      h_castor_hf_diff_5        = new TH1D((add + string("_h_castor_hf_diff_5")).c_str(),"",100,0,10000);
      h_castor_hf_diff_10       = new TH1D((add + string("_h_castor_hf_diff_10")).c_str(),"",100,0,10000);
      h_castor_gap_hf           = new TH1D((add + string("_h_castor_gap_hf")).c_str(),"",100,0,50);
      h_castor_nogap_hf         = new TH1D((add + string("_h_castor_nogap_hf")).c_str(),"",100,0,50);
      h_eff                     = new TH1D((add + string("_h_eff")).c_str(),"",9,-0.5,8.5);

      h_hf_cut_single           = new TH1D((add + string("_h_hf_cut_single")).c_str(),"",21,-0.125,5.125);
      h_hf_cut_single_noise     = new TH1D((add + string("_h_hf_cut_single_noise")).c_str(),"",21,-0.125,5.125);
      h_hf_cut_double           = new TH1D((add + string("_h_hf_cut_double")).c_str(),"",21,-0.125,5.125);
      h_hf_cut_double_noise     = new TH1D((add + string("_h_hf_cut_double_noise")).c_str(),"",21,-0.125,5.125);

      h_hf_fsc_p                = new TH2D((add + string("_h_hf_fsc_p")).c_str(),"",200,0,200,200,0,80000);
      h_hf_fsc_m                = new TH2D((add + string("_h_hf_fsc_m")).c_str(),"",200,0,200,200,0,80000);

      h_hf_hits_coll_lumi       = new TProfile((add + string("_h_hf_hits_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_minus_lumi      = new TProfile((add + string("_h_hf+_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_plus_lumi       = new TProfile((add + string("_h_hf-_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_noise_lumi      = new TProfile((add + string("_h_hf_hits_noise_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_coll_lumi       = new TProfile((add + string("_h_hf_totE_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_minus_lumi      = new TProfile((add + string("_h_hf+_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_plus_lumi       = new TProfile((add + string("_h_hf-_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_noise_lumi      = new TProfile((add + string("_h_hf_totE_noise_lumi")).c_str(),"",2000,0,2000);
      h_lumi                    = new TH1D((add + string("_h_lumi")).c_str(),"",2000,0,2000);
      h_lumi_3GeV               = new TH1D((add + string("_h_lumi_3GeV")).c_str(),"",2000,0,2000);

      h_perf_hf_rechits_single_3gev        = new TH1D((add + string("_h_perf_hf_rechits_single_3gev")).c_str(),"",100,0,1000);
      h_perf_hf_rechits_double_1dot5gev    = new TH1D((add + string("_h_perf_hf_rechits_double_1dot5gev")).c_str(),"",100,0,1000);
      h_perf_hf_totE_single_3gev           = new TH1D((add + string("_h_perf_hf_totE_single_3gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_double_1dot5gev       = new TH1D((add + string("_h_perf_hf_totE_double_1dot5gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_eta_single_3gev       = new TH1D((add + string("_h_perf_hf_totE_eta_single_3gev")).c_str(),"",100,-5.2,5.2);
      h_perf_hf_totE_eta_double_1dot5gev   = new TH1D((add + string("_h_perf_hf_totE_eta_double_1dot5gev")).c_str(),"",100,-5.2,5.2);


      h_eff->GetXaxis()->SetBinLabel(1,"MinBias");
      h_eff->GetXaxis()->SetBinLabel(2,"ZeroBias");
      h_eff->GetXaxis()->SetBinLabel(3,"Noise");
      h_eff->GetXaxis()->SetBinLabel(4,"Beam Gas");
      h_eff->GetXaxis()->SetBinLabel(5,"HF double >2GeV");
      h_eff->GetXaxis()->SetBinLabel(6,"HF double >3GeV");
      h_eff->GetXaxis()->SetBinLabel(7,"HF single >3GeV");
      h_eff->GetXaxis()->SetBinLabel(8,"HF double >2GeV (noise or bg)");
      h_eff->GetXaxis()->SetBinLabel(9,"HF double >3GeV (noise or bg)");


      for(int iEvent=0; iEvent<tree->GetEntries(); iEvent++)
        {
          if(iEvent==MAXEVT) break;
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << (MAXEVT>0?MAXEVT:tree->GetEntries()) << endl;
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
          double sum_cas_e = 0;
          for (vector<RecHitCASTOR>::const_iterator it = event->CASTOR.begin(); it < event->CASTOR.end(); ++it)
            {//break;
              if(it->GetModuleId() < 3)
                sum_cas_e_em += it->Energy;
              else if(it->GetModuleId() < 5)
                sum_cas_e_had += it->Energy;
            }
          sum_cas_e = sum_cas_e_had + sum_cas_e_em;

          int hf_n = event->HFtowers.size();
          int hf_zero_count = ForwardRecord::nMaxHFMRecHits - hf_n;
          double hf_double_energy_max = -1;
          double hf_single_energy_max = -1;
          double hf_m_energy_max = -1;
          double hf_p_energy_max = -1;
          double hf_p_energy = 0;
          double hf_m_energy = 0;
          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              //cout << it->Eta << " " << it->Energy << endl;
              if(it->Eta > 0. && it->Energy > hf_p_energy_max)
                hf_p_energy_max = it->Energy;
              if(it->Eta <= 0. && it->Energy > hf_m_energy_max)
                hf_m_energy_max = it->Energy;

              if(it->Eta > 0.)
                hf_p_energy += it->Energy;
              else
                hf_m_energy += it->Energy;

            }
          double hf_pm_energy = hf_p_energy + hf_m_energy;
          //cout << "rechit: " << hf_pm_energy << endl;

          hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
          hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
          //Filling HISTOS
          if(coll)                                                  h_zero_count_zb_coll->Fill(hf_zero_count);
          if(noise || beam_gas)                                     h_zero_count_zb_no_coll->Fill(hf_zero_count);

          if(coll)                                                  h_hf_hits_coll->Fill(hf_double_energy_max);
          if(coll)                                                  h_hf_hits_plus->Fill(hf_p_energy_max);
          if(coll)                                                  h_hf_hits_minus->Fill(hf_m_energy_max);
          if(noise)                                                 h_hf_hits_noise->Fill(hf_double_energy_max);
          if(beam_gas)                                              h_hf_hits_beam_gas->Fill(hf_double_energy_max);

          if(coll && hf_double_energy_max < 3)                      h_castor_hf_diff_3->Fill(sum_cas_e);
          if(coll && hf_double_energy_max < 5)                      h_castor_hf_diff_5->Fill(sum_cas_e);
          if(coll && hf_double_energy_max < 10)                     h_castor_hf_diff_10->Fill(sum_cas_e);

          if(coll && sum_cas_e <= 700)                              h_castor_gap_hf->Fill(hf_single_energy_max);
          if(coll && sum_cas_e >  700)                              h_castor_nogap_hf->Fill(hf_single_energy_max);

          if(min_bias)                                              h_eff->Fill(0);
          if(coll)                                                  h_eff->Fill(1);
          if(noise)                                                 h_eff->Fill(2);
          if(beam_gas)                                              h_eff->Fill(3);
          if(coll && hf_double_energy_max > 2)                      h_eff->Fill(4);
          if(coll && hf_double_energy_max > 3)                      h_eff->Fill(5);
          if(coll && hf_single_energy_max > 3)                      h_eff->Fill(6);
          if((noise || beam_gas) && hf_double_energy_max > 2)       h_eff->Fill(7);
          if((noise || beam_gas) && hf_double_energy_max > 3)       h_eff->Fill(8);

          for (double cut=0; cut < 5.5; cut+=0.25)
            {
              if(coll && hf_double_energy_max > cut)                h_hf_cut_double->Fill(cut);
              if((noise || beam_gas) && hf_double_energy_max > cut) h_hf_cut_double_noise->Fill(cut);
              if(coll && hf_single_energy_max > cut)                h_hf_cut_single->Fill(cut);
              if((noise || beam_gas) && hf_single_energy_max > cut) h_hf_cut_single_noise->Fill(cut);
            }

          if(coll)                                                  h_hf_fsc_p->Fill(hf_p_energy_max,fsc_sum_plus);
          if(coll)                                                  h_hf_fsc_m->Fill(hf_m_energy_max,fsc_sum_minus);

          if(coll && hf_double_energy_max > 1.5)                    h_hf_hits_coll_lumi->Fill(event->lumiNb,hf_double_energy_max);
          if(coll && hf_double_energy_max > 1.5)                    h_hf_hits_minus_lumi->Fill(event->lumiNb,hf_m_energy_max);
          if(coll && hf_double_energy_max > 1.5)                    h_hf_hits_plus_lumi->Fill(event->lumiNb,hf_p_energy_max);
          if((noise || beam_gas))                                   h_hf_hits_noise_lumi->Fill(event->lumiNb,hf_double_energy_max);
          if(coll && hf_double_energy_max > 1.5)                    h_hf_totE_coll_lumi->Fill(event->lumiNb,hf_pm_energy);
          if(coll && hf_double_energy_max > 1.5)                    h_hf_totE_minus_lumi->Fill(event->lumiNb,hf_m_energy);
          if(coll && hf_double_energy_max > 1.5)                    h_hf_totE_plus_lumi->Fill(event->lumiNb,hf_p_energy);
          if((noise || beam_gas))                                   h_hf_totE_noise_lumi->Fill(event->lumiNb,hf_pm_energy);
          if(coll && hf_double_energy_max > 1.5)                    h_lumi_3GeV->Fill(event->lumiNb);
          if(coll)                                                  h_lumi->Fill(event->lumiNb);

          if(coll && hf_single_energy_max > 3)                      h_perf_hf_rechits_single_3gev->Fill(hf_n);
          if(coll && hf_double_energy_max > 1.5)                    h_perf_hf_rechits_double_1dot5gev->Fill(hf_n);

          if(coll && hf_single_energy_max > 3)                      h_perf_hf_totE_single_3gev->Fill(hf_pm_energy);
          if(coll && hf_double_energy_max > 1.5)                    h_perf_hf_totE_double_1dot5gev->Fill(hf_pm_energy);

          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              if(coll && hf_double_energy_max > 1.5)
                h_perf_hf_totE_eta_double_1dot5gev->Fill(it->Eta,it->Energy);
              if(coll && hf_single_energy_max > 3)
                h_perf_hf_totE_eta_single_3gev->Fill(it->Eta,it->Energy);
            }

        }

      //******************************************AFTER EVENT LOOP*******************************************
      double integ_lumi = 1854.344875;// nb^-1
      double n_total = double(tree->GetEntries());
      double n_zb = h_eff->GetBinContent(1) * (MAXEVT>0?n_total/double(MAXEVT):1.) * 40.;
      double n_mb = h_eff->GetBinContent(4) * (MAXEVT>0?n_total/double(MAXEVT):1.) * 40.;

      cout << endl << "Cross section ZB: " << n_zb/integ_lumi << " --- cross section MB: " << n_mb/integ_lumi << endl;
      h_perf_hf_totE_eta_double_1dot5gev->Scale(1./double(MAXEVT>0?MAXEVT:n_total));
      h_perf_hf_totE_eta_single_3gev->Scale(1./double(MAXEVT>0?MAXEVT:n_total));

    }

  //********************************************AFTER FILE LOOP************************************************

  out_file->Write();
  out_file->Save();

  return 0;
}
