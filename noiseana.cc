#define MAXEVT -1000000

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TTree.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <utility>

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
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/DataNoise/*.root"); sample_name.push_back("noise"); sample_type.push_back(DATA);

  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos_noise.root","RECREATE");


  TH1D* h_hf_hits_noise;
  TH1D* h_hf_hits_beamgas;

  TH2D* h_noise_tracks_hf;
  TH2D* h_beamgas_tracks_hf;

  TH1D* h_hf_cut_single_noise;
  TH1D* h_hf_cut_single_beamgas;
  TH1D* h_hf_cut_double_noise;
  TH1D* h_hf_cut_double_beamgas;


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

      tree->SetBranchStatus("*", 1);
//       tree->SetBranchStatus("event", 1);
//       tree->SetBranchStatus("genProcessID", 1);
//       tree->SetBranchStatus("GEN*", 1);
//       tree->SetBranchStatus("*HFtowers*", 1);
//       tree->SetBranchStatus("CASTOR*", 1);
//       tree->SetBranchStatus("*lumi*", 1);
//       tree->SetBranchStatus("*Lumi*", 1);

//       tree->SetBranchStatus("HLT_PAZeroBias_v1",1);
//       tree->SetBranchStatus("HLT_PAL1Tech53_MB_SingleTrack_v1",1);
//       tree->SetBranchStatus("HLT_PARandom_v1",1);
//       tree->SetBranchStatus("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",1);
//       tree->SetBranchStatus("L1Tech_BPTX_plus_AND_NOT_minus.v0_DecisionBeforeMask",1);
//       tree->SetBranchStatus("L1Tech_BPTX_minus_AND_not_plus.v0_DecisionBeforeMask",1);

      //________________________________________________________________________________


      AnalysisEvent* event = 0;
      tree->SetBranchAddress("event", &event);

      float fsc_sum_minus = 0;
      float fsc_sum_plus = 0;

      int zero_bias;
      int zero_bias_prescale_L1;
      int zero_bias_prescale_HLT;
      int min_bias;
      int random;
      int random_prescale_HLT;
      int bptx_p_m;
      int bptx_p_nm;
      int bptx_np_m;
      int bptx_quiet;
      tree->SetBranchAddress("L1_ZeroBias_algPrescale",&zero_bias_prescale_L1);
      tree->SetBranchAddress("HLT_PAZeroBias_v1_Prescl",&zero_bias_prescale_HLT);
      tree->SetBranchAddress("HLT_PAZeroBias_v1",&zero_bias);
      tree->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1",&min_bias);
      tree->SetBranchAddress("HLT_PARandom_v1",&random);
      tree->SetBranchAddress("HLT_PARandom_v1_Prescl",&random_prescale_HLT);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",&bptx_p_m);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_DecisionBeforeMask",&bptx_p_nm);
      tree->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_DecisionBeforeMask",&bptx_np_m);
      tree->SetBranchAddress("L1Tech_BPTX_quiet.v0_DecisionBeforeMask",&bptx_quiet);

      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];
      h_hf_hits_noise           = new TH1D((add + string("_h_hf_hits_noise")).c_str(),"",200,0,20);
      h_hf_hits_beamgas        = new TH1D((add + string("_h_hf_hits_beamgas")).c_str(),"",200,0,20);

      h_hf_cut_single_noise     = new TH1D((add + string("_h_hf_cut_single_noise")).c_str(),"",101,-0.05,10.05);
      h_hf_cut_single_beamgas   = new TH1D((add + string("_h_hf_cut_single_beamgas")).c_str(),"",101,-0.05,10.05);
      h_hf_cut_double_noise     = new TH1D((add + string("_h_hf_cut_double_noise")).c_str(),"",101,-0.05,10.05);
      h_hf_cut_double_beamgas   = new TH1D((add + string("_h_hf_cut_double_beamgas")).c_str(),"",101,-0.05,10.05);

      h_noise_tracks_hf   = new TH2D((add + string("_h_noise_tracks_hf")).c_str(),"",200,0,20,20,0,20);
      h_beamgas_tracks_hf   = new TH2D((add + string("_h_beamgas_tracks_hf")).c_str(),"",200,0,20,20,0,20);

      for (int i=-41; i<=-30; i++)
        {
          ostringstream hf_str, hf_cut_str;
          hf_str  << add << "h_hf_rings_single_" << i;
          hf_cut_str  << add << "h_hf_rings_cut_single_" << i;
          TH1D* h_hf_rings_single = new TH1D(hf_str.str().c_str(),hf_str.str().c_str(),2000,0,200);
          TH1D* h_hf_rings_cut_single = new TH1D(hf_cut_str.str().c_str(),hf_cut_str.str().c_str(),101,-0.05,10.05);
        }
      for (int i=30; i<=41; i++)
        {
          ostringstream hf_str, hf_cut_str;
          hf_str  << add << "h_hf_rings_single_" << i;
          hf_cut_str  << add << "h_hf_rings_cut_single_" << i;
          TH1D* h_hf_rings_single = new TH1D(hf_str.str().c_str(),hf_str.str().c_str(),2000,0,200);
          TH1D* h_hf_rings_cut_single = new TH1D(hf_cut_str.str().c_str(),hf_cut_str.str().c_str(),101,-0.05,10.05);
        }

      set<int> unpaired;
      unpaired.insert(465);
      unpaired.insert(473);
      unpaired.insert(479);
      unpaired.insert(482);
      unpaired.insert(487);
      unpaired.insert(490);
      unpaired.insert(496);
      unpaired.insert(499);
      unpaired.insert(504);
      unpaired.insert(507);
      unpaired.insert(513);
      unpaired.insert(516);
      unpaired.insert(521);
      unpaired.insert(524);
      unpaired.insert(530);
      unpaired.insert(533);
      unpaired.insert(538);
      unpaired.insert(541);
      unpaired.insert(547);
      unpaired.insert(550);
      unpaired.insert(555);
      unpaired.insert(558);
      unpaired.insert(564);
      unpaired.insert(567);
      unpaired.insert(572);
      unpaired.insert(575);
      unpaired.insert(581);
      unpaired.insert(584);
      unpaired.insert(589);
      unpaired.insert(592);
      unpaired.insert(598);
      unpaired.insert(601);
      unpaired.insert(606);
      unpaired.insert(609);
      unpaired.insert(615);
      unpaired.insert(618);
      unpaired.insert(623);
      unpaired.insert(626);
      unpaired.insert(632);
      unpaired.insert(635);
      unpaired.insert(640);
      unpaired.insert(643);
      unpaired.insert(649);
      unpaired.insert(652);
      unpaired.insert(657);
      unpaired.insert(660);
      unpaired.insert(666);
      unpaired.insert(674);
      unpaired.insert(697);
      unpaired.insert(705);
      unpaired.insert(711);
      unpaired.insert(714);
      unpaired.insert(719);
      unpaired.insert(722);
      unpaired.insert(728);
      unpaired.insert(731);
      unpaired.insert(736);
      unpaired.insert(739);
      unpaired.insert(745);
      unpaired.insert(748);
      unpaired.insert(753);
      unpaired.insert(756);
      unpaired.insert(762);
      unpaired.insert(765);
      unpaired.insert(770);
      unpaired.insert(773);
      unpaired.insert(779);
      unpaired.insert(782);
      unpaired.insert(787);
      unpaired.insert(790);
      unpaired.insert(796);
      unpaired.insert(799);
      unpaired.insert(804);
      unpaired.insert(807);
      unpaired.insert(813);
      unpaired.insert(816);
      unpaired.insert(821);
      unpaired.insert(824);
      unpaired.insert(830);
      unpaired.insert(833);
      unpaired.insert(838);
      unpaired.insert(841);
      unpaired.insert(847);
      unpaired.insert(855);

      for(int iEvent=0; iEvent<tree->GetEntries(); iEvent++)
        {
          if(iEvent==MAXEVT) break;
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << (MAXEVT>0?MAXEVT:tree->GetEntries()) << endl;
          tree->GetEntry(iEvent);
          //          if(event->runNb != 210885)
          //continue;

          bool noise=0, beamgas=0;

          beamgas      = (bptx_np_m || bptx_p_nm); // only single beam
          noise         = bptx_quiet;// && !bptx_np_m && !bptx_p_nm; //not both and not single beam

           if(unpaired.count(event->bxNb))
            {
              noise = 0;
              beamgas = 1;
            }

          if(!noise && !beamgas) //not intersted
            continue;


          const int prescale      = random_prescale_HLT;
          const double lumiPerLS     = event->instLuminosity * event->instLuminosityCorr * 1e6;
          const double lumiPerTime     = lumiPerLS / 23.31;
          const double evtWeight  = lumiPerTime?double(prescale) / lumiPerTime:0.;

          //---------------------------------------------CASTOR
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





          //---------------------------------------------HF
          int hf_n = event->HFtowers.size();
          int hf_zero_count = ForwardRecord::nMaxHFMRecHits - hf_n;
          double hf_double_energy_max = 0;
          double hf_single_energy_max = 0;
          double hf_m_energy_max = 0;
          double hf_p_energy_max = 0;
          double hf_p_energy = 0;
          double hf_m_energy = 0;
          map<int,double> hf_ring_map_single_max;
          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              const int Ieta = it->Eta > 0?it->IetaAbs:-it->IetaAbs;
              if(hf_ring_map_single_max.count(Ieta) == 0)
                hf_ring_map_single_max.insert(pair<int,double>(Ieta,it->Energy));
              else if (hf_ring_map_single_max[Ieta] < it->Energy)
                hf_ring_map_single_max[Ieta]=it->Energy;


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


          //Booking run histograms

          ostringstream run_str;
          run_str << event->runNb;
          TH1D* h_run_events_single = NULL;
          TH1D* h_run_events_double = NULL;
          TH1D* h_run_events = NULL;
          h_run_events_single = (TH1D*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events_single")).c_str()));
          h_run_events_double = (TH1D*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events_double")).c_str()));
          h_run_events = (TH1D*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events")).c_str()));
          if(h_run_events_single == NULL)
            {
              h_run_events_single = new TH1D((add+run_str.str()+string("_h_run_events_single")).c_str(),run_str.str().c_str(),2000,0,2e5);
              h_run_events_double = new TH1D((add+run_str.str()+string("_h_run_events_double")).c_str(),run_str.str().c_str(),2000,0,2e5);
              h_run_events = new TH1D((add+run_str.str()+string("_h_run_events")).c_str(),run_str.str().c_str(),2000,0,2e5);
            }

          //HF Rings

          int ieta = -41;
          for (int i=0; i<24; i++)
            {

              double hot = hf_ring_map_single_max.count(ieta) ? hf_ring_map_single_max[ieta] : 0;

              //map<int,double>::const_iterator it = hf_ring_map_single_max.begin();
              //while(noise && it != hf_ring_map_single_max.end())
              ostringstream hf_str, hf_cut_str;
              hf_str  << add << "h_hf_rings_single_" << ieta;
              hf_cut_str  << add << "h_hf_rings_cut_single_" << ieta;
              TH1D* h_hf_rings_single = NULL;
              TH1D* h_hf_rings_cut_single = NULL;
              h_hf_rings_single = (TH1D*)(out_file->Get((add+string("/")+hf_str.str()).c_str()));
              h_hf_rings_cut_single = (TH1D*)(out_file->Get((add+string("/")+hf_cut_str.str()).c_str()));
              if(h_hf_rings_single == NULL || h_hf_rings_cut_single == NULL)
                {
                  cerr << "no histo " << ieta << endl;
                  return 1;
                }
              h_hf_rings_single->Fill(hot);

              for (double cut=0; cut <= 10; cut+=0.1)
                {
                  if(hot >= cut) h_hf_rings_cut_single->Fill(cut);
                }

              //++it;
              if(ieta==-30) ieta=30;
              else ieta++;
            }


          //---------------------------------------------Filling HISTOS
          if(noise)                                                h_hf_hits_noise->Fill(hf_double_energy_max);
          if(beamgas)                                              h_hf_hits_beamgas->Fill(hf_double_energy_max);

          for (double cut=0; cut <= 10; cut+=0.1)
            {
              if((noise || beamgas) && hf_double_energy_max >= cut) h_hf_cut_double_noise->Fill(cut);
              if((beamgas || beamgas) && hf_double_energy_max >= cut) h_hf_cut_double_beamgas->Fill(cut);
              if((noise || beamgas) && hf_single_energy_max >= cut) h_hf_cut_single_noise->Fill(cut);
              if((beamgas || beamgas) && hf_single_energy_max >= cut) h_hf_cut_single_beamgas->Fill(cut);
            }

          if(noise)                                                 h_noise_tracks_hf->Fill(hf_single_energy_max,event->Tracks.size());
          if(beamgas)                                               h_beamgas_tracks_hf->Fill(hf_single_energy_max,event->Tracks.size());

          if((noise || beamgas) && hf_single_energy_max > 3)        h_run_events_single->Fill(lumiPerTime,evtWeight);
          if((noise || beamgas) && hf_double_energy_max > 1.5)      h_run_events_double->Fill(lumiPerTime,evtWeight);
          if((noise || beamgas) )                                   h_run_events->Fill(lumiPerTime,evtWeight);
                                                                        

        }

      //******************************************AFTER EVENT LOOP*******************************************
      double n_total = double(tree->GetEntries());
      cout << n_total << " events done" << endl;
    }

  //********************************************AFTER FILE LOOP************************************************

  out_file->Write();
  out_file->Save();

  return 0;
}
