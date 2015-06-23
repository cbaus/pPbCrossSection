#define MAXEVT -10000

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

#include "CastorCorrFactorpp2015.h"


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

void BinLogX(TH1* h);

int main()
{
  TH1::SetDefaultSumw2();
  //style();

  //*************************************************************INPUT***********************************************************
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/data_247324.root"); sample_name.push_back("data247324"); sample_type.push_back(DATA);

  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos_noise.root","RECREATE");


  TH1D* h_hf_hits_noise;
  TH1D* h_hf_hits_beamgas;
  TH1D* h_castor_sumE_noise;
  TH1D* h_castor_sumE_beam_plus;
  TH1D* h_castor_sumE_beam_minus;

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

      tree->SetCacheSize(50000000);
      tree->AddBranchToCache("*");

      if (tree->GetEntries() == 0) {
        cout << "No events found in file \"" << sample_fname[sample] << "\"" << endl;
        return 0;
      }

      //________________________________________________________________________________

      tree->SetBranchStatus("*", 1);
      tree->SetBranchStatus("EB*",0);
      tree->SetBranchStatus("EE*",0);
      tree->SetBranchStatus("HBHE*",0);
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
      int random;
      int random_prescale_HLT;
      int bptx_p_m;
      int bptx_p;
      int bptx_m;
      int bptx_quiet;
      tree->SetBranchAddress("L1_ZeroBias_algPrescale",&zero_bias_prescale_L1);
      tree->SetBranchAddress("HLT_ZeroBias_part0_v1_Prescl",&zero_bias_prescale_HLT);
      tree->SetBranchAddress("HLT_ZeroBias_part0_v1",&zero_bias);
      tree->SetBranchAddress("HLT_Random_v2",&random);
      tree->SetBranchAddress("HLT_Random_v2_Prescl",&random_prescale_HLT);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",&bptx_p_m);
      tree->SetBranchAddress("L1Tech_BPTX_minus.v0_DecisionBeforeMask",&bptx_m);
      tree->SetBranchAddress("L1Tech_BPTX_plus.v0_DecisionBeforeMask",&bptx_p);
      tree->SetBranchAddress("L1Tech_BPTX_quiet.v0_DecisionBeforeMask",&bptx_quiet);

      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];
      h_hf_hits_noise           = new TH1D((add + string("_h_hf_hits_noise")).c_str(),"",50,log10(0.5),log10(50));
      h_hf_hits_beamgas        = new TH1D((add + string("_h_hf_hits_beamgas")).c_str(),"",50,log10(0.5),log10(50));
      BinLogX(h_hf_hits_noise);
      BinLogX(h_hf_hits_beamgas);

      h_castor_sumE_noise           = new TH1D((add + string("_h_castor_sumE_noise")).c_str(),"",50,log10(0.01),log10(100));
      h_castor_sumE_beam_minus        = new TH1D((add + string("_h_castor_sumE_beam_minus")).c_str(),"",50,log10(0.01),log10(100));
      h_castor_sumE_beam_plus        = new TH1D((add + string("_h_castor_sumE_beam_plus")).c_str(),"",50,log10(0.01),log10(100));
      BinLogX(h_castor_sumE_noise);
      BinLogX(h_castor_sumE_beam_minus);
      BinLogX(h_castor_sumE_beam_plus);

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

      for(int iEvent=0; iEvent<tree->GetEntries(); iEvent++)
        {
          if(iEvent==MAXEVT) break;
          if(iEvent % 1000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << (MAXEVT>0?MAXEVT:tree->GetEntries()) << endl;
          tree->GetEntry(iEvent);
          //          if(event->runNb != 210885)
          //continue;

          bool noise=0, beamgas=0;

          beamgas      = (bptx_m != bptx_p); // only single beam
          noise         = !bptx_p && !bptx_m; //not both and not single beam

          if(!noise && !beamgas) //not intersted
            continue;


          const int prescale      = random_prescale_HLT;
          const double lumiPerLS     = event->instLuminosity * event->instLuminosityCorr * 1e6;
          const double lumiPerTime     = lumiPerLS / 23.31;
          const double evtWeight  = lumiPerTime?double(prescale) / lumiPerTime:0.;

          //---------------------------------------------CASTOR
          double sum_CAS_E_em = 0;
          double sum_CAS_E_had = 0;
          vector<double> sum_CAS_E_mod(16,0.);

          for (vector<RecHitCASTOR>::const_iterator it = event->CASTOR.begin(); it < event->CASTOR.end(); ++it)
            {//break;
              double RecHitGeV =  it->Energy;
              const int sec = it->GetSectorId();
              const int mod = it->GetModuleId();

              if(sample_type[sample] == DATA)
                {
                  //RecHitGeV -= castor::channelPedMean[sec][mod];
                  #warning "remove for 3.8 T"
                  RecHitGeV *= castor::corr0Tto38T[sec-1][mod-1];
                  RecHitGeV *= castor::channelGainQE[sec][mod];
                  RecHitGeV *= castor::absCasEscaleFactor;
                }
              RecHitGeV *= (int)castor::channelQuality[sec][mod];

              sum_CAS_E_mod[it->GetModuleId()] += RecHitGeV;

              if(it->GetModuleId() < 2)
                sum_CAS_E_em += RecHitGeV;
              else if(it->GetModuleId() < 5)
                sum_CAS_E_had += RecHitGeV;
            }

          const double sum_CAS_E = sum_CAS_E_had + sum_CAS_E_em;
          const bool castor_tag = sum_CAS_E > 5.6; //esstimated by sebastian for data only!


          //---------------------------------------------HF
          int hf_n = event->HFtowers.size();
          int hf_zero_count = ForwardRecord::nMaxHFMRecHits - hf_n;
          double hf_double_energy_max = 0; // if no tower it is 0 for cut later
          double hf_single_energy_max = 0;
          double hf_m_energy_max = 0;
          double hf_p_energy_max = 0;
          double hf_p_energy = 0;
          double hf_m_energy = 0;
          map<int,double> hf_ring_map_single_max;
          //Tower loop
          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              const int Ieta = it->Eta > 0?it->IetaAbs:-it->IetaAbs;
              //make new entry or if exists and new is larger overwrite
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
          if(noise)                                                h_hf_hits_noise->Fill(hf_single_energy_max);
          if(beamgas)                                              h_hf_hits_beamgas->Fill(hf_single_energy_max);
          if(noise)                                                h_castor_sumE_noise->Fill(sum_CAS_E);
          if(beamgas && bptx_m)                                    h_castor_sumE_beam_minus->Fill(sum_CAS_E);
          if(beamgas && bptx_p)                                    h_castor_sumE_beam_plus->Fill(sum_CAS_E);

          for (double cut=0; cut <= 10; cut+=0.1)
            {
              if((noise) && hf_double_energy_max >= cut) h_hf_cut_double_noise->Fill(cut);
              if((beamgas) && hf_double_energy_max >= cut) h_hf_cut_double_beamgas->Fill(cut);
              if((noise) && hf_single_energy_max >= cut) h_hf_cut_single_noise->Fill(cut);
              if((beamgas) && hf_single_energy_max >= cut) h_hf_cut_single_beamgas->Fill(cut);
            }

          if(noise)                                                 h_noise_tracks_hf->Fill(hf_single_energy_max,event->Tracks.size());
          if(beamgas)                                               h_beamgas_tracks_hf->Fill(hf_single_energy_max,event->Tracks.size());

          if((noise) && hf_single_energy_max > 4)        h_run_events_single->Fill(lumiPerTime,evtWeight);
          if((noise) && hf_double_energy_max > 3)      h_run_events_double->Fill(lumiPerTime,evtWeight);
          if((noise) )                                   h_run_events->Fill(lumiPerTime,evtWeight);


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

void BinLogX(TH1* h)
{

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);

  }
  axis->Set(bins, new_bins);
  delete new_bins;
}
