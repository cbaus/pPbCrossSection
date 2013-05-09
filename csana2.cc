#define MAXEVT -50000

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
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210614/*_*.root"); sample_name.push_back("data210614"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210837/*_*.root"); sample_name.push_back("data210837"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210885/*_*.root"); sample_name.push_back("data210885"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data211607/*_*.root"); sample_name.push_back("data211607"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("Epos"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos_SL/*.root"); sample_name.push_back("Epos_SL"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Hijing/*.root"); sample_name.push_back("Hijing"); sample_type.push_back(MC);
  sample_fname.push_back("/afs/cern.ch/work/c/cbaus/private/trees/treeMC.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightDPMjet_v2/treeMC.root"); sample_name.push_back("Starlight_DPMJet");  sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlighDPMjet/treeMC.root"); sample_name.push_back("Starlight_DPMJet_Old");  sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightPythia/treeMC.root"); sample_name.push_back("Starlight_Pythia");  sample_type.push_back(MC);

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
  TH1D* h_mc_rapidity;
  TH1D* h_mc_eta_e;
  TH1D* h_mc_eff;

  TH2D* h_hf_fsc_p;
  TH2D* h_hf_fsc_m;
  
  TH2D* h_random_trig_tracks_hf;

  TH1D* h_hf_cut_single;
  TH1D* h_hf_cut_single_noise;
  TH1D* h_hf_cut_single_beamgas;
  TH1D* h_hf_cut_double;
  TH1D* h_hf_cut_double_noise;
  TH1D* h_hf_cut_double_beamgas;

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

  TH1D* h_mc_diffraction_single;
  TH1D* h_mc_diffraction_double;
  TH1D* h_mc_diffraction_SD;
  TH1D* h_mc_diffraction_DD;
  TH1D* h_mc_diffraction_CD;
  TH1D* h_mc_diffraction_ND;
  TH1D* h_mc_diffraction_all;
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
      int bptx_p_m;
      int bptx_p_nm;
      int bptx_np_m;
      int bptx_quiet;
      tree->SetBranchAddress("L1_ZeroBias_algPrescale",&zero_bias_prescale_L1);
      tree->SetBranchAddress("HLT_PAZeroBias_v1_Prescl",&zero_bias_prescale_HLT);
      tree->SetBranchAddress("HLT_PAZeroBias_v1",&zero_bias);
      tree->SetBranchAddress("HLT_PAL1Tech53_MB_SingleTrack_v1",&min_bias);
      tree->SetBranchAddress("HLT_PARandom_v1",&random);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_DecisionBeforeMask",&bptx_p_m);
      tree->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_DecisionBeforeMask",&bptx_p_nm);
      tree->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_DecisionBeforeMask",&bptx_np_m);
      tree->SetBranchAddress("L1Tech_BPTX_quiet.v0_DecisionBeforeMask",&bptx_quiet);

      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];
      h_zero_count_zb_coll      = new TH1D((add + string("_h_zero_count_zb_coll")).c_str(),"",100,728,1728);
      h_zero_count_zb_no_coll   = new TH1D((add + string("_h_zero_count_zb_no_coll")).c_str(),"",100,728,1728);
      h_hf_hits_coll            = new TH1D((add + string("_h_hf_hits_coll")).c_str(),"",200,0,200);
      h_hf_hits_plus            = new TH1D((add + string("_h_hfp_hits_coll")).c_str(),"",200,0,200);
      h_hf_hits_minus           = new TH1D((add + string("_h_hfm_hits_coll")).c_str(),"",200,0,200);
      h_hf_hits_noise           = new TH1D((add + string("_h_hf_hits_noise")).c_str(),"",200,0,20);
      h_hf_hits_beam_gas        = new TH1D((add + string("_h_hf_hits_beam_gas")).c_str(),"",200,0,20);
      h_castor_hf_diff_3        = new TH1D((add + string("_h_castor_hf_diff_3")).c_str(),"",100,0,10000);
      h_castor_hf_diff_5        = new TH1D((add + string("_h_castor_hf_diff_5")).c_str(),"",100,0,10000);
      h_castor_hf_diff_10       = new TH1D((add + string("_h_castor_hf_diff_10")).c_str(),"",100,0,10000);
      h_castor_gap_hf           = new TH1D((add + string("_h_castor_gap_hf")).c_str(),"",100,0,50);
      h_castor_nogap_hf         = new TH1D((add + string("_h_castor_nogap_hf")).c_str(),"",100,0,50);

      h_hf_cut_single           = new TH1D((add + string("_h_hf_cut_single")).c_str(),"",51,-0.05,5.05);
      h_hf_cut_single_noise     = new TH1D((add + string("_h_hf_cut_single_noise")).c_str(),"",51,-0.05,5.05);
      h_hf_cut_single_beamgas   = new TH1D((add + string("_h_hf_cut_single_beamgas")).c_str(),"",51,-0.05,5.05);
      h_hf_cut_double           = new TH1D((add + string("_h_hf_cut_double")).c_str(),"",51,-0.05,5.05);
      h_hf_cut_double_noise     = new TH1D((add + string("_h_hf_cut_double_noise")).c_str(),"",51,-0.05,5.05);
      h_hf_cut_double_beamgas   = new TH1D((add + string("_h_hf_cut_double_beamgas")).c_str(),"",51,-0.05,5.05);


      h_hf_fsc_p                = new TH2D((add + string("_h_hf_fsc_p")).c_str(),"",200,0,200,200,0,80000);
      h_hf_fsc_m                = new TH2D((add + string("_h_hf_fsc_m")).c_str(),"",200,0,200,200,0,80000);

      h_random_trig_tracks_hf   = new TH2D((add + string("_h_random_trig_tracks_hf")).c_str(),"",200,0,20,20,0,20);

      h_hf_hits_coll_lumi       = new TProfile((add + string("_h_hf_hits_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_minus_lumi      = new TProfile((add + string("_h_hfp_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_plus_lumi       = new TProfile((add + string("_h_hfm_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_noise_lumi      = new TProfile((add + string("_h_hf_hits_noise_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_coll_lumi       = new TProfile((add + string("_h_hf_totE_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_minus_lumi      = new TProfile((add + string("_h_hfp_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_plus_lumi       = new TProfile((add + string("_h_hfm_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_noise_lumi      = new TProfile((add + string("_h_hf_totE_noise_lumi")).c_str(),"",2000,0,2000);
      h_lumi                    = new TH1D((add + string("_h_lumi")).c_str(),"",2000,0,2000);
      h_lumi_3GeV               = new TH1D((add + string("_h_lumi_3GeV")).c_str(),"",2000,0,2000);

      h_perf_hf_rechits_single_3gev        = new TH1D((add + string("_h_perf_hf_rechits_single_3gev")).c_str(),"",100,0,1000);
      h_perf_hf_rechits_double_1dot5gev    = new TH1D((add + string("_h_perf_hf_rechits_double_1dot5gev")).c_str(),"",100,0,1000);
      h_perf_hf_totE_single_3gev           = new TH1D((add + string("_h_perf_hf_totE_single_3gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_double_1dot5gev       = new TH1D((add + string("_h_perf_hf_totE_double_1dot5gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_eta_single_3gev       = new TH1D((add + string("_h_perf_hf_totE_eta_single_3gev")).c_str(),"",100,-5.2,5.2);
      h_perf_hf_totE_eta_double_1dot5gev   = new TH1D((add + string("_h_perf_hf_totE_eta_double_1dot5gev")).c_str(),"",100,-5.2,5.2);

      h_mc_diffraction_single = new TH1D((add + string("_h_mc_diffraction_single")).c_str(),"",100,-9,0);
      h_mc_diffraction_double = new TH1D((add + string("_h_mc_diffraction_double")).c_str(),"",100,-9,0);
      h_mc_diffraction_SD     = new TH1D((add + string("_h_mc_diffraction_SD")).c_str(),"",100,-9,0);
      h_mc_diffraction_DD     = new TH1D((add + string("_h_mc_diffraction_DD")).c_str(),"",100,-9,0);
      h_mc_diffraction_CD     = new TH1D((add + string("_h_mc_diffraction_CD")).c_str(),"",100,-9,0);
      h_mc_diffraction_ND     = new TH1D((add + string("_h_mc_diffraction_ND")).c_str(),"",100,-9,0);
      h_mc_diffraction_all    = new TH1D((add + string("_h_mc_diffraction_all")).c_str(),"",100,-9,0);

      h_mc_rapidity           = new TH1D((add + string("_h_mc_rapidity")).c_str(),"",100,-12,12);
      h_mc_eta_e         = new TH1D((add + string("_h_mc_eta_e")).c_str(),"",100,-12,12);
      h_mc_eff                = new TH1D((add + string("_h_mc_eff")).c_str(),"",9,-0.5,8.5);
      h_mc_eff->GetXaxis()->SetBinLabel(1,"All");
      h_mc_eff->GetXaxis()->SetBinLabel(2,"HF single > 3 GeV");
      h_mc_eff->GetXaxis()->SetBinLabel(3,"HF single > 4 GeV");
      h_mc_eff->GetXaxis()->SetBinLabel(4,"HF single > 5 GeV");
      h_mc_eff->GetXaxis()->SetBinLabel(5,"HF double > 1.5 GeV");
      h_mc_eff->GetXaxis()->SetBinLabel(6,"HF double > 3 GeV");
      h_mc_eff->GetXaxis()->SetBinLabel(7,"nTrack #geq 1");
      h_mc_eff->GetXaxis()->SetBinLabel(8,"CASTOR E_{tot} > 12.5 GeV");

      double n_total = double(tree->GetEntries());
      if(MAXEVT<n_total && MAXEVT>0)
        n_total = double(MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);
          //          if(event->runNb != 210885)
          //continue;

          bool coll=0, noise=0, beam_gas=0;

          coll          = zero_bias && bptx_p_m; //double beam
          beam_gas      = (bptx_np_m || bptx_p_nm); // only single beam
          noise         = bptx_quiet;// && !bptx_np_m && !bptx_p_nm; //not both and not single beam

          if(sample_type[sample] == MC)
            {
              beam_gas = 0;
              noise = 0;
              coll = 1;
              min_bias = 1;
            }

          if(!coll && !noise && !beam_gas && !min_bias) //not intersted
            continue;



          //---------------------------------------------CASTOR
          double sum_cas_e_em = 0;
          double sum_cas_e_had = 0;
          
          for (vector<RecHitCASTOR>::const_iterator it = event->CASTOR.begin(); it < event->CASTOR.end(); ++it)
            {//break;
              if(it->GetModuleId() < 2)
                sum_cas_e_em += it->Energy;
              else if(it->GetModuleId() < 5)
                sum_cas_e_had += it->Energy;
            }
          
          const double sum_cas_e = sum_cas_e_had + sum_cas_e_em;
          const bool castor_tag = sum_cas_e > 12.5;





          //---------------------------------------------HF
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
          bool hf_double_tag = (hf_m_energy_max > 1.5) && hf_p_energy_max > 1.5;

          //---------------------------------------------GEN Particles
          const double s = 5020*5020;
          double m_x=0, xi_x=0;
          double m_y=0, xi_y=0;
          double rapGap=0;
          bool SD = event->genProcessID == 103 || event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD && !DD && !CD;
          if(sample_type[sample] == MC)
            {
              map<double,GenParticle*> rapidityMassMap;
              if(event->GEN.size() == 0)
                {
                  cerr << endl << " Empty event... skipping. (" << iEvent<< ")" << endl;
                  continue;
                }
              for (vector<GenParticle>::iterator it = event->GEN.begin(); it < event->GEN.end(); ++it)
                {
                  //cout << distance(event->GEN.begin(),it) << " " << it->Id << " " << it->GetEnergy() << endl;
                  if (it->Status != 1)
                    continue;
                                          
                  if (it->Id > 1e9) //skip fragments
                    continue;
                  double Rapidity= it->GetRapidity();
                  h_mc_rapidity->Fill(Rapidity);
                  h_mc_eta_e->Fill(it->GetEta(),it->GetEnergy());
                  rapidityMassMap.insert(pair<double,GenParticle*>(it->GetRapidity(),&(*it)));
                }
              map<double,GenParticle*>::const_iterator thisIsIt; //start of m_x
              for (map<double,GenParticle*>::const_iterator it = rapidityMassMap.begin(); it != rapidityMassMap.end(); ++it)
                {
                  if (it ==  rapidityMassMap.begin())
                    continue;
                  map<double,GenParticle*>::const_iterator it1=--it;
                  ++it;
                  double thisRapGap = fabs(it->second->GetRapidity() - it1->second->GetRapidity());
                  //cout << it->second->GetRapidity() << " " <<  it1->second->GetRapidity();
                  if(!TMath::Finite(thisRapGap) || TMath::IsNaN(thisRapGap))
                    continue;
                    //cerr << "   !!! ! PARTICLES WITH INFINITE RAP DETECTED " << distance(event->GEN.begin(),it) << " " << it->GetRapidity() << " " << it->Id << endl;

                  if (thisRapGap > rapGap)
                    {
                      //cout << " blib";
                      rapGap = thisRapGap;
                      thisIsIt = it;
                    }
                  //cout << endl;
                }

              map<double,GenParticle*>::const_iterator it = rapidityMassMap.begin();
              TLorentzVector vecX(it->second->Px,it->second->Py,it->second->Pz,it->second->GetEnergy());
              ++it;
              for (; it != thisIsIt; ++it)
                {
                  TLorentzVector vec(it->second->Px,it->second->Py,it->second->Pz,it->second->GetEnergy());
                  vecX += vec;
                }
              m_x = vecX.M();

              map<double,GenParticle*>::const_iterator it2 = thisIsIt;
              TLorentzVector vecY(it2->second->Px,it2->second->Py,it2->second->Pz,it2->second->GetEnergy());
              ++it2;
              for (; it2 != rapidityMassMap.end(); ++it2)
                {
                  TLorentzVector vec(it2->second->Px,it2->second->Py,it2->second->Pz,it2->second->GetEnergy());
                  vecY += vec;
                }
              m_y = vecY.M();
              if (m_x < m_y)
                {
                  double help = m_y;
                  m_y = m_x;
                  m_x = help;
                }
              xi_x = m_x*m_x / s;
              xi_y = m_y*m_y / s;
              //cout << "----------- !! M_X=" << m_x << "    M_Y=" << m_y << endl;
            }

          //Counting events
          const int prescale      = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS     = event->instLuminosity * event->instLuminosityCorr * 1e6;
          const double evtWeight  = lumiPerLS?double(prescale) / lumiPerLS:0.;


          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //Booking run histograms

          ostringstream run_str;
          run_str << event->runNb;
          TH1D* h_run_events_single = NULL;
          TH1D* h_run_events_double = NULL;
          TH1D* h_run_events_lumi = NULL;
          h_run_events_single = (TH1D*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events_single_")).c_str()));
          h_run_events_double = (TH1D*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events_double_")).c_str()));
          h_run_events_lumi = (TProfile*)(out_file->Get((add+string("/")+add+run_str.str()+string("_h_run_events_lumi_")).c_str()));
          if(h_run_events_single == NULL)
            {
              h_run_events_single = new TH1D((add+run_str.str()+string("_h_run_events_single_")).c_str(),run_str.str().c_str(),2000,0,2000);
              h_run_events_double = new TH1D((add+run_str.str()+string("_h_run_events_double_")).c_str(),run_str.str().c_str(),2000,0,2000);
              h_run_events_lumi = new TProfile((add+run_str.str()+string("_h_run_events_lumi_")).c_str(),run_str.str().c_str(),2000,0,2000);
            }


          //---------------------------------------------Filling HISTOS
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

          if(coll)                                                  h_mc_eff->Fill(0);
          if(coll && hf_single_energy_max > 3)                      h_mc_eff->Fill(1);
          if(coll && hf_single_energy_max > 4)                      h_mc_eff->Fill(2);
          if(coll && hf_single_energy_max > 5)                      h_mc_eff->Fill(3);
          if(coll && hf_double_energy_max > 1.5)                    h_mc_eff->Fill(4);
          if(coll && hf_double_energy_max > 3)                      h_mc_eff->Fill(5);
          if(coll && event->Tracks.size()>=1)                       h_mc_eff->Fill(6);
          if(coll && castor_tag)                                    h_mc_eff->Fill(7);

          for (double cut=0; cut <= 5; cut+=0.1)
            {
              if(coll && hf_double_energy_max > cut)                h_hf_cut_double->Fill(cut);
              if((noise || beam_gas) && hf_double_energy_max > cut) h_hf_cut_double_noise->Fill(cut);
              if((beam_gas || beam_gas) && hf_double_energy_max > cut) h_hf_cut_double_beamgas->Fill(cut);
              if(coll && hf_single_energy_max > cut)                h_hf_cut_single->Fill(cut);
              if((noise || beam_gas) && hf_single_energy_max > cut) h_hf_cut_single_noise->Fill(cut);
              if((beam_gas || beam_gas) && hf_single_energy_max > cut) h_hf_cut_single_beamgas->Fill(cut);
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

          if(coll && hf_single_energy_max > 3)                      h_mc_diffraction_single->Fill(log10(xi_x));
          if(coll && hf_double_tag)                                 h_mc_diffraction_double->Fill(log10(xi_x));
          if(coll && SD)                                            h_mc_diffraction_SD->Fill(log10(xi_x));
          if(coll && DD)                                            h_mc_diffraction_DD->Fill(log10(xi_x));
          if(coll && CD)                                            h_mc_diffraction_CD->Fill(log10(xi_x));
          if(coll && ND)                                            h_mc_diffraction_ND->Fill(log10(xi_x));
          if(coll)                                                  h_mc_diffraction_all->Fill(log10(xi_x));

          if(coll && hf_single_energy_max > 3)                      h_run_events_single->Fill(event->lumiNb,evtWeight);
          if(coll && hf_double_energy_max > 1.5)                    h_run_events_double->Fill(event->lumiNb,evtWeight);
                                                                    h_run_events_lumi->Fill(event->lumiNb,lumiPerLS);

          if(noise)                                                 h_random_trig_tracks_hf->Fill(hf_single_energy_max,event->Tracks.size());

        }

      //******************************************AFTER EVENT LOOP*******************************************
      double integ_lumi = 1854.344875;// nb^-1
      double n_zb = h_mc_eff->GetBinContent(1) * n_total;
      double n_mb = h_mc_eff->GetBinContent(4) * n_total;

      cout << endl << "Cross section ZB: " << n_zb/integ_lumi << " --- cross section MB: " << n_mb/integ_lumi << endl;
      h_perf_hf_totE_eta_double_1dot5gev->Scale(1./n_total);
      h_perf_hf_totE_eta_single_3gev->Scale(1./n_total);

      h_mc_rapidity->Scale(1./n_total/h_mc_rapidity->GetBinWidth(1));
      h_mc_eta_e   ->Scale(1./n_total/h_mc_eta_e->GetBinWidth(1));

    }

  //********************************************AFTER SAMPLE LOOP************************************************

  out_file->Write();
  out_file->Save();

  return 0;
}
