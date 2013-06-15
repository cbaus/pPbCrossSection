#define _MAXEVT 50000
#define _SkipHFRings 1

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

enum e_type
  {
    MC = 0,
    DATA
  };

using namespace std;

int IetaToRing(int ring);
void BinLogX(TH1* h);


int main()
{
  TH1::SetDefaultSumw2();

  vector<string> sample_fname;
  vector<string> sample_name;
  vector<e_type> sample_type;

  const int array_size = 9;
  double array_cuts_single[array_size] = {5.4, 6, 6.6, 7,   8, 9,   10, 15,  20};
  double array_cuts_double[array_size] = {1.8, 2, 2.2, 2.5, 3, 3.5, 4,  4.5, 5,};
  vector<double> cut_energies_single(array_cuts_single,array_cuts_single+array_size);
  vector<double> cut_energies_double(array_cuts_double,array_cuts_double+array_size);

  //*************************************************************INPUT***********************************************************
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210614/*_*.root"); sample_name.push_back("data210614"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210885/*_*.root"); sample_name.push_back("data210885"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data210998/*_*.root"); sample_name.push_back("data210998"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211000/*.root"); sample_name.push_back("data211000"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211001/*.root"); sample_name.push_back("data211001"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211032/*.root"); sample_name.push_back("data211032"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211256/*.root"); sample_name.push_back("data211256"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211371/*.root"); sample_name.push_back("data211371"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211390/*.root"); sample_name.push_back("data211390"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211460/*.root"); sample_name.push_back("data211460"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211532/*.root"); sample_name.push_back("data211532"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211538/*.root"); sample_name.push_back("data211538"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data211607/*_*.root"); sample_name.push_back("data211607"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("Epos"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("EposSDWeight2"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos_SL/*.root"); sample_name.push_back("Epos_SL"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Hijing/*.root"); sample_name.push_back("Hijing"); sample_type.push_back(MC);
  // sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightDPMjet_v2/treeMC.root"); sample_name.push_back("Starlight_DPMJet");  sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightPythia/treeMC.root"); sample_name.push_back("Starlight_Pythia");  sample_type.push_back(MC);

  TFile f("plots/hf_cuts_noise.root");
  TVectorD* hf_m_cuts_light_noise = NULL;
  TVectorD* hf_p_cuts_light_noise = NULL;
  TVectorD* hf_m_cuts_medium_noise = NULL;
  TVectorD* hf_p_cuts_medium_noise = NULL;
  TVectorD* hf_m_cuts_heavy_noise = NULL;
  TVectorD* hf_p_cuts_heavy_noise = NULL;
  hf_m_cuts_light_noise=(TVectorD*)f.Get("hf_m_cuts_light_noise");
  hf_p_cuts_light_noise=(TVectorD*)f.Get("hf_p_cuts_light_noise");
  hf_m_cuts_medium_noise=(TVectorD*)f.Get("hf_m_cuts_medium_noise");
  hf_p_cuts_medium_noise=(TVectorD*)f.Get("hf_p_cuts_medium_noise");
  hf_m_cuts_heavy_noise=(TVectorD*)f.Get("hf_m_cuts_heavy_noise");
  hf_p_cuts_heavy_noise=(TVectorD*)f.Get("hf_p_cuts_heavy_noise");
  if(!hf_m_cuts_light_noise || !hf_p_cuts_light_noise || !hf_m_cuts_medium_noise || !hf_p_cuts_medium_noise || !hf_m_cuts_heavy_noise || !hf_p_cuts_heavy_noise) {cerr << "error" << endl; return 1;}
  cout << "Heavy cuts:" << endl;
  hf_p_cuts_heavy_noise->Print();

  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos.root","RECREATE");

  TH1D* h_zero_count_zb_coll;
  TH1D* h_zero_count_zb_no_coll;
  TH1D* h_hf_hits_coll_single;
  TH1D* h_hf_hits_coll_double;
  TH1D* h_hf_hits_noise_single;
  TH1D* h_hf_hits_noise_double;
  TH1D* h_hf_hits_noise_plus;
  TH1D* h_hf_hits_noise_minus;
  TH1D* h_hf_hits_plus;
  TH1D* h_hf_hits_minus;
  TH1D* h_castor_hf_diff_3;
  TH1D* h_castor_hf_diff_5;
  TH1D* h_castor_hf_diff_10;
  TH1D* h_castor_gap_hf;
  TH1D* h_castor_nogap_hf;
  TH1D* h_mc_rapidity;
  TH1D* h_mc_eta_e;
  TH1D* h_mc_eff;
  TH1D* h_mc_effgain_single;
  TH1D* h_mc_effgain_double;

  TH2D* h_random_trig_tracks_hf;

  vector<TH1D*> v_run_events_single(array_size,NULL);
  vector<TH1D*> v_run_events_double(array_size,NULL);

  TH1D* h_run_events_lumi;

  TH1D* h_hf_cut_single;
  TH1D* h_hf_cut_single_noise;
  TH1D* h_hf_cut_double;
  TH1D* h_hf_cut_double_noise;
  TH1D* h_hf_new_cut_single;
  TH1D* h_hf_new_cut_double;

  TH1D* h_hf_noise_all_lumi;
  TH1D* h_hf_noise_selected_single_lumi;
  TH1D* h_hf_noise_selected_double_lumi;

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
  TH1D* h_mc_diffraction_SD1;
  TH1D* h_mc_diffraction_SD2;
  TH1D* h_mc_diffraction_DD;
  TH1D* h_mc_diffraction_CD;
  TH1D* h_mc_diffraction_ND;
  TH1D* h_mc_diffraction_all;
  TH1D* h_mc_lrg_single;
  TH1D* h_mc_lrg_double;
  TH1D* h_mc_lrg_SD1;
  TH1D* h_mc_lrg_SD2;
  TH1D* h_mc_lrg_DD;
  TH1D* h_mc_lrg_CD;
  TH1D* h_mc_lrg_ND;
  TH1D* h_mc_lrg_all;
  TH2D* h_mc_lrg_xi;
  TH2D* h_mc_lrg_xiy;
  TH2D* h_mc_xix_xiy;
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

//       tree->SetBranchStatus("*", 1);
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
      if(sample_type[sample] == DATA)
        {
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
        }
      //________________________________________________________________________________

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());
      string add = sample_name[sample];

      ostringstream run_str;
      run_str << event->runNb;

      h_zero_count_zb_coll        = new TH1D((add + string("_h_zero_count_zb_coll")).c_str(),"",100,728,1728);
      h_zero_count_zb_no_coll     = new TH1D((add + string("_h_zero_count_zb_no_coll")).c_str(),"",100,728,1728);
      h_hf_hits_coll_single       = new TH1D((add + string("_h_hf_hits_coll_single")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_coll_double       = new TH1D((add + string("_h_hf_hits_coll_double")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_noise_single      = new TH1D((add + string("_h_hf_hits_noise_single")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_noise_double      = new TH1D((add + string("_h_hf_hits_noise_double")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_noise_plus      = new TH1D((add + string("_h_hfp_hits_noise")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_noise_minus      = new TH1D((add + string("_h_hfm_hits_noise")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_plus              = new TH1D((add + string("_h_hfp_hits_coll")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_minus             = new TH1D((add + string("_h_hfm_hits_coll")).c_str(),"",100,log10(0.5),log10(500.));

      BinLogX(h_hf_hits_coll_single);
      BinLogX(h_hf_hits_coll_double);
      BinLogX(h_hf_hits_noise_single);
      BinLogX(h_hf_hits_noise_double);
      BinLogX(h_hf_hits_noise_plus);
      BinLogX(h_hf_hits_noise_minus);
      BinLogX(h_hf_hits_plus);
      BinLogX(h_hf_hits_minus);

      h_castor_hf_diff_3          = new TH1D((add + string("_h_castor_hf_diff_3")).c_str(),"",100,0,10000);
      h_castor_hf_diff_5          = new TH1D((add + string("_h_castor_hf_diff_5")).c_str(),"",100,0,10000);
      h_castor_hf_diff_10         = new TH1D((add + string("_h_castor_hf_diff_10")).c_str(),"",100,0,10000);
      h_castor_gap_hf             = new TH1D((add + string("_h_castor_gap_hf")).c_str(),"",100,0,50);
      h_castor_nogap_hf           = new TH1D((add + string("_h_castor_nogap_hf")).c_str(),"",100,0,50);
      if(sample_type[sample] == DATA)
        {
          for (int i = 0; i<array_size; i++)
            {
              ostringstream ss_single; ss_single << add <<"_h_run_events_single_" << cut_energies_single[i];
              ostringstream ss_double; ss_double << add <<"_h_run_events_double_" << cut_energies_double[i];
              v_run_events_single[i] = new TH1D(ss_single.str().c_str(),run_str.str().c_str(),2000,0,2000);
              v_run_events_double[i] = new TH1D(ss_double.str().c_str(),run_str.str().c_str(),2000,0,2000);
            }

          h_run_events_lumi = new TProfile((add+string("_h_run_events_lumi")).c_str(),run_str.str().c_str(),2000,0,2000,"s");


          h_random_trig_tracks_hf     = new TH2D((add + string("_h_random_trig_tracks_hf")).c_str(),"",200,0,20,20,0,20);

          h_hf_noise_all_lumi              = new TH1D((add + string("_h_hf_noise_all_lumi")).c_str(),"",2000,0,2000);
          h_hf_noise_selected_single_lumi  = new TH1D((add + string("_h_hf_noise_selected_single_lumi")).c_str(),"",2000,0,2000);
          h_hf_noise_selected_double_lumi  = new TH1D((add + string("_h_hf_noise_selected_double_lumi")).c_str(),"",2000,0,2000);
        }

      h_hf_cut_single             = new TH1D((add + string("_h_hf_cut_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_single_noise       = new TH1D((add + string("_h_hf_cut_single_noise")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_double             = new TH1D((add + string("_h_hf_cut_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_double_noise       = new TH1D((add + string("_h_hf_cut_double_noise")).c_str(),"",201,-0.05,20.05);

      h_hf_new_cut_single         = new TH1D((add + string("_h_hf_new_cut_single")).c_str(),"",4,-0.5,3.5);
      h_hf_new_cut_double         = new TH1D((add + string("_h_hf_new_cut_double")).c_str(),"",4,-0.5,3.5);


      h_hf_hits_coll_lumi         = new TProfile((add + string("_h_hf_hits_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_minus_lumi        = new TProfile((add + string("_h_hfp_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_plus_lumi         = new TProfile((add + string("_h_hfm_hits_lumi")).c_str(),"",2000,0,2000);
      h_hf_hits_noise_lumi        = new TProfile((add + string("_h_hf_hits_noise_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_coll_lumi         = new TProfile((add + string("_h_hf_totE_coll_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_minus_lumi        = new TProfile((add + string("_h_hfp_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_plus_lumi         = new TProfile((add + string("_h_hfm_totE_lumi")).c_str(),"",2000,0,2000);
      h_hf_totE_noise_lumi        = new TProfile((add + string("_h_hf_totE_noise_lumi")).c_str(),"",2000,0,2000);
      h_lumi                      = new TH1D((add + string("_h_lumi")).c_str(),"",2000,0,2000);
      h_lumi_3GeV                 = new TH1D((add + string("_h_lumi_3GeV")).c_str(),"",2000,0,2000);

      h_perf_hf_rechits_single_3gev        = new TH1D((add + string("_h_perf_hf_rechits_single_3gev")).c_str(),"",100,0,1000);
      h_perf_hf_rechits_double_1dot5gev    = new TH1D((add + string("_h_perf_hf_rechits_double_1dot5gev")).c_str(),"",100,0,1000);
      h_perf_hf_totE_single_3gev           = new TH1D((add + string("_h_perf_hf_totE_single_3gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_double_1dot5gev       = new TH1D((add + string("_h_perf_hf_totE_double_1dot5gev")).c_str(),"",500,0,10000);
      h_perf_hf_totE_eta_single_3gev       = new TH1D((add + string("_h_perf_hf_totE_eta_single_3gev")).c_str(),"",100,-5.2,5.2);
      h_perf_hf_totE_eta_double_1dot5gev   = new TH1D((add + string("_h_perf_hf_totE_eta_double_1dot5gev")).c_str(),"",100,-5.2,5.2);
      if(sample_type[sample] == MC)
        {
          h_mc_diffraction_single = new TH1D((add + string("_h_mc_diffraction_single")).c_str(),"",100,-9,2);
          h_mc_diffraction_double = new TH1D((add + string("_h_mc_diffraction_double")).c_str(),"",100,-9,2);
          h_mc_diffraction_SD1    = new TH1D((add + string("_h_mc_diffraction_SD1")).c_str(),"",100,-9,2);
          h_mc_diffraction_SD2    = new TH1D((add + string("_h_mc_diffraction_SD2")).c_str(),"",100,-9,2);
          h_mc_diffraction_DD     = new TH1D((add + string("_h_mc_diffraction_DD")).c_str(),"",100,-9,2);
          h_mc_diffraction_CD     = new TH1D((add + string("_h_mc_diffraction_CD")).c_str(),"",100,-9,2);
          h_mc_diffraction_ND     = new TH1D((add + string("_h_mc_diffraction_ND")).c_str(),"",100,-9,2);
          h_mc_diffraction_all    = new TH1D((add + string("_h_mc_diffraction_all")).c_str(),"",100,-9,2);
          h_mc_lrg_single         = new TH1D((add + string("_h_mc_lrg_single")).c_str(),"",500,-1,20);
          h_mc_lrg_double         = new TH1D((add + string("_h_mc_lrg_double")).c_str(),"",500,-1,20);
          h_mc_lrg_SD1            = new TH1D((add + string("_h_mc_lrg_SD1")).c_str(),"",500,-1,20);
          h_mc_lrg_SD2            = new TH1D((add + string("_h_mc_lrg_SD2")).c_str(),"",500,-1,20);
          h_mc_lrg_CD             = new TH1D((add + string("_h_mc_lrg_CD")).c_str(),"",500,-1,20);
          h_mc_lrg_DD             = new TH1D((add + string("_h_mc_lrg_DD")).c_str(),"",500,-1,20);
          h_mc_lrg_ND             = new TH1D((add + string("_h_mc_lrg_ND")).c_str(),"",500,-1,20);
          h_mc_lrg_all            = new TH1D((add + string("_h_mc_lrg_all")).c_str(),"",500,-1,20);
          h_mc_lrg_xi             = new TH2D((add + string("_h_mc_lrg_xi")).c_str(),"",200,-9.5,2,500,-1,20);
          h_mc_lrg_xiy             = new TH2D((add + string("_h_mc_lrg_xiy")).c_str(),"",200,-9.5,2,500,-1,20);
          h_mc_xix_xiy             = new TH2D((add + string("_h_mc_xix_xiy")).c_str(),"",200,-9.5,2,200,-9.5,2);

          h_mc_rapidity           = new TH1D((add + string("_h_mc_rapidity")).c_str(),"",100,-12,12);
          h_mc_eta_e              = new TH1D((add + string("_h_mc_eta_e")).c_str(),"",100,-12,12);
          h_mc_eff                = new TH1D((add + string("_h_mc_eff")).c_str(),"",8,-0.5,7.5);
          h_mc_eff->GetXaxis()->SetBinLabel(1,"All");
          h_mc_eff->GetXaxis()->SetBinLabel(2,"HF single > 3 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(3,"HF single > 4.5 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(4,"HF single > 6 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(5,"HF double > 2 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(6,"HF double > 3 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(7,"no. of tracks #geq 1");
          h_mc_eff->GetXaxis()->SetBinLabel(8,"CASTOR E_{tot} > 12.5 GeV");
          h_mc_effgain_single     = new TH1D((add + string("_h_mc_effgain_single")).c_str(),"",7,-0.5,6.5);
          h_mc_effgain_single->GetXaxis()->SetBinLabel(1,"All");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(2,"HF single > 3 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(3,"HF single > 4.5 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(4,"HF double > 2 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(5,"HF double > 3 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(6,"no. of tracks #geq 1");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(7,"CASTOR E_{tot} > 12.5 GeV");
          h_mc_effgain_double     = new TH1D((add + string("_h_mc_effgain_double")).c_str(),"",6,-0.5,5.5);
          h_mc_effgain_double->GetXaxis()->SetBinLabel(1,"All");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(2,"HF single > 3 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(3,"HF single > 4.5 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(4,"HF single > 6 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(5,"no. of tracks  #geq 1");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(6,"CASTOR E_{tot} > 12.5 GeV");
        }

      double n_total = double(tree->GetEntries());
      if(_MAXEVT<n_total && _MAXEVT>0)
        n_total = double(_MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);
          //          if(event->runNb != 210885)
          //continue;

          bool coll=0, noise=0;

          if(sample_type[sample] == DATA)
            {
              coll          = zero_bias && bptx_p_m; //double beam
              noise         = random && bptx_quiet;// && !bptx_np_m && !bptx_p_nm; //not both and not single beam
            }
          else if(sample_type[sample] == MC)
            {
              noise = 0;
              coll = 1;
              min_bias = 1;
            }

          if(!coll && !noise && !min_bias) //not intersted
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
          double hf_double_energy_max = 0;
          double hf_single_energy_max = 0;
          double hf_m_energy_max = 0;
          double hf_p_energy_max = 0;
          double hf_p_energy = 0;
          double hf_m_energy = 0;
          bool hf_m_light_cut = false;
          bool hf_m_medium_cut = false;
          bool hf_m_heavy_cut = false;
          bool hf_p_light_cut = false;
          bool hf_p_medium_cut = false;
          bool hf_p_heavy_cut = false;
          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              if(_SkipHFRings && it->IetaAbs == 41)
                continue;
              const int Ieta = it->Eta > 0?it->IetaAbs:-it->IetaAbs;
              //cout << it->Eta << " " << it->Energy << endl;
              if(it->Eta > 0. && it->Energy > hf_p_energy_max)
                hf_p_energy_max = it->Energy;
              if(it->Eta <= 0. && it->Energy > hf_m_energy_max)
                hf_m_energy_max = it->Energy;

              if(it->Eta > 0.)
                hf_p_energy += it->Energy;
              else
                hf_m_energy += it->Energy;

              if(Ieta < 0 && it->Energy >= (*hf_m_cuts_light_noise)[IetaToRing(Ieta)])  hf_m_light_cut = true;
              if(Ieta < 0 && it->Energy >= (*hf_m_cuts_medium_noise)[IetaToRing(Ieta)]) hf_m_medium_cut = true;
              if(Ieta < 0 && it->Energy >= (*hf_m_cuts_heavy_noise)[IetaToRing(Ieta)])  hf_m_heavy_cut = true;
              if(Ieta > 0 && it->Energy >= (*hf_p_cuts_light_noise)[IetaToRing(Ieta)])  hf_p_light_cut = true;
              if(Ieta > 0 && it->Energy >= (*hf_p_cuts_medium_noise)[IetaToRing(Ieta)]) hf_p_medium_cut = true;
              if(Ieta > 0 && it->Energy >= (*hf_p_cuts_heavy_noise)[IetaToRing(Ieta)])  hf_p_heavy_cut = true;
            }
          double hf_pm_energy = hf_p_energy + hf_m_energy;
          //cout << "rechit: " << hf_pm_energy << endl;

          hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
          hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
          bool hf_single_tag = hf_single_energy_max >= 6;
          bool hf_double_tag = hf_double_energy_max >= 2;

          bool hf_light_double_tag = hf_p_light_cut && hf_m_light_cut;
          bool hf_medium_double_tag = hf_p_medium_cut && hf_m_medium_cut;
          bool hf_heavy_double_tag = hf_p_heavy_cut && hf_m_heavy_cut;
          bool hf_light_single_tag = hf_p_light_cut || hf_m_light_cut;
          bool hf_medium_single_tag = hf_p_medium_cut || hf_m_medium_cut;
          bool hf_heavy_single_tag = hf_p_heavy_cut || hf_m_heavy_cut;




          //---------------------------------------------GEN Particles
          const double s = 5020*5020;
          double m_m=0, m_x=0, xi_x=0;
          double m_p=0, m_y=0, xi_y=0;
          double rapGap=-1;
          bool SD1 = event->genProcessID == 103; //lead dissociates
          bool SD2 = event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD1 && !SD2 && !DD && !CD;
          if(sample_type[sample] == MC)
            {
              multimap<double,GenParticle*> rapidityMassMap;
              if(event->GEN.size() == 0)
                {
                  cerr << endl << " Empty event... skipping. (" << iEvent<< ")" << endl;
                  continue;
                }

              for (vector<GenParticle>::iterator it = event->GEN.begin(); it < event->GEN.end(); ++it)
                {
                  if (it->Status != 1)
                    continue;

                  if (it->Id > 1e9) //skip fragments
                    continue;

                  if(abs(it->Id) == 12 || abs(it->Id) == 14 || abs(it->Id) == 16 || abs(it->Id) == 13)
                    continue; // skip muon + neutrinos

                  double Rapidity= it->GetRapidity();
                  h_mc_rapidity->Fill(Rapidity);
                  h_mc_eta_e->Fill(it->GetEta(),it->GetEnergy());
                  rapidityMassMap.insert(pair<double,GenParticle*>(it->GetRapidity(),&(*it)));
                }
              assert(rapidityMassMap.size() >= 2); //Assume 2 or more final usable particles

              multimap<double,GenParticle*>::const_iterator thisIsIt; //start of m_x
              int n=0;
              for (multimap<double,GenParticle*>::const_iterator it = rapidityMassMap.begin(); it != rapidityMassMap.end(); ++it)
                {
                  //if(SD2) cout << n++ << " " << it->second->Id << ": E=" << it->second->GetEnergy() << "  --  eta=" << it->second->GetEta() << "  --  y=" << it->second->GetRapidity() << endl;
                  if (it == rapidityMassMap.begin())
                    continue;
                  multimap<double,GenParticle*>::const_iterator it1=--it;
                  ++it;
                  double thisRapGap = fabs(it->second->GetRapidity() - it1->second->GetRapidity());
                  //if(SD2) cout << it->second->GetRapidity() << " " << it1->second->GetRapidity();
                  if(!TMath::Finite(thisRapGap) || TMath::IsNaN(thisRapGap))
                    continue;
                  //cerr << "   !!! ! PARTICLES WITH INFINITE RAP DETECTED " << distance(event->GEN.begin(),it) << " " << it->GetRapidity() << " " << it->Id << endl;

                  if (thisRapGap > rapGap)
                    {
                      //if(SD2) cout << " blib";
                      rapGap = thisRapGap;
                      thisIsIt = it;
                    }
                  //if(SD2) cout << endl;
                }
              TLorentzVector vecM(0,0,0,0);
              multimap<double,GenParticle*>::const_iterator it = rapidityMassMap.begin();
              for (; it != thisIsIt; ++it)
                {
                  if(it->second->GetEta() > 1e9) //skip fragments only for mass calculation
                    continue;
                  TLorentzVector vec(it->second->Px,it->second->Py,it->second->Pz,it->second->GetEnergy());
                  vecM += vec;
                }
              m_m = vecM.M();

              TLorentzVector vecP(0,0,0,0);
              multimap<double,GenParticle*>::const_iterator it2 = thisIsIt;
              for (; it2 != rapidityMassMap.end(); ++it2)
                {
                  if(it->second->GetEta() > 1e9) //skip fragments only for mass calculation
                    continue;
                  TLorentzVector vec(it2->second->Px,it2->second->Py,it2->second->Pz,it2->second->GetEnergy());
                  vecP += vec;
                }
              m_p = vecP.M();
              //if (SD2 && rapGap< 2) cout << "HERE HERE HERE " << m_x << " " << m_y << endl; //WARNING

              m_x = TMath::Max(m_m,m_p);
              m_y = TMath::Min(m_m,m_p);

              xi_x = m_x*m_x / s;
              xi_y = m_y*m_y / s;
              //cout << "----------- !! M_X=" << m_x << "    M_Y=" << m_y << endl;
            }
          //Counting events
          const int prescale        = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS    = event->instLuminosity * event->instLuminosityCorr * 1e6;
          double evtWeight    = 1;
          if(sample_type[sample]==DATA)
            evtWeight = double(prescale);
          else if(sample_name[sample]=="EposSDWeight2" && (SD1 || SD2))
            evtWeight *= 2.;
          const double noiseWeight = random_prescale_HLT;


          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //---------------------------------------------Filling HISTOS
          if(coll)                                                  h_zero_count_zb_coll->Fill(hf_zero_count,evtWeight);
          if(noise)                                                 h_zero_count_zb_no_coll->Fill(hf_zero_count,noiseWeight);

          if(coll)                                                  h_hf_hits_coll_single->Fill(hf_single_energy_max);
          if(coll)                                                  h_hf_hits_coll_double->Fill(hf_double_energy_max);
          if(noise)                                                 h_hf_hits_noise_single->Fill(hf_single_energy_max);
          if(noise)                                                 h_hf_hits_noise_double->Fill(hf_double_energy_max);
          if(noise)                                                 h_hf_hits_noise_plus->Fill(hf_p_energy_max);
          if(noise)                                                 h_hf_hits_noise_minus->Fill(hf_m_energy_max);
          if(coll)                                                  h_hf_hits_plus->Fill(hf_p_energy_max);
          if(coll)                                                  h_hf_hits_minus->Fill(hf_m_energy_max);

          if(coll && hf_double_energy_max < 3)                      h_castor_hf_diff_3->Fill(sum_cas_e,evtWeight);
          if(coll && hf_double_energy_max < 5)                      h_castor_hf_diff_5->Fill(sum_cas_e,evtWeight);
          if(coll && hf_double_energy_max < 10)                     h_castor_hf_diff_10->Fill(sum_cas_e,evtWeight);

          if(coll && sum_cas_e <= 700)                              h_castor_gap_hf->Fill(hf_single_energy_max,evtWeight);
          if(coll && sum_cas_e >  700)                              h_castor_nogap_hf->Fill(hf_single_energy_max,evtWeight);

          if(coll && hf_double_tag)                                 h_hf_hits_coll_lumi->Fill(event->lumiNb,hf_double_energy_max,evtWeight);
          if(coll && hf_double_tag)                                 h_hf_hits_minus_lumi->Fill(event->lumiNb,hf_m_energy_max,evtWeight);
          if(coll && hf_double_tag)                                 h_hf_hits_plus_lumi->Fill(event->lumiNb,hf_p_energy_max,evtWeight);
          if(noise)                                                 h_hf_hits_noise_lumi->Fill(event->lumiNb,hf_double_energy_max,noiseWeight);
          if(coll && hf_double_tag)                                 h_hf_totE_coll_lumi->Fill(event->lumiNb,hf_pm_energy,evtWeight);
          if(coll && hf_double_tag)                                 h_hf_totE_minus_lumi->Fill(event->lumiNb,hf_m_energy,evtWeight);
          if(coll && hf_double_tag)                                 h_hf_totE_plus_lumi->Fill(event->lumiNb,hf_p_energy,evtWeight);
          if(noise)                                                 h_hf_totE_noise_lumi->Fill(event->lumiNb,hf_pm_energy,noiseWeight);
          if(coll && hf_double_tag)                                 h_lumi_3GeV->Fill(event->lumiNb,evtWeight);
          if(coll)                                                  h_lumi->Fill(event->lumiNb,evtWeight);

          if(coll && hf_single_tag)                                 h_perf_hf_rechits_single_3gev->Fill(hf_n,evtWeight);
          if(coll && hf_double_tag)                                 h_perf_hf_rechits_double_1dot5gev->Fill(hf_n,evtWeight);

          if(coll && hf_single_tag)                                 h_perf_hf_totE_single_3gev->Fill(hf_pm_energy,evtWeight);
          if(coll && hf_double_tag)                                 h_perf_hf_totE_double_1dot5gev->Fill(hf_pm_energy,evtWeight);

          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              if(coll && hf_double_tag)
                h_perf_hf_totE_eta_double_1dot5gev->Fill(it->Eta,it->Energy);
              if(coll && hf_single_tag)
                h_perf_hf_totE_eta_single_3gev->Fill(it->Eta,it->Energy);
            } //no event weight. one would need to count weighted number of events

          if(coll)                                                  h_hf_new_cut_single->Fill(0.,evtWeight);
          if(coll && hf_light_single_tag)                           h_hf_new_cut_single->Fill(1.,evtWeight);
          if(coll && hf_medium_single_tag)                          h_hf_new_cut_single->Fill(2.,evtWeight);
          if(coll && hf_heavy_single_tag)                           h_hf_new_cut_single->Fill(3.,evtWeight);
          if(coll)                                                  h_hf_new_cut_double->Fill(0.,evtWeight);
          if(coll && hf_light_double_tag)                           h_hf_new_cut_double->Fill(1.,evtWeight);
          if(coll && hf_medium_double_tag)                          h_hf_new_cut_double->Fill(2.,evtWeight);
          if(coll && hf_heavy_double_tag)                           h_hf_new_cut_double->Fill(3.,evtWeight);

          if(sample_type[sample]==DATA)
            {

          for (double cut=0; cut < 20.1; cut+=0.1)
            {
              if(coll && hf_double_energy_max >= cut)               h_hf_cut_double->Fill(cut,evtWeight);
              if(noise && hf_double_energy_max >= cut)              h_hf_cut_double_noise->Fill(cut,noiseWeight);
              if(coll && hf_single_energy_max >= cut)               h_hf_cut_single->Fill(cut,evtWeight);
              if(noise && hf_single_energy_max >= cut)              h_hf_cut_single_noise->Fill(cut,noiseWeight);
            }


              for (int i = 0; i<array_size; i++)
                {
                  //create new tags in loop
                  bool hf_single_tag2 = hf_single_energy_max >= cut_energies_single[i];
                  bool hf_double_tag2 = hf_double_energy_max >= cut_energies_double[i];
                  if(coll && hf_single_tag2)                        v_run_events_single[i]->Fill(event->lumiNb,evtWeight);
                  if(coll && hf_double_tag2)                        v_run_events_double[i]->Fill(event->lumiNb,evtWeight);
                }

              if(coll)                                              h_run_events_lumi->Fill(event->lumiNb,lumiPerLS);

              if(noise)                                             h_hf_noise_all_lumi            ->Fill(event->lumiNb,noiseWeight);
              if(noise && hf_single_tag)                            h_hf_noise_selected_single_lumi->Fill(event->lumiNb,noiseWeight);
              if(noise && hf_double_tag)                            h_hf_noise_selected_double_lumi->Fill(event->lumiNb,noiseWeight);
            }


          if(sample_type[sample] == MC)
            {
              for (double cut=0; cut < 20.1; cut+=0.1)
                {
                  if(coll && hf_double_energy_max >= cut)               h_hf_cut_double->Fill(cut,evtWeight);
                  if(coll && hf_single_energy_max >= cut)               h_hf_cut_single->Fill(cut,evtWeight);
                }

              if(coll)                                                  h_mc_eff->Fill(0.,evtWeight);
              if(coll && hf_single_energy_max > 3)                      h_mc_eff->Fill(1.,evtWeight);
              if(coll && hf_single_energy_max > 4.5)                    h_mc_eff->Fill(2.,evtWeight);
              if(coll && hf_single_energy_max > 6)                      h_mc_eff->Fill(3.,evtWeight);
              if(coll && hf_double_energy_max > 2)                      h_mc_eff->Fill(4.,evtWeight);
              if(coll && hf_double_energy_max > 3)                      h_mc_eff->Fill(5.,evtWeight);
              if(coll && event->Tracks.size()>=1)                       h_mc_eff->Fill(6.,evtWeight);
              if(coll && castor_tag)                                    h_mc_eff->Fill(7.,evtWeight);

              if(!(hf_single_energy_max > 6))                                      h_mc_effgain_single->Fill(0.,evtWeight);
              if(!(hf_single_energy_max > 6) && hf_single_energy_max > 3)          h_mc_effgain_single->Fill(1.,evtWeight);
              if(!(hf_single_energy_max > 6) && hf_single_energy_max > 4.5)        h_mc_effgain_single->Fill(2.,evtWeight);
              if(!(hf_single_energy_max > 6) && hf_double_energy_max > 2)          h_mc_effgain_single->Fill(3.,evtWeight);
              if(!(hf_single_energy_max > 6) && hf_double_energy_max > 3)          h_mc_effgain_single->Fill(4.,evtWeight);
              if(!(hf_single_energy_max > 6) && event->Tracks.size()>=1)           h_mc_effgain_single->Fill(5.,evtWeight);
              if(!(hf_single_energy_max > 6) && castor_tag)                        h_mc_effgain_single->Fill(6.,evtWeight);

              if(!(hf_double_energy_max > 2))                                      h_mc_effgain_double->Fill(0.,evtWeight);
              if(!(hf_double_energy_max > 2) && hf_single_energy_max > 3)          h_mc_effgain_double->Fill(1.,evtWeight);
              if(!(hf_double_energy_max > 2) && hf_single_energy_max > 4.5)        h_mc_effgain_double->Fill(2.,evtWeight);
              if(!(hf_double_energy_max > 2) && hf_single_energy_max > 6)          h_mc_effgain_double->Fill(3.,evtWeight);
              if(!(hf_double_energy_max > 2) && event->Tracks.size()>=1)           h_mc_effgain_double->Fill(4.,evtWeight);
              if(!(hf_double_energy_max > 2) && castor_tag)                        h_mc_effgain_double->Fill(5.,evtWeight);


              if(coll && hf_single_tag)                                 h_mc_diffraction_single->Fill(log10(xi_x),evtWeight);
              if(coll && hf_double_tag)                                 h_mc_diffraction_double->Fill(log10(xi_x),evtWeight);
              if(coll && SD1)                                           h_mc_diffraction_SD1->Fill(log10(xi_x),evtWeight);
              if(coll && SD2)                                           h_mc_diffraction_SD2->Fill(log10(xi_x),evtWeight);
              if(coll && DD)                                            h_mc_diffraction_DD->Fill(log10(xi_x),evtWeight);
              if(coll && CD)                                            h_mc_diffraction_CD->Fill(log10(xi_x),evtWeight);
              if(coll && ND)                                            h_mc_diffraction_ND->Fill(log10(xi_x),evtWeight);
              if(coll)                                                  h_mc_diffraction_all->Fill(log10(xi_x),evtWeight);

              if(coll && hf_single_tag)                                 h_mc_lrg_single->Fill(rapGap,evtWeight);
              if(coll && hf_double_tag)                                 h_mc_lrg_double->Fill(rapGap,evtWeight);
              if(coll && SD1)                                           h_mc_lrg_SD1->Fill(rapGap,evtWeight);
              if(coll && SD2)                                           h_mc_lrg_SD2->Fill(rapGap,evtWeight);
              if(coll && DD)                                            h_mc_lrg_DD->Fill(rapGap,evtWeight);
              if(coll && CD)                                            h_mc_lrg_CD->Fill(rapGap,evtWeight);
              if(coll && ND)                                            h_mc_lrg_ND->Fill(rapGap,evtWeight);
              if(coll)                                                  h_mc_lrg_all->Fill(rapGap,evtWeight);
              if(coll)                                                  h_mc_lrg_xi->Fill(log10(xi_x),rapGap);
              if(coll)                                                  h_mc_lrg_xiy->Fill(log10(xi_y),rapGap);
              if(coll)                                                  h_mc_xix_xiy->Fill(log10(xi_x),log10(xi_y));
            }


        }

      //******************************************AFTER EVENT LOOP*******************************************

      cout << endl << "--- Processed " << n_total << " events." << endl << endl;
      h_perf_hf_totE_eta_double_1dot5gev->Scale(1./n_total);
      h_perf_hf_totE_eta_single_3gev->Scale(1./n_total);

      if(sample_type[sample] == MC)
        {
          h_mc_rapidity->Scale(1./n_total/h_mc_rapidity->GetBinWidth(1));
          h_mc_eta_e   ->Scale(1./n_total/h_mc_eta_e->GetBinWidth(1));
        }

    }

  //********************************************AFTER SAMPLE LOOP************************************************

  out_file->Write();
  out_file->Save();

  return 0;
}

int IetaToRing(int ring)
///converts ieta to numbers from 0 to 24
{
  if(ring < 0)
    return ring + 41;
  else
    return ring - 30;
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
