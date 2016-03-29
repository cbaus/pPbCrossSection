#define _MAXEVT 50000
#define _SkipHFRings 1 //skip 41 and 29 as suggested by HCAL DPG
#define _HFEnergyScale 1.0 //1.0 //0.8
#define _HFEnergyCalibration 0 //0 or 1 (rescale MC) or 2 this does not scale MC but data according to raddam from lev

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
#include <iomanip>
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

int IetaToRing(int ieta);
void BinLogX(TH1* h);


int main()
{
  TH1::SetDefaultSumw2();

  vector<string> sample_fname;
  vector<string> sample_name;
  vector<e_type> sample_type;


  vector<double> cut_energies_single;
  vector<double> cut_energies_double;
  for(double x=1; x<=6; x+=0.25)
    {
      cut_energies_single.push_back(x*2);
      cut_energies_double.push_back(x);
    }

  //*************************************************************INPUT***********************************************************
  sample_fname.push_back("/tmp/cbaus/data259163.root"); sample_name.push_back("data259163"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/data_247324.root"); sample_name.push_back("data247324"); sample_type.push_back(DATA);

  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/pythiaz2star.root"); sample_name.push_back("PythiaZ2Star"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/pythiamonash.root"); sample_name.push_back("PythiaMonash"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/pythiambr.root"); sample_name.push_back("PythiaMBR"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/epos.root"); sample_name.push_back("Epos"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_diffraction/cbaus/pp13TeV/inel_cross_section/qgsjetii.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);


#if _HFEnergyCalibration == 1
  TFile calibfile("plots/hf_calibration_data.root");
  TVectorD* hf_calibration = NULL;
  hf_calibration=(TVectorD*)calibfile.Get("hf_calibration");
  hf_calibration->Print();
#endif

  vector<double> c_lev_m,c_lev_m_e,c_lev_p,c_lev_p_e;
  //TGraphErrors* gr_lev_p = new TGraphErrors(c_lev_m.size(),&c_lev_m.front());
  c_lev_m.push_back(1.07); //-29 //12 //innermost
  c_lev_m.push_back(0.97);
  c_lev_m.push_back(0.95);
  c_lev_m.push_back(0.99);
  c_lev_m.push_back(0.96);
  c_lev_m.push_back(0.91);
  c_lev_m.push_back(0.92);
  c_lev_m.push_back(0.86);
  c_lev_m.push_back(0.80);
  c_lev_m.push_back(0.72);
  c_lev_m.push_back(0.69);
  c_lev_m.push_back(0.83);
  c_lev_m.push_back(0.73); //-41 //0

  c_lev_p.push_back(1.01);
  c_lev_p.push_back(0.94);
  c_lev_p.push_back(0.91);
  c_lev_p.push_back(0.89);
  c_lev_p.push_back(0.87);
  c_lev_p.push_back(0.92);
  c_lev_p.push_back(0.84);
  c_lev_p.push_back(0.85);
  c_lev_p.push_back(0.83);
  c_lev_p.push_back(0.67);
  c_lev_p.push_back(0.61);
  c_lev_p.push_back(0.75);
  c_lev_p.push_back(0.66);

  c_lev_m_e.push_back(0.23);
  c_lev_m_e.push_back(0.04);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.03);
  c_lev_m_e.push_back(0.04);
  c_lev_m_e.push_back(0.06);
  c_lev_m_e.push_back(0.10);
  c_lev_m_e.push_back(0.09);

  c_lev_p_e.push_back(0.21);
  c_lev_p_e.push_back(0.04);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.03);
  c_lev_p_e.push_back(0.07);
  c_lev_p_e.push_back(0.09);
  c_lev_p_e.push_back(0.08);

  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos_new.root","RECREATE");

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
  TH1D* h_hf_hitdistr_coll;
  TH1D* h_hf_hitdistr_coll_p;
  TH1D* h_hf_hitdistr_coll_m;
  TH1D* h_hf_hitdistr_noise;
  TH1D* h_hf_hitdistr_noise_p;
  TH1D* h_hf_hitdistr_noise_m;
  TH1D* h_hf_triggered_ring_single;
  TH1D* h_hf_triggered_ring_double;
  TH1D* h_hf_triggered_ring_noise_single;
  TH1D* h_hf_triggered_ring_noise_double;
  TH1D* h_castor_hf_diff_3;
  TH1D* h_castor_hf_diff_5;
  TH1D* h_castor_hf_diff_10;
  TH1D* h_castor_gap_hf;
  TH1D* h_castor_nogap_hf;
  TH1D* h_mc_rapidity;
  TH1D* h_mc_eta_e;
  TH1D* h_mc_eta_e_SD1;
  TH1D* h_mc_eta_e_SD2;
  TH1D* h_mc_eta_e_CD;
  TH1D* h_mc_eta_e_DD;
  TH1D* h_mc_eta_e_ND;
  TH1D* h_mc_eff;
  TH1D* h_mc_effgain_single;
  TH1D* h_mc_effgain_double;

  TH2D* h_random_trig_tracks_hf;

  vector<TH1D*> v_run_events_single(cut_energies_single.size(),NULL);
  vector<TH1D*> v_run_events_double( cut_energies_double.size(),NULL);

  vector<TH1D*> v_hf_hits_rings(24,NULL);
  vector<TH1D*> v_hf_hits_rings_noise(24,NULL);

  TH1D* h_run_events_lumi;

  TH1D* h_hf_cut_single;
  TH1D* h_hf_cut_single_noise;
  TH1D* h_hf_cut_double;
  TH1D* h_hf_cut_double_noise;

  TH1D* h_hf_cut_ND_single;
  TH1D* h_hf_cut_SD1_single;
  TH1D* h_hf_cut_SD2_single;
  TH1D* h_hf_cut_DD_single;
  TH1D* h_hf_cut_CD_single;

  TH1D* h_hf_cut_ND_double;
  TH1D* h_hf_cut_SD1_double;
  TH1D* h_hf_cut_SD2_double;
  TH1D* h_hf_cut_DD_double;
  TH1D* h_hf_cut_CD_double;

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

  TH1D* h_perf_no_of_towers_single;
  TH1D* h_perf_no_of_towers_double;
  TH1D* h_perf_no_of_towers_aboveth_single;
  TH1D* h_perf_no_of_towers_aboveth_double;
  TH1D* h_perf_hf_totE_single_3gev;
  TH1D* h_perf_hf_totE_double_1dot5gev;
  TH1D* h_perf_hf_totE_ZBSingleTrack;
  TH1D* h_perf_hf_totE_ZBSingleTrack_noise;
  TH1D* h_perf_hf_totE_eta_single_3gev;
  TH1D* h_perf_hf_totE_eta_double_1dot5gev;
  TH1D* h_perf_hf_totE_eta_lev_p;
  TH1D* h_perf_hf_totE_eta_lev_m;
  TH1I* h_perf_hf_totE_eta_lev_n_p;
  TH1I* h_perf_hf_totE_eta_lev_n_m;

  TH1D* h_mc_diffraction_single;
  TH1D* h_mc_diffraction_double;
  TH1D* h_mc_diffraction_SD1;
  TH1D* h_mc_diffraction_SD2;
  TH1D* h_mc_diffraction_DD;
  TH1D* h_mc_diffraction_CD;
  TH1D* h_mc_diffraction_ND;
  TH1D* h_mc_diffraction_all;

  TH1D* h_mc_p_single;
  TH1D* h_mc_p_double;
  TH1D* h_mc_p_single_sel;
  TH1D* h_mc_p_double_sel;

  TH1D* h_mc_p_single_cut;
  TH1D* h_mc_p_double_cut;
  TH1D* h_mc_p_single_bg;
  TH1D* h_mc_p_double_bg;

  TH1D* h_mc_pt_single;
  TH1D* h_mc_pt_double;
  TH1D* h_mc_pt_single_sel;
  TH1D* h_mc_pt_double_sel;

  TH1D* h_mc_xi_single_cut;
  TH1D* h_mc_xi_double_cut;
  TH1D* h_mc_xi_single_bg;
  TH1D* h_mc_xi_double_bg;

  TH2D* h_mc_p_Ehf_correlation_single;
  TH2D* h_mc_p_Ehf_correlation_double;
  TH2D* h_mc_pt_Ehf_correlation_single;
  TH2D* h_mc_pt_Ehf_correlation_double;

  TH1D* h_mc_xisd_single;
  TH1D* h_mc_xisd_double;
  TH1D* h_mc_xisd_SD1;
  TH1D* h_mc_xisd_SD2;
  TH1D* h_mc_xisd_DD;
  TH1D* h_mc_xisd_CD;
  TH1D* h_mc_xisd_ND;
  TH1D* h_mc_xisd_all;

  TH1D* h_mc_xisd_plus_single;
  TH1D* h_mc_xisd_plus_double;
  TH1D* h_mc_xisd_plus_SD1;
  TH1D* h_mc_xisd_plus_SD2;
  TH1D* h_mc_xisd_plus_DD;
  TH1D* h_mc_xisd_plus_CD;
  TH1D* h_mc_xisd_plus_ND;
  TH1D* h_mc_xisd_plus_all;

  TH1D* h_mc_xisd_minus_single;
  TH1D* h_mc_xisd_minus_double;
  TH1D* h_mc_xisd_minus_SD1;
  TH1D* h_mc_xisd_minus_SD2;
  TH1D* h_mc_xisd_minus_DD;
  TH1D* h_mc_xisd_minus_CD;
  TH1D* h_mc_xisd_minus_ND;
  TH1D* h_mc_xisd_minus_all;

  TH1D* h_mc_diff_p_single_SD1;
  TH1D* h_mc_diff_p_single_SD2;
  TH1D* h_mc_diff_p_single_DD;
  TH1D* h_mc_diff_p_single_CD;
  TH1D* h_mc_diff_p_single_ND;
  TH1D* h_mc_diff_p_single_all;
  TH1D* h_mc_diff_p_double_SD1;
  TH1D* h_mc_diff_p_double_SD2;
  TH1D* h_mc_diff_p_double_DD;
  TH1D* h_mc_diff_p_double_CD;
  TH1D* h_mc_diff_p_double_ND;
  TH1D* h_mc_diff_p_double_all;

  TH1D* h_mc_diff_pt_single_SD1;
  TH1D* h_mc_diff_pt_single_SD2;
  TH1D* h_mc_diff_pt_single_DD;
  TH1D* h_mc_diff_pt_single_CD;
  TH1D* h_mc_diff_pt_single_ND;
  TH1D* h_mc_diff_pt_single_all;
  TH1D* h_mc_diff_pt_double_SD1;
  TH1D* h_mc_diff_pt_double_SD2;
  TH1D* h_mc_diff_pt_double_DD;
  TH1D* h_mc_diff_pt_double_CD;
  TH1D* h_mc_diff_pt_double_ND;
  TH1D* h_mc_diff_pt_double_all;

  TH1D* h_mc_diff_e_single_SD1;
  TH1D* h_mc_diff_e_single_SD2;
  TH1D* h_mc_diff_e_single_DD;
  TH1D* h_mc_diff_e_single_CD;
  TH1D* h_mc_diff_e_single_ND;
  TH1D* h_mc_diff_e_single_all;
  TH1D* h_mc_diff_e_double_SD1;
  TH1D* h_mc_diff_e_double_SD2;
  TH1D* h_mc_diff_e_double_DD;
  TH1D* h_mc_diff_e_double_CD;
  TH1D* h_mc_diff_e_double_ND;
  TH1D* h_mc_diff_e_double_all;

  TH2D* h_mc_lrg_xi;
  TH2D* h_mc_lrg_xiy;
  TH2D* h_mc_xix_xiy;
  TH2D* h_mc_xix_xiy_sdsel;
  TH2D* h_mc_xix_xiy_ndsel;
  TH2D* h_mc_unfold;

  vector<double> deta_lev, eta_lev_m,eta_lev_p;
  //eta_lev_p.push_back(2.865);
  eta_lev_p.push_back(2.975);
  eta_lev_p.push_back(3.15);
  eta_lev_p.push_back(3.325);
  eta_lev_p.push_back(3.5);
  eta_lev_p.push_back(3.675);
  eta_lev_p.push_back(3.85);
  eta_lev_p.push_back(4.025);
  eta_lev_p.push_back(4.2);
  eta_lev_p.push_back(4.375);
  eta_lev_p.push_back(4.55);
  eta_lev_p.push_back(4.725);
  eta_lev_p.push_back(4.9);
  //eta_lev_p.push_back(5.2);

  const int neta_lev = eta_lev_p.size()-1; //-1 because of last bin low edge

  for (int j=eta_lev_p.size()-1; j>=0; j--)
    {
      eta_lev_m.push_back(-eta_lev_p[j]);
    }
  for (int j=0; j<int(eta_lev_p.size()); j++)
    {
      cout << j << " " << eta_lev_p[j] << " " << eta_lev_m[j] << endl;
    }



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

      double n_total = double(tree->GetEntries());
      int n_castor_tag = 0;
      int n_hf_single_tag = 0;
      int n_hf_double_tag = 0;
      if (n_total == 0.) {
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

      tree->SetBranchStatus("EB*",0);
      tree->SetBranchStatus("EE*",0);
      tree->SetBranchStatus("HBHE*",0);

      //________________________________________________________________________________


      AnalysisEvent* event = 0;
      tree->SetBranchAddress("event", &event);

      int zero_bias;
      int zero_bias_prescale_L1;
      int zero_bias_prescale_HLT;
      int random;
      int random_prescale_HLT;
      int bptx_p_m;
      int bptx_p_nm;
      int bptx_np_m;
      int bptx_quiet;
      if(sample_type[sample] == DATA)
        {
          tree->SetBranchAddress("L1_ZeroBias_algPrescale",&zero_bias_prescale_L1);
          if (sample_name[sample] == "data247324")
            {
              tree->SetBranchAddress("HLT_ZeroBias_part0_v1_Prescl",&zero_bias_prescale_HLT);
              tree->SetBranchAddress("HLT_ZeroBias_part0_v1",&zero_bias);
              tree->SetBranchAddress("HLT_Random_v2",&random);
              tree->SetBranchAddress("HLT_Random_v2_Prescl",&random_prescale_HLT);
            }
          else
            {
              tree->SetBranchAddress("HLT_ZeroBias_v1_Prescl",&zero_bias_prescale_HLT);
              tree->SetBranchAddress("HLT_ZeroBias_v1",&zero_bias);
              tree->SetBranchAddress("HLT_Random_v1",&random);
              tree->SetBranchAddress("HLT_Random_v1_Prescl",&random_prescale_HLT);
            }
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
      h_hf_hits_coll_single       = new TH1D((add + string("_h_hf_hits_coll_single")).c_str(),"",40,log10(0.5),log10(500.));
      h_hf_hits_coll_double       = new TH1D((add + string("_h_hf_hits_coll_double")).c_str(),"",40,log10(0.5),log10(500.));
      h_hf_hits_plus              = new TH1D((add + string("_h_hfp_hits_coll")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hits_minus             = new TH1D((add + string("_h_hfm_hits_coll")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hitdistr_coll          = new TH1D((add + string("_h_hf_hitdistr_coll")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hitdistr_coll_p        = new TH1D((add + string("_h_hf_hitdistr_coll_p")).c_str(),"",100,log10(0.5),log10(500.));
      h_hf_hitdistr_coll_m        = new TH1D((add + string("_h_hf_hitdistr_coll_m")).c_str(),"",100,log10(0.5),log10(500.));

      BinLogX(h_hf_hits_coll_single);
      BinLogX(h_hf_hits_coll_double);
      BinLogX(h_hf_hits_minus);
      BinLogX(h_hf_hits_plus);
      BinLogX(h_hf_hitdistr_coll);
      BinLogX(h_hf_hitdistr_coll_p);
      BinLogX(h_hf_hitdistr_coll_m);

      h_hf_triggered_ring_single  = new TH1D((add + string("_h_hf_triggered_ring_single")).c_str(),"",83,-41.5,41.5);
      h_hf_triggered_ring_double  = new TH1D((add + string("_h_hf_triggered_ring_double")).c_str(),"",83,-41.5,41.5);
      h_hf_triggered_ring_noise_single  = new TH1D((add + string("_h_hf_triggered_ring_noise_single")).c_str(),"",83,-41.5,41.5);
      h_hf_triggered_ring_noise_double  = new TH1D((add + string("_h_hf_triggered_ring_noise_double")).c_str(),"",83,-41.5,41.5);

      h_castor_hf_diff_3          = new TH1D((add + string("_h_castor_hf_diff_3")).c_str(),"",100,0,10000);
      h_castor_hf_diff_5          = new TH1D((add + string("_h_castor_hf_diff_5")).c_str(),"",100,0,10000);
      h_castor_hf_diff_10         = new TH1D((add + string("_h_castor_hf_diff_10")).c_str(),"",100,0,10000);
      h_castor_gap_hf             = new TH1D((add + string("_h_castor_gap_hf")).c_str(),"",100,0,50);
      h_castor_nogap_hf           = new TH1D((add + string("_h_castor_nogap_hf")).c_str(),"",100,0,50);

      for (int i = 0; i<24; i++)
        {
          ostringstream ss_single; ss_single << add <<"_h_hf_hits_rings_" << i;
          v_hf_hits_rings[i] = new TH1D(ss_single.str().c_str(),run_str.str().c_str(),100,log10(0.5),log10(500.));
          ostringstream ss_noise; ss_noise << add <<"_h_hf_hits_rings_noise_" << i;
          v_hf_hits_rings_noise[i] = new TH1D(ss_noise.str().c_str(),run_str.str().c_str(),100,log10(0.5),log10(500.));

          BinLogX(v_hf_hits_rings[i]);
          BinLogX(v_hf_hits_rings_noise[i]);
        }

      if(sample_type[sample] == DATA)
        {

          h_hf_hits_noise_single      = new TH1D((add + string("_h_hf_hits_noise_single")).c_str(),"",40,log10(0.5),log10(500.));
          h_hf_hits_noise_double      = new TH1D((add + string("_h_hf_hits_noise_double")).c_str(),"",40,log10(0.5),log10(500.));
          h_hf_hits_noise_plus        = new TH1D((add + string("_h_hfp_hits_noise")).c_str(),"",100,log10(0.5),log10(500.));
          h_hf_hits_noise_minus       = new TH1D((add + string("_h_hfm_hits_noise")).c_str(),"",100,log10(0.5),log10(500.));
          h_hf_hitdistr_noise         = new TH1D((add + string("_h_hf_hitdistr_noise")).c_str(),"",501,-50.,50.);
          h_hf_hitdistr_noise_p       = new TH1D((add + string("_h_hf_hitdistr_noise_p")).c_str(),"",501,-50.,50.);
          h_hf_hitdistr_noise_m       = new TH1D((add + string("_h_hf_hitdistr_noise_m")).c_str(),"",501,-50.,50.);
          BinLogX(h_hf_hits_noise_single);
          BinLogX(h_hf_hits_noise_double);
          BinLogX(h_hf_hits_noise_plus);
          BinLogX(h_hf_hits_noise_minus);
          // BinLogX(h_hf_hitdistr_noise);
          // BinLogX(h_hf_hitdistr_noise_p);
          // BinLogX(h_hf_hitdistr_noise_m);

          for (int i = 0; i<int(cut_energies_single.size()); i++)
            {
              ostringstream ss_single; ss_single << add <<"_h_run_events_single_" << cut_energies_single[i];
              ostringstream ss_double; ss_double << add <<"_h_run_events_double_" << cut_energies_double[i];
              v_run_events_single[i] = new TH1D(ss_single.str().c_str(),run_str.str().c_str(),2000,0,2000);
              v_run_events_double[i] = new TH1D(ss_double.str().c_str(),run_str.str().c_str(),2000,0,2000);
            }

          h_run_events_lumi = new TProfile((add+string("_h_run_events_lumi")).c_str(),run_str.str().c_str(),2000,0,2000,"s");


          h_random_trig_tracks_hf          = new TH2D((add + string("_h_random_trig_tracks_hf")).c_str(),"",200,0,20,20,0,20);

          h_hf_noise_all_lumi              = new TH1D((add + string("_h_hf_noise_all_lumi")).c_str(),"",2000,0,2000);
          h_hf_noise_selected_single_lumi  = new TH1D((add + string("_h_hf_noise_selected_single_lumi")).c_str(),"",2000,0,2000);
          h_hf_noise_selected_double_lumi  = new TH1D((add + string("_h_hf_noise_selected_double_lumi")).c_str(),"",2000,0,2000);
        }

      h_hf_cut_single             = new TH1D((add + string("_h_hf_cut_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_single_noise       = new TH1D((add + string("_h_hf_cut_single_noise")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_double             = new TH1D((add + string("_h_hf_cut_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_double_noise       = new TH1D((add + string("_h_hf_cut_double_noise")).c_str(),"",201,-0.05,20.05);

      h_hf_cut_ND_single             = new TH1D((add + string("_h_hf_cut_ND_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_SD1_single             = new TH1D((add + string("_h_hf_cut_SD1_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_SD2_single             = new TH1D((add + string("_h_hf_cut_SD2_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_CD_single             = new TH1D((add + string("_h_hf_cut_CD_single")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_DD_single             = new TH1D((add + string("_h_hf_cut_DD_single")).c_str(),"",201,-0.05,20.05);

      h_hf_cut_ND_double             = new TH1D((add + string("_h_hf_cut_ND_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_SD1_double             = new TH1D((add + string("_h_hf_cut_SD1_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_SD2_double             = new TH1D((add + string("_h_hf_cut_SD2_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_CD_double             = new TH1D((add + string("_h_hf_cut_CD_double")).c_str(),"",201,-0.05,20.05);
      h_hf_cut_DD_double             = new TH1D((add + string("_h_hf_cut_DD_double")).c_str(),"",201,-0.05,20.05);

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

      h_perf_no_of_towers_single           = new TH1D((add + string("_h_perf_no_of_towers_single")).c_str(),"",100,0,1000);
      h_perf_no_of_towers_double           = new TH1D((add + string("_h_perf_no_of_towers_double")).c_str(),"",100,0,1000);
      h_perf_no_of_towers_aboveth_single   = new TH1D((add + string("_h_perf_no_of_towers_aboveth_single")).c_str(),"",100,0,1000);
      h_perf_no_of_towers_aboveth_double   = new TH1D((add + string("_h_perf_no_of_towers_aboveth_double")).c_str(),"",100,0,1000);
      h_perf_hf_totE_single_3gev           = new TH1D((add + string("_h_perf_hf_totE_single_3gev")).c_str(),"",500,0,10);
      h_perf_hf_totE_double_1dot5gev       = new TH1D((add + string("_h_perf_hf_totE_double_1dot5gev")).c_str(),"",500,0,10);
      h_perf_hf_totE_ZBSingleTrack         = new TH1D((add + string("_h_perf_hf_totE_ZBSingleTrack")).c_str(),"",500,0,10);
      h_perf_hf_totE_ZBSingleTrack_noise   = new TH1D((add + string("_h_perf_hf_totE_ZBSingleTrack_noise")).c_str(),"",500,0,10);
      h_perf_hf_totE_eta_single_3gev       = new TH1D((add + string("_h_perf_hf_totE_eta_single_3gev")).c_str(),"",100,-5.2,5.2);
      h_perf_hf_totE_eta_double_1dot5gev   = new TH1D((add + string("_h_perf_hf_totE_eta_double_1dot5gev")).c_str(),"",100,-5.2,5.2);
      h_perf_hf_totE_eta_lev_m     = new TH1D((add + string("_h_perf_hf_totE_eta_lev_m")).c_str(),"",neta_lev,eta_lev_m.data());
      h_perf_hf_totE_eta_lev_p     = new TH1D((add + string("_h_perf_hf_totE_eta_lev_p")).c_str(),"",neta_lev,eta_lev_p.data());
      h_perf_hf_totE_eta_lev_n_m   = new TH1I((add + string("_h_perf_hf_totE_eta_lev_n_m")).c_str(),"",neta_lev,eta_lev_m.data());
      h_perf_hf_totE_eta_lev_n_p   = new TH1I((add + string("_h_perf_hf_totE_eta_lev_n_p")).c_str(),"",neta_lev,eta_lev_p.data());
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

          h_mc_xisd_single         = new TH1D((add + string("_h_mc_xisd_single")).c_str(),"",50,-9.5,2);
          h_mc_xisd_double         = new TH1D((add + string("_h_mc_xisd_double")).c_str(),"",50,-9.5,2);
          h_mc_xisd_SD1            = new TH1D((add + string("_h_mc_xisd_SD1")).c_str(),"",50,-9.5,2);
          h_mc_xisd_SD2            = new TH1D((add + string("_h_mc_xisd_SD2")).c_str(),"",50,-9.5,2);
          h_mc_xisd_CD             = new TH1D((add + string("_h_mc_xisd_CD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_DD             = new TH1D((add + string("_h_mc_xisd_DD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_ND             = new TH1D((add + string("_h_mc_xisd_ND")).c_str(),"",50,-9.5,2);
          h_mc_xisd_all            = new TH1D((add + string("_h_mc_xisd_all")).c_str(),"",50,-9.5,2);

          h_mc_xisd_plus_single    = new TH1D((add + string("_h_mc_xisd_plus_single")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_double    = new TH1D((add + string("_h_mc_xisd_plus_double")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_SD1       = new TH1D((add + string("_h_mc_xisd_plus_SD1")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_SD2       = new TH1D((add + string("_h_mc_xisd_plus_SD2")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_CD        = new TH1D((add + string("_h_mc_xisd_plus_CD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_DD        = new TH1D((add + string("_h_mc_xisd_plus_DD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_ND        = new TH1D((add + string("_h_mc_xisd_plus_ND")).c_str(),"",50,-9.5,2);
          h_mc_xisd_plus_all       = new TH1D((add + string("_h_mc_xisd_plus_all")).c_str(),"",50,-9.5,2);

          h_mc_xisd_minus_single   = new TH1D((add + string("_h_mc_xisd_minus_single")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_double   = new TH1D((add + string("_h_mc_xisd_minus_double")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_SD1      = new TH1D((add + string("_h_mc_xisd_minus_SD1")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_SD2      = new TH1D((add + string("_h_mc_xisd_minus_SD2")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_CD       = new TH1D((add + string("_h_mc_xisd_minus_CD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_DD       = new TH1D((add + string("_h_mc_xisd_minus_DD")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_ND       = new TH1D((add + string("_h_mc_xisd_minus_ND")).c_str(),"",50,-9.5,2);
          h_mc_xisd_minus_all      = new TH1D((add + string("_h_mc_xisd_minus_all")).c_str(),"",50,-9.5,2);

          h_mc_diff_p_single_SD1  = new TH1D((add + string("_h_mc_diff_p_single_SD1")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_single_SD2  = new TH1D((add + string("_h_mc_diff_p_single_SD2")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_single_CD   = new TH1D((add + string("_h_mc_diff_p_single_CD")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_single_DD   = new TH1D((add + string("_h_mc_diff_p_single_DD")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_single_ND   = new TH1D((add + string("_h_mc_diff_p_single_ND")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_single_all  = new TH1D((add + string("_h_mc_diff_p_single_all")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_SD1  = new TH1D((add + string("_h_mc_diff_p_double_SD1")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_SD2  = new TH1D((add + string("_h_mc_diff_p_double_SD2")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_CD   = new TH1D((add + string("_h_mc_diff_p_double_CD")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_DD   = new TH1D((add + string("_h_mc_diff_p_double_DD")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_ND   = new TH1D((add + string("_h_mc_diff_p_double_ND")).c_str(),"",101,-0.5,100.5);
          h_mc_diff_p_double_all  = new TH1D((add + string("_h_mc_diff_p_double_all")).c_str(),"",101,-0.5,100.5);

          h_mc_diff_pt_single_SD1  = new TH1D((add + string("_h_mc_diff_pt_single_SD1")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_single_SD2  = new TH1D((add + string("_h_mc_diff_pt_single_SD2")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_single_CD   = new TH1D((add + string("_h_mc_diff_pt_single_CD")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_single_DD   = new TH1D((add + string("_h_mc_diff_pt_single_DD")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_single_ND   = new TH1D((add + string("_h_mc_diff_pt_single_ND")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_single_all  = new TH1D((add + string("_h_mc_diff_pt_single_all")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_SD1  = new TH1D((add + string("_h_mc_diff_pt_double_SD1")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_SD2  = new TH1D((add + string("_h_mc_diff_pt_double_SD2")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_CD   = new TH1D((add + string("_h_mc_diff_pt_double_CD")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_DD   = new TH1D((add + string("_h_mc_diff_pt_double_DD")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_ND   = new TH1D((add + string("_h_mc_diff_pt_double_ND")).c_str(),"",81,-0.025,4.025);
          h_mc_diff_pt_double_all  = new TH1D((add + string("_h_mc_diff_pt_double_all")).c_str(),"",81,-0.025,4.025);

          h_mc_lrg_xi             = new TH2D((add + string("_h_mc_lrg_xi")).c_str(),"",200,-9.5,2,500,-1,20);
          h_mc_lrg_xiy            = new TH2D((add + string("_h_mc_lrg_xiy")).c_str(),"",200,-9.5,2,500,-1,20);
          h_mc_xix_xiy            = new TH2D((add + string("_h_mc_xix_xiy")).c_str(),"",200,-9.5,2,200,-9.5,2);
          h_mc_xix_xiy_sdsel      = new TH2D((add + string("_h_mc_xix_xiy_sdsel")).c_str(),"",200,-9.5,2,200,-9.5,2);
          h_mc_xix_xiy_ndsel      = new TH2D((add + string("_h_mc_xix_xiy_ndsel")).c_str(),"",200,-9.5,2,200,-9.5,2);
          h_mc_unfold             = new TH2D((add + string("_h_mc_unfold")).c_str(),"",15,0,200,15,0,200);

          h_mc_diff_e_single_SD1  = new TH1D((add + string("_h_mc_diff_e_single_SD1")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_single_SD2  = new TH1D((add + string("_h_mc_diff_e_single_SD2")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_single_CD   = new TH1D((add + string("_h_mc_diff_e_single_CD")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_single_DD   = new TH1D((add + string("_h_mc_diff_e_single_DD")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_single_ND   = new TH1D((add + string("_h_mc_diff_e_single_ND")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_single_all  = new TH1D((add + string("_h_mc_diff_e_single_all")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_SD1  = new TH1D((add + string("_h_mc_diff_e_double_SD1")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_SD2  = new TH1D((add + string("_h_mc_diff_e_double_SD2")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_CD   = new TH1D((add + string("_h_mc_diff_e_double_CD")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_DD   = new TH1D((add + string("_h_mc_diff_e_double_DD")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_ND   = new TH1D((add + string("_h_mc_diff_e_double_ND")).c_str(),"",50,log10(0.5),log10(700));
          h_mc_diff_e_double_all  = new TH1D((add + string("_h_mc_diff_e_double_all")).c_str(),"",50,log10(0.5),log10(700));
          BinLogX(h_mc_diff_e_single_SD1);
          BinLogX(h_mc_diff_e_single_SD2);
          BinLogX(h_mc_diff_e_single_CD);
          BinLogX(h_mc_diff_e_single_DD);
          BinLogX(h_mc_diff_e_single_ND);
          BinLogX(h_mc_diff_e_single_all);
          BinLogX(h_mc_diff_e_double_SD1);
          BinLogX(h_mc_diff_e_double_SD2);
          BinLogX(h_mc_diff_e_double_CD);
          BinLogX(h_mc_diff_e_double_DD);
          BinLogX(h_mc_diff_e_double_ND);
          BinLogX(h_mc_diff_e_double_all);

	  h_mc_p_single        = new TH1D((add + string("_h_mc_p_single")).c_str(),"",50,log10(1e-3),log10(5e2));
	  h_mc_p_double        = new TH1D((add + string("_h_mc_p_double")).c_str(),"",50,log10(1e-3),log10(5e2));
	  h_mc_p_single_sel    = new TH1D((add + string("_h_mc_p_single_sel")).c_str(),"",50,log10(1e-3),log10(5e2));
	  h_mc_p_double_sel    = new TH1D((add + string("_h_mc_p_double_sel")).c_str(),"",50,log10(1e-3),log10(5e2));
	  h_mc_p_single_cut = new TH1D((add + string("_h_mc_p_single_cut")).c_str(),"",101,-0.5,100.5);
	  h_mc_p_double_cut = new TH1D((add + string("_h_mc_p_double_cut")).c_str(),"",101,-0.5,100.5);
	  h_mc_p_single_bg = new TH1D((add + string("_h_mc_p_single_bg")).c_str(),"",101,-0.5,100.5);
	  h_mc_p_double_bg = new TH1D((add + string("_h_mc_p_double_bg")).c_str(),"",101,-0.5,100.5);
          BinLogX(h_mc_p_single);
          BinLogX(h_mc_p_double);
          BinLogX(h_mc_p_single_sel);
          BinLogX(h_mc_p_double_sel);

	  h_mc_pt_single        = new TH1D((add + string("_h_mc_pt_single")).c_str(),"",50,log10(1e-3),log10(1e2));
	  h_mc_pt_double        = new TH1D((add + string("_h_mc_pt_double")).c_str(),"",50,log10(1e-3),log10(1e2));
	  h_mc_pt_single_sel    = new TH1D((add + string("_h_mc_pt_single_sel")).c_str(),"",50,log10(1e-3),log10(1e2));
	  h_mc_pt_double_sel    = new TH1D((add + string("_h_mc_pt_double_sel")).c_str(),"",50,log10(1e-3),log10(1e2));
	  h_mc_xi_single_cut = new TH1D((add + string("_h_mc_xi_single_cut")).c_str(),"",101,-10.05,0.05);
	  h_mc_xi_double_cut = new TH1D((add + string("_h_mc_xi_double_cut")).c_str(),"",101,-10.05,0.05);
	  h_mc_xi_single_bg = new TH1D((add + string("_h_mc_xi_single_bg")).c_str(),"",101,-10.05,0.05);
	  h_mc_xi_double_bg = new TH1D((add + string("_h_mc_xi_double_bg")).c_str(),"",101,-10.05,0.05);
          BinLogX(h_mc_pt_single);
          BinLogX(h_mc_pt_double);
          BinLogX(h_mc_pt_single_sel);
          BinLogX(h_mc_pt_double_sel);

          h_mc_pt_Ehf_correlation_single = new TH2D((add + string("_h_mc_pt_Ehf_correlation_single")).c_str(),"",100,0,9,100,0,400);
          h_mc_pt_Ehf_correlation_double = new TH2D((add + string("_h_mc_pt_Ehf_correlation_double")).c_str(),"",100,0,7,100,0,300);
          h_mc_p_Ehf_correlation_single  = new TH2D((add + string("_h_mc_p_Ehf_correlation_single")).c_str(),"",100,0,800,100,0,400);
          h_mc_p_Ehf_correlation_double  = new TH2D((add + string("_h_mc_p_Ehf_correlation_double")).c_str(),"",100,0,500,100,0,300);

          h_mc_rapidity           = new TH1D((add + string("_h_mc_rapidity")).c_str(),"",100,-12,12);
          h_mc_eta_e              = new TH1D((add + string("_h_mc_eta_e")).c_str(),"",100,-12,12);
          h_mc_eta_e_SD1          = new TH1D((add + string("_h_mc_eta_e_SD1")).c_str(),"",100,-12,12);
          h_mc_eta_e_SD2          = new TH1D((add + string("_h_mc_eta_e_SD2")).c_str(),"",100,-12,12);
          h_mc_eta_e_CD           = new TH1D((add + string("_h_mc_eta_e_CD")).c_str(),"",100,-12,12);
          h_mc_eta_e_DD           = new TH1D((add + string("_h_mc_eta_e_DD")).c_str(),"",100,-12,12);
          h_mc_eta_e_ND           = new TH1D((add + string("_h_mc_eta_e_ND")).c_str(),"",100,-12,12);

          h_mc_eff                = new TH1D((add + string("_h_mc_eff")).c_str(),"",9,-0.5,8.5);
          h_mc_eff->GetXaxis()->SetBinLabel(1,"All");
          h_mc_eff->GetXaxis()->SetBinLabel(2,"HF single > 6 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(3,"HF single > 8 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(4,"HF single > 10 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(5,"HF double > 3 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(6,"HF double > 4 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(7,"HF double > 5 GeV");
          h_mc_eff->GetXaxis()->SetBinLabel(8,"no. of tracks #geq 1");
          h_mc_eff->GetXaxis()->SetBinLabel(9,"CASTOR E_{tot} > 12.5 GeV");
          h_mc_effgain_single     = new TH1D((add + string("_h_mc_effgain_single")).c_str(),"",7,-0.5,6.5);
          h_mc_effgain_single->GetXaxis()->SetBinLabel(1,"All");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(2,"HF single > 6 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(3,"HF double > 3 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(4,"HF double > 4 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(5,"HF double > 5 GeV");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(6,"no. of tracks #geq 1");
          h_mc_effgain_single->GetXaxis()->SetBinLabel(7,"CASTOR E_{tot} > 12.5 GeV");
          h_mc_effgain_double     = new TH1D((add + string("_h_mc_effgain_double")).c_str(),"",7,-0.5,6.5);
          h_mc_effgain_double->GetXaxis()->SetBinLabel(1,"All");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(2,"HF single > 6 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(3,"HF single > 8 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(4,"HF single > 10 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(5,"HF double > 3 GeV");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(6,"no. of tracks  #geq 1");
          h_mc_effgain_double->GetXaxis()->SetBinLabel(7,"CASTOR E_{tot} > 12.5 GeV");
        }

      if(_MAXEVT<n_total && _MAXEVT>0)
        n_total = double(_MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 100 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);
          //          if(event->runNb != 210885)
          //continue;

          bool coll=0, noise=0;

          if(sample_type[sample] == DATA)
            {
              coll          = zero_bias && bptx_p_m; //double beam
              noise         = random && !bptx_p_m;// && !bptx_np_m && !bptx_p_nm; //not both and not single beam
            }
          else if(sample_type[sample] == MC)
            {
              noise = 0;
              coll = 1;
            }

          if(!coll && !noise) //not intersted
            continue;

          if (sample_type[sample] == DATA && event->lumiNb <= 88 && event->runNb == 247324)
            continue; //hack for wrong json



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
          int hf_n_aboveth = 0;
          int hf_n = event->HFtowers.size();
          int hf_zero_count = ForwardRecord::nMaxHFMRecHits - hf_n;
          double hf_double_energy_max = 0;
          double hf_single_energy_max = 0;
          vector<double> v_hf_rings_energy_max(24,0.);

          double hf_m_energy_max = 0;
          double hf_p_energy_max = 0;
          double hf_p_energy = 0;
          double hf_m_energy = 0;
          int triggered_ring = 0;
          int triggered_ring_m = 0;
          int triggered_ring_p = 0;
          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              if(_SkipHFRings && (it->IetaAbs == 41 || it->IetaAbs == 29))
                continue;
              const double eta = it->Eta;
              const int Ieta = it->Eta > 0?it->IetaAbs:-it->IetaAbs;
              double tower_e = it->Energy * _HFEnergyScale;
              double c_lev = eta<0?c_lev_m[-Ieta-29]:c_lev_p[Ieta-29];
#if _HFEnergyCalibration == 1
              if(sample_type[sample]==MC)
                {
                  //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
                  tower_e /= (*hf_calibration)[IetaToRing(Ieta)];
                }
#endif
#if _HFEnergyCalibration == 2
              if(sample_type[sample]==DATA)
                {
                  //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
                  tower_e /= c_lev;//(*hf_calibration)[IetaToRing(Ieta)];
                }
#endif
              //cout << it->IetaAbs << " " << eta << " " << tower_e << endl;

              if(eta < 0 && tower_e > v_hf_rings_energy_max[IetaToRing(Ieta)])
                v_hf_rings_energy_max[IetaToRing(Ieta)] = tower_e;
              if(eta > 0 && tower_e > v_hf_rings_energy_max[IetaToRing(Ieta)])
                v_hf_rings_energy_max[IetaToRing(Ieta)] = tower_e;

              if(eta > 0. && tower_e > hf_p_energy_max)
                {
                  hf_p_energy_max = tower_e;
                  triggered_ring_p = Ieta;
                }
              if(eta <= 0. && tower_e > hf_m_energy_max)
                {
                  hf_m_energy_max = tower_e;
                  triggered_ring_m = Ieta;
                }

              if(eta > 0.)
                hf_p_energy += tower_e;
              else
                hf_m_energy += tower_e;

              if(tower_e > 3.)
                hf_n_aboveth++;

              if(coll)               h_hf_hitdistr_coll->   Fill(tower_e);
              if(coll && eta <= 0.)  h_hf_hitdistr_coll_m-> Fill(tower_e);
              if(coll && eta > 0.)   h_hf_hitdistr_coll_p-> Fill(tower_e);
              if(noise)              h_hf_hitdistr_noise->  Fill(tower_e);
              if(noise && eta <= 0.) h_hf_hitdistr_noise_m->Fill(tower_e);
              if(noise && eta > 0.)  h_hf_hitdistr_noise_p->Fill(tower_e);
            }
          double hf_pm_energy = hf_p_energy + hf_m_energy;
          //cout << "rechit: " << hf_pm_energy << endl;

          triggered_ring = hf_m_energy_max>hf_p_energy_max?triggered_ring_m:triggered_ring_p;

          hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
          hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
          bool hf_single_tag = hf_single_energy_max >= 5;
          bool hf_double_tag = hf_double_energy_max >= 3;



          //---------------------------------------------GEN Particles
          const double s = 5020*5020;
          double xi_p=-10;
          double xi_m=-10;
          double rapGap=-1;
          double ymin = -5.2;
          double ymax = 5.2;
          double gen_HF_E=0;
          double gen_HF_pt_double=0;
          double gen_HF_pt_single=0;
          double gen_HF_p_double=0;
          double gen_HF_p_single=0;
          bool SD1 = event->genProcessID == 103; //lead dissociates
          bool SD2 = event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD1 && !SD2 && !DD && !CD;
          double xi_sd = 0;

          //Counting events
          const int prescale        = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS    = event->instLuminosity * event->instLuminosityCorr * 1e6;
          double evtWeight    = 1;

          if(sample_type[sample] == MC)
            {
	      double gen_HF_pt_minus = 0;
	      double gen_HF_pt_plus = 0;
	      double gen_HF_p_minus = 0;
	      double gen_HF_p_plus = 0;
              multimap<double,GenParticle*> rapidityMap;
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

                  const double Rapidity= it->GetRapidity();
                  rapidityMap. insert(pair<double,GenParticle*>(Rapidity,&(*it)));

                  if(abs(it->Id) == 12 || abs(it->Id) == 14 || abs(it->Id) == 16 || abs(it->Id) == 13)
                    continue; // skip muon + neutrinos

                  const double eta=it->GetEta();
                  const double energy=it->GetEnergy();
                  const double pt=it->GetTransverseMomentum();
                  const double p=it->GetMomentum();

                  if(3 < fabs(eta) && fabs(eta) < 5 && energy > gen_HF_E)
                    gen_HF_E = energy;

                  if(-5 < eta && eta < -3 && pt > gen_HF_pt_minus)
                    gen_HF_pt_minus = pt;
                  if(-5 < eta && eta < -3 && p > gen_HF_p_minus)
                    gen_HF_p_minus = p;

                  if(+3 < eta && eta < +5 && pt > gen_HF_pt_plus)
                    gen_HF_pt_plus = pt;
                  if(+3 < eta && eta < +5 && p > gen_HF_p_plus)
                    gen_HF_p_plus = p;

                  if (pt > 0.2 && Rapidity > ymin && Rapidity < ymax)
		    h_mc_rapidity->Fill(Rapidity);

                  h_mc_eta_e->Fill(eta,evtWeight);
                  if(SD1)
                    h_mc_eta_e_SD1->Fill(eta,evtWeight);
                  if(SD2)
                    h_mc_eta_e_SD2->Fill(eta,evtWeight);
                  if(CD)
                    h_mc_eta_e_CD->Fill(eta,evtWeight);
                  if(DD)
                    h_mc_eta_e_DD->Fill(eta,evtWeight);
                  if(ND)
                    h_mc_eta_e_ND->Fill(eta,evtWeight);
		}
	      gen_HF_pt_single = TMath::Max(gen_HF_pt_minus,gen_HF_pt_plus);
	      gen_HF_pt_double = TMath::Min(gen_HF_pt_minus,gen_HF_pt_plus);

	      gen_HF_p_single = TMath::Max(gen_HF_p_minus,gen_HF_p_plus);
	      gen_HF_p_double = TMath::Min(gen_HF_p_minus,gen_HF_p_plus);

              //find lrg
              if(rapidityMap.size() >= 1) //Assume 1 or more final usable particles
                {
                  multimap<double,GenParticle*>::const_iterator ItIsAtRapGap = rapidityMap.begin(); //will later be start of m_x
                  for (multimap<double,GenParticle*>::const_iterator it = rapidityMap.begin(); it != rapidityMap.end(); ++it)
                    {
                      if (it == rapidityMap.begin()) //skip first event in order to compare to previous particle eta
                        continue;
                      double currRapGap = -1;
                      multimap<double,GenParticle*>::const_iterator itPrev=it;
                      --itPrev;
                      //eta not used
                      currRapGap = fabs(it->second->GetRapidity() - itPrev->second->GetRapidity());

                      if(!TMath::Finite(currRapGap) || TMath::IsNaN(currRapGap))
                        continue;

                      //cout << itPrev->second->momentum().eta() << itPrev->second->pdg_id() << " " << " deltaY=" << rapGap;
                      if (currRapGap > rapGap)
                        {
                          rapGap = currRapGap;
                          ItIsAtRapGap = it;
                          //cout << " (updated rapgap) ";
                        }
                      //cout << endl;
                    }

                  //split into m_x and m_y system
                  TLorentzVector vecM(0,0,0,0);
                  multimap<double,GenParticle*>::const_iterator it = rapidityMap.begin();
                  for (; it != ItIsAtRapGap; ++it)
                    {
                      const TLorentzVector vec(it->second->Px,it->second->Py,it->second->Pz,it->second->GetEnergy());
                      vecM += vec;
                      //cout << "vecM " << it->second->second->pdg_id() << " " << vec4.px() << " " << vec4.py() << " " << vec4.pz() << " " << vec4.e() << " invM=" << vec.M() << endl;
                    }

                  TLorentzVector vecP(0,0,0,0);
                  multimap<double,GenParticle*>::const_iterator it2 = ItIsAtRapGap;
                  for (; it2 != rapidityMap.end(); ++it2)
                    {
                      const TLorentzVector vec(it2->second->Px,it2->second->Py,it2->second->Pz,it2->second->GetEnergy());
                      vecP += vec;
                      //cout << "vecP " << it2->second->pdg_id() << " " << vec4.px() << " " << vec4.py() << " " << vec4.pz() << " " << vec4.e() << " invM=" << vec.M() << endl;
                    }
                  const double sqrts = 13000.;
                  const double mass_m = vecM.M();
                  const double mass_p = vecP.M();

                  //const double mass_x = TMath::Max(mass_m, mass_p);
                  //const double mass_y = TMath::Min(mass_m, mass_p);
                  //cout << "mx=" << mass_x << " my=" << mass_y << endl;
                  //cout << "log10(mx)=" << log10(mass_x) << " log10(my)=" << log10(mass_y) << endl;

                  xi_p = mass_p*mass_p / pow(sqrts,2);
                  xi_m = mass_m*mass_m / pow(sqrts,2);
                  xi_sd = max(xi_p, xi_m);

                  //const double xi_x = mass_x*mass_x / pow(sqrts,2);
                  //const double xi_y = mass_y*mass_y / pow(sqrts,2);

                  // const bool isSD = log10(mass_m) < 0.5;
                  // const bool isDD = 0.5 < log10(mass_m) && log10(mass_m) < 1.1;

                  // const double deltaRapGap = -log(pow(mass_p,2)*pow(mass_m,2)/(pow(sqrts,2)*pow(0.93827,2))); //deltarap > 3 === -log(xi) but no idea why this is identity and the m1 here comes directly from roberts code
                }
            }
          if(sample_type[sample]==DATA)
            evtWeight = double(prescale);
          const double noiseWeight = random_prescale_HLT;

          if (castor_tag) n_castor_tag++;
          if (hf_single_tag) n_hf_single_tag++;
          if (hf_double_tag) n_hf_double_tag++;

          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //---------------------------------------------Filling HISTOS
          if(coll)                                                  h_zero_count_zb_coll->Fill(hf_zero_count,evtWeight);
          if(noise)                                                 h_zero_count_zb_no_coll->Fill(hf_zero_count,noiseWeight);

          for (int i = 0; i<24; i++)
            {
              if(coll)                                              v_hf_hits_rings[i]      ->Fill(v_hf_rings_energy_max[i]);
              if(noise)                                             v_hf_hits_rings_noise[i]->Fill(v_hf_rings_energy_max[i]);
            }

          if(coll)                                                  h_hf_hits_coll_single->Fill(hf_single_energy_max);
          if(coll)                                                  h_hf_hits_coll_double->Fill(hf_double_energy_max);

          if(coll)                                                  h_hf_hits_plus->Fill(hf_p_energy_max);
          if(coll)                                                  h_hf_hits_minus->Fill(hf_m_energy_max);

          if(coll && hf_double_energy_max < 3)                      h_castor_hf_diff_3->Fill(sum_CAS_E,evtWeight);
          if(coll && hf_double_energy_max < 5)                      h_castor_hf_diff_5->Fill(sum_CAS_E,evtWeight);
          if(coll && hf_double_energy_max < 10)                     h_castor_hf_diff_10->Fill(sum_CAS_E,evtWeight);

          if(coll && sum_CAS_E <= 700)                              h_castor_gap_hf->Fill(hf_single_energy_max,evtWeight);
          if(coll && sum_CAS_E >  700)                              h_castor_nogap_hf->Fill(hf_single_energy_max,evtWeight);

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

          if(coll && hf_single_tag)                                 h_perf_no_of_towers_single->Fill(hf_n,evtWeight);
          if(coll && hf_double_tag)                                 h_perf_no_of_towers_double->Fill(hf_n,evtWeight);
          if(coll && hf_single_tag)                                 h_perf_no_of_towers_aboveth_single->Fill(hf_n_aboveth,evtWeight);
          if(coll && hf_double_tag)                                 h_perf_no_of_towers_aboveth_double->Fill(hf_n_aboveth,evtWeight);

          if(coll && hf_single_tag)                                 h_perf_hf_totE_single_3gev->Fill(hf_pm_energy/1000.,evtWeight);
          if(coll && hf_double_tag)                                 h_perf_hf_totE_double_1dot5gev->Fill(hf_pm_energy/1000.,evtWeight);
          if(coll && event->Tracks.size()>=1)                       h_perf_hf_totE_ZBSingleTrack->Fill(hf_pm_energy/1000.,evtWeight);
          if(noise && event->Tracks.size()>=1)                      h_perf_hf_totE_ZBSingleTrack_noise->Fill(hf_pm_energy/1000.,evtWeight);

          if(coll && hf_single_tag)                                 h_hf_triggered_ring_single->Fill(triggered_ring);
          if(coll && hf_double_tag)                                 h_hf_triggered_ring_double->Fill(triggered_ring_m);
          if(coll && hf_double_tag)                                 h_hf_triggered_ring_double->Fill(triggered_ring_p);
          if(noise && hf_single_tag)                                 h_hf_triggered_ring_noise_single->Fill(triggered_ring);
          if(noise && hf_double_tag)                                 h_hf_triggered_ring_noise_double->Fill(triggered_ring_m);
          if(noise && hf_double_tag)                                 h_hf_triggered_ring_noise_double->Fill(triggered_ring_p);

          for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
            {
              if(_SkipHFRings && (it->IetaAbs == 41 || it->IetaAbs == 29))
                continue;

              if(coll && hf_double_tag)
                {
                  h_perf_hf_totE_eta_double_1dot5gev->Fill(it->Eta,it->Energy);
                }
              if(coll && hf_single_tag)
                h_perf_hf_totE_eta_single_3gev->Fill(it->Eta,it->Energy);
              if(coll && hf_single_tag)// event->Tracks.size()>=1)
                {
                  const double eta = it->Eta;
                  const int Ieta = eta>0?it->IetaAbs:-it->IetaAbs;
                  const double c_lev = eta<0?c_lev_m[-Ieta-29]:c_lev_p[Ieta-29];
                  double tower_e = it->Energy;
#if _HFEnergyCalibration == 1
                  if(sample_type[sample]==MC)
                    {
                      //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
                      tower_e /= (*hf_calibration)[IetaToRing(Ieta)];
                      //cout << Ieta << " " << eta << " " << c_lev << endl;
                    }
#endif
#if _HFEnergyCalibration == 2
                  if(sample_type[sample]==DATA)
                    {
                      //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
                      tower_e /= c_lev;
                      //cout << Ieta << " " << eta << " " << c_lev << endl;
                    }
#endif
                  h_perf_hf_totE_eta_lev_m->Fill(it->Eta,tower_e);
                  h_perf_hf_totE_eta_lev_p->Fill(it->Eta,tower_e);
                  h_perf_hf_totE_eta_lev_n_m->Fill(it->Eta);
                  h_perf_hf_totE_eta_lev_n_p->Fill(it->Eta);
                }
            } //no event weight. one would need to count weighted number of events

          if(sample_type[sample]==DATA)
            {

              if(noise)                                                 h_hf_hits_noise_plus->Fill(hf_p_energy_max);
              if(noise)                                                 h_hf_hits_noise_minus->Fill(hf_m_energy_max);
              if(noise)                                                 h_hf_hits_noise_single->Fill(hf_single_energy_max);
              if(noise)                                                 h_hf_hits_noise_double->Fill(hf_double_energy_max);

              for (double cut=0; cut < 20.1; cut+=0.1)
                {
                  if(coll && hf_double_energy_max >= cut)               h_hf_cut_double->Fill(cut,evtWeight);
                  if(noise && hf_double_energy_max >= cut)              h_hf_cut_double_noise->Fill(cut,noiseWeight);
                  if(coll && hf_single_energy_max >= cut)               h_hf_cut_single->Fill(cut,evtWeight);
                  if(noise && hf_single_energy_max >= cut)              h_hf_cut_single_noise->Fill(cut,noiseWeight);
                }


              for (int i = 0; i<int(cut_energies_single.size()); i++)
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

                  if(coll && ND && hf_single_energy_max >= cut)         h_hf_cut_ND_single->Fill(cut,evtWeight);
                  if(coll && SD1&& hf_single_energy_max >= cut)         h_hf_cut_SD1_single->Fill(cut,evtWeight);
                  if(coll && SD2&& hf_single_energy_max >= cut)         h_hf_cut_SD2_single->Fill(cut,evtWeight);
                  if(coll && CD && hf_single_energy_max >= cut)         h_hf_cut_CD_single->Fill(cut,evtWeight);
                  if(coll && DD && hf_single_energy_max >= cut)         h_hf_cut_DD_single->Fill(cut,evtWeight);

                  if(coll && ND && hf_double_energy_max >= cut)         h_hf_cut_ND_double->Fill(cut,evtWeight);
                  if(coll && SD1&& hf_double_energy_max >= cut)         h_hf_cut_SD1_double->Fill(cut,evtWeight);
                  if(coll && SD2&& hf_double_energy_max >= cut)         h_hf_cut_SD2_double->Fill(cut,evtWeight);
                  if(coll && CD && hf_double_energy_max >= cut)         h_hf_cut_CD_double->Fill(cut,evtWeight);
                  if(coll && DD && hf_double_energy_max >= cut)         h_hf_cut_DD_double->Fill(cut,evtWeight);

                }

              if(coll)                                                  h_mc_eff->Fill(0.,evtWeight);
              if(coll && hf_single_energy_max >= 6)                      h_mc_eff->Fill(1.,evtWeight);
              if(coll && hf_single_energy_max >= 8)                      h_mc_eff->Fill(2.,evtWeight);
              if(coll && hf_single_energy_max >= 10)                     h_mc_eff->Fill(3.,evtWeight);
              if(coll && hf_double_energy_max >= 3)                      h_mc_eff->Fill(4.,evtWeight);
              if(coll && hf_double_energy_max >= 4)                     h_mc_eff->Fill(5.,evtWeight);
              if(coll && hf_double_energy_max >= 5)                      h_mc_eff->Fill(6.,evtWeight);
              if(coll && event->Tracks.size()>=1)                       h_mc_eff->Fill(7.,evtWeight);
              if(coll && castor_tag)                                    h_mc_eff->Fill(8.,evtWeight);

              if(!(hf_single_energy_max > 8))                                     h_mc_effgain_single->Fill(0.,evtWeight);
              if(!(hf_single_energy_max > 8) && hf_single_energy_max > 6)         h_mc_effgain_single->Fill(1.,evtWeight);
              if(!(hf_single_energy_max > 8) && hf_double_energy_max > 3)         h_mc_effgain_single->Fill(2.,evtWeight);
              if(!(hf_single_energy_max > 8) && hf_double_energy_max > 4)         h_mc_effgain_single->Fill(3.,evtWeight);
              if(!(hf_single_energy_max > 8) && hf_double_energy_max > 5)         h_mc_effgain_single->Fill(4.,evtWeight);
              if(!(hf_single_energy_max > 8) && event->Tracks.size()>=1)          h_mc_effgain_single->Fill(5.,evtWeight);
              if(!(hf_single_energy_max > 8) && castor_tag)                       h_mc_effgain_single->Fill(6.,evtWeight);

              if(!(hf_double_energy_max > 4))                                      h_mc_effgain_double->Fill(0.,evtWeight);
              if(!(hf_double_energy_max > 4) && hf_single_energy_max > 6)          h_mc_effgain_double->Fill(1.,evtWeight);
              if(!(hf_double_energy_max > 4) && hf_single_energy_max > 8)          h_mc_effgain_double->Fill(2.,evtWeight);
              if(!(hf_double_energy_max > 4) && hf_single_energy_max > 10)         h_mc_effgain_double->Fill(3.,evtWeight);
              if(!(hf_double_energy_max > 4) && hf_double_energy_max > 2)          h_mc_effgain_double->Fill(4.,evtWeight);
              if(!(hf_double_energy_max > 4) && event->Tracks.size()>=1)           h_mc_effgain_double->Fill(5.,evtWeight);
              if(!(hf_double_energy_max > 4) && castor_tag)                        h_mc_effgain_double->Fill(6.,evtWeight);

              if(coll && hf_single_tag)                                 h_mc_diffraction_single->Fill(log10(xi_m),evtWeight);
              if(coll && hf_double_tag)                                 h_mc_diffraction_double->Fill(log10(xi_m),evtWeight);
              if(coll && SD1)                                           h_mc_diffraction_SD1->Fill(log10(xi_m),evtWeight);
              if(coll && SD2)                                           h_mc_diffraction_SD2->Fill(log10(xi_m),evtWeight);
              if(coll && DD)                                            h_mc_diffraction_DD->Fill(log10(xi_m),evtWeight);
              if(coll && CD)                                            h_mc_diffraction_CD->Fill(log10(xi_m),evtWeight);
              if(coll && ND)                                            h_mc_diffraction_ND->Fill(log10(xi_m),evtWeight);
              if(coll)                                                  h_mc_diffraction_all->Fill(log10(xi_m),evtWeight);

              if(coll && SD1)                                           h_mc_diff_p_single_SD1->Fill(gen_HF_p_single,evtWeight);
              if(coll && SD2)                                           h_mc_diff_p_single_SD2->Fill(gen_HF_p_single,evtWeight);
              if(coll && DD)                                            h_mc_diff_p_single_DD->Fill(gen_HF_p_single,evtWeight);
              if(coll && CD)                                            h_mc_diff_p_single_CD->Fill(gen_HF_p_single,evtWeight);
              if(coll && ND)                                            h_mc_diff_p_single_ND->Fill(gen_HF_p_single,evtWeight);
              if(coll)                                                  h_mc_diff_p_single_all->Fill(gen_HF_p_single,evtWeight);
              if(coll && SD1)                                           h_mc_diff_p_double_SD1->Fill(gen_HF_p_double,evtWeight);
              if(coll && SD2)                                           h_mc_diff_p_double_SD2->Fill(gen_HF_p_double,evtWeight);
              if(coll && DD)                                            h_mc_diff_p_double_DD->Fill(gen_HF_p_double,evtWeight);
              if(coll && CD)                                            h_mc_diff_p_double_CD->Fill(gen_HF_p_double,evtWeight);
              if(coll && ND)                                            h_mc_diff_p_double_ND->Fill(gen_HF_p_double,evtWeight);
              if(coll)                                                  h_mc_diff_p_double_all->Fill(gen_HF_p_double,evtWeight);

              if(coll && SD1)                                           h_mc_diff_pt_single_SD1->Fill(gen_HF_pt_single,evtWeight);
              if(coll && SD2)                                           h_mc_diff_pt_single_SD2->Fill(gen_HF_pt_single,evtWeight);
              if(coll && DD)                                            h_mc_diff_pt_single_DD->Fill(gen_HF_pt_single,evtWeight);
              if(coll && CD)                                            h_mc_diff_pt_single_CD->Fill(gen_HF_pt_single,evtWeight);
              if(coll && ND)                                            h_mc_diff_pt_single_ND->Fill(gen_HF_pt_single,evtWeight);
              if(coll)                                                  h_mc_diff_pt_single_all->Fill(gen_HF_pt_single,evtWeight);
              if(coll && SD1)                                           h_mc_diff_pt_double_SD1->Fill(gen_HF_pt_double,evtWeight);
              if(coll && SD2)                                           h_mc_diff_pt_double_SD2->Fill(gen_HF_pt_double,evtWeight);
              if(coll && DD)                                            h_mc_diff_pt_double_DD->Fill(gen_HF_pt_double,evtWeight);
              if(coll && CD)                                            h_mc_diff_pt_double_CD->Fill(gen_HF_pt_double,evtWeight);
              if(coll && ND)                                            h_mc_diff_pt_double_ND->Fill(gen_HF_pt_double,evtWeight);
              if(coll)                                                  h_mc_diff_pt_double_all->Fill(gen_HF_pt_double,evtWeight);

              if(coll && SD1)                                           h_mc_diff_e_single_SD1->Fill(hf_single_energy_max,evtWeight);
              if(coll && SD2)                                           h_mc_diff_e_single_SD2->Fill(hf_single_energy_max,evtWeight);
              if(coll && DD)                                            h_mc_diff_e_single_DD->Fill(hf_single_energy_max,evtWeight);
              if(coll && CD)                                            h_mc_diff_e_single_CD->Fill(hf_single_energy_max,evtWeight);
              if(coll && ND)                                            h_mc_diff_e_single_ND->Fill(hf_single_energy_max,evtWeight);
              if(coll)                                                  h_mc_diff_e_single_all->Fill(hf_single_energy_max,evtWeight);
              if(coll && SD1)                                           h_mc_diff_e_double_SD1->Fill(hf_double_energy_max,evtWeight);
              if(coll && SD2)                                           h_mc_diff_e_double_SD2->Fill(hf_double_energy_max,evtWeight);
              if(coll && DD)                                            h_mc_diff_e_double_DD->Fill(hf_double_energy_max,evtWeight);
              if(coll && CD)                                            h_mc_diff_e_double_CD->Fill(hf_double_energy_max,evtWeight);
              if(coll && ND)                                            h_mc_diff_e_double_ND->Fill(hf_double_energy_max,evtWeight);
              if(coll)                                                  h_mc_diff_e_double_all->Fill(hf_double_energy_max,evtWeight);

              if(coll && hf_single_tag)                                 h_mc_xisd_single->Fill(log10(xi_sd),evtWeight);
              if(coll && hf_double_tag)                                 h_mc_xisd_double->Fill(log10(xi_sd),evtWeight);
              if(coll && SD1)                                           h_mc_xisd_SD1->Fill(log10(xi_sd),evtWeight);
              if(coll && SD2)                                           h_mc_xisd_SD2->Fill(log10(xi_sd),evtWeight);
              if(coll && DD)                                            h_mc_xisd_DD->Fill(log10(xi_sd),evtWeight);
              if(coll && CD)                                            h_mc_xisd_CD->Fill(log10(xi_sd),evtWeight);
              if(coll && ND)                                            h_mc_xisd_ND->Fill(log10(xi_sd),evtWeight);
              if(coll)                                                  h_mc_xisd_all->Fill(log10(xi_sd),evtWeight);

              if(coll && hf_single_tag)                                 h_mc_xisd_minus_single->Fill(log10(xi_m),evtWeight);
              if(coll && hf_double_tag)                                 h_mc_xisd_minus_double->Fill(log10(xi_m),evtWeight);
              if(coll && SD1)                                           h_mc_xisd_minus_SD1->Fill(log10(xi_m),evtWeight);
              if(coll && SD2)                                           h_mc_xisd_minus_SD2->Fill(log10(xi_m),evtWeight);
              if(coll && DD)                                            h_mc_xisd_minus_DD->Fill(log10(xi_m),evtWeight);
              if(coll && CD)                                            h_mc_xisd_minus_CD->Fill(log10(xi_m),evtWeight);
              if(coll && ND)                                            h_mc_xisd_minus_ND->Fill(log10(xi_m),evtWeight);
              if(coll)                                                  h_mc_xisd_minus_all->Fill(log10(xi_m),evtWeight);

              if(coll && hf_single_tag)                                 h_mc_xisd_plus_single->Fill(log10(xi_p),evtWeight);
              if(coll && hf_double_tag)                                 h_mc_xisd_plus_double->Fill(log10(xi_p),evtWeight);
              if(coll && SD1)                                           h_mc_xisd_plus_SD1->Fill(log10(xi_p),evtWeight);
              if(coll && SD2)                                           h_mc_xisd_plus_SD2->Fill(log10(xi_p),evtWeight);
              if(coll && DD)                                            h_mc_xisd_plus_DD->Fill(log10(xi_p),evtWeight);
              if(coll && CD)                                            h_mc_xisd_plus_CD->Fill(log10(xi_p),evtWeight);
              if(coll && ND)                                            h_mc_xisd_plus_ND->Fill(log10(xi_p),evtWeight);
              if(coll)                                                  h_mc_xisd_plus_all->Fill(log10(xi_p),evtWeight);

              if(coll)                                                  h_mc_lrg_xi->Fill(log10(xi_sd),rapGap);
              if(coll)                                                  h_mc_lrg_xiy->Fill(log10(xi_p),rapGap);
              if(coll)                                                  h_mc_xix_xiy->Fill(log10(xi_m),log10(xi_p));
              if(coll && hf_single_tag)                                 h_mc_xix_xiy_sdsel->Fill(log10(xi_m),log10(xi_p));
              if(coll
                 && hf_double_tag
                 && event->nVertex >= 1
                 && event->vertexIsFake == 0
                 && event->Tracks.size() >= 3)                          h_mc_xix_xiy_ndsel->Fill(log10(xi_m),log10(xi_p));

              if(coll)                                                  h_mc_p_single->Fill(gen_HF_p_single,evtWeight);
              if(coll)                                                  h_mc_p_double->Fill(gen_HF_p_double,evtWeight);
              if(coll && hf_single_tag)                                 h_mc_p_single_sel->Fill(gen_HF_p_single,evtWeight);
              if(coll && hf_double_tag)                                 h_mc_p_double_sel->Fill(gen_HF_p_double,evtWeight);

              if(coll)                                                  h_mc_pt_single->Fill(gen_HF_pt_single,evtWeight);
              if(coll)                                                  h_mc_pt_double->Fill(gen_HF_pt_double,evtWeight);
              if(coll && hf_single_tag)                                 h_mc_pt_single_sel->Fill(gen_HF_pt_single,evtWeight);
              if(coll && hf_double_tag)                                 h_mc_pt_double_sel->Fill(gen_HF_pt_double,evtWeight);

              if(coll)                                                  h_mc_pt_Ehf_correlation_single->Fill(log10(xi_sd),hf_single_energy_max);
              if(coll)                                                  h_mc_pt_Ehf_correlation_double->Fill(log10(xi_sd),hf_double_energy_max);
              if(coll)                                                  h_mc_p_Ehf_correlation_single->Fill(gen_HF_p_single,hf_single_energy_max);
              if(coll)                                                  h_mc_p_Ehf_correlation_double->Fill(gen_HF_p_double,hf_double_energy_max);




              for (double cut=0; cut < 100; cut+=1.)
                { // bin 1 contains all events
                  if(coll && hf_single_tag && gen_HF_p_single >= cut)  h_mc_p_single_cut->Fill(cut,evtWeight);
                  if(coll && hf_double_tag && gen_HF_p_double >= cut)  h_mc_p_double_cut->Fill(cut,evtWeight);
                  if(coll && !hf_single_tag && gen_HF_p_single >= cut) h_mc_p_single_bg->Fill(cut,evtWeight);
                  if(coll && !hf_double_tag && gen_HF_p_double >= cut) h_mc_p_double_bg->Fill(cut,evtWeight);
                }

              for (double cut=-10; cut < 0; cut+=0.1)
                { // bin 1 contains all events
                  if(coll && hf_single_tag && log10(xi_sd) >= cut)  h_mc_xi_single_cut->Fill(cut,evtWeight);
                  if(coll && hf_double_tag && log10(xi_sd) >= cut)  h_mc_xi_double_cut->Fill(cut,evtWeight);
                  if(coll && !hf_single_tag && log10(xi_sd) >= cut) h_mc_xi_single_bg->Fill(cut,evtWeight);
                  if(coll && !hf_double_tag && log10(xi_sd) >= cut) h_mc_xi_double_bg->Fill(cut,evtWeight);
                }

              if(coll) h_mc_unfold->Fill(hf_single_energy_max,gen_HF_E);
            }


        }

      //******************************************AFTER EVENT LOOP*******************************************

      cout << endl << "--- Processed " << n_total << " events." << endl << endl;
      h_perf_hf_totE_eta_double_1dot5gev->Scale(1./n_total);
      h_perf_hf_totE_eta_single_3gev->Scale(1./n_total);

      h_hf_hitdistr_coll->Scale(1./n_total);
      h_hf_hitdistr_coll_m->Scale(1./n_total);
      h_hf_hitdistr_coll_p->Scale(1./n_total);

      if(sample_type[sample] == DATA)
        {
          h_hf_hitdistr_noise->Scale(1./n_total);
          h_hf_hitdistr_noise_p->Scale(1./n_total);
          h_hf_hitdistr_noise_m->Scale(1./n_total);
        }


      if(sample_type[sample] == MC)
        {
          h_mc_rapidity   ->Scale(1./n_total/h_mc_rapidity->GetBinWidth(1));
          h_mc_eta_e      ->Scale(1./n_total/h_mc_eta_e->GetBinWidth(1));
          h_mc_eta_e_SD1  ->Scale(1./n_total/h_mc_eta_e_SD1->GetBinWidth(1));
          h_mc_eta_e_SD2  ->Scale(1./n_total/h_mc_eta_e_SD2->GetBinWidth(1));
          h_mc_eta_e_CD   ->Scale(1./n_total/h_mc_eta_e_CD->GetBinWidth(1));
          h_mc_eta_e_DD   ->Scale(1./n_total/h_mc_eta_e_DD->GetBinWidth(1));
          h_mc_eta_e_ND   ->Scale(1./n_total/h_mc_eta_e_ND->GetBinWidth(1));
        }

      for (int j=0;j<=neta_lev;j++)//bin_size reduced for skipped rings
        {
          int bin=j+1;
          if(n_hf_single_tag)
            h_perf_hf_totE_eta_lev_m->SetBinContent(bin,h_perf_hf_totE_eta_lev_m->GetBinContent(bin)/(eta_lev_m[j+1]-eta_lev_m[j])/double(n_hf_single_tag));
          else
            h_perf_hf_totE_eta_lev_m->SetBinContent(bin,0);

          if(n_hf_single_tag)
            h_perf_hf_totE_eta_lev_p->SetBinContent(bin,h_perf_hf_totE_eta_lev_p->GetBinContent(bin)/(eta_lev_p[j+1]-eta_lev_p[j])/double(n_hf_single_tag));
          else
            h_perf_hf_totE_eta_lev_p->SetBinContent(bin,0);
        }

    }

  //********************************************AFTER SAMPLE LOOP************************************************

  out_file->Write();
  out_file->Save();
  out_file->Close();

  return 0;
}

int IetaToRing(int ieta)
///converts ieta to numbers from 0 to 23
{
  //rings range from [-41,29] and [29,41]
  //29 is skipped in trees
  //41 is skipped usually because of HCAL prescription
  if(ieta < 0)
    return ieta + 41;
  else
    return ieta + 12 - 30;
}

int RingToIeta(int ring)
///converts rings from 0 to 23 to ieta
{
  //rings range from [-41,29] and [29,41]
  //29 is skipped in trees
  //41 is skipped usually because of HCAL prescription
  if(ring < 12)
    return ring - 41;
  else
    return ring - 12 + 30;
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