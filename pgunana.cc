#define _MAXEVT -10000
#define _SkipHFRings 1 //skip 41 and 29 as suggested by HCAL DPG
#define _HFEnergyScale 1.0 //1.0 //0.8
#define _HFEnergyCalibration 1 //0 or 1 (rescale MC) or 2 this does not scale MC but data according to raddam from lev

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
  for(double x=1; x<=10; x+=0.5)
    {
      cut_energies_single.push_back(x*2);
      cut_energies_double.push_back(x);
    }

  //*************************************************************INPUT***********************************************************
  sample_fname.push_back("/afs/cern.ch/user/c/cbaus/work/public/pgun_5_3_8_HI_patch2_tag_v26/treeMC.root"); sample_name.push_back("pgun"); sample_type.push_back(MC);


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

  TFile* out_file = new TFile("histos_deleteme.root","RECREATE");

  TH2D* h_hfreco_vs_hfgen;
  TH2D* h_hfreco_vs_hfgen_sum;


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

      // h_hf_hitdistr_coll_m       = new TH1D((add + string("_h_hf_hitdistr_coll_m")).c_str(),"",100,log10(0.5),log10(500.));

      // BinLogX(h_hf_hits_coll_single);

      h_hfreco_vs_hfgen    = new TH2D((add + string("_h_hfreco_vs_hfgen")).c_str(),"",40,0,100,40,0,100);
      h_hfreco_vs_hfgen_sum= new TH2D((add + string("_h_hfreco_vs_hfgen_sum")).c_str(),"",40,0,100,40,0,100);

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
               //double c_lev = eta<0?c_lev_m[-Ieta-29]:c_lev_p[Ieta-29];
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

             }
           double hf_pm_energy = hf_p_energy + hf_m_energy;
           //cout << "rechit: " << hf_pm_energy << endl;

           triggered_ring = hf_m_energy_max>hf_p_energy_max?triggered_ring_m:triggered_ring_p;

           hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
           hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
           bool hf_single_tag = hf_single_energy_max >= 8;
           bool hf_double_tag = hf_double_energy_max >= 4;



          //---------------------------------------------GEN Particles
          const double s = 5020*5020;
          double m_m=0, m_x=0, xi_x=0;
          double m_p=0, m_y=0, xi_y=0;
          double rapGap=-1;
          double gen_HF_E=0;
          double gen_eta=-99;
          double gen_HF_p=0;
          bool SD1 = event->genProcessID == 103; //lead dissociates
          bool SD2 = event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD1 && !SD2 && !DD && !CD;

          //Counting events
          const int prescale        = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS    = event->instLuminosity * event->instLuminosityCorr * 1e6;
          double evtWeight    = 1;

          if(sample_type[sample] == MC)
            {
              multimap<double,GenParticle*> phiMassMap;
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

                  double eta=it->GetEta();
                  double energy=it->GetEnergy();
                  double p=it->GetMomentum();

                  if(3 < fabs(eta) && fabs(eta) < 5 && energy > gen_HF_E)
                    gen_HF_E = energy;
                  if(3 < fabs(eta) && fabs(eta) < 5 && p > gen_HF_p)
                    gen_HF_p = p;
                  gen_eta=eta;

                  double phi= it->GetPhi();
                  phiMassMap.insert(pair<double,GenParticle*>(phi,&(*it)));
                }

            }

          if(sample_type[sample]==DATA)
            evtWeight = double(prescale);
          else if(sample_name[sample]=="EposDiffWeight150" && (SD1 || SD2 || DD || CD))
            evtWeight *= 1.50;
          else if(sample_name[sample]=="EposDiffWeight200" && (SD1 || SD2 || DD || CD))
            evtWeight *= 2.00;
          else if(sample_name[sample]=="EposDiffWeight299" && (SD1 || SD2 || DD || CD))
            evtWeight *= 2.99;
          else if(sample_name[sample]=="EposDiffWeightOpt" && (SD1 || SD2 || DD || CD))
            evtWeight *= 1.12;
          else if(sample_name[sample]=="QGSJetIIDiffWeight150" && (SD1 || SD2 || DD || CD))
            evtWeight *= 1.50;
          else if(sample_name[sample]=="QGSJetIIDiffWeight200" && (SD1 || SD2 || DD || CD))
            evtWeight *= 2.00;
          else if(sample_name[sample]=="QGSJetIIDiffWeight452" && (SD1 || SD2 || DD || CD))
            evtWeight *= 4.52;
          else if(sample_name[sample]=="QGSJetIIDiffWeightOpt" && (SD1 || SD2 || DD || CD))
            evtWeight *= 1.50;
          const double noiseWeight = random_prescale_HLT;


          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //---------------------------------------------Filling HISTOS

          if(3.5 < gen_eta && gen_eta < 4.5)
            {
              h_hfreco_vs_hfgen->Fill(gen_HF_p,hf_single_energy_max);
              h_hfreco_vs_hfgen_sum->Fill(gen_HF_p,hf_p_energy);
            }
        }

      //******************************************AFTER EVENT LOOP*******************************************

      cout << endl << "--- Processed " << n_total << " events." << endl << endl;

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
