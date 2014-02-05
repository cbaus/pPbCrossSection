#define _MAXEVT 20000
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
#include <iterator>     // std::distance
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
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210614/*_*.root"); sample_name.push_back("data210614"); sample_type.push_back(DATA);
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210885/*_*.root"); sample_name.push_back("data210885"); sample_type.push_back(DATA);
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
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("EposDiffWeight150"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("EposDiffWeight200"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("EposDiffWeight299"); sample_type.push_back(MC);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("EposDiffWeightOpt"); sample_type.push_back(MC);
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos_SL/*.root"); sample_name.push_back("Epos_SL"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Hijing/*.root"); sample_name.push_back("Hijing"); sample_type.push_back(MC);
  //sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);
  // sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetIIDiffWeight150"); sample_type.push_back(MC);
  // sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetIIDiffWeight200"); sample_type.push_back(MC);
  // sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetIIDiffWeight452"); sample_type.push_back(MC);
  // sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetIIDiffWeightOpt"); sample_type.push_back(MC);
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightDPMjet_v2/treeMC.root"); sample_name.push_back("Starlight_DPMJet");  sample_type.push_back(MC);
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightPythia/treeMC.root"); sample_name.push_back("Starlight_Pythia");  sample_type.push_back(MC);


  //**************************************************************OUTPUT*********************************************************

  TFile* out_file = new TFile("histos_deleteme.root","RECREATE");

  TH1D* h_mc_castor_neutron_energy;
  TH1D* h_mc_castor_neutron_energy_isolated;
  TH1D* h_mc_castor_kaon_energy_isolated;
  TH1D* h_mc_castor_lambda_energy_isolated;
  TH1D* h_castor_neutron_energy;


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

      h_mc_castor_neutron_energy          = new TH1D((add + string("_h_mc_castor_neutron_energy")).c_str(),"",100,0,4000);
      h_mc_castor_neutron_energy_isolated = new TH1D((add + string("_h_mc_castor_neutron_energy_isolated")).c_str(),"",100,0,4000);
      h_mc_castor_kaon_energy_isolated    = new TH1D((add + string("_h_mc_castor_kaon_energy_isolated")).c_str(),"",100,0,4000);
      h_mc_castor_lambda_energy_isolated  = new TH1D((add + string("_h_mc_castor_lambda_energy_isolated")).c_str(),"",100,0,4000);
      h_castor_neutron_energy             = new TH1D((add + string("_h_castor_neutron_energy")).c_str(),"",100,0,4000);

      double n_total = double(tree->GetEntries());
      if(_MAXEVT<n_total && _MAXEVT>0)
        n_total = double(_MAXEVT);

      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 1000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);

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

          //Counting events
          const int prescale        = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS    = event->instLuminosity * event->instLuminosityCorr * 1e6;
          double evtWeight    = 1;
          const double s = 5020*5020;
          bool SD1 = event->genProcessID == 103; //lead dissociates
          bool SD2 = event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD1 && !SD2 && !DD && !CD;
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
	  //           int hf_n = event->HFtowers.size();
	  //           int hf_zero_count = ForwardRecord::nMaxHFMRecHits - hf_n;
	  //           double hf_double_energy_max = 0;
	  //           double hf_single_energy_max = 0;
	  //           vector<double> v_hf_rings_energy_max(24,0.);

	  //           double hf_m_energy_max = 0;
	  //           double hf_p_energy_max = 0;
	  //           double hf_p_energy = 0;
	  //           double hf_m_energy = 0;
	  //           int triggered_ring = 0;
	  //           int triggered_ring_m = 0;
	  //           int triggered_ring_p = 0;
	  //           for (vector<TowerHF>::const_iterator it = event->HFtowers.begin(); it < event->HFtowers.end(); ++it)
	  //             {
	  //               if(_SkipHFRings && (it->IetaAbs == 41 || it->IetaAbs == 29))
	  //                 continue;
	  //               const double eta = it->Eta;
	  //               const int Ieta = it->Eta > 0?it->IetaAbs:-it->IetaAbs;
	  //               double tower_e = it->Energy * _HFEnergyScale;
	  //               double c_lev = eta<0?c_lev_m[-Ieta-29]:c_lev_p[Ieta-29];
	  // #if _HFEnergyCalibration == 1
	  //               if(sample_type[sample]==MC)
	  //                 {
	  //                   //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
	  //                   tower_e /= (*hf_calibration)[IetaToRing(Ieta)];
	  //                 }
	  // #endif#if _HFEnergyCalibration == 2
	  //               if(sample_type[sample]==DATA)
	  //                 {
	  //                   //cout << "Ieta: " << Ieta << "tower e: " << tower_e << " -> " << tower_e/c_lev << endl;
	  //                   tower_e /= c_lev;//(*hf_calibration)[IetaToRing(Ieta)];
	  //                 }
	  // #endif
	  //               //cout << it->IetaAbs << " " << eta << " " << tower_e << endl;

	  //               if(eta < 0 && tower_e > v_hf_rings_energy_max[IetaToRing(Ieta)])
	  //                 v_hf_rings_energy_max[IetaToRing(Ieta)] = tower_e;
	  //               if(eta > 0 && tower_e > v_hf_rings_energy_max[IetaToRing(Ieta)])
	  //                 v_hf_rings_energy_max[IetaToRing(Ieta)] = tower_e;

	  //               if(eta > 0. && tower_e > hf_p_energy_max)
	  //                 {
	  //                   hf_p_energy_max = tower_e;
	  //                   triggered_ring_p = Ieta;
	  //                 }
	  //               if(eta <= 0. && tower_e > hf_m_energy_max)
	  //                 {
	  //                   hf_m_energy_max = tower_e;
	  //                   triggered_ring_m = Ieta;
	  //                 }

	  //               if(eta > 0.)
	  //                 hf_p_energy += tower_e;
	  //               else
	  //                 hf_m_energy += tower_e;

	  //               if(coll)               h_hf_hitdistr_coll->   Fill(tower_e);
	  //               if(coll && eta <= 0.)  h_hf_hitdistr_coll_m-> Fill(tower_e);
	  //               if(coll && eta > 0.)   h_hf_hitdistr_coll_p-> Fill(tower_e);
	  //               if(noise)              h_hf_hitdistr_noise->  Fill(tower_e);
	  //               if(noise && eta <= 0.) h_hf_hitdistr_noise_m->Fill(tower_e);
	  //               if(noise && eta > 0.)  h_hf_hitdistr_noise_p->Fill(tower_e);
	  //             }
	  //           double hf_pm_energy = hf_p_energy + hf_m_energy;
	  //           //cout << "rechit: " << hf_pm_energy << endl;

	  //           triggered_ring = hf_m_energy_max>hf_p_energy_max?triggered_ring_m:triggered_ring_p;

	  //           hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
	  //           hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
	  //           bool hf_single_tag = hf_single_energy_max >= 8;
	  //           bool hf_double_tag = hf_double_energy_max >= 4;



          //---------------------------------------------GEN Particles

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
                  double phi= it->GetPhi();

                  if (!(-6.6 < eta && eta < -5.2)) //castor acc
                    continue;

                  phiMassMap.insert(pair<double,GenParticle*>(phi,&(*it)));
                }

	      if(phiMassMap.size() == 0)
		continue;
	      else if(phiMassMap.size() == 1)
		{
		  multimap<double,GenParticle*>::iterator it = phiMassMap.begin();
		  if(abs(it->second->Id) == 2112) h_mc_castor_neutron_energy->Fill(it->second->GetEnergy());
		  if(abs(it->second->Id) == 2112) h_mc_castor_neutron_energy_isolated->Fill(it->second->GetEnergy());
		  if(abs(it->second->Id) ==  130) h_mc_castor_kaon_energy_isolated->Fill(it->second->GetEnergy());
		  if(abs(it->second->Id) ==  3122) h_mc_castor_lambda_energy_isolated->Fill(it->second->GetEnergy());
		}
	      else
		{
		  assert(phiMassMap.size() >= 2); //Assume 2 or more final usable particles

		  for (multimap<double,GenParticle*>::iterator it = phiMassMap.begin(); it != phiMassMap.end(); ++it)
		    {
		  
		      //if(SD2) cout << n++ << " " << it->second->Id << ": E=" << it->second->GetEnergy() << "  --  eta=" << it->second->GetEta() << "  --  y=" << it->second->GetPhi() << endl;

		      multimap<double,GenParticle*>::const_iterator itback;
		      multimap<double,GenParticle*>::const_iterator itforth;

		      if (it == phiMassMap.begin())
			{
			  itback=std::next(phiMassMap.end(),-1);
			  itforth=std::next(it,1);
			}
		      else if (it == std::next(phiMassMap.end(),-1))
			{
			  itback=std::next(it,-1);
			  itforth=phiMassMap.begin();
			}
		      else
			{
			  itback=std::next(it,-1);
			  itforth=std::next(it,1);
			}

		      const double angleDiffForth = fmod(fabs(it->second->GetPhi() - itforth->second->GetPhi()), TMath::Pi()/2);
		      const double angleDiffBack = fmod(fabs(it->second->GetPhi() - itback->second->GetPhi()), TMath::Pi()/2);
		      const double angleDiff = TMath::Min(angleDiffForth,angleDiffBack)/TMath::DegToRad();
		      const bool isolated = angleDiff > 8;
		      //if(abs(it->second->Id) == 2112 && angleDiff > 8.) cout << distance(phiMassMap.begin(),it) << " " << itback->second->GetPhi() << " "  << it->second->GetPhi() << " " << itforth->second->GetPhi() << " " << angleDiff << endl;
		      if(!TMath::Finite(angleDiff) || TMath::IsNaN(angleDiff))
			continue;

		      if(abs(it->second->Id) == 2112)             h_mc_castor_neutron_energy->Fill(it->second->GetEnergy(),evtWeight);
		      if(abs(it->second->Id) == 2112 && isolated) h_mc_castor_neutron_energy_isolated->Fill(it->second->GetEnergy(),evtWeight);
		      if(abs(it->second->Id) ==  130 && isolated) h_mc_castor_kaon_energy_isolated->Fill(it->second->GetEnergy(),evtWeight);
		      if(abs(it->second->Id) ==  3122 && isolated) h_mc_castor_lambda_energy_isolated->Fill(it->second->GetEnergy(),evtWeight);


		    }
		}
	    }


          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //---------------------------------------------Filling HISTOS
	}
	 

      //******************************************AFTER EVENT LOOP*******************************************

      cout << endl << "--- Processed " << n_total << " events." << endl << endl;
      if(sample_type[sample] == DATA)
        {
	  h_castor_neutron_energy;
        }


      if(sample_type[sample] == MC)
        {
          h_mc_castor_neutron_energy         ->Scale(1./n_total/h_mc_castor_neutron_energy->GetBinWidth(1));
	  h_mc_castor_neutron_energy_isolated->Scale(1./n_total/h_mc_castor_neutron_energy_isolated->GetBinWidth(1));
	  h_mc_castor_kaon_energy_isolated->Scale(1./n_total/h_mc_castor_kaon_energy_isolated->GetBinWidth(1));
	  h_mc_castor_lambda_energy_isolated->Scale(1./n_total/h_mc_castor_lambda_energy_isolated->GetBinWidth(1));
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
