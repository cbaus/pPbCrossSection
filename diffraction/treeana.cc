#define _MAXEVT 5e6
#define _SkipHFRings 1 //skip 41 and 29 as suggested by HCAL DPG
#define _HFEnergyScale 1.0 //1.0 //0.8
#define _HFEnergyCalibration 2 //0 or 1 (rescale MC) or 2 this does not scale MC but data according to raddam from lev
#define _SimulateCastorNoise 1

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

#include "CastorTreeVariables.h"
#include "ParticleInfo.h"

#include "CastorCorrFactorpPb2013.h"

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
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data210998/*_*.root"); sample_name.push_back("data210998"); sample_type.push_back(DATA);
  //sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211000/*.root"); sample_name.push_back("data211000"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211001/*.root"); sample_name.push_back("data211001"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211032/*.root"); sample_name.push_back("data211032"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211256/*.root"); sample_name.push_back("data211256"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211371/*.root"); sample_name.push_back("data211371"); sample_type.push_back(DATA);
  //sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211390/*.root"); sample_name.push_back("data211390"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211460/*.root"); sample_name.push_back("data211460"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211538/*.root"); sample_name.push_back("data211538"); sample_type.push_back(DATA);
  //sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data211607/*_*.root"); sample_name.push_back("data211607"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos_SL/*.root"); sample_name.push_back("Epos_SL"); sample_type.push_back(MC);


//_______________________ THESE ARE STUDIED. FORGET ABOVE
  // sample_fname.push_back("root://eoscms//eos/cms/store/caf/user/cbaus/pPb2013/trees/Data211532/*.root"); sample_name.push_back("data211532"); sample_type.push_back(DATA);
  // sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Data210885/*_*.root"); sample_name.push_back("data210885"); sample_type.push_back(DATA);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Epos/*.root"); sample_name.push_back("Epos"); sample_type.push_back(MC);
  sample_fname.push_back("/afs/cern.ch/work/c/cbaus/public/castortree/pPb_QGSJetII/treeMC.root"); sample_name.push_back("QGSJetII"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/Hijing/*.root"); sample_name.push_back("Hijing"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/DPMJet/treeMC.root"); sample_name.push_back("DPMJet"); sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightDPMjet_v2/treeMC.root"); sample_name.push_back("Starlight_DPMJet");  sample_type.push_back(MC);
  sample_fname.push_back("root://eoscms//eos/cms/store/group/phys_heavyions/cbaus/trees/StarlightPythia/treeMC.root"); sample_name.push_back("Starlight_Pythia");  sample_type.push_back(MC);


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

  TFile noise_file("histos_noise.root","READ");
  TH1D* h_noise_input = (TH1D*)(noise_file.Get("data210885/cas_with_rapgap_noise_bxexcl"));
  //noise_file.Close();


  TH1D* h_mc_rapidity;
  TH1D* h_mc_eta_e;
  TH1D* h_mc_eta_e_SD1;
  TH1D* h_mc_eta_e_SD2;
  TH1D* h_mc_eta_e_CD;
  TH1D* h_mc_eta_e_DD;
  TH1D* h_mc_eta_e_ND;

  TH1D* h_mc_cas_plus_with_rapgap_dy0;
  TH1D* h_mc_cas_minus_with_rapgap_dy0;
  TH1D* h_mc_cas_plus_with_rapgap_dy5;
  TH1D* h_mc_cas_minus_with_rapgap_dy5;
  TH1D* h_mc_cas_plus_with_rapgap_dy10;
  TH1D* h_mc_cas_minus_with_rapgap_dy10;

  TH1D* h_gen_cas_plus_with_rapgap_dy0;
  TH1D* h_gen_cas_minus_with_rapgap_dy0;
  TH1D* h_gen_cas_plus_with_rapgap_dy1;
  TH1D* h_gen_cas_minus_with_rapgap_dy1;
  TH1D* h_gen_cas_plus_with_rapgap_dy2;
  TH1D* h_gen_cas_minus_with_rapgap_dy2;
  TH1D* h_gen_cas_plus_with_rapgap_dy3;
  TH1D* h_gen_cas_minus_with_rapgap_dy3;

  TH1D* h_cas_with_rapgap_zb;
  TH1D* h_cas_with_rapgap_noise;
  TH1D* h_cas_with_rapgap_noise_bxexcl;
  TH1D* h_cas_with_rapgap_noise_unpaired;
  TH1D* h_cas_with_noise_bx;
  TH1D* h_cas_with_noise_bx8;
  TH1D* h_cas_with_noise_deltabx;
  TH1D* h_hf_with_noise_bx;
  TH1D* h_hf_with_noise_deltabx;
  TH1D* h_cas_with_rapgap_dy0;
  TH1D* h_cas_with_rapgap_dy1;
  TH1D* h_cas_with_rapgap_dy2;
  TH1D* h_cas_with_rapgap_dy3;

  TH1D* h_cas_z_deltabx_m10;
  TH1D* h_cas_z_deltabx_m2;
  TH1D* h_cas_z_deltabx_m1;
  TH1D* h_cas_z_deltabx_p0;
  TH1D* h_cas_z_deltabx_p1;

  //****************************************************************sth for noise*********************************************************
  set<int> unpaired,paired; //true at least for 210885 and 211532
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

  paired.insert(1);
  paired.insert(9);
  paired.insert(18);
  paired.insert(26);
  paired.insert(35);
  paired.insert(43);
  paired.insert(52);
  paired.insert(60);
  paired.insert(69);
  paired.insert(77);
  paired.insert(86);
  paired.insert(94);
  paired.insert(103);
  paired.insert(111);
  paired.insert(120);
  paired.insert(128);
  paired.insert(137);
  paired.insert(145);
  paired.insert(154);
  paired.insert(162);
  paired.insert(171);
  paired.insert(179);
  paired.insert(188);
  paired.insert(196);
  paired.insert(233);
  paired.insert(241);
  paired.insert(250);
  paired.insert(258);
  paired.insert(267);
  paired.insert(275);
  paired.insert(284);
  paired.insert(292);
  paired.insert(301);
  paired.insert(309);
  paired.insert(318);
  paired.insert(326);
  paired.insert(335);
  paired.insert(343);
  paired.insert(352);
  paired.insert(360);
  paired.insert(369);
  paired.insert(377);
  paired.insert(386);
  paired.insert(394);
  paired.insert(403);
  paired.insert(411);
  paired.insert(420);
  paired.insert(428);
  paired.insert(892);
  paired.insert(900);
  paired.insert(909);
  paired.insert(917);
  paired.insert(926);
  paired.insert(934);
  paired.insert(943);
  paired.insert(951);
  paired.insert(960);
  paired.insert(968);
  paired.insert(977);
  paired.insert(985);
  paired.insert(994);
  paired.insert(1002);
  paired.insert(1011);
  paired.insert(1019);
  paired.insert(1028);
  paired.insert(1036);
  paired.insert(1045);
  paired.insert(1053);
  paired.insert(1062);
  paired.insert(1070);
  paired.insert(1079);
  paired.insert(1087);
  paired.insert(1124);
  paired.insert(1132);
  paired.insert(1141);
  paired.insert(1149);
  paired.insert(1158);
  paired.insert(1166);
  paired.insert(1175);
  paired.insert(1183);
  paired.insert(1192);
  paired.insert(1200);
  paired.insert(1209);
  paired.insert(1217);
  paired.insert(1226);
  paired.insert(1234);
  paired.insert(1243);
  paired.insert(1251);
  paired.insert(1260);
  paired.insert(1268);
  paired.insert(1277);
  paired.insert(1285);
  paired.insert(1294);
  paired.insert(1302);
  paired.insert(1311);
  paired.insert(1319);
  paired.insert(1356);
  paired.insert(1364);
  paired.insert(1373);
  paired.insert(1381);
  paired.insert(1390);
  paired.insert(1398);
  paired.insert(1407);
  paired.insert(1415);
  paired.insert(1424);
  paired.insert(1432);
  paired.insert(1441);
  paired.insert(1449);
  paired.insert(1458);
  paired.insert(1466);
  paired.insert(1475);
  paired.insert(1483);
  paired.insert(1492);
  paired.insert(1500);
  paired.insert(1509);
  paired.insert(1517);
  paired.insert(1526);
  paired.insert(1534);
  paired.insert(1543);
  paired.insert(1551);
  paired.insert(1588);
  paired.insert(1596);
  paired.insert(1605);
  paired.insert(1613);
  paired.insert(1622);
  paired.insert(1630);
  paired.insert(1639);
  paired.insert(1647);
  paired.insert(1656);
  paired.insert(1664);
  paired.insert(1673);
  paired.insert(1681);
  paired.insert(1690);
  paired.insert(1698);
  paired.insert(1707);
  paired.insert(1715);
  paired.insert(1724);
  paired.insert(1732);
  paired.insert(1783);
  paired.insert(1791);
  paired.insert(1800);
  paired.insert(1808);
  paired.insert(1817);
  paired.insert(1825);
  paired.insert(1834);
  paired.insert(1842);
  paired.insert(1851);
  paired.insert(1859);
  paired.insert(1868);
  paired.insert(1876);
  paired.insert(1885);
  paired.insert(1893);
  paired.insert(1902);
  paired.insert(1910);
  paired.insert(1919);
  paired.insert(1927);
  paired.insert(1936);
  paired.insert(1944);
  paired.insert(1953);
  paired.insert(1961);
  paired.insert(1970);
  paired.insert(1978);
  paired.insert(2015);
  paired.insert(2023);
  paired.insert(2032);
  paired.insert(2040);
  paired.insert(2049);
  paired.insert(2057);
  paired.insert(2066);
  paired.insert(2074);
  paired.insert(2083);
  paired.insert(2091);
  paired.insert(2100);
  paired.insert(2108);
  paired.insert(2117);
  paired.insert(2125);
  paired.insert(2134);
  paired.insert(2142);
  paired.insert(2151);
  paired.insert(2159);
  paired.insert(2168);
  paired.insert(2176);
  paired.insert(2185);
  paired.insert(2193);
  paired.insert(2202);
  paired.insert(2210);
  paired.insert(2247);
  paired.insert(2255);
  paired.insert(2264);
  paired.insert(2272);
  paired.insert(2281);
  paired.insert(2289);
  paired.insert(2298);
  paired.insert(2306);
  paired.insert(2315);
  paired.insert(2323);
  paired.insert(2332);
  paired.insert(2340);
  paired.insert(2349);
  paired.insert(2357);
  paired.insert(2366);
  paired.insert(2374);
  paired.insert(2383);
  paired.insert(2391);
  paired.insert(2400);
  paired.insert(2408);
  paired.insert(2417);
  paired.insert(2425);
  paired.insert(2434);
  paired.insert(2442);
  paired.insert(2479);
  paired.insert(2487);
  paired.insert(2496);
  paired.insert(2504);
  paired.insert(2513);
  paired.insert(2521);
  paired.insert(2530);
  paired.insert(2538);
  paired.insert(2547);
  paired.insert(2555);
  paired.insert(2564);
  paired.insert(2572);
  paired.insert(2581);
  paired.insert(2589);
  paired.insert(2598);
  paired.insert(2606);
  paired.insert(2615);
  paired.insert(2623);
  paired.insert(2674);
  paired.insert(2682);
  paired.insert(2691);
  paired.insert(2699);
  paired.insert(2708);
  paired.insert(2716);
  paired.insert(2725);
  paired.insert(2733);
  paired.insert(2742);
  paired.insert(2750);
  paired.insert(2759);
  paired.insert(2767);
  paired.insert(2776);
  paired.insert(2784);
  paired.insert(2793);
  paired.insert(2801);
  paired.insert(2810);
  paired.insert(2818);
  paired.insert(2827);
  paired.insert(2835);
  paired.insert(2844);
  paired.insert(2852);
  paired.insert(2861);
  paired.insert(2869);
  paired.insert(2906);
  paired.insert(2914);
  paired.insert(2923);
  paired.insert(2931);
  paired.insert(2940);
  paired.insert(2948);
  paired.insert(2957);
  paired.insert(2965);
  paired.insert(2974);
  paired.insert(2982);
  paired.insert(2991);
  paired.insert(2999);
  paired.insert(3008);
  paired.insert(3016);
  paired.insert(3025);
  paired.insert(3033);
  paired.insert(3042);
  paired.insert(3050);
  paired.insert(3059);
  paired.insert(3067);
  paired.insert(3104);
  paired.insert(3112);
  paired.insert(3121);
  paired.insert(3129);
  paired.insert(3138);
  paired.insert(3146);
  paired.insert(3155);
  paired.insert(3163);
  paired.insert(3172);
  paired.insert(3180);
  paired.insert(3189);
  paired.insert(3197);
  paired.insert(3206);
  paired.insert(3214);
  paired.insert(3223);
  paired.insert(3231);
  paired.insert(3240);
  paired.insert(3248);
  paired.insert(3257);
  paired.insert(3265);
  paired.insert(3274);
  paired.insert(3282);
  paired.insert(3291);
  paired.insert(3299);



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

      out_file->mkdir(sample_name[sample].c_str());
      out_file->cd(sample_name[sample].c_str());

      if(sample_type[sample] == MC)
        {
          h_mc_rapidity           = new TH1D((string("mc_rapidity")).c_str(),"",40,-12,12);
          h_mc_eta_e              = new TH1D((string("mc_eta_e")).c_str(),"",40,-12,12);
          h_mc_eta_e_SD1          = new TH1D((string("mc_eta_e_SD1")).c_str(),"",40,-12,12);
          h_mc_eta_e_SD2          = new TH1D((string("mc_eta_e_SD2")).c_str(),"",40,-12,12);
          h_mc_eta_e_CD           = new TH1D((string("mc_eta_e_CD")).c_str(),"",40,-12,12);
          h_mc_eta_e_DD           = new TH1D((string("mc_eta_e_DD")).c_str(),"",40,-12,12);
          h_mc_eta_e_ND           = new TH1D((string("mc_eta_e_ND")).c_str(),"",40,-12,12);

          h_mc_cas_plus_with_rapgap_dy0 = new TH1D((string("mc_cas_plus_with_rapgap_dy0")).c_str(),"",40,log10(10),log10(12000));
          h_mc_cas_minus_with_rapgap_dy0 = new TH1D((string("mc_cas_minus_with_rapgap_dy0")).c_str(),"",40,log10(10),log10(4000));
          h_mc_cas_plus_with_rapgap_dy5 = new TH1D((string("mc_cas_plus_with_rapgap_dy5")).c_str(),"",40,log10(10),log10(12000));
          h_mc_cas_minus_with_rapgap_dy5 = new TH1D((string("mc_cas_minus_with_rapgap_dy5")).c_str(),"",40,log10(10),log10(4000));
          h_mc_cas_plus_with_rapgap_dy10 = new TH1D((string("mc_cas_plus_with_rapgap_dy10")).c_str(),"",40,log10(10),log10(12000));
          h_mc_cas_minus_with_rapgap_dy10 = new TH1D((string("mc_cas_minus_with_rapgap_dy10")).c_str(),"",40,log10(10),log10(4000));
          BinLogX(h_mc_cas_plus_with_rapgap_dy0);
          BinLogX(h_mc_cas_minus_with_rapgap_dy0);
          BinLogX(h_mc_cas_plus_with_rapgap_dy5);
          BinLogX(h_mc_cas_minus_with_rapgap_dy5);
          BinLogX(h_mc_cas_plus_with_rapgap_dy10);
          BinLogX(h_mc_cas_minus_with_rapgap_dy10);

          h_gen_cas_plus_with_rapgap_dy0 = new TH1D((string("gen_cas_plus_with_rapgap_dy0")).c_str(),"",40,log10(10),log10(12000));
          h_gen_cas_minus_with_rapgap_dy0 = new TH1D((string("gen_cas_minus_with_rapgap_dy0")).c_str(),"",40,log10(10),log10(4000));
          h_gen_cas_plus_with_rapgap_dy1 = new TH1D((string("gen_cas_plus_with_rapgap_dy1")).c_str(),"",40,log10(10),log10(12000));
          h_gen_cas_minus_with_rapgap_dy1 = new TH1D((string("gen_cas_minus_with_rapgap_dy1")).c_str(),"",40,log10(10),log10(4000));
          h_gen_cas_plus_with_rapgap_dy2 = new TH1D((string("gen_cas_plus_with_rapgap_dy2")).c_str(),"",40,log10(10),log10(12000));
          h_gen_cas_minus_with_rapgap_dy2 = new TH1D((string("gen_cas_minus_with_rapgap_dy2")).c_str(),"",40,log10(10),log10(4000));
          h_gen_cas_plus_with_rapgap_dy3 = new TH1D((string("gen_cas_plus_with_rapgap_dy3")).c_str(),"",40,log10(10),log10(12000));
          h_gen_cas_minus_with_rapgap_dy3 = new TH1D((string("gen_cas_minus_with_rapgap_dy3")).c_str(),"",40,log10(10),log10(4000));
          BinLogX(h_gen_cas_plus_with_rapgap_dy0);
          BinLogX(h_gen_cas_minus_with_rapgap_dy0);
          BinLogX(h_gen_cas_plus_with_rapgap_dy1);
          BinLogX(h_gen_cas_minus_with_rapgap_dy1);
          BinLogX(h_gen_cas_plus_with_rapgap_dy2);
          BinLogX(h_gen_cas_minus_with_rapgap_dy2);
          BinLogX(h_gen_cas_plus_with_rapgap_dy3);
          BinLogX(h_gen_cas_minus_with_rapgap_dy3);
        }

      h_cas_with_rapgap_dy0 = new TH1D((string("cas_with_rapgap_dy0")).c_str(),"",40,log10(0.5),log10(7000));
      h_cas_with_rapgap_dy1 = new TH1D((string("cas_with_rapgap_dy1")).c_str(),"",40,log10(0.5),log10(7000));
      h_cas_with_rapgap_dy2 = new TH1D((string("cas_with_rapgap_dy2")).c_str(),"",40,log10(0.5),log10(7000));
      h_cas_with_rapgap_dy3 = new TH1D((string("cas_with_rapgap_dy3")).c_str(),"",40,log10(0.5),log10(7000));
      BinLogX(h_cas_with_rapgap_dy0);
      BinLogX(h_cas_with_rapgap_dy1);
      BinLogX(h_cas_with_rapgap_dy2);
      BinLogX(h_cas_with_rapgap_dy3);


      h_cas_with_rapgap_zb = new TH1D((string("cas_with_rapgap_zb")).c_str(),"",40,log10(0.5),log10(7000));
      BinLogX(h_cas_with_rapgap_zb);
      if(sample_type[sample] == DATA)
        {
	  h_cas_with_rapgap_noise = new TH1D((string("cas_with_rapgap_noise")).c_str(),"",40,log10(0.5),log10(7000));
	  h_cas_with_rapgap_noise_bxexcl = new TH1D((string("cas_with_rapgap_noise_bxexcl")).c_str(),"",40,log10(0.5),log10(7000));
	  h_cas_with_rapgap_noise_unpaired = new TH1D((string("cas_with_rapgap_noise_unpaired")).c_str(),"",40,log10(0.5),log10(7000));
	  h_cas_with_noise_bx = new TProfile((string("cas_with_noise_bx")).c_str(),"",3701,-0.5,3700.5);
          h_cas_with_noise_bx8 = new TH1D((string("cas_with_noise_bx8")).c_str(),"",100,-10,90);
	  h_hf_with_noise_bx  = new TProfile((string("hf_with_noise_bx")).c_str(),"",3701,-0.5,3700.5);
	  h_cas_with_noise_deltabx = new TProfile((string("cas_with_noise_deltabx")).c_str(),"",201,-100.5,100.5);
	  h_hf_with_noise_deltabx  = new TProfile((string("hf_with_noise_detlabx")).c_str(),"",201,-100.5,100.5);
          BinLogX(h_cas_with_rapgap_noise);
          BinLogX(h_cas_with_rapgap_noise_bxexcl);
          BinLogX(h_cas_with_rapgap_noise_unpaired);
	}
      h_cas_z_deltabx_m10 = new TProfile((string("cas_z_deltabx_m10")).c_str(),"",16,-0.5,15.5);
      h_cas_z_deltabx_m2 = new TProfile((string("cas_z_deltabx_m2")).c_str(),"",16,-0.5,15.5);
      h_cas_z_deltabx_m1 = new TProfile((string("cas_z_deltabx_m1")).c_str(),"",16,-0.5,15.5);
      h_cas_z_deltabx_p0 = new TProfile((string("cas_z_deltabx_p0")).c_str(),"",16,-0.5,15.5);
      h_cas_z_deltabx_p1 = new TProfile((string("cas_z_deltabx_p1")).c_str(),"",16,-0.5,15.5);

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

      ostringstream run_str;
      run_str << event->runNb;

      double n_total = double(tree->GetEntries());
      int maxevts = sample_type[sample] == MC ? _MAXEVT/10 : _MAXEVT;
      if(maxevts<n_total && maxevts>0)
        n_total = double(maxevts);


      for(int iEvent=0; iEvent<n_total; iEvent++)
        {
          if(iEvent % 10000 == 0) cout << sample+1 << "/" << sample_name.size() << " -- " << sample_name[sample].c_str() << " -- Entry: " << iEvent << " / " << n_total << endl;
          tree->GetEntry(iEvent);
          bool coll=0, noise=0;

          if(sample_type[sample] == DATA)
            {
              coll          = zero_bias && bptx_p_m; //double beam
              noise         = random;// && bptx_quiet;// && !bptx_np_m && !bptx_p_nm; //not both and not single beam
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
                  RecHitGeV *= castor::channelGainQE[sec][mod];
                  RecHitGeV *= castor::absEscaleFactor;
                }
              RecHitGeV *= (int)castor::channelQuality[sec][mod];

              sum_CAS_E_mod[it->GetModuleId()] += RecHitGeV;

              if(it->GetModuleId() < 2)
                sum_CAS_E_em += RecHitGeV;
              else if(it->GetModuleId() < 5)
                sum_CAS_E_had += RecHitGeV;
            }

          double sum_CAS_E = sum_CAS_E_had + sum_CAS_E_em;
          if(sample_type[sample] == MC)
            sum_CAS_E += h_noise_input->GetRandom();


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
            }

          double hf_pm_energy = hf_p_energy + hf_m_energy;
          //cout << "rechit: " << hf_pm_energy << endl;

          triggered_ring = hf_m_energy_max>hf_p_energy_max?triggered_ring_m:triggered_ring_p;

          hf_double_energy_max = TMath::Min(hf_m_energy_max,hf_p_energy_max);
          hf_single_energy_max = TMath::Max(hf_m_energy_max,hf_p_energy_max);
	  double sum_HF_E = hf_p_energy + hf_m_energy;
          bool hf_single_tag = hf_single_energy_max >= 8;
          bool hf_double_tag = hf_double_energy_max >= 4;



          //---------------------------------------------GEN Particles
          const double s = 5020*5020;
          double m_m=0, m_x=0, xi_x=0;
          double m_p=0, m_y=0, xi_y=0;
          double rapGapPlus=-1;
          double rapGapMinus=-1;
          double rapGap=-1;
          double ymin = -5.2;
          double ymax = 5.2;
          double gen_HF_p_double=0;
          double gen_HF_p_single=0;
          bool SD1 = event->genProcessID == 103; //lead dissociates
          bool SD2 = event->genProcessID == 104;
          bool DD = event->genProcessID == 105;
          bool CD = event->genProcessID == 102;
          bool ND = !SD1 && !SD2 && !DD && !CD;
          double gen_CAS_E_plus  = 0;
          double gen_CAS_E_minus = 0;
          double gen_HF_p_minus = -1;
          double gen_HF_p_plus  = -1;
          double gen_central_tracks  = 0;
          int nTracks;

          //Counting events
          const int prescale        = zero_bias_prescale_L1*zero_bias_prescale_HLT;
          const double lumiPerLS    = event->instLuminosity * event->instLuminosityCorr * 1e6;
          double evtWeight    = 1;

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

                  if (abs(it->Id) > 1e9) //skip fragments
                    continue;

                  if(abs(it->Id) == 12 || abs(it->Id) == 14 || abs(it->Id) == 16 || abs(it->Id) == 13)
                    continue; // skip muon + neutrinos

                  const double eta=it->GetEta();
                  const double energy=it->GetEnergy();
                  const double pt=it->GetTransverseMomentum();
                  const double p=it->GetMomentum();
                  const double charge=it->GetCharge();

                  if(5.2 < eta && eta < 6.6)
                    gen_CAS_E_plus += energy;
                  if(-6.6 < eta && eta < -5.2)
                    gen_CAS_E_minus += energy;

                  if(-5 < eta && eta < -3 && p > gen_HF_p_minus)
                    gen_HF_p_minus = p;
                  if(+3 < eta && eta < +5 && p > gen_HF_p_plus)
                    gen_HF_p_plus = p;
                  if(fabs(eta) < 2.5 && pt > 0.2 && charge != 0)
                    gen_central_tracks++;

                  double Rapidity= eta;//it->GetRapidity();
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

		  if (eta > ymin && eta < ymax)
		    {
		      if (pt > 0.2)  rapidityMassMap.insert(pair<double,GenParticle*>(eta,&(*it)));
		    }
		}
	      gen_HF_p_single = TMath::Max(gen_HF_p_minus,gen_HF_p_plus);
	      gen_HF_p_double = TMath::Min(gen_HF_p_minus,gen_HF_p_plus);

              if(rapidityMassMap.size() >= 1) //Assume 1 or more final usable particles
                {
		  multimap<double,GenParticle*>::const_iterator thisIsIt; //start of m_x
		  int n=0;

		  rapGapMinus = fabs(rapidityMassMap.begin()->second->GetEta() - ymin);
		  multimap<double,GenParticle*>::const_iterator it = rapidityMassMap.end();
		  --it;
		  rapGapPlus = fabs(it->second->GetEta() - ymax);
		  rapGap = TMath::Max(rapGapPlus,rapGapMinus);
		}
            }

	  //calculating delta bx
	  int deltabx = 1e9;
	  int evtBx = event->bxNb;
	  bool isUnpaired = unpaired.count(evtBx);
	  bool isPaired = paired.count(evtBx);
	  bool isNoBeam = bptx_quiet;

          if (isNoBeam && isPaired)
            cerr << "Warning: Paired and no beam:" << evtBx << endl;
          if (isNoBeam && isUnpaired)
            cerr << "Warning: Unpaired and no beam:" << evtBx << endl;

	  if(noise)
	    {
	      int prev = -1;
	      int curr = -1;
	      for (auto &bx : paired)
		{
		  curr = bx;
		  //cout << iEvent << " " << evtBx << " " << curr << endl;
		  if(abs(evtBx-prev) < abs(evtBx-curr)) //check if previous from list was closer
		    break;
		  else
		    prev = curr;
		}
	      deltabx = evtBx - prev; //bx of event - closest filled bx
	      if (deltabx != 0 && bptx_quiet==0 && !isUnpaired)  cout << "delta=" << deltabx << " bptx:" << bptx_quiet << " unpaired:" << isUnpaired << endl;
	    }

          nTracks =  event->Tracks.size();

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

	  if(evtWeight <= 0)
	    cerr << "Event weight = " << evtWeight << endl;
          //cout << prescale << " " << event->instLuminosity << " " <<  event->instLuminosityCorr << endl;

          //---------------------------------------------Filling HISTOS
	  if(sample_type[sample] == MC)
            {
              if (rapGapPlus == -1.)
                rapGapPlus = 11; //WARNING
              if (rapGapPlus >= 00.) h_mc_cas_plus_with_rapgap_dy0->Fill(gen_CAS_E_plus,evtWeight);
              if (rapGapPlus >= 00.) h_mc_cas_minus_with_rapgap_dy0->Fill(gen_CAS_E_minus,evtWeight);
              if (rapGapPlus >= 05.) h_mc_cas_plus_with_rapgap_dy5->Fill(gen_CAS_E_plus,evtWeight);
              if (rapGapPlus >= 05.) h_mc_cas_minus_with_rapgap_dy5->Fill(gen_CAS_E_minus,evtWeight);
              if (rapGapPlus >= 10.) h_mc_cas_plus_with_rapgap_dy10->Fill(gen_CAS_E_plus,evtWeight);
              if (rapGapPlus >= 10.) h_mc_cas_minus_with_rapgap_dy10->Fill(gen_CAS_E_minus,evtWeight);
            }
          bool trig_hf_p, trig_track, trig_hf_m;
          trig_hf_p = hf_p_energy_max > 8.;
          trig_hf_m = hf_m_energy_max > 8.;
          trig_track = nTracks >= 1;

          bool gen_trig_hf_p, gen_trig_track, gen_trig_hf_m;
          gen_trig_hf_p = gen_HF_p_plus > 21.3;
          gen_trig_hf_m = gen_HF_p_minus > 21.3;
          gen_trig_track = gen_central_tracks >= 1;


          if (coll &&  trig_hf_p &&  trig_track &&  trig_hf_m) h_cas_with_rapgap_dy0->Fill(sum_CAS_E,evtWeight);
          if (coll &&  trig_hf_p &&  trig_track && !trig_hf_m) h_cas_with_rapgap_dy1->Fill(sum_CAS_E,evtWeight);
          if (coll &&  trig_hf_p && !trig_track && !trig_hf_m) h_cas_with_rapgap_dy2->Fill(sum_CAS_E,evtWeight);
          if (coll && !trig_hf_p && !trig_track && !trig_hf_m) h_cas_with_rapgap_dy3->Fill(sum_CAS_E,evtWeight);

          if(sample_type[sample] == MC)
            {
              if (gen_trig_hf_p &&  gen_trig_track &&  gen_trig_hf_m) h_gen_cas_minus_with_rapgap_dy0->Fill(gen_CAS_E_minus,evtWeight);
              if (gen_trig_hf_p &&  gen_trig_track && !gen_trig_hf_m) h_gen_cas_minus_with_rapgap_dy1->Fill(gen_CAS_E_minus,evtWeight);
              if (gen_trig_hf_p && !gen_trig_track && !gen_trig_hf_m) h_gen_cas_minus_with_rapgap_dy2->Fill(gen_CAS_E_minus,evtWeight);
              if (!gen_trig_hf_p && !gen_trig_track && !gen_trig_hf_m) h_gen_cas_minus_with_rapgap_dy3->Fill(gen_CAS_E_minus,evtWeight);

              if (gen_trig_hf_p &&  gen_trig_track &&  gen_trig_hf_m) h_gen_cas_plus_with_rapgap_dy0->Fill(gen_CAS_E_plus,evtWeight);
              if (gen_trig_hf_p &&  gen_trig_track && !gen_trig_hf_m) h_gen_cas_plus_with_rapgap_dy1->Fill(gen_CAS_E_plus,evtWeight);
              if (gen_trig_hf_p && !gen_trig_track && !gen_trig_hf_m) h_gen_cas_plus_with_rapgap_dy2->Fill(gen_CAS_E_plus,evtWeight);
              if (!gen_trig_hf_p && !gen_trig_track && !gen_trig_hf_m) h_gen_cas_plus_with_rapgap_dy3->Fill(gen_CAS_E_plus,evtWeight);
            }

          if (coll) h_cas_with_rapgap_zb->Fill(sum_CAS_E);
          if (noise && isNoBeam) h_cas_with_rapgap_noise->Fill(sum_CAS_E);
          if (noise && isNoBeam && (deltabx > 15 || deltabx < -5)) h_cas_with_rapgap_noise_bxexcl->Fill(sum_CAS_E);
          if (noise && isUnpaired) h_cas_with_rapgap_noise_unpaired->Fill(sum_CAS_E);
          if (noise) h_cas_with_noise_bx->Fill(evtBx,sum_CAS_E);
          if (noise && evtBx == 8) h_cas_with_noise_bx8->Fill(evtBx);
          if (noise) h_cas_with_noise_deltabx->Fill(deltabx,sum_CAS_E);
          if (noise) h_hf_with_noise_bx->Fill(evtBx,sum_HF_E);
          if (noise) h_hf_with_noise_deltabx->Fill(deltabx,sum_HF_E);

          if (noise && isNoBeam && trig_track)
            cout << "Warning fake reconstructed tracks. Event " << iEvent << endl;

          if (noise)
            {
              TH1D* tofill = 0;
              if(deltabx < -2 && deltabx > -11)
                tofill = h_cas_z_deltabx_m10;
              if(deltabx == -2)
                tofill = h_cas_z_deltabx_m2;
              if(deltabx == -1)
                tofill = h_cas_z_deltabx_m1;
              if(deltabx == 0)
                tofill = h_cas_z_deltabx_p0;
              if(deltabx == 1)
                tofill = h_cas_z_deltabx_p1;
              if(tofill)
                for(int k=0; k<int(sum_CAS_E_mod.size()); ++k)
                  tofill->Fill(k,sum_CAS_E_mod[k]);
            }
        }

      //******************************************AFTER EVENT LOOP*******************************************

      cout << endl << "--- Processed " << n_total << " events." << endl << endl;

      if(sample_type[sample] == DATA)
        {
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
