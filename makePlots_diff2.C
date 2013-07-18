#include <TColor.h>
#include <iomanip>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TMarker.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"


#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void makePlots_diff2()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> type;
  type.push_back(string("single"));
  type.push_back(string("double"));

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  name.push_back(string("HIJING 1.383"));
  name.push_back(string("QGSJetII-04"));

  cout << " & SD [$\\%$] & DD [$\\%$] & CD [$\\%$] & ND [$\\%$] & Sum [$\\%$] & "
       << "Ratio_{MC} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "Ratio_{Data} $\\frac{\\text{double-arm}}{\\text{single-arm}}$"
       << "\\\\\\hline" << endl;


  for(int i=0; i<int(list.size()); i++)
    {
      cout << i+1 << "/" << int(list.size()) << endl;
      TFile* file = TFile::Open("histos_mc.root");
      file->cd();

      TH1D* sd1_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_SD1")).c_str());
      TH1D* sd2_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_SD2")).c_str());
      TH1D* dd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_DD")).c_str());
      TH1D* cd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_CD")).c_str());
      TH1D* nd_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_ND")).c_str());
      TH1D* all_single=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_single_all")).c_str());

      TH1D* sd1_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_SD1")).c_str());
      TH1D* sd2_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_SD2")).c_str());
      TH1D* dd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_DD")).c_str());
      TH1D* cd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_CD")).c_str());
      TH1D* nd_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_ND")).c_str());
      TH1D* all_double=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_e_double_all")).c_str());

      int nbins = all_single->GetNbinsX();

      int nbinsel_single = all_single->FindBin(8.);
      double n_sd1_single = sd1_single->Integral(0,nbins+1);
      double n_sd2_single = sd2_single->Integral(0,nbins+1);
      double n_dd_single = dd_single->Integral(0,nbins+1);
      double n_cd_single = cd_single->Integral(0,nbins+1);
      double n_nd_single = nd_single->Integral(0,nbins+1);
      double n_all_single = all_single->Integral(0,nbins+1);

      double n_sel_sd1_single = sd1_single->Integral(nbinsel_single,nbins+1);
      double n_sel_sd2_single = sd2_single->Integral(nbinsel_single,nbins+1);
      double n_sel_dd_single = dd_single->Integral(nbinsel_single,nbins+1);
      double n_sel_cd_single = cd_single->Integral(nbinsel_single,nbins+1);
      double n_sel_nd_single = nd_single->Integral(nbinsel_single,nbins+1);
      double n_sel_all_single = all_single->Integral(nbinsel_single,nbins+1);


      int nbinsel_double = all_double->FindBin(2.5);
      double n_sd1_double = sd1_double->Integral(0,nbins+1);
      double n_sd2_double = sd2_double->Integral(0,nbins+1);
      double n_dd_double = dd_double->Integral(0,nbins+1);
      double n_cd_double = cd_double->Integral(0,nbins+1);
      double n_nd_double = nd_double->Integral(0,nbins+1);
      double n_all_double = all_double->Integral(0,nbins+1);

      double n_sel_sd1_double = sd1_double->Integral(nbinsel_double,nbins+1);
      double n_sel_sd2_double = sd2_double->Integral(nbinsel_double,nbins+1);
      double n_sel_dd_double = dd_double->Integral(nbinsel_double,nbins+1);
      double n_sel_cd_double = cd_double->Integral(nbinsel_double,nbins+1);
      double n_sel_nd_double = nd_double->Integral(nbinsel_double,nbins+1);
      double n_sel_all_double = all_double->Integral(nbinsel_double,nbins+1);

      cout << fixed << setprecision(1) //SAME FOR BOTH: EXISTS ONCE
           << "No Selection & "
           << (n_sd1_single +n_sd2_single)/n_all_single*100. << " & "
           << n_dd_single/n_all_single*100. << " & "
           << n_cd_single/n_all_single*100. << " & "
           << n_nd_single/n_all_single*100. << " & "
           << "100" << "\\\\" << endl;

      cout << fixed << setprecision(1)
           << "$E_{\\text{HF}}>8$ GeV (single-arm) & "
           << (n_sel_sd1_single + n_sel_sd2_single)/n_all_single*100. << " & "
           << n_sel_dd_single/n_all_single*100. << " & "
           << n_sel_cd_single/n_all_single*100. << " & "
           << n_sel_nd_single/n_all_single*100. << " & "
           << n_sel_all_single/n_all_single*100.<< " & "
           << "\\multirow{2}{*}{" << (n_sel_all_double/n_all_double)/(n_sel_all_single/n_all_single)*100 << "}"
           << "\\multirow{2}{*}{" << 1.97752/2.01004*100 << "}"
           << "\\\\" << endl;
      cout << fixed << setprecision(1)
           << "$E_{\\text{HF}}>2.5$ GeV (double-arm) & "
           << (n_sel_sd1_double + n_sel_sd2_double)/n_all_double*100. << " & "
           << n_sel_dd_double/n_all_double*100. << " & "
           << n_sel_cd_double/n_all_double*100. << " & "
           << n_sel_nd_double/n_all_double*100. << " & "
           << n_sel_all_double/n_all_double*100. << " & "
           << " " //multicolumn
           << "\\\\\\hline" << endl;

      ///////////////DRAWING/////////
      for(int cur=0; cur<int(type.size()); cur++)
        {
          TH1D* sd1 = type[cur]=="single"?sd1_single:sd1_double;
          TH1D* sd2 = type[cur]=="single"?sd2_single:sd2_double;
          TH1D* dd = type[cur]=="single"?dd_single:dd_double;
          TH1D* cd = type[cur]=="single"?cd_single:cd_double;
          TH1D* nd = type[cur]=="single"?nd_single:nd_double;
          TH1D* all = type[cur]=="single"?all_single:all_double;

          sd1->Scale(1./double(all->GetEntries()));
          sd2->Scale(1./double(all->GetEntries()));
          dd->Scale(1./double(all->GetEntries()));
          cd->Scale(1./double(all->GetEntries()));
          nd->Scale(1./double(all->GetEntries()));
          all->Scale(1./double(all->GetEntries()));

          sd1->SetLineWidth(2.5);
          sd2->SetLineWidth(2.5);
          dd->SetLineWidth(2.5);
          cd->SetLineWidth(2.5);
          nd->SetLineWidth(2.5);
          all->SetLineWidth(2.5);

          cd->SetMarkerColor(TColor::GetColor(120,255,120));
          dd->SetMarkerColor(TColor::GetColor(255,90,90));
          sd2->SetMarkerColor(TColor::GetColor(60,60,255));
          sd1->SetMarkerColor(TColor::GetColor(00,00,150));
          nd->SetMarkerColor(kGray);

          sd1->SetLineColor(sd1->GetMarkerColor());
          sd2->SetLineColor(sd2->GetMarkerColor());
          dd->SetLineColor(dd->GetMarkerColor());
          cd->SetLineColor(cd->GetMarkerColor());
          nd->SetLineColor(nd->GetMarkerColor());

          sd1->SetFillColor(sd1->GetMarkerColor());
          sd2->SetFillColor(sd2->GetMarkerColor());
          dd->SetFillColor(dd->GetMarkerColor());
          cd->SetFillColor(cd->GetMarkerColor());
          nd->SetFillColor(nd->GetMarkerColor());

          ostringstream txt;
          txt.str(""); txt << "SD1 (Underflow: " << fixed << setprecision(1) << sd1->GetBinContent(0)*1e3 << "x10^{-3})"; sd1->SetTitle(txt.str().c_str());
          txt.str(""); txt << "SD2 (Underflow: " << fixed << setprecision(1) << sd2->GetBinContent(0)*1e3 << "x10^{-3})"; sd2->SetTitle(txt.str().c_str());
          txt.str(""); txt << "DD (Underflow: " << fixed << setprecision(1) << dd->GetBinContent(0)*1e3 << "x10^{-3})"; dd->SetTitle(txt.str().c_str());
          txt.str(""); txt << "CD (Underflow: " << fixed << setprecision(1) << cd->GetBinContent(0)*1e3 << "x10^{-3})"; cd->SetTitle(txt.str().c_str());
          txt.str(""); txt << "ND (Underflow: " << fixed << setprecision(1) << nd->GetBinContent(0)*1e3 << "x10^{-3})"; nd->SetTitle(txt.str().c_str());
          txt.str(""); txt << "all (Underflow: " << fixed << setprecision(1) << all->GetBinContent(0)*1e3 << "x10^{-3})"; all->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "SD1"; sd1->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "SD2"; sd2->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "DD"; dd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "CD"; cd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "ND"; nd->SetTitle(txt.str().c_str());
          // txt.str(""); txt << "all)"; all->SetTitle(txt.str().c_str());

          TCanvas* can1 = new TCanvas;
          THStack* hs = new THStack("hs","EPOS");
          hs->Add(sd1);
          hs->Add(sd2);
          hs->Add(dd);
          hs->Add(cd);
          hs->Add(nd);
          hs->Draw("HIST");

          hs->SetTitle(";E_{HF} [GeV];events (normalised)");
          hs->GetXaxis()->SetTitleOffset(hs->GetXaxis()->GetTitleOffset()*1.15);
          hs->SetMinimum(5e-5);
          hs->SetMaximum(10.); //off by factor 1.3???
          can1->SetLogy();
          can1->SetLogx();

          TLegend* leg1 = can1->BuildLegend(0.25,0.63,0.45,0.93);
          SetLegAtt(leg1);
          leg1->Draw();
          CMSText(0,0,1,name[i],type[cur]=="single"?"single-arm selection":"double-arm selection");

          TLine* line = new TLine(type[cur]=="single"?8:2.5,5e-5,type[cur]=="single"?8:2.5,1e-2);
          line->SetLineWidth(2);
          line->SetLineStyle(2);
          line->Draw("SAME");

          can1->SaveAs((string("plots/diff_energy_")+type[cur]+string("_")+list[i]+string(".pdf")).c_str());
        }
    }
}
