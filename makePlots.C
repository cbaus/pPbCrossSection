#include <iomanip>

TVectorD corr_fac_em(2);
TVectorD corr_fac_noise(2);
TVectorD corr_fac_mc(2);
TVectorD corr_fac_all(2);

const double fac_error = 1.0; 
const double cut_value_single = 6.*fac_error;
const double cut_value_double = 2.*fac_error;

void makePlots()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> type; type.push_back("single"); type.push_back("double");

  for(int n=0; n<int(type.size()); n++)
    {

  TFile* file = TFile::Open("histos_old.root");
  TFile* file2 = TFile::Open("histos_mc.root");
  TFile* file3 = TFile::Open("plots/hf_cuts_noise.root");

  TH1D* a=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n]).c_str());
  TH1D* a2=(TH1D*)file->Get((string("data210885/data210885_h_hf_cut_") + type[n] + string("_noise")).c_str());
  TH1D* b=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_cut_") + type[n]).c_str());
  TH1D* b2=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_new_cut_")+ type[n]).c_str());
  TH1D* c=(TH1D*)file->Get((string("Epos/Epos_h_hf_cut_")+ type[n]).c_str());
  TH1D* c2=(TH1D*)file->Get((string("Epos/Epos_h_hf_new_cut_")+ type[n]).c_str());
  TH1D* d=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_cut_") + type[n]).c_str());
  TH1D* d2=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_new_cut_")+ type[n]).c_str());
  TH1D* e=(TH1D*)file->Get((string("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_") + type[n]).c_str());
  TH1D* f=(TH1D*)file->Get((string("Starlight_Pythia/Starlight_Pythia_h_hf_cut_") + type[n]).c_str());
  TVectorD* hf_cuts_equivalent = (TVectorD*)file3->Get("hf_cuts_equivalent");

  a->Scale(1./double(a->GetBinContent(1)));
  a2->Scale(1./double(a2->GetBinContent(1)));
  b->Scale(1./double(b->GetBinContent(1)));
  b2->Scale(1./double(b2->GetBinContent(1)));
  c->Scale(1./double(c->GetBinContent(1)));
  c2->Scale(1./double(c2->GetBinContent(1)));
  d->Scale(1./double(d->GetBinContent(1)));
  d2->Scale(1./double(d2->GetBinContent(1)));
  e->Scale(1./double(e->GetBinContent(1)));
  cerr << "warning" << endl;
//   e->Scale(1./double(e->GetBinContent(1))/ 2100. * 195.);
  f->Scale(1./double(f->GetBinContent(1))/ 195. * 122.);

  a->SetLineWidth(3);
  a2->SetLineWidth(3);
  b->SetLineWidth(3);
  c->SetLineWidth(3);
  d->SetLineWidth(3);
  e->SetLineWidth(3);
  f->SetLineWidth(3);
  a2->SetLineColor(kCyan-2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+2);
  f->SetLineColor(kMagenta-5);
  b->SetLineStyle(9);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");
  f->SetTitle("SL+Pythia");

  b->GetYaxis()->SetRangeUser(0.9,1.001);
  b->GetXaxis()->SetTitle("cut value / GeV");
  b->GetYaxis()->SetTitle("efficiency #epsilon");


  TGraphErrors* b3 = new TGraphErrors(3);
  TGraphErrors* c3 = new TGraphErrors(3);
  TGraphErrors* d3 = new TGraphErrors(3);
  for (int j=0; j<3; j++)
    {
      b3->SetPoint(j,(*hf_cuts_equivalent)[j],b2->GetBinContent(j+2));
      c3->SetPoint(j,(*hf_cuts_equivalent)[j],c2->GetBinContent(j+2));
      d3->SetPoint(j,(*hf_cuts_equivalent)[j],d2->GetBinContent(j+2));
    }
  b3->SetMarkerColor(b->GetLineColor());
  c3->SetMarkerColor(c->GetLineColor());
  d3->SetMarkerColor(d->GetLineColor());
  b3->SetMarkerSize(1.5);
  c3->SetMarkerSize(1.5);
  d3->SetMarkerSize(1.5);

  TCanvas* can1 = new TCanvas;
  b->Draw("HIST L");
  c->Draw("HIST L SAME");
  d->Draw("HIST L SAME");
  TLegend* leg = can1->BuildLegend();
  b3->Draw("P");
  c3->Draw("P");
  d3->Draw("P");
  leg->SetFillColor(kWhite);
  leg->SetX1(0.65);
  leg->SetX2(0.85);
  leg->SetY1(0.2);
  leg->SetY2(0.4);
  leg->Draw();

  a->GetYaxis()->SetRangeUser(0,1.01);
  a->GetXaxis()->SetTitle("cut value / GeV");
  a->GetYaxis()->SetTitle("event fraction");

  TH1D* copy = new TH1D;
  a->Copy(*copy);  TH1D* a1 = copy;
  a1->SetMarkerColor(kRed);
  a1->SetName("data_copy");
  a1->SetTitle("corrected data");

  TH1D* copy2 = new TH1D;
  a->Copy(*copy2);
  TH1D* h_corr = copy2;
  h_corr->SetName("correctionfactor");
  h_corr->SetTitle("correction factor;cut value / GeV;correction factor");

  int startPlot=0;
  for(int i=1; i<=a1->GetNbinsX(); i++)
    {
      const double f_em = 0.5* (e->GetBinContent(i) + f->GetBinContent(i)) * 195./2100.;
      const double f_noise = a2->GetBinContent(i)/a2->GetBinContent(1);
      const double f_mc = (c->GetBinContent(i)/c->GetBinContent(1) + d->GetBinContent(i)/d->GetBinContent(1))/2.;
      const double n_sel_zb = a->GetBinContent(i);
      const double n_zb = 1;

      double n_inel = 0;
      double corr = -1;
      if(f_noise<1) //if eff can be calculated
        {
          n_inel=(n_sel_zb/f_mc-n_zb*f_noise)/(1-f_noise) * (1.-f_em);
          if(n_inel)
            corr=n_sel_zb/n_inel;
          a1->SetBinContent(i,n_inel);
          h_corr->SetBinContent(i,corr);
        }
      else
        startPlot = i+1;

      if((i-1)%10 == 0)
      cout
        << setprecision(4)
        << endl << i << "(" << a1->GetBinCenter(i) << ")"
        << endl << "f_mc= " << f_mc
        << endl << "f_em= " << f_em
        << endl << "f_noise= " << f_noise
        << endl << "n_sel_zb= " << n_sel_zb
        << endl << "n_inel= " << n_inel
        << endl << "corr factor=" << corr << endl << endl;

      if(i==a->FindBin(cut_value_single) && type[n]==string("single"))
        {
          corr_fac_em[0] = f_em;
          corr_fac_mc[0] = f_mc;
          corr_fac_noise[0] = f_noise;
          corr_fac_all[0] = corr;
        }
      if(i==a->FindBin(cut_value_double) && type[n]==string("double"))
        {
          corr_fac_em[1] = f_em;
          corr_fac_mc[1] = f_mc;
          corr_fac_noise[1] = f_noise;
          corr_fac_all[1] = corr;

          cout << i << "EM error: " << fabs(e->GetBinContent(i) - f->GetBinContent(i)) * 195./2100. / f_em << endl;
          cout << "MC error: " << fabs((c->GetBinContent(i)/c->GetBinContent(1) - d->GetBinContent(i)/d->GetBinContent(1))/2.) / f_mc << endl;
        }
    }
  a1->GetXaxis()->SetRange(startPlot,a1->GetNbinsX());
  h_corr->GetXaxis()->SetRange(startPlot,h_corr->GetNbinsX());
  a1->SetLineColor(kRed);
  CMSPreliminary();
  can1->SaveAs((string("plots/full_p_space_eff_")+type[n]+string(".eps")).c_str());
  can1->SaveAs((string("plots/full_p_space_eff_")+type[n]+string(".pdf")).c_str());
  can1->SaveAs((string("plots/full_p_space_eff_")+type[n]+string(".png")).c_str());

  TCanvas* can2 = new TCanvas;
  a2->Draw("HIST l");
  //a1->Draw("HIST l SAME");
  //a2->Draw("HIST l SAME");
  e->Draw("HIST l SAME");
  f->Draw("HIST l SAME");
  TLegend* leg2 = can2->BuildLegend();
  leg2->SetX1(0.45);
  leg2->SetX2(0.80);
  leg2->SetY1(0.65);
  leg2->SetY2(0.85);
  leg2->SetFillColor(kWhite);
  leg2->Draw();
  CMSPreliminary();
  can2->SaveAs((string("plots/starlight_eff_")+type[n]+string(".eps")).c_str());
  can2->SaveAs((string("plots/starlight_eff_")+type[n]+string(".pdf")).c_str());
  can2->SaveAs((string("plots/starlight_eff_")+type[n]+string(".png")).c_str());

//   TCanvas* can3 = new TCanvas;
//   h_corr->GetYaxis()->SetRangeUser(0.5,3.);
//   h_corr->Draw();
//   TF1* line = new TF1("line","1",0,10);
//   line->Draw("SAME");
//   CMSPreliminary();

  //file->Close();
  //file2->Close();
    }

  TFile outfile("plots/corr_factors.root","recreate");
  corr_fac_em.Write("corr_fac_em");
  corr_fac_noise.Write("corr_fac_noise");
  corr_fac_mc.Write("corr_fac_mc");
  corr_fac_all.Write("corr_fac_all");
  outfile.Close();
}
