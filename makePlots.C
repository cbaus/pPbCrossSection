void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, TH1D* e, TH1D* f, string type);

void makePlots()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_old.root");
  TFile* file2 = TFile::Open("histos.root");
  TH1D* a=file->Get("data/data_h_hf_cut_single");
  TH1D* a2=file->Get("data/data_h_hf_cut_single_noise");
  TH1D* b=file->Get("Hijing/Hijing_h_hf_cut_single");
  TH1D* c=file->Get("Epos/Epos_h_hf_cut_single");
  TH1D* d=file->Get("QGSJetII/QGSJetII_h_hf_cut_single");
  TH1D* e=file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_single");

  Show(a,a2,b,c,d,e,string("single"));

  TH1D* aa=file->Get("data/data_h_hf_cut_double");
  TH1D* aa2=file->Get("data/data_h_hf_cut_double_noise");
  TH1D* bb=file->Get("Hijing/Hijing_h_hf_cut_double");
  TH1D* cc=file->Get("Epos/Epos_h_hf_cut_double");
  TH1D* dd=file->Get("QGSJetII/QGSJetII_h_hf_cut_double");
  TH1D* ee=file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_double");

  Show(aa,aa2,bb,cc,dd,ee,string("double"));
}

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c, TH1D* d, TH1D* e, string type)
{

  a->Scale(1./double(a->GetBinContent(1)));
  a2->Scale(1./double(a2->GetBinContent(1)));
  b->Scale(1./double(b->GetBinContent(1)));
  c->Scale(1./double(c->GetBinContent(1)));
  d->Scale(1./double(d->GetBinContent(1)));
  e->Scale(1./double(e->GetBinContent(1)));

  a->SetLineWidth(2);
  a2->SetLineWidth(2);
  b->SetLineWidth(2);
  c->SetLineWidth(2);
  d->SetLineWidth(2);
  e->SetLineWidth(2);
  a2->SetLineColor(kCyan-2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+2);
  b->SetLineStyle(9);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");

  b->GetYaxis()->SetRangeUser(0.9,1.001);
  b->GetXaxis()->SetTitle("cut value / GeV");
  b->GetYaxis()->SetTitle("efficiency #epsilon");

  TCanvas* c1 = new TCanvas;
  b->Draw("HIST L");
  c->Draw("HIST L SAME");
  d->Draw("HIST L SAME");
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->SetX1(0.3);
  leg->SetX2(0.5);
  leg->SetY1(0.3);
  leg->SetY2(0.5);
  leg->Draw();
  
  a->GetYaxis()->SetRangeUser(0,1.01);
  a->GetXaxis()->SetTitle("cut value / GeV");
  a->GetYaxis()->SetTitle("event fraction");

  TH1D* copy = new TH1D;
  a->Copy(*copy);
  TH1D* a1 = copy;
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
      const double f_em = e->GetBinContent(i)/e->GetBinContent(1) / 2100. * 195.;
      const double f_noise = a2->GetBinContent(i)/a2->GetBinContent(1);
      const double f_mc = (b->GetBinContent(i)/b->GetBinContent(1) + c->GetBinContent(i)/c->GetBinContent(1) + d->GetBinContent(i)/d->GetBinContent(1))/3.;
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
    

      cout 
        << endl << i << "(" << a1->GetBinCenter(i) << ")"
        << endl << "f_mc= " << f_mc
        << endl << "f_em= " << f_em
        << endl << "f_noise= " << f_noise
        << endl << "n_sel_zb= " << n_sel_zb
        << endl << "n_inel= " << n_inel
        << endl << "corr factor=" << corr << endl << endl;
    }
  a1->GetXaxis()->SetRange(startPlot,a1->GetNbinsX());
  h_corr->GetXaxis()->SetRange(startPlot,h_corr->GetNbinsX());
  a1->SetLineColor(kRed);
  c1->SaveAs((string("plots/full_p_space_eff_")+type+string(".eps")).c_str());
  c1->SaveAs((string("plots/full_p_space_eff_")+type+string(".pdf")).c_str());
  c1->SaveAs((string("plots/full_p_space_eff_")+type+string(".png")).c_str());

  TCanvas* c2 = new TCanvas;
  a->Draw("HIST l");
  a1->Draw("HIST l SAME");
  a2->Draw("HIST l SAME");
  e->Draw("HIST l SAME");
  TLegend* leg2 = c2->BuildLegend();
  leg2->SetX1(0.45);
  leg2->SetX2(0.80);
  leg2->SetY1(0.65);
  leg2->SetY2(0.85);
  leg2->SetFillColor(kWhite);
  leg2->Draw();
  c2->SaveAs((string("plots/data_eff_")+type+string(".eps")).c_str());
  c2->SaveAs((string("plots/data_eff_")+type+string(".pdf")).c_str());
  c2->SaveAs((string("plots/data_eff_")+type+string(".png")).c_str());

  TCanvas* c3 = new TCanvas;
  h_corr->GetYaxis()->SetRangeUser(0.5,3.);
  h_corr->Draw();
  TF1* line = new TF1("line","1",0,5);
  line->Draw("SAME");
  c3->SaveAs((string("plots/corr_factor_")+type+string(".eps")).c_str());
  c3->SaveAs((string("plots/corr_factor_")+type+string(".pdf")).c_str());
  c3->SaveAs((string("plots/corr_factor_")+type+string(".png")).c_str());

}
