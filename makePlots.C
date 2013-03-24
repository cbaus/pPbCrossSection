void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, TH1D* e);

void makePlots()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");
  TH1D* a=file->Get("data/data_h_hf_cut_single");
  TH1D* a2=file->Get("data/data_h_hf_cut_single_noise");
  TH1D* b=file->Get("Hijing/Hijing_h_hf_cut_single");
  TH1D* c=file->Get("Epos/Epos_h_hf_cut_single");
  TH1D* d=file->Get("Epos/Epos_h_hf_cut_single");//QGSJetII/QGSJetII_h_hf_cut_single");

  Show(a,a2,b,c,d);

  TH1D* aa=file->Get("data/data_h_hf_cut_double");
  TH1D* aa2=file->Get("data/data_h_hf_cut_double_noise");
  TH1D* bb=file->Get("Hijing/Hijing_h_hf_cut_double");
  TH1D* cc=file->Get("Epos/Epos_h_hf_cut_double");
  TH1D* dd=file->Get("Epos/Epos_h_hf_cut_double");//QGSJetII/QGSJetII_h_hf_cut_double");

  Show(aa,aa2,bb,cc,dd);
}

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c, TH1D* d)
{

  a->Scale(1./double(a->GetBinContent(1)));
  a2->Scale(1./double(a2->GetBinContent(1)));
  b->Scale(1./double(b->GetBinContent(1)));
  c->Scale(1./double(c->GetBinContent(1)));
  d->Scale(1./double(d->GetBinContent(1)));

  a->SetLineWidth(2);
  a2->SetLineWidth(2);
  b->SetLineWidth(2);
  c->SetLineWidth(2);
  d->SetLineWidth(2);
  a2->SetLineColor(kCyan-2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");

  b->GetYaxis()->SetRangeUser(0.9,1.01);
  b->GetXaxis()->SetTitle("cut value / GeV");
  b->GetYaxis()->SetTitle("efficiency #epsilon");

  TCanvas* c1 = new TCanvas;
  b->Draw("HIST L");
  c->Draw("HIST L SAME");
  d->Draw("HIST L SAME");
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();
  
  a->GetYaxis()->SetRangeUser(0,1.01);
  a->GetXaxis()->SetTitle("cut value / GeV");
  a->GetYaxis()->SetTitle("efficiency #epsilon");

  TH1D* copy = new TH1D;
  a->Copy(*copy);
  TH1D* a1 = copy;
  a1->SetMarkerColor(kRed);
  a1->SetName("data_copy");
  a1->SetTitle("data (corrected by Epos eff.)");
  for(int i=2; i<=a1->GetNbinsX(); i++)
    {
      const double f_noise = a2->GetBinContent(i)/a2->GetBinContent(1);
      const double f_zb = a->GetBinContent(i)/a->GetBinContent(1);
      const double f_mc = b->GetBinContent(i)/b->GetBinContent(1);
      double tmp;
      if(f_zb-f_mc==0)
        tmp = 1e9;
      else
        tmp=(f_noise-f_zb)/(f_zb-f_mc);
      const double f_prime = tmp;
      const double y = a2->GetBinContent(1);
      const double x_plus_y = a->GetBinContent(1);
      a1->SetBinContent(i,f_prime/(f_prime+1.)*x_plus_y);

      cout 
        << endl << i
        << endl << "f_mc= " << f_mc
        << endl << "f_zb= " << f_zb
        << endl << "f_noise= " << f_noise
        << endl << "f_prime*y= " << f_prime*y
        << endl << "f_prime/(f_prime+1)*x_plus_y= " << f_prime/(f_prime+1.)*x_plus_y << endl << endl;
    }
  a1->GetXaxis()->SetRange(4,a1->GetNbinsX());
  a1->SetLineColor(kRed);

  TCanvas* c2 = new TCanvas;
  a->Draw("HIST l");
  a1->Draw("HIST l SAME");
  a2->Draw("HIST l SAME");
  TLegend* leg2 = c2->BuildLegend();
  leg2->SetFillColor(kWhite);
  leg2->Draw();

}
