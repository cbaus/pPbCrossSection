void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d);

void makePlots6()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");
  TH1D* a=file->Get("data/data_h_perf_hf_totE_single_3gev");
  TH1D* b=file->Get("Hijing/Hijing_h_perf_hf_totE_single_3gev");
  TH1D* c=file->Get("Epos/Epos_h_perf_hf_totE_single_3gev");
  TH1D* d=c;//file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_single_3gev");

  Show(a,b,c,d);

  TFile* file = TFile::Open("histos.root");
  TH1D* e=file->Get("data/data_h_perf_hf_totE_double_1dot5gev");
  TH1D* f=file->Get("Hijing/Hijing_h_perf_hf_totE_double_1dot5gev");
  TH1D* g=file->Get("Epos/Epos_h_perf_hf_totE_double_1dot5gev");
  TH1D* h=g;//file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_double_1dot5gev");

  Show(e,f,g,h);
}

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d)
{
  const int normbin1 = a->FindBin(1000);
  const int normbin2 = a->FindBin(1500);
  a->Scale(1./double(a->GetEntries()));
  b->Scale(a->Integral(normbin1,normbin2)/b->Integral(normbin1,normbin2));
  c->Scale(a->Integral(normbin1,normbin2)/c->Integral(normbin1,normbin2));
  d->Scale(a->Integral(normbin1,normbin2)/d->Integral(normbin1,normbin2));

  a->SetMarkerSize(1.2);
  a->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  d->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);

  a->SetLineColor(kBlack);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);

  a->SetTitle("zero bias");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");

  a->GetXaxis()->SetRangeUser(0,5000);
  a->GetXaxis()->SetTitle("total HF recHit E / GeV");
  a->GetYaxis()->SetTitle("events (normalised)");

  TCanvas* c1 = new TCanvas;
  a->Draw("");
  b->Draw("HIST SAME");
  c->Draw("HIST SAME");
  d->Draw("HIST SAME");
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();

  c1->SetLogy();  
}
