void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type);

void makePlots5()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=file->Get("data/data_h_perf_hf_totE_eta_single_3gev");
  TH1D* b=file->Get("Hijing/Hijing_h_perf_hf_totE_eta_single_3gev");
  TH1D* c=file->Get("Epos/Epos_h_perf_hf_totE_eta_single_3gev");
  TH1D* d=file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_single_3gev");

  Show(a,b,c,d,"single");

  TFile* file = TFile::Open("histos_old.root");
  TH1D* e=file->Get("data/data_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* f=file->Get("Hijing/Hijing_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* g=file->Get("Epos/Epos_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* h=file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_double_1dot5gev");

  Show(e,f,g,h,"double");
}

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type)
{
  a->Scale(1./double(a->Integral()));
  b->Scale(1./double(b->Integral()));
  c->Scale(1./double(c->Integral()));
  d->Scale(1./double(d->Integral()));

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

  //a->GetYaxis()->SetRangeUser(1e-6,1.01);
  a->GetXaxis()->SetTitle("#eta");
  a->GetYaxis()->SetTitle("E / a.u.");

  TCanvas* c1 = new TCanvas;
  a->Draw("");
  b->Draw("HIST SAME");
  c->Draw("HIST SAME");
  d->Draw("HIST SAME");
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SaveAs((string("plots/hf_perf_2_")+type+string(".eps")).c_str());
  c1->SaveAs((string("plots/hf_perf_2_")+type+string(".pdf")).c_str());
  c1->SaveAs((string("plots/hf_perf_2_")+type+string(".png")).c_str());
  
}
