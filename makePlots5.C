void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type);

void makePlots5()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_old.root");
  TFile* file2 = TFile::Open("histos_mc.root");
  TH1D* a= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_single_3gev");
  TH1D* b= (TH1D*)file2->Get("Hijing/Hijing_h_perf_hf_totE_eta_single_3gev");
  TH1D* c= (TH1D*)file2->Get("Epos/Epos_h_perf_hf_totE_eta_single_3gev");
  TH1D* d= (TH1D*)file2->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_single_3gev");

  Show(a,b,c,d,"single");

  TFile* file = TFile::Open("histos_old.root");
  TFile* file2 = TFile::Open("histos_mc.root");
  TH1D* e= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* f= (TH1D*)file2->Get("Hijing/Hijing_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* g= (TH1D*)file2->Get("Epos/Epos_h_perf_hf_totE_eta_double_1dot5gev");
  TH1D* h= (TH1D*)file2->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_double_1dot5gev");

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

  double maximumy = 1.2 * TMath::Max(TMath::Max(a->GetMaximum(),b->GetMaximum()),TMath::Max(c->GetMaximum(),d->GetMaximum()));
  a->GetYaxis()->SetRangeUser(0,maximumy);
  a->GetXaxis()->SetTitle("#eta");
  a->GetYaxis()->SetTitle("E / a.u.");

  TCanvas* c1 = new TCanvas;
  a->Draw("");
  b->Draw("HIST SAME");
  c->Draw("HIST SAME");
  d->Draw("HIST SAME");
  TLegend* leg = new TLegend(0.3,0.7,0.7,0.9);
  leg->AddEntry(a,"","P");
  leg->AddEntry(b,"","L");
  leg->AddEntry(c,"","L");
  leg->AddEntry(d,"","L");
  SetLegAtt(leg);
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SaveAs((string("plots/hf_perf_2_")+type+string(".pdf")).c_str());

}
