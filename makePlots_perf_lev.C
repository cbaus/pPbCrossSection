///HF plots E over eta

void Show(TH1D* a,TH1D* b,TH1D* c, TH1D* epossl, TH1D* d, string type);

void makePlots_perf_lev()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();


  TFile* file = TFile::Open("histos_nohfcalib.root");
  TFile* file2 = TFile::Open("histos_nohfcalib.root");
  TH1D* e= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_lev_m");
  TH1D* f= (TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_eta_lev_m");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_m");
  TH1D* epossl= (TH1D*)file->Get("Epos_SL/Epos_SL_h_perf_hf_totE_eta_lev_m");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_m");

  Show(e,f,g,epossl,h,"lev_minus");

  TFile* file = TFile::Open("histos_nohfcalib.root");
  TFile* file2 = TFile::Open("histos_nohfcalib.root");
  TH1D* e= (TH1D*)file->Get("data210885/data210885_h_perf_hf_totE_eta_lev_p");
  TH1D* f= (TH1D*)file->Get("Hijing/Hijing_h_perf_hf_totE_eta_lev_p");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_p");
  TH1D* epossl= (TH1D*)file->Get("Epos_SL/Epos_SL_h_perf_hf_totE_eta_lev_p");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_p");

  Show(e,f,g,epossl,h,"lev_plus");
}

void Show(TH1D* a, TH1D* b, TH1D* c, TH1D* epossl, TH1D* d, string type)
{
//   a->Scale(1./double(a->Integral()));
//   b->Scale(1./double(b->Integral()));
//   c->Scale(1./double(c->Integral()));
//   epossl->Scale(1./double(c->Integral()));
//   d->Scale(1./double(d->Integral()));

  a->SetMarkerSize(1.2);
  a->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  epossl->SetLineWidth(2.5);
  d->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  epossl->SetMarkerColor(kCyan);
  d->SetMarkerColor(kGreen+2);

  a->SetLineColor(kBlack);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  epossl->SetLineColor(kCyan);
  d->SetLineColor(kGreen+2);

  a->SetTitle("zero bias");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  epossl->SetTitle("EPOS (SL)");
  d->SetTitle("QGSJetII");

  double maximumy = 1.2 * TMath::Max(TMath::Max(a->GetMaximum(),b->GetMaximum()),TMath::Max(c->GetMaximum(),d->GetMaximum()));
  a->GetYaxis()->SetRangeUser(0,60);
  a->GetXaxis()->SetTitle("#eta");
  a->GetYaxis()->SetTitle("E [GeV]");

  TCanvas* c1 = new TCanvas;
  a->Draw("HIST P");
  b->Draw("HIST L SAME");
  c->Draw("HIST L SAME");
  epossl->Draw("HIST L SAME");
  d->Draw("HIST L SAME");
  TLegend* leg = new TLegend(0.3,0.7,0.7,0.9);
  leg->AddEntry(a,"","P");
  leg->AddEntry(b,"","L");
  leg->AddEntry(c,"","L");
  leg->AddEntry(epossl,"","L");
  leg->AddEntry(d,"","L");
  SetLegAtt(leg);
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SaveAs((string("plots/hf_perf_2_")+type+string(".pdf")).c_str());

}
