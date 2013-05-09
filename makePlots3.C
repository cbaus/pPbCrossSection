{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file2 = TFile::Open("histos.root");
  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=file->Get("data/data_h_hf-_hits_min");
  TH1D* a2=file->Get("data/data_h_hf_hits_noise");
  TH1D* b=file->Get("Hijing/Hijing_h_hf-_hits_min");
  TH1D* c=file->Get("Epos/Epos_h_hf-_hits_min");
  TH1D* d=file->Get("QGSJetII/QGSJetII_h_hf-_hits_min");
  TH1D* e=file2->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfm_hits_coll");
  TH1D* f=file2->Get("Starlight_Pythia/Starlight_Pythia_h_hfm_hits_coll");

  const int normbin = a->FindBin(10);
  const int normendbin = a->GetNbinsX();
  double eposscale=double(c->GetEntries())/double(e->GetEntries());
  double eposscale2=double(c->GetEntries())/double(f->GetEntries());
  a->Scale(1./double(a->GetEntries()));
  a2->Scale(1./double(a2->GetEntries()));
  b->Scale(a->Integral(normbin,normendbin)/b->Integral(normbin,normendbin));
  double eposnorm=a->Integral(normbin,normendbin)/c->Integral(normbin,normendbin);
  c->Scale(eposnorm);
  d->Scale(a->Integral(normbin,normendbin)/d->Integral(normbin,normendbin));
  e->Scale(eposscale*195./2100.*eposnorm);
  f->Scale(eposscale2*122./2100.*eposnorm);

  a->SetMarkerSize(1.25);
  a->SetLineWidth(2.5);
  a2->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  d->SetLineWidth(2.5);
  e->SetLineWidth(2.5);
  f->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  a2->SetMarkerColor(kCyan+2);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);
  e->SetMarkerColor(kMagenta+4);
  f->SetMarkerColor(kMagenta-5);

  b->SetMarkerStyle(4);
  c->SetMarkerStyle(4);
  d->SetMarkerStyle(4);
  e->SetMarkerStyle(34);
  f->SetMarkerStyle(34);

  a->SetLineColor(kBlack);
  a2->SetLineColor(kCyan+2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+4);
  f->SetLineColor(kMagenta-5);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");
  f->SetTitle("SL+Pythia");

  a->GetYaxis()->SetRangeUser(1e-7,1.01);
  a->GetXaxis()->SetTitle("hottest HF tower E / GeV");
  a->GetYaxis()->SetTitle("events (normalised)");

//   double fac = 1.25;
//   TH1D* bs = new TH1D ("asd","asd",a->GetNbinsX(),0,200*fac);
//   for (int i=0; i<a->GetNbinsX()+1;i++)
//     {
//       bs->SetBinContent(i,a->GetBinContent(i));
//     }
//   a=bs;

  TCanvas* c1 = new TCanvas;
  a->Draw("P");
  a2->Draw("L SAME");
  b->Draw("L SAME");
  c->Draw("L SAME");
  d->Draw("L SAME");
  e->Draw("L SAME");
  f->Draw("L SAME");
  TLegend* leg = c1->BuildLegend(0.4,0.67,0.78,0.88);
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SetLogy();

  TCanvas* c2 = new TCanvas;
  TH1D* acopy = a.Clone("acopy");
  acopy->GetXaxis()->SetRangeUser(0,10);
  acopy->Draw("P");
  a2->Draw("L SAME");
  b->Draw("L SAME");
  c->Draw("L SAME");
  d->Draw("L SAME");
  e->Draw("L SAME");
  f->Draw("L SAME");
  TLegend* leg = c1->BuildLegend(0.4,0.67,0.78,0.88);
  leg->SetFillColor(kWhite);
  leg->Draw();
  c2->SetLogy();
  
}
