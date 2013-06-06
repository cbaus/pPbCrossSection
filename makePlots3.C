void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string prefix);

void makePlots3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  {
  TFile* file2 = TFile::Open("histos_old.root");
  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfp_hits_coll");
  TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_hits_noise");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfp_hits_coll");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfp_hits_coll");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfp_hits_coll");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfp_hits_coll");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfp_hits_coll");
  Show(a,a2,b,c,d,e,f,"p");
  }
  {
  TFile* file2 = TFile::Open("histos_old.root");
  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfm_hits_coll");
  TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_hits_noise");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfm_hits_coll");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfm_hits_coll");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfm_hits_coll");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfm_hits_coll");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfm_hits_coll");
  Show(a,a2,b,c,d,e,f,"m");
  }
}
void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f, string prefix)
  {
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
  e->SetMarkerSize(1.25);
  f->SetMarkerSize(1.25);

  // a->SetLineWidth(2);
  // a2->SetLineWidth(1.5);
  // b->SetLineWidth(2);
  // c->SetLineWidth(2);
  // d->SetLineWidth(2);
  // e->SetLineWidth(2);
  // f->SetLineWidth(2);


  a->SetMarkerColor(kBlack);
  a2->SetMarkerColor(kCyan+2);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);
  e->SetMarkerColor(kMagenta+3);
  f->SetMarkerColor(kMagenta-5);

  a2->SetMarkerStyle(34);
  b->SetMarkerStyle(4);
  c->SetMarkerStyle(4);
  d->SetMarkerStyle(4);
  e->SetMarkerStyle(22);
  f->SetMarkerStyle(23);

  a->SetLineColor(kBlack);
  a2->SetLineColor(kCyan+2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+3);
  f->SetLineColor(kMagenta-5);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");
  f->SetTitle("SL+Pythia");

  a->GetXaxis()->SetLimits(0.5,200);
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
  TLegend* leg = c1->BuildLegend(0.6,0.67,0.78,0.88);
  SetLegAtt(leg);
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  CMSPreliminary();
  c1->SaveAs((string("plots/hf_") + prefix + string("_signal")+string(".pdf")).c_str());

  // TCanvas* c2 = new TCanvas;
  // TH1D* acopy = a->Clone("acopy");
  // acopy->GetXaxis()->SetRangeUser(0,10);
  // acopy->Draw("P");
  // a2->Draw("L SAME");
  // b->Draw("L SAME");
  // c->Draw("L SAME");
  // d->Draw("L SAME");
  // e->Draw("L SAME");
  // f->Draw("L SAME");
  // TLegend* leg = c1->BuildLegend(0.4,0.67,0.78,0.88);
  // leg->SetFillColor(kWhite);
  // leg->Draw();
  // c2->SetLogy();
  // CMSPreliminary();
  // c2->SaveAs((string("plots/hf_p_signal_zoom_")+string(".eps")).c_str());
  // c2->SaveAs((string("plots/hf_p_signal_zoom_")+string(".pdf")).c_str());
  // c2->SaveAs((string("plots/hf_p_signal_zoom_")+string(".png")).c_str());

}
