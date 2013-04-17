{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file2 = TFile::Open("histos_old.root");
  TFile* file = TFile::Open("histos_old.root");
  TH1D* a=file2->Get("data/data_h_hf_hits_coll");
  TH1D* a2=file2->Get("data/data_h_hf_hits_noise");
  TH1D* b=file->Get("Hijing/Hijing_h_hf_hits_coll");
  TH1D* c=file->Get("Epos/Epos_h_hf_hits_coll");
  TH1D* d=file->Get("QGSJetII/QGSJetII_h_hf_hits_coll");
  TH1D* e=file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll");

  const int normbin = a->FindBin(20);
  double eposscale=double(c->GetEntries())/double(e->GetEntries());
  a->Scale(1./double(a->GetEntries()));
  a2->Scale(1./double(a2->GetEntries()));
  b->Scale(a->GetBinContent(normbin)/b->GetBinContent(normbin));
  double eposnorm=a->GetBinContent(normbin)/c->GetBinContent(normbin);
  c->Scale(eposnorm);
  d->Scale(a->GetBinContent(normbin)/d->GetBinContent(normbin));
  e->Scale(eposscale*195./2100.*eposnorm);
  cout << " ada ds asd asd asd !!!! " << c->GetEntries() << " " << e->GetEntries() << " " << eposscale*195./2100. <<  "" << eposscale*eposnorm << endl;

  a->SetMarkerSize(1.2);
  a->SetLineWidth(2.5);
  a2->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  d->SetLineWidth(2.5);
  e->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  a2->SetMarkerColor(kCyan+2);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);
  e->SetMarkerColor(kMagenta+2);

  a->SetLineColor(kBlack);
  a2->SetLineColor(kCyan+2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+2);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");

  a->GetYaxis()->SetRangeUser(1e-7,1.01);
  a->GetXaxis()->SetTitle("hottest cell E / GeV");
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
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SetLogy();
  
}
