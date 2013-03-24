{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");
  TH1D* a=file->Get("data/data_h_hf_hits_coll");
  TH1D* a2=file->Get("data/data_h_hf_hits_noise");
  TH1D* b=file->Get("Hijing/Hijing_h_hf_hits_coll");
  TH1D* c=file->Get("Epos/Epos_h_hf_hits_coll");
  TH1D* d=c;//file->Get("QGSJetII/QGSJetII_h_hf_hits_coll");

  const int normbin = a->FindBin(20);
  a->Scale(1./double(a->GetEntries()));
  a2->Scale(1./double(a2->GetEntries()));
  b->Scale(a->GetBinContent(normbin)/b->GetBinContent(normbin));
  c->Scale(a->GetBinContent(normbin)/c->GetBinContent(normbin));
  d->Scale(a->GetBinContent(normbin)/d->GetBinContent(normbin));

  a->SetMarkerSize(1.2);
  a->SetLineWidth(2.5);
  a2->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);
  d->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  a2->SetMarkerColor(kCyan+2);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);

  a->SetLineColor(kBlack);
  a2->SetLineColor(kCyan+2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");

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
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SetLogy();
  
}
