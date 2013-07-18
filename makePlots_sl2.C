{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_SL.root");
  TH1D* a=file->Get("data/data_h_hf_hits_coll");
  TH1D* b=file->Get("Epos/Epos_h_hf_hits_coll");
  TH1D* c=file->Get("Epos_SL/Epos_SL_h_hf_hits_coll");

  const int normbin = a->FindBin(10);
  const int normendbin = a->GetNbinsX();
  a->Scale(1./double(a->GetEntries()));
  b->Scale(a->Integral(normbin,normendbin)/b->Integral(normbin,normendbin));
  c->Scale(a->Integral(normbin,normendbin)/c->Integral(normbin,normendbin));

  a->SetMarkerSize(1.25);
  a->SetLineWidth(2.5);
  b->SetLineWidth(2.5);
  c->SetLineWidth(2.5);


  a->SetMarkerColor(kBlack);
  b->SetMarkerColor(kBlue);
  c->SetMarkerColor(kRed);

  a->SetLineColor(kBlack);
  b->SetLineColor(kBlue);
  c->SetLineColor(kRed);

  a->SetTitle("zero bias");
  b->SetTitle("EPOS (GFlash)");
  c->SetTitle("EPOS (SL)");

  a->GetYaxis()->SetRangeUser(1e-7,1.01);
  a->GetXaxis()->SetTitle("hottest HF rec hit E / GeV");
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
  b->Draw("L SAME");
  c->Draw("L SAME");
  TLegend* leg = c1->BuildLegend(0.4,0.67,0.78,0.88);
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->SetLogy();
  c1->SaveAs((string("plots/showerlibrary_2")+string(".eps")).c_str());
  c1->SaveAs((string("plots/showerlibrary_2")+string(".pdf")).c_str());
  c1->SaveAs((string("plots/showerlibrary_2")+string(".png")).c_str());

//   TCanvas* c2 = new TCanvas;
//   TH1D* acopy = a.Clone("acopy");
//   acopy->GetXaxis()->SetRangeUser(0,10);
//   acopy->Draw("P");
//   a2->Draw("L SAME");
//   b->Draw("L SAME");
//   c->Draw("L SAME");
//   d->Draw("L SAME");
//   e->Draw("L SAME");
//   TLegend* leg = c1->BuildLegend(0.4,0.67,0.78,0.88);
//   leg->SetFillColor(kWhite);
//   leg->Draw();
//   c2->SetLogy();
  
}
