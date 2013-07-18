{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");
  TH1D* a1=file->Get("data/data_h_hf+_totE_lumi");
  TH1D* a2=file->Get("data/data_h_hf-_totE_lumi");

  a1->Rebin(50);
  a2->Rebin(50);

  a1->SetLineWidth(2);
  a2->SetLineWidth(2);

  a1->SetMarkerColor(kRed);
  a2->SetMarkerColor(kBlue);

  a1->SetLineColor(kRed);
  a2->SetLineColor(kBlue);

  a1->SetTitle("HF+");
  a2->SetTitle("HF-");

  //b->GetYaxis()->SetRangeUser(0.9,1.01);
  a2->GetXaxis()->SetTitle("lumi section");
  a2->GetYaxis()->SetTitle("1/N_{tot} #sumE / GeV");

  TCanvas* c1 = new TCanvas;
  a1->Fit("pol1");
  a2->Fit("pol1");
  a2->Draw("PROFS");
  a1->Draw("PROFS SAME");
  TLegend* leg = c1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->Draw();



}
