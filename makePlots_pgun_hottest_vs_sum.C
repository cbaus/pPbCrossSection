{
  gStyle->SetOptStat(0);

  TFile* file = TFile::Open("histos_deleteme.root");
  TH2D* hottest=file->Get("pgun/pgun_h_hfreco_vs_hfgen");
  TH2D* sum=file->Get("pgun/pgun_h_hfreco_vs_hfgen_sum");

  TProfile* hottestprof = hottest->ProfileX("_pfx", 1, -1, "s");
  TProfile* sumprof = sum->ProfileX("_pfx_sum", 1, -1, "s");

  TLine* line = new TLine(0,0,100,100);
  line->SetLineWidth(1.5);

  hottest->SetTitle("#pi^{+} particle gun 3.5<#eta<4.5; p_{HF,gen} [GeV/c]; E_{HF} [GeV]");
  sum->SetTitle("#pi^{+} particle gun 3.5<#eta<4.5; p_{HF,gen} [GeV/c]; #sumE_{HF,reco} [GeV]");

  hottestprof->SetMarkerSize(1.5);
  sumprof->SetMarkerSize(1.5);
  hottestprof->SetLineWidth(1.5);
  sumprof->SetLineWidth(1.5);
  hottestprof->SetMarkerColor(kRed);
  sumprof->SetMarkerColor(kRed);
  hottestprof->SetLineColor(kRed);
  sumprof->SetLineColor(kRed);

  TCanvas* c1 = new TCanvas;
  hottest->Draw("COLZ");
  c1->SetLogz();
  hottestprof->Draw("SAME");
  line->Draw("SAME");
  c1->SaveAs("plots/pgun_ehf_vs_genehf.pdf");
  c1->SaveAs("plots/pgun_ehf_vs_genehf.png");

  TCanvas* c2 = new TCanvas;
  sum->Draw("COLZ");
  c2->SetLogz();
  sumprof->Draw("SAME");
  line->Draw("SAME");
  c2->SaveAs("plots/pgun_sumehf_vs_genehf.pdf");
  c2->SaveAs("plots/pgun_sumehf_vs_genehf.png");

  TCanvas* c3 = new TCanvas;
  hottestprof->Divide(sumprof);
  hottestprof->Draw();
}
