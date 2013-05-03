void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, TH1D* e, TH1D* f);

void makePlots_rechits1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();


  TFile* file = TFile::Open("histos_old.root");
  TFile* file2 = TFile::Open("histos_rechits.root");
  TH1D* a=file->Get("data/data_h_hf_cut_single_noise");
  TH1D* a2=file2->Get("data/data_h_hf_cut_single_noise");
  TH1D* b=file->Get("Epos/Epos_h_hf_cut_single");
  TH1D* b2=file2->Get("Epos/Epos_h_hf_cut_single");
  TH1D* c=file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_single");
  TH1D* c2=file2->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_cut_single");

  Show(a,a2,b,b2,c,c2);
}

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* b2, TH1D* c, TH1D* c2)
{

  a->Scale(1./double(a->GetBinContent(1)));
  a2->Scale(1./double(a2->GetBinContent(1)));
  b->Scale(1./double(b->GetBinContent(1)));
  b2->Scale(1./double(b2->GetBinContent(1)));
  c->Scale(1./double(c->GetBinContent(1)));
  c2->Scale(1./double(c2->GetBinContent(1)));

  a->SetLineWidth(2);
  a2->SetLineWidth(2);
  b->SetLineWidth(2);
  b2->SetLineWidth(2);
  c->SetLineWidth(2);
  c2->SetLineWidth(2);
  a->SetLineColor(kCyan-2);
  a2->SetLineColor(kCyan-2);
  b->SetLineColor(kRed);
  b2->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  c2->SetLineColor(kBlue);
  a2->SetLineStyle(7);
  b2->SetLineStyle(7);
  c2->SetLineStyle(7);

  a->SetTitle("noise (Towers)");
  a2->SetTitle("noise (RecHits)");
  b->SetTitle("EPOS (Towers)");
  b2->SetTitle("EPOS (RecHits)");
  c->SetTitle("SL+DPMJet (Towers)");
  c2->SetTitle("SL+DPMJet (RecHits)");

  b->GetYaxis()->SetRangeUser(0.9,1.001);
  b->GetXaxis()->SetTitle("cut value / GeV");
  b->GetYaxis()->SetTitle("efficiency #epsilon");

  TCanvas* can1 = new TCanvas;
  b->Draw("HIST L");
  b2->Draw("HIST L SAME");
  TLegend* leg = can1->BuildLegend();
  leg->SetFillColor(kWhite);
  leg->SetX1(0.3);
  leg->SetX2(0.7);
  leg->SetY1(0.3);
  leg->SetY2(0.5);
  leg->Draw();
  
  a->GetYaxis()->SetRangeUser(0,1.01);
  a->GetXaxis()->SetTitle("cut value / GeV");
  a->GetYaxis()->SetTitle("event fraction");
  can1->SaveAs((string("plots/rechits_1")+string(".eps")).c_str());
  can1->SaveAs((string("plots/rechits_1")+string(".pdf")).c_str());
  can1->SaveAs((string("plots/rechits_1")+string(".png")).c_str());


  TCanvas* can2 = new TCanvas;
  a->Draw("HIST l");
  a2->Draw("HIST l SAME");
  c->Draw("HIST l SAME");
  c2->Draw("HIST l SAME");
  TLegend* leg2 = can2->BuildLegend();
  leg2->SetX1(0.45);
  leg2->SetX2(0.80);
  leg2->SetY1(0.65);
  leg2->SetY2(0.85);
  leg2->SetFillColor(kWhite);
  leg2->Draw();
  can2->SaveAs((string("plots/rechits_2")+string(".pdf")).c_str());
  can2->SaveAs((string("plots/rechits_2")+string(".png")).c_str());
  can2->SaveAs((string("plots/rechits_2")+string(".eps")).c_str());

}
