void Show(TH1D* a,TH1D* b,TH1D* c=0,TH1D* d=0, TH1D* e=0, TH1D* f=0);

void makePlots_sl1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();


  TFile* file = TFile::Open("histos_SL.root");
  TH1D* a=file->Get("Epos_SL/Epos_SL_h_hf_cut_single");
  TH1D* b=file->Get("Epos_SL/Epos_SL_h_hf_cut_double");
  TH1D* a2=file->Get("Epos/Epos_h_hf_cut_single");
  TH1D* b2=file->Get("Epos/Epos_h_hf_cut_double");

  Show(a,a2,b,b2);
}

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* b2, TH1D* c, TH1D* c2)
{

  a->Scale(1./double(a->GetBinContent(1)));
  a2->Scale(1./double(a2->GetBinContent(1)));
  b->Scale(1./double(b->GetBinContent(1)));
  b2->Scale(1./double(b2->GetBinContent(1)));

  a->SetLineWidth(2);
  a2->SetLineWidth(2);
  b->SetLineWidth(2);
  b2->SetLineWidth(2);
  a->SetLineColor(kRed);
  a2->SetLineColor(kRed);
  b->SetLineColor(kBlue);
  b2->SetLineColor(kBlue);
  a2->SetLineStyle(7);
  b2->SetLineStyle(7);

  a->SetTitle("single-arm EPOS (GFlash)");
  a2->SetTitle("single-arm EPOS (SL)");
  b->SetTitle("double-arm EPOS (GFlash)");
  b2->SetTitle("double-arm EPOS (SL)");

  a->GetYaxis()->SetRangeUser(0.9,1.001);
  a->GetXaxis()->SetTitle("cut value / GeV");
  a->GetYaxis()->SetTitle("efficiency #epsilon");


  TCanvas* can1 = new TCanvas;
  a->Draw("HIST l");
  a2->Draw("HIST l SAME");
  b->Draw("HIST l SAME");
  b2->Draw("HIST l SAME");
  TLegend* leg2 = can1->BuildLegend();
  leg2->SetX1(0.23);
  leg2->SetX2(0.58);
  leg2->SetY1(0.23);
  leg2->SetY2(0.44);
  leg2->SetFillColor(kWhite);
  leg2->Draw();
  can1->SaveAs((string("plots/showerlibrary_1")+string(".eps")).c_str());
  can1->SaveAs((string("plots/showerlibrary_1")+string(".pdf")).c_str());
  can1->SaveAs((string("plots/showerlibrary_1")+string(".png")).c_str());

}
