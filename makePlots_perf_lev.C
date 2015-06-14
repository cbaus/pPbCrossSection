///HF plots E over eta

void Show(TH1D* a,TH1D* b,TH1D* c, TH1D* epossl, TH1D* d, TH1D* e, TH1D* f, string type);

void makePlots_perf_lev()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();


  TFile* file = TFile::Open("histos.root");
  TFile* file2 = TFile::Open("histos.root");
  TH1D* data= (TH1D*)file->Get("data247324/data247324_h_perf_hf_totE_eta_lev_m");
  TH1D* mc1= (TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_hf_totE_eta_lev_m");
  TH1D* mc2= (TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_hf_totE_eta_lev_m");
  TH1D* mc3= (TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_hf_totE_eta_lev_m");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_m");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_m");
  // TH1D* i=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_perf_hf_totE_eta_lev_m");
  // TH1D* j=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_perf_hf_totE_eta_lev_m");

  Show(data,mc1,mc2,mc3,h,0,0,"lev_minus");

  TFile* file = TFile::Open("histos.root");
  TFile* file2 = TFile::Open("histos.root");
  TH1D* data= (TH1D*)file->Get("data247324/data247324_h_perf_hf_totE_eta_lev_p");
  TH1D* mc1= (TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_hf_totE_eta_lev_p");
  TH1D* mc2= (TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_hf_totE_eta_lev_p");
  TH1D* mc2= (TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_hf_totE_eta_lev_p");
  TH1D* g= (TH1D*)file->Get("Epos/Epos_h_perf_hf_totE_eta_lev_p");
  TH1D* epossl= (TH1D*)file->Get("Epos_SL/Epos_SL_h_perf_hf_totE_eta_lev_p");
  TH1D* h= (TH1D*)file->Get("QGSJetII/QGSJetII_h_perf_hf_totE_eta_lev_p");
  // TH1D* i=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_perf_hf_totE_eta_lev_p");
  // TH1D* j=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_perf_hf_totE_eta_lev_p");

  Show(data,mc1,mc2,mc3,h,0,0,"lev_plus");
}

void Show(TH1D* data, TH1D* mc1, TH1D* mc2, TH1D* mc3, TH1D* mc4, TH1D* sl1, TH1D* sl2, string type)
{
//   data->Scale(1./double(data->Integral()));
//   mc1->Scale(1./double(mc1->Integral()));
//   mc2->Scale(1./double(mc2->Integral()));
//   mc3->Scale(1./double(mc2->Integral()));
//   mc4->Scale(1./double(mc4->Integral()));
  // sl1->Scale(195./2130.);
  // sl2->Scale(122./2130.);
  // sl1->SetBit(TH1::kIsAverage);
  // sl2->SetBit(TH1::kIsAverage);
  // if(sl1->Add(sl2) == kFALSE) {new TCanvas; sl1->Draw(); new TCanvas; sl2->Draw(); cerr << "add failed" << endl; return;}
  // TH1D* sl = sl1;

  data->SetMarkerSize(1.2);
  data->SetLineWidth(2.5);
  mc1->SetLineWidth(2.5);
  if(mc2) mc2->SetLineWidth(2.5);
  if(mc3) mc3->SetLineWidth(2.5);
  if(mc4) mc4->SetLineWidth(2.5);
  // sl->SetLineWidth(2.5);


  data->SetMarkerColor(kBlack);
  mc1->SetMarkerColor(kRed);
  if(mc2) mc2->SetMarkerColor(kBlue);
  if(mc3) mc3->SetMarkerColor(kCyan);
  if(mc4) mc4->SetMarkerColor(kGreen+2);
  // sl->SetMarkerColor(kOrange-9);
  // sl->SetFillColor(sl->GetMarkerColor());

  data->SetLineColor(kBlack);
  mc1->SetLineColor(kRed);
  if(mc2) mc2->SetLineColor(kBlue);
  if(mc3) mc3->SetLineColor(kCyan);
  if(mc4) mc4->SetLineColor(kGreen+2);
  // sl->SetLineColor(kOrange-9);

  data->SetTitle("unbiased trigger + CASTOR E_{tot}>5.6 GeV");
  mc1->SetTitle("Pythia6 Z2*");
  if(mc2) mc2->SetTitle("Pythia8 Monash Tune");
  if(mc3) mc3->SetTitle("EPOS (SL)");
  if(mc4) mc4->SetTitle("QGSJetII");
  // sl->SetTitle("#gamma-p (STARLIGHT+DPMJET/PYTHIA)");

  double maximumy = 1.2 * TMath::Max(data->GetMaximum(),mc1->GetMaximum());
  if (mc2&&mc4) maximumy = 1.2 * TMath::Max(TMath::Max(data->GetMaximum(),mc1->GetMaximum()),TMath::Max(mc2->GetMaximum(),mc4->GetMaximum()));
  data->GetYaxis()->SetRangeUser(0,100);
  data->GetXaxis()->SetTitle("#eta");
  data->GetYaxis()->SetTitle("1/N dE_{HF,tot}/d#eta [GeV]");

  TCanvas* c1 = new TCanvas;
  data->Draw("HIST P");
  mc1->Draw("HIST L SAME");
  if(mc2) mc2->Draw("HIST L SAME");
  // sl->Draw("HIST F SAME");
  if(mc3) mc3->Draw("HIST L SAME");
  if(mc4) mc4->Draw("HIST L SAME");
  data->Draw("SAME AXIS");
  TLegend* leg = new TLegend(0.3,0.7,0.7,0.9);
  leg->AddEntry(data,"","P");
  leg->AddEntry(mc1,"","L");
  if(mc2) leg->AddEntry(mc2,"","L");
  if(mc3) leg->AddEntry(mc3,"","L");
  if(mc4) leg->AddEntry(mc4,"","L");
  // leg->AddEntry(sl,"#gammap (not stacked)","F");
  SetLegAtt(leg);
  leg->SetFillColor(kWhite);
  leg->Draw();
  //c1->SaveAs((string("plots/hf_perf_2_")+type+string(".pdf")).c_str());

}

TH1D* merge(TH1D* bg1, TH1D* signal)
{
  TList *list = new TList;
  list->Add(bg1);
  list->Add(signal);
  ostringstream newname; newname << signal->GetName() << "_merged";
  TH1D* newsignal = (TH1D*)signal->Clone(newname.str().c_str());
  newsignal->Reset();
  newsignal->Merge(list);
  return newsignal;
}
