///HF distribution of number of rechits

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, string type);

void makePlots_perf_hf_norechits()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");
  TH1D* data=(TH1D*)file->Get("data247324/data247324_h_perf_no_of_towers_single");
  TH1D* mc1=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_no_of_towers_single");
  TH1D* mc2=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_no_of_towers_single");
  TH1D* mc3=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_no_of_towers_single");

  Show(data,mc1,mc2,mc3,"single");

  TFile* file = TFile::Open("histos.root");
  TH1D* data=(TH1D*)file->Get("data247324/data247324_h_perf_no_of_towers_double");
  TH1D* mc1=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_no_of_towers_double");
  TH1D* mc2=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_no_of_towers_double");
  TH1D* mc3=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_no_of_towers_double");

  Show(data,mc1,mc2,mc3,"double");

  TFile* file = TFile::Open("histos.root");
  TH1D* data=(TH1D*)file->Get("data247324/data247324_h_perf_no_of_towers_aboveth_single");
  TH1D* mc1=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_perf_no_of_towers_aboveth_single");
  TH1D* mc2=(TH1D*)file->Get("PythiaMonash/PythiaMonash_h_perf_no_of_towers_aboveth_single");
  TH1D* mc3=(TH1D*)file->Get("PythiaMBR/PythiaMBR_h_perf_no_of_towers_aboveth_single");

  Show(data,mc1,mc2,mc3,"singleth");
}

void Show(TH1D* data,TH1D* mc1,TH1D* mc2,TH1D* mc3, string type)
{
  // const int normbin = data->FindBin(300);
  // const int normendbin = data->GetNbinsX();
  // data->Scale(1./double(data->Integral()));
  // if(mc1) mc1->Scale(data->Integral(normbin,normendbin)/mc1->Integral(normbin,normendbin));
  // if(mc2) mc2->Scale(data->Integral(normbin,normendbin)/mc2->Integral(normbin,normendbin));
  // if(mc3) mc3->Scale(data->Integral(normbin,normendbin)/mc3->Integral(normbin,normendbin));

  data->Scale(1./data->Integral());
  if (mc1) mc1->Scale(1./mc1->Integral());
  if (mc2) mc2->Scale(1./mc2->Integral());
  if (mc3) mc3->Scale(1./mc3->Integral());


  data->SetMarkerSize(1.2);
  data->SetLineWidth(2.5);
  if(mc1) mc1->SetLineWidth(2.5);
  if(mc2) mc2->SetLineWidth(2.5);
  if(mc3) mc3->SetLineWidth(2.5);


  data->SetMarkerColor(kBlack);
  if(mc1) mc1->SetMarkerColor(kRed);
  if(mc2) mc2->SetMarkerColor(kBlue);
  if(mc3) mc3->SetMarkerColor(kGreen+2);

  data->SetLineColor(kBlack);
  if(mc1) mc1->SetLineColor(kRed);
  if(mc2) mc2->SetLineColor(kBlue);
  if(mc3) mc3->SetLineColor(kGreen+2);

  data->SetTitle("unbias trigger");
  if(mc1) mc1->SetTitle("Pythia6 Z2*");
  if(mc2) mc2->SetTitle("Pythia8 Monash Tune");
  if(mc3) mc3->SetTitle("Pythia8 MBR Tune");

  data->GetYaxis()->SetRangeUser(1e-6,1e2);
  data->GetXaxis()->SetRangeUser(-10,500);
  data->GetXaxis()->SetNdivisions(504);
  data->GetXaxis()->SetTitle("no. of towers");
  data->GetYaxis()->SetTitle("N/N_{tot}");

  TCanvas* c1 = new TCanvas;
  data->Draw("");
  if(mc1) mc1->Draw("HIST SAME");
  if(mc2) mc2->Draw("HIST SAME");
  if(mc3) mc3->Draw("HIST SAME");
  TLegend* leg = c1->BuildLegend(0.22,0.72,0.62,0.92);
  SetLegAtt(leg);
  leg->Draw();

  if (type=="single")
    CMSText(3,0,1,"Single-arm selection");
  if (type=="double")
    CMSText(3,0,1,"Double-arm selection");
  if (type=="singleth")
    CMSText(3,0,1,"Single-arm selection","E_{tower}>0.5 GeV");


  c1->SetLogy();
  c1->SaveAs((string("plots/hf_perf_hf_norechits_")+type+string(".eps")).c_str());
  c1->SaveAs((string("plots/hf_perf_hf_norechits_")+type+string(".pdf")).c_str());
}
