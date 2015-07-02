#include <sstream>
void makePlots_noise_beamgas()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TCanvas* can1 = new TCanvas;

      TFile* file = TFile::Open("histos_noise.root");
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "data247324/data247324_h_hf_hits_noise";
      cout << filename_ss.str() << endl;
      TH1D* noise=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data247324/data247324_h_hf_hits_beamgas";
      cout << filename_ss.str() << endl;
      TH1D* beamgas=file->Get(filename_ss.str().c_str());

      double integral = noise->Integral();
      //integral = 1;
      noise->Scale(1./noise->GetBinWidth(1)/noise->Integral());
      beamgas->Scale(1./beamgas->GetBinWidth(1)/beamgas->Integral());
      // int rebin = 50;
      // noise->Rebin(rebin);
      // noise->Scale(1./double(rebin));

      // beamgas->Rebin(rebin);
      // beamgas->Scale(1./double(rebin));


      noise->SetMarkerColor(kBlack);
      beamgas->SetMarkerColor(kRed);
      noise->SetMarkerStyle(21);
      beamgas->SetMarkerStyle(25);

      noise->SetLineWidth(3);
      beamgas->SetLineWidth(2);
      noise->SetLineColor(noise->GetMarkerColor());
      beamgas->SetLineColor(beamgas->GetMarkerColor());
      noise->SetLineStyle(7);
      beamgas->SetLineStyle(1);

      // noise->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);
      // beamgas->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);

      noise->GetYaxis()->SetRangeUser(1e-5,1.3);
      beamgas->GetYaxis()->SetRangeUser(1e-5,1.3);

      double overflow_noise = noise->GetBinContent(noise->GetNbinsX());
      double overflow_beamgas = beamgas->GetBinContent(beamgas->GetNbinsX());

      ostringstream ss_title_noise; ss_title_noise << "Noise (overflow: " << scientific << setprecision(0) << overflow_noise << ");#it{E}_{HF} [GeV];event fraction";
      ostringstream ss_title_beamgas; ss_title_beamgas << "Beam gas (overflow: " << scientific << setprecision(0) << overflow_beamgas << ");#it{E}_{HF} [GeV];event fraction";

      noise->SetTitle(ss_title_noise.str().c_str());
      beamgas->SetTitle(ss_title_beamgas.str().c_str());

      can1->cd();
      noise->Draw("HIST L");
      beamgas->Draw("SAME");

      can1->SetLogy();
      TLegend* leg1 = can1->BuildLegend(0.45,0.6,0.85,0.7);
  SetLegAtt(leg1);
  leg1->SetFillColor(kWhite);
  leg1->Draw();

  CMSText(1,0,2,"Single-arm selection");
  can1->SaveAs((string("plots/noise_beamgas")+string(".eps")).c_str());
}
