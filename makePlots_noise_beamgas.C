void makePlots_noise_beamgas()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TCanvas* can1 = new TCanvas;

      TFile* file = TFile::Open("histos_noise2.root");
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "noise/noise_h_hf_hits_noise";
      cout << filename_ss.str() << endl;         
      TH1D* noise=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "noise/noise_h_hf_hits_beamgas";
      cout << filename_ss.str() << endl;         
      TH1D* beamgas=file->Get(filename_ss.str().c_str());
      
      double integral = noise->GetBinContent(1);
      //integral = 1;
      noise->Scale(1./integral);
      beamgas->Scale(1./integral / 84. * (3563.-296.));
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

      noise->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);
      beamgas->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);

      noise->GetYaxis()->SetRangeUser(1e-5,1.3);
      beamgas->GetYaxis()->SetRangeUser(1e-5,1.3);

      noise->SetTitle("noise;E_{HF} [GeV];event fraction");
      beamgas->SetTitle("beam gas;E_{HF} [GeV];event fraction");

      can1->cd();
      noise->Draw("HIST L");
      beamgas->Draw("SAME");

      can1->SetLogy();
  TLegend* leg1 = can1->BuildLegend(0.55,0.6,0.95,0.7);
  SetLegAtt(leg1);
  leg1->SetFillColor(kWhite);
  leg1->Draw();

  CMSText(1,0,1,"single-arm selection");
  can1->SaveAs((string("plots/noise_beamgas")+string(".pdf")).c_str());
}

