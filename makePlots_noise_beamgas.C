#include <sstream>
//void Show(string, TH1D*, TH1D*, TH1D*, string, string);
void Show(string prefix, TH1D* noise, TH1D* beamgas1=0, string title1="", TH1D* beamgas2=0, string title2="");
void makePlots_noise_beamgas()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  TFile* file = TFile::Open("histos_noise.root");
  ostringstream filename_ss;
  filename_ss.str(""); filename_ss << "data259163/data259163_h_hf_hits_noise";
  cout << filename_ss.str() << endl;
  TH1D* noise=(TH1D*)file->Get(filename_ss.str().c_str());
  filename_ss.str(""); filename_ss << "data259163/data259163_h_hf_hits_beamgas";
  TH1D* beamgas=(TH1D*)file->Get(filename_ss.str().c_str());
  filename_ss.str(""); filename_ss << "noisesim/noisesim_h_hf_hits_noise";
  cout << filename_ss.str() << endl;
  TH1D* noisesim=(TH1D*)file->Get(filename_ss.str().c_str());

  //Show("hf", noise, beamgas, "BPTX+XOR-", noisesim, "In simulation");
  Show("hf", noise, noisesim, "In simulation");

  filename_ss.str(""); filename_ss << "data259163/data259163_h_castor_sumE_noise";
  cout << filename_ss.str() << endl;
  TH1D* noise=(TH1D*)file->Get(filename_ss.str().c_str());
  filename_ss.str(""); filename_ss << "data259163/data259163_h_castor_sumE_beam_plus";
  cout << filename_ss.str() << endl;
  TH1D* plus=(TH1D*)file->Get(filename_ss.str().c_str());
  filename_ss.str(""); filename_ss << "data259163/data259163_h_castor_sumE_beam_minus";
  cout << filename_ss.str() << endl;
  TH1D* minus=(TH1D*)file->Get(filename_ss.str().c_str());

  filename_ss.str(""); filename_ss << "noisesim/noisesim_h_castor_sumE_noise";
  cout << filename_ss.str() << endl;
  TH1D* castorsim=(TH1D*)file->Get(filename_ss.str().c_str());


  //Show("castor", noise, plus, "BPTX+ only", minus, "BPTX- only");
  Show("castor", noise, castorsim, "In simulation");

}
void Show(string prefix, TH1D* noise, TH1D* beamgas1, string title1, TH1D* beamgas2, string title2)
{
  TCanvas* can1 = new TCanvas;

      double integral = noise->Integral();
      //integral = 1;
      noise->Scale(1./noise->Integral(0,noise->GetNbinsX()));
      if(beamgas1) beamgas1->Scale(1./beamgas1->Integral(0,beamgas1->GetNbinsX()));
      if(beamgas2) beamgas2->Scale(1./beamgas2->Integral(0,beamgas2->GetNbinsX()));
      // int rebin = 50;
      // noise->Rebin(rebin);
      // noise->Scale(1./double(rebin));

      // beamgas1->Rebin(rebin);
      // beamgas1->Scale(1./double(rebin));


      noise->SetMarkerColor(kBlack);
      if(beamgas1) beamgas1->SetMarkerColor(kRed);
      if(beamgas2) beamgas2->SetMarkerColor(kBlue);
      noise->SetMarkerStyle(21);
      if(beamgas1) beamgas1->SetMarkerStyle(25);
      if(beamgas2) beamgas2->SetMarkerStyle(24);

      noise->SetLineWidth(3);
      if(beamgas1) beamgas1->SetLineWidth(2);
      if(beamgas2) beamgas2->SetLineWidth(2);
      noise->SetLineColor(noise->GetMarkerColor());
      if(beamgas1) beamgas1->SetLineColor(beamgas1->GetMarkerColor());
      if(beamgas2) beamgas2->SetLineColor(beamgas2->GetMarkerColor());
      noise->SetLineStyle(7);
      if(beamgas1) beamgas1->SetLineStyle(1);
      if(beamgas2) beamgas2->SetLineStyle(1);

      // noise->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);
      // beamgas1->GetXaxis()->SetRange(2,noise->GetNbinsX()/2);

      noise->GetYaxis()->SetRangeUser(1e-6,10);
      // if(beamgas1) beamgas1->GetYaxis()->SetRangeUser(1e-6,1.3);
      // if(beamgas2) beamgas2->GetYaxis()->SetRangeUser(1e-6,1.3);

      double overflow_noise = noise->GetBinContent(noise->GetNbinsX());
      double overflow_beamgas1 = 0; if(beamgas1) overflow_beamgas1 = beamgas1->GetBinContent(beamgas1->GetNbinsX());
      double overflow_beamgas2 = 0; if(beamgas2) overflow_beamgas2 = beamgas2->GetBinContent(beamgas2->GetNbinsX());

      ostringstream ss_title_noise;
      ostringstream ss_title_beamgas1;
      ostringstream ss_title_beamgas2;

      {
        ss_title_noise << "No BPTX";
        if (overflow_noise)
          ss_title_noise << " (overflow: " << scientific << setprecision(0) << overflow_noise << ")";
        ss_title_noise << ";#it{E} [GeV];event fraction";
      }

      if(beamgas1)
        {
          ss_title_beamgas1 << title1;
          if (overflow_beamgas1)
            ss_title_beamgas1 << " (overflow: " << scientific << setprecision(0) << overflow_beamgas1 << ")";
          ss_title_beamgas1 << ";#it{E} [GeV];event fraction";
        }
      if(beamgas2)
        {
          ss_title_beamgas2 << title1;
          if (overflow_beamgas2)
            ss_title_beamgas2 << " (overflow: " << scientific << setprecision(0) << overflow_beamgas2 << ")";
          ss_title_beamgas2 << ";#it{E} [GeV];event fraction";
        }

      noise->SetTitle(ss_title_noise.str().c_str());
      if(beamgas1) beamgas1->SetTitle(ss_title_beamgas1.str().c_str());
      if(beamgas2) beamgas2->SetTitle(ss_title_beamgas2.str().c_str());

      can1->cd();
      noise->Draw("HIST L");
      if(beamgas1) beamgas1->Draw("SAME");
      if(beamgas2) beamgas2->Draw("SAME");

      can1->SetLogy();
      // can1->SetLogx();
      TLegend* leg1 = can1->BuildLegend(0.22,0.82,0.57,0.93);
  SetLegAtt(leg1);
  leg1->SetFillColor(kWhite);
  leg1->Draw();

  CMSText(3,prefix=="hf"?0:0,prefix=="hf"?1:1,prefix=="hf"?"HF hottest tower (single-arm)":"CASTOR (sum first 5 mod)","Random Trigger");
  can1->SaveAs((string("plots/noise_beamgas_") + prefix +string(".eps")).c_str());
}
