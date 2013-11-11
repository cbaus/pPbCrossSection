
void makePlots_noise3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TCanvas* can1 = new TCanvas;
  TCanvas* can2 = new TCanvas;

  vector<int> run_num;
  run_num.push_back(210676);
  run_num.push_back(210737);
  run_num.push_back(210738);
  run_num.push_back(210759);
  run_num.push_back(210818);
  run_num.push_back(210837);
  run_num.push_back(210855);
  run_num.push_back(210885);
  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos_noise2.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "noise" << "/noise" << runname << "_h_run_events_single";
      cout << filename_ss.str() << endl;
      TH1D* a=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "noise" << "/noise" << runname << "_h_run_events_double";
      TH1D* a2=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "noise" << "/noise" << runname << "_h_run_events";
      TH1D* atot=file->Get(filename_ss.str().c_str());

      a->Divide(atot);

      a2->Divide(atot);

      // int rebin = 50;
      // a->Rebin(rebin);
      // a->Scale(1./double(rebin)/a->GetBinWidth(1));

      // a2->Rebin(rebin);
      // a2->Scale(1./double(rebin)/a2->GetBinWidth(1));

      a->SetLineColor(run+1);
      a2->SetLineColor(run+1);

      a->SetMarkerColor(run+1);
      a2->SetMarkerColor(run+1);

      a->GetXaxis()->SetRangeUser(0,100000);
      a2->GetXaxis()->SetRangeUser(0,100000);

      runname_ss << ";L_{i} / (s b)^{-1};d#frac{N_{sel}}{N_{Noise}} / dL_{i}";
      a->SetTitle(runname_ss.str().c_str());
      a2->SetTitle(runname_ss.str().c_str());

      can1->cd();
      a->Draw(run==0?"":"SAME");

      can2->cd();
      a2->Draw(run==0?"":"SAME");
    }

  can1->cd();
  TLegend* leg1 = can1->BuildLegend(0.25,0.5,0.35,0.88);
  leg1->SetFillColor(kWhite);
  leg1->Draw();

  can2->cd();
  leg1->Draw();
  CMSPreliminary();
  can1->SaveAs((string("plots/noise_lumi_double")+string(".eps")).c_str());
  can1->SaveAs((string("plots/noise_lumi_double")+string(".pdf")).c_str());
  can1->SaveAs((string("plots/noise_lumi_double")+string(".png")).c_str());
}

