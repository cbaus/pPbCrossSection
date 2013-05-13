void makePlots_cs1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<int> run_num;
  run_num.push_back(210614);
  run_num.push_back(210837);
  run_num.push_back(210885);
  run_num.push_back(211607);

  //Histogram with all runs
  TGraphErrors* h_run_single = new TGraphErrors(run_num.size()); h_run_single.SetName("h_run_single");
  TGraphErrors* h_run_double = new TGraphErrors(run_num.size()); h_run_single.SetName("h_run_double");

  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << runname << "_h_run_events_single_";
      TH1D* a=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << runname << "_h_run_events_double_";
      TH1D* a2=file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << runname << "_h_run_events_lumi_";
      TH1D* alumi=file->Get(filename_ss.str().c_str());

      TFile f("plots/corr_factors.root");
      TVectorD* corr_fac_em = NULL;
      TVectorD* corr_fac_noise = NULL;
      TVectorD* corr_fac_mc = NULL;
      TVectorD* corr_fac_all = NULL;
      corr_fac_em = (TVectorD*)f.Get("corr_fac_em");
      corr_fac_noise = (TVectorD*)f.Get("corr_fac_noise");
      corr_fac_mc = (TVectorD*)f.Get("corr_fac_mc");
      corr_fac_all = (TVectorD*)f.Get("corr_fac_all");
      if(!corr_fac_all || !corr_fac_noise || !corr_fac_mc || !corr_fac_em) {cerr << "error" << endl; return;}
      cout << "All correction:" << endl;
      corr_fac_all->Print();

      //PILEUP!!
      for(int i=0; i<=alumi->GetNbinsX();i++)
        {
          double lumiPerLS=alumi->GetBinContent(i);
          const double lambda = lumiPerLS*2.1/11246./23.31/296.;
          const double p1 = (*corr_fac_mc)[0];
          double pileupcorr = 1;
          for (int k=2; k<100; k++)
            {
              double add =  p1/(1.-pow(1.-p1,k)) * k * TMath::Poisson(k,lambda);;
              pileupcorr += add;
              if(add < 0.0001)
                break;
            }
          //cout << "lambda=" << lambda << " " <<pileupcorr << endl;
          a->SetBinContent(i,a->GetBinContent(i)*pileupcorr);
        }
      a->Scale(1./(*corr_fac_all)[0]);
      a2->Scale(1./(*corr_fac_all)[1]);

      a->SetLineWidth(2);
      a2->SetLineWidth(2);

      a->SetTitle("single-arm;lumisection;#sigma_{inel} / b");
      a2->SetTitle("double-arm;lumisection;#sigma_{inel} / b");

      TCanvas* c1 = new TCanvas;
      a->Draw();
      TFitResultPtr fit_single = a->Fit("pol0","QS");
      h_run_single->SetPoint(run,run_num[run],fit_single->Parameter(0));
      h_run_single->SetPointError(run,0,fit_single->ParError(0));

      TCanvas* c2 = new TCanvas;
      a2->Draw("");
      TFitResultPtr fit_double = a2->Fit("pol0","QS");
      h_run_double->SetPoint(run,run_num[run],fit_double->Parameter(0));
      h_run_double->SetPointError(run,0,fit_double->ParError(0));
      file->Close();
    }
  h_run_double->SetLineColor(kRed);
  h_run_double->SetMarkerColor(kRed);

  TCanvas* c3 = new TCanvas;
  h_run_single->Draw("AP");
  h_run_double->Draw("P");
}
