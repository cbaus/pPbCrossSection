
void makePlots_cs1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<int> run_num;
  //run_num.push_back(210614);
  run_num.push_back(210885);
//   run_num.push_back(211000);
//   run_num.push_back(211001);
//   run_num.push_back(211032);
//   run_num.push_back(211256);
//   run_num.push_back(211371);
//   run_num.push_back(211390);
//   run_num.push_back(211532);
//   run_num.push_back(211538);
//   run_num.push_back(211607);

  //Histogram with all runs
  TH1D* h_run_single = new TH1D("h_run_single","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_single");
  TH1D* h_run_double = new TH1D("h_run_double","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_double");

  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos.root");
      ostringstream runname_ss;
      runname_ss << run_num[run];
      string runname = runname_ss.str();
      ostringstream filename_ss;
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_single_";
      TH1D* h_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_double_";
      TH1D* h_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_run_events_lumi_";
      TH1D* h_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single_noise";
      TH1D* h_noise_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double_noise";
      TH1D* h_noise_double=(TH1D*)file->Get(filename_ss.str().c_str());


      TFile f("plots/corr_factors.root");
      TVectorD* f_em = NULL;
      TVectorD* f_mc = NULL;
      TVectorD* f_all = NULL;
      f_em = (TVectorD*)f.Get("corr_fac_em");
      f_mc = (TVectorD*)f.Get("corr_fac_mc");
      f_all = (TVectorD*)f.Get("corr_fac_all");
      if(!f_all || !f_mc || !f_em) {cerr << "error" << endl; return;}
      cout << "All correction:" << endl;
      f_all->Print();

      //Noise
      const double f_noise_single = h_noise_single->GetBinContent(31)/h_noise_single->GetBinContent(1);
      const double f_noise_double = h_noise_double->GetBinContent(21)/h_noise_double->GetBinContent(1);

      //PILEUP!!
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
          double lumiPerLS=h_lumi->GetBinContent(i);
          const double lambda = lumiPerLS*2.1/11246./23.31/296.;
          const double p1_single = (*f_mc)[0];
          const double p1_double = (*f_mc)[1];
          double f_pileup_single = 1;
          double f_pileup_double = 1;
          for (int k=2; k<100; k++)
            {
              double add_single =  p1_single/(1.-pow(1.-p1_single,k)) * k * TMath::Poisson(k,lambda);;
              double add_double =  p1_double/(1.-pow(1.-p1_double,k)) * k * TMath::Poisson(k,lambda);;
              f_pileup_single += add_single;
              f_pileup_double += add_double;
              if(add_single < 0.0001 && add_double < 0.0001)
                break;
            }
          //cout << "lambda=" << lambda << " " <<f_pileup << endl;
          double n_cut_single = h_single->GetBinContent(i);
          double n_noise_single = f_noise_single * n_cut_single;
          double n_em_single = (*f_em)[0] * n_cut_single;
          h_single->SetBinContent(i, (n_cut_single/(*f_mc)[0]) - n_noise_single - n_em_single) / ((1./f_pileup_single) - f_noise_single);
          double n_cut_double = h_double->GetBinContent(i);
          double n_noise_double = f_noise_double * n_cut_double;
          double n_em_double = (*f_em)[1] * n_cut_double;
          if(i<100)
            {
              cout << n_cut_single<< " " << n_cut_single/(*f_mc)[0]<< " " << n_noise_single<< " " << n_em_single<< " " << ((1./f_pileup_single) - f_noise_single) << endl;
              cout << n_cut_double<< " " << n_cut_double/(*f_mc)[1]<< " " << n_noise_double<< " " << n_em_double<< " " << ((1./f_pileup_double) - f_noise_double) << endl;
            }
          h_double->SetBinContent(i, (n_cut_double/(*f_mc)[1]) - n_noise_double - n_em_double) / ((1./f_pileup_double) - f_noise_double);
        }


      
      //h_single->Scale(1./(*f_all)[0]);
      //h_double->Scale(1./(*f_all)[1]);
      
//       int rebin = 50;
//       h_single->Rebin(rebin);
//       h_single->Scale(1./double(rebin));
      
//       h_double->Rebin(rebin);
//       h_double->Scale(1./double(rebin));

      h_single->GetXaxis()->SetRangeUser(0,1200);
      h_double->GetXaxis()->SetRangeUser(0,1200);

      h_single->SetLineWidth(2);
      h_double->SetLineWidth(2);

      h_single->SetTitle("single-arm;lumisection;#sigma_{inel} / b");
      h_double->SetTitle("double-arm;lumisection;#sigma_{inel} / b");

      //TCanvas* c1 = new TCanvas;
      //h_single->Draw();
      TFitResultPtr fit_single = h_single->Fit("pol0","QS");
      h_run_single->SetBinContent(run+1,run,fit_single->Parameter(0));
      h_run_single->SetBinError(run+1,0,fit_single->ParError(0));

      //TCanvas* c2 = new TCanvas;
      //h_double->Draw("");
      TFitResultPtr fit_double = h_double->Fit("pol0","QS");
      h_run_double->SetBinContent(run+1,run,fit_double->Parameter(0));
      h_run_double->SetBinError(run+1,0,fit_double->ParError(0));
      h_run_double->GetXaxis()->SetBinLabel(run+1,runname.c_str());
      //file->Close();
      cout << "Single: " << fit_single->Parameter(0) << " --- Double: " << fit_double->Parameter(0) << endl;
    }
  h_run_double->SetLineColor(kRed);
  h_run_double->SetMarkerColor(kRed);
  h_run_double->SetTitle("double-arm;run number;#sigma_{inel} / b");

  h_run_single->SetTitle("single-arm");

  h_run_double->GetYaxis()->SetRangeUser(2.0,2.5);


  TCanvas* can3 = new TCanvas;
  h_run_double->Draw("P");
  h_run_single->Draw("PSAME");
  TLegend* leg3 = can3->BuildLegend(0.65,0.75,0.8,0.85);
  leg3->SetFillColor(kWhite);
  leg3->Draw();

}
