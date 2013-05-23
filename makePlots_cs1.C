
const double cut_value_single = 6.6.;
const double cut_value_double = 2.2.;

void makePlots_cs1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<int> run_num;
  run_num.push_back(210614);
  run_num.push_back(210885);
  run_num.push_back(211000);
  run_num.push_back(211001);
  run_num.push_back(211032);
  run_num.push_back(211256);
  run_num.push_back(211371);
  run_num.push_back(211390);
  run_num.push_back(211532);
  run_num.push_back(211538);
  run_num.push_back(211607);

  //Histogram with all runs
  TH1D* h_run_single = new TH1D("h_run_single","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_single");
  TH1D* h_run_double = new TH1D("h_run_double","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_double");
  TH1D* h_run_single_Pbp = new TH1D("h_run_single_Pbp","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_single");
  TH1D* h_run_double_Pbp = new TH1D("h_run_double_Pbp","",run_num.size(),-0.5,run_num.size()-0.5);// h_run_single->SetName("h_run_double");

  for (int run=0; run<int(run_num.size()); run++)
    {
      cout << endl << " Processing ... run: " << run_num[run] << endl << endl;
      TFile* file = TFile::Open("histos_error_p.root");
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
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_single";
      TH1D* h_zb_single=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_cut_double";
      TH1D* h_zb_double=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_all_lumi";
      TH1D* h_hf_noise_all_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_selected_single_lumi";
      TH1D* h_hf_noise_selected_single_lumi=(TH1D*)file->Get(filename_ss.str().c_str());
      filename_ss.str(""); filename_ss << "data" << runname << "/data" << runname << "_h_hf_noise_selected_double_lumi";
      TH1D* h_hf_noise_selected_double_lumi=(TH1D*)file->Get(filename_ss.str().c_str());


      TFile f("plots/corr_factors.root");
      TVectorD* f_em = NULL;
      TVectorD* f_mc = NULL;
      TVectorD* f_all = NULL;
      f_em = (TVectorD*)f.Get("corr_fac_em");
      f_mc = (TVectorD*)f.Get("corr_fac_mc");
      f_all = (TVectorD*)f.Get("corr_fac_all");
      if(!f_all || !f_mc || !f_em) {cerr << "error" << endl; return;}
      cout << "MC Eff. Correction:" << endl;
      f_mc->Print();

      //Noise
      double f_noise_single = h_noise_single->GetBinContent(h_noise_single->FindBin(cut_value_single))/h_noise_single->GetBinContent(1);
      double f_noise_double = h_noise_double->GetBinContent(h_noise_double->FindBin(cut_value_double))/h_noise_double->GetBinContent(1);
      h_hf_noise_selected_single_lumi->Divide(h_hf_noise_all_lumi);
      h_hf_noise_selected_single_lumi->Rebin(20);
      h_hf_noise_selected_single_lumi->Scale(1./20.);
      h_hf_noise_selected_double_lumi->Divide(h_hf_noise_all_lumi);
      h_hf_noise_selected_double_lumi->Rebin(20);
      h_hf_noise_selected_double_lumi->Scale(1./20.);

      //ZB
      const double f_zb_single = h_zb_single->GetBinContent(1)/h_zb_single->GetBinContent(h_zb_single->FindBin(cut_value_single));
      const double f_zb_double = h_zb_double->GetBinContent(1)/h_zb_double->GetBinContent(h_zb_double->FindBin(cut_value_double));

      //Loop over LumiSections
      for(int i=0; i<=h_lumi->GetNbinsX();i++)
        {
      //Noise per LumiSection
//           cout << "f_noise=" << f_noise_single;
//           f_noise_single = h_hf_noise_selected_single_lumi->GetBinContent(h_hf_noise_selected_single_lumi->FindBin(h_lumi->GetBinCenter(i)));
//           cout << " f_noise=" << f_noise_single << endl;
//           f_noise_double = h_hf_noise_selected_double_lumi->GetBinContent(h_hf_noise_selected_double_lumi->FindBin(h_lumi->GetBinCenter(i)));
      //PILEUP!!
          const double lumiPerLS=h_lumi->GetBinContent(i);
          const double lambda = lumiPerLS*2.1/11246./23.31/296.;
          const double p1_single = (*f_mc)[0];
          const double p1_double = (*f_mc)[1];
          const double f_pileup_single = 1;
          const double f_pileup_double = 1;
          for (int k=2; k<100; k++)
            {
              double add_single =  p1_single/(1.-pow(1.-p1_single,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
              double add_double =  p1_double/(1.-pow(1.-p1_double,k)) * k * TMath::Poisson(k,lambda)/TMath::Poisson(1,lambda);
              f_pileup_single += add_single;
              f_pileup_double += add_double;
              if(add_single < 0.0001 && add_double < 0.0001)
                break;
            }
          //cout << "lambda=" << lambda << " " << f_pileup_single << endl;

      //Other Corrections
          if(lumiPerLS>0.)
            {
          const double n_cut_single = h_single->GetBinContent(i) / lumiPerLS;
          const double n_zb = (f_zb_single*n_cut_single);
          const double n_noise_single = f_noise_single * n_zb;
          const double n_em_single = (*f_em)[0] * n_cut_single;
          const double n_cut_double = h_double->GetBinContent(i) / lumiPerLS;
          const double n_zb_double = (f_zb_double*n_cut_double);
          const double n_noise_double = f_noise_double * (f_zb_double*n_cut_double);
          const double n_em_double = (*f_em)[1] * n_cut_double;
          
          if(i<200)
            {
              cout << n_cut_single<< " " << n_cut_single/(*f_mc)[0]<< " " << n_noise_single<< " " << n_em_single<< " " << ((1./f_pileup_single) - f_noise_single) << " " << f_pileup_single << endl;
              //cout << n_zb_double << " " << n_cut_double<< " " << n_cut_double/(*f_mc)[1]<< " " << n_noise_double<< " " << n_em_double<< " " << ((1./f_pileup_double) - f_noise_double) << endl;
            }
          if( ((1./f_pileup_double) - f_noise_double) != 0 && ((1./f_pileup_single) - f_noise_single) != 0)
            {
              h_single->SetBinContent(i, (n_cut_single/(*f_mc)[0]) - n_noise_single - n_em_single) / ((1./f_pileup_single) - f_noise_single);
              h_double->SetBinContent(i, (n_cut_double/(*f_mc)[1]) - n_noise_double - n_em_double) / ((1./f_pileup_double) - f_noise_double);
              h_single->SetBinError(i, h_single->GetBinError(i) / lumiPerLS);
              h_double->SetBinError(i, h_double->GetBinError(i) / lumiPerLS);
            }
          else cerr << "div 0 in bin: " << i << endl;
            }
          else if (lumiPerLS<0.) cerr << "lumi neg: " << i << endl;
        }


      
      //h_single->Scale(1./(*f_all)[0]);
      //h_double->Scale(1./(*f_all)[1]);
      


      h_double->SetMarkerColor(kRed);
      h_double->SetLineColor(kRed);
      
      h_single->SetLineWidth(2);
      h_double->SetLineWidth(2);

      h_single->SetTitle("single-arm;lumisection;#sigma_{inel} / b");
      h_double->SetTitle("double-arm;lumisection;#sigma_{inel} / b");

      if(run_num[run] == 210885)
        TCanvas* c1 = new TCanvas;

      TFitResultPtr fit_single = h_single->Fit("pol1","QS");
      TFitResultPtr fit_double = h_double->Fit("pol1","QS");

//       int rebin = 1;
//       h_single->Rebin(rebin);
//       h_single->Scale(1./double(rebin));
      
//       h_double->Rebin(rebin);
//       h_double->Scale(1./double(rebin));

      h_single->GetXaxis()->SetRangeUser(0,1200);
      h_double->GetXaxis()->SetRangeUser(0,1200);
      h_single->GetYaxis()->SetRangeUser(2.0,2.8);
      h_double->GetYaxis()->SetRangeUser(2.0,2.8);

      if(run_num[run] == 210885)
        {
          h_single->Draw();
          h_double->Draw("SAME");
          CMSPreliminary();
          DataText(false,true,"run #210885");
          c1->SaveAs((string("plots/CS_run")+string(".eps")).c_str());
          c1->SaveAs((string("plots/CS_run")+string(".pdf")).c_str());
          c1->SaveAs((string("plots/CS_run")+string(".png")).c_str());
        }

      if(run_num[run] <= 211256) //pPb or Pbp?
        {
          h_run_single->SetBinContent(run+1,run,fit_single->Parameter(0));
          h_run_single->SetBinError(run+1,0,fit_single->ParError(0));
          
          h_run_double->SetBinContent(run+1,run,fit_double->Parameter(0));
          h_run_double->SetBinError(run+1,0,fit_double->ParError(0));
        }
      else
        {
          h_run_single_Pbp->SetBinContent(run+1,run,fit_single->Parameter(0));
          h_run_single_Pbp->SetBinError(run+1,0,fit_single->ParError(0));
          
          h_run_double_Pbp->SetBinContent(run+1,run,fit_double->Parameter(0));
          h_run_double_Pbp->SetBinError(run+1,0,fit_double->ParError(0));
        }
      h_run_double->GetXaxis()->SetBinLabel(run+1,runname.c_str());
      cout << "Single: " << fit_single->Parameter(0) << " --- Double: " << fit_double->Parameter(0) << endl;

      if(run_num[run] == 210885)
        c1 = new TCanvas; //please use a new one to paint your fit crap
    }
  h_run_double->GetYaxis()->SetRangeUser(1.8,2.5);
  h_run_double->SetLineColor(kRed);
  h_run_double->SetMarkerColor(kRed);
  h_run_double->SetTitle("double-arm (pPb);run number;#sigma_{inel} / b");

  h_run_single->SetLineColor(kRed-1);
  h_run_single->SetMarkerColor(kRed-1);
  h_run_single->SetMarkerStyle(34);
  h_run_single->SetTitle("single-arm (pPb)");

  h_run_double_Pbp->SetLineColor(kBlue);
  h_run_double_Pbp->SetMarkerColor(kBlue);
  h_run_double_Pbp->SetTitle("double-arm (Pbp);run number;#sigma_{inel} / b");

  h_run_single_Pbp->SetLineColor(kBlue-1);
  h_run_single_Pbp->SetMarkerColor(kBlue-1);
  h_run_single_Pbp->SetMarkerStyle(34);
  h_run_single_Pbp->SetTitle("single-arm (Pbp)");



  TCanvas* can3 = new TCanvas;
  h_run_double->Draw("P");
  h_run_single->Draw("PSAME");
  h_run_double_Pbp->Draw("PSAME");
  h_run_single_Pbp->Draw("PSAME");
  TLegend* leg3 = can3->BuildLegend(0.6,0.7,0.85,0.9);
  leg3->SetFillColor(kWhite);
  leg3->Draw();
  CMSPreliminary();
  DataText(true,true);
  can3->SaveAs((string("plots/CS_runs")+string(".eps")).c_str());
  can3->SaveAs((string("plots/CS_runs")+string(".pdf")).c_str());
  can3->SaveAs((string("plots/CS_runs")+string(".png")).c_str());

  TCanvas* can4 = new TCanvas;
  h_run_double->Draw("P");
  h_run_double_Pbp->Draw("PSAME");
  TLegend* leg4 = can3->BuildLegend(0.6,0.7,0.85,0.9);
  leg4->SetFillColor(kWhite);
  leg4->Draw();
  CMSPreliminary();
  DataText(true,true);
  can4->SaveAs((string("plots/CS_runs_pas")+string(".eps")).c_str());
  can4->SaveAs((string("plots/CS_runs_pas")+string(".pdf")).c_str());
  can4->SaveAs((string("plots/CS_runs_pas")+string(".png")).c_str());


}
