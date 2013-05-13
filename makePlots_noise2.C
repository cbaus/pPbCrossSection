TVectorD hf_m_cuts_light_noise(12);
TVectorD hf_p_cuts_light_noise(12);
TVectorD hf_m_cuts_medium_noise(12);
TVectorD hf_p_cuts_medium_noise(12);
TVectorD hf_m_cuts_heavy_noise(12);
TVectorD hf_p_cuts_heavy_noise(12);

void makePlots_noise2()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_noise.root");
  TCanvas* c1 = new TCanvas;
  TCanvas* c0;
  c1->Divide(6,4);

  vector<double> x;
  vector<double> xerr;
  vector<double> y;
  vector<double> yerr;
  vector<double> yeff;

  string add = "noise";
  int count=0;
  for (int i=-41; i<=41;count++)
    {
      cout << " Processing ... " << i << " count: " << count << endl;
      ostringstream hf_str;
      hf_str  << add << "h_hf_rings_cut_single_" << i;
      TH1D* h_hf_rings_cut_single = NULL;
      h_hf_rings_cut_single = (TH1D*)(file->Get((add+string("/")+hf_str.str()).c_str()));

      if(h_hf_rings_cut_single != NULL)
        {
          c1->cd(i<0?count+1:count+13)->SetLogy();
          h_hf_rings_cut_single->Scale(1./double(h_hf_rings_cut_single->GetBinContent(1)));

          TVectorD* tvecp;
          TVectorD* tvecm;
          double cut_light_value = 1.-pow(1.-0.06,1./24.);
          double cut_medium_value = 1.-pow(1.-0.03,1./24.);
          double cut_heavy_value = 1.-pow(1.-0.01,1./24.);
          double cut_value;
          for(int n=0; n<3; n++)
            {
              if(n==0)
                {
                  cut_value = cut_light_value;
                  tvecm = hf_m_cuts_light_noise;
                  tvecp = hf_p_cuts_light_noise;
                }
              else if(n==1)
                {
                  cut_value = cut_medium_value;
                  tvecm = hf_m_cuts_medium_noise;
                  tvecp = hf_p_cuts_medium_noise;
                }
              else if(n==2)
                {
                  cut_value = cut_heavy_value;
                  tvecm = hf_m_cuts_heavy_noise;
                  tvecp = hf_p_cuts_heavy_noise;
                }

              cout << "Cut Value = " << cut_value << endl;

              TF1* fitfunc = new TF1("a","expo");
              TFitResultPtr fit = h_hf_rings_cut_single->Fit(fitfunc,"QS","",2,4);
              if (int(fit)!=0)
                {
                  cerr << "  ! FIT BROKEN" << endl;
                  goto 100;
                }
              double c = fit->Parameter(0);
              double slope = fit->Parameter(1);
              double inverse = (log(cut_value)-c)/slope;
              cout << "Ring " << i << " = " << inverse << endl;
              int j=0;
              double inversenew=inverse;
              do
                {
                  inverse=inversenew;
                  fit = h_hf_rings_cut_single->Fit(fitfunc,"QS","",inverse-1,inverse+1);
                  c = fit->Parameter(0);
                  slope = fit->Parameter(1);
                  inversenew = (log(cut_value)-c)/slope;
                  cout << "Ring " << i << " = " << inversenew << endl;
                } while(fabs(inverse-inversenew)>0.05 && j++<10);

            100:
              h_hf_rings_cut_single->SetTitle("; cut value / GeV;event fraction");
              h_hf_rings_cut_single->GetXaxis()->SetRangeUser(0.,10.);
              h_hf_rings_cut_single->GetYaxis()->SetRangeUser(0.0001,1.01);
              h_hf_rings_cut_single->Draw();

              if(i<0)
                (*tvecm)[count] = inverse;
              else
                (*tvecp)[count] = inverse;
            }
        }
      else
        {
          cout << "skipping ... " << i << endl;
        }
      if(i==-30)
        {
          i=30;
          count = -1;
        }
      else
        i++;
    }


  TFile f("plots/hf_cuts_noise.root","recreate");
  hf_m_cuts_light_noise.Write("hf_m_cuts_light_noise");
  hf_p_cuts_light_noise.Write("hf_p_cuts_light_noise");
  hf_m_cuts_medium_noise.Write("hf_m_cuts_medium_noise");
  hf_p_cuts_medium_noise.Write("hf_p_cuts_medium_noise");
  hf_m_cuts_heavy_noise.Write("hf_m_cuts_heavy_noise");
  hf_p_cuts_heavy_noise.Write("hf_p_cuts_heavy_noise");

  f.Close();

}
