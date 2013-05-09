void makePlots_cs1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos.root");

  TH1D* a=file->Get("data2/data2210885_h_run_events_single_");
  TH1D* a2=file->Get("data2/data2210885_h_run_events_double_");
  TH1D* alumi=file->Get("data2/data2210885_h_run_events_lumi_");

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
      cout << "lambda=" << lambda << " " <<pileupcorr << endl;
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
  a->Fit("pol0");

  TCanvas* c2 = new TCanvas;
  a2->Draw("");
  a2->Fit("pol0");
}
