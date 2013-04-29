void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void makePlots_eff()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> list;
  list.push_back(string("Epos"));
  //list.push_back(string("Hijing"));
  //list.push_back(string("QGSJetII"));

  
  for(int i=0; i<int(list.size()); i++)
    {
      cout << i << " " << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");

      cout << string(list[i]+string("/")+list[i]) << endl;
      TH1D* a=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_eff")).c_str());


      a->Scale(1./double(a->GetBinContent(1)));

      a->SetMarkerSize(1.2);
      a->SetMarkerStyle(21);
      
      a->GetYaxis()->SetRangeUser(0.9,1.01);

      a->SetTitle(";;selection efficiency");

      TCanvas* c1 = new TCanvas;
      a->Draw("HIST P");
      c1->SaveAs((string("plots/mc_eff_")+list[i]+string(".png")).c_str());
      c1->SaveAs((string("plots/mc_eff_")+list[i]+string(".eps")).c_str());
      c1->SaveAs((string("plots/mc_eff_")+list[i]+string(".pdf")).c_str());
  
    }
}
