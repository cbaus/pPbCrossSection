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
      TFile* file = TFile::Open("histos_mc.root");

      TH1D* h_eff=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_eff")).c_str());
      TH1D* h_effgain_single=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_effgain_single")).c_str());
      TH1D* h_effgain_double=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_effgain_double")).c_str());


      h_effgain_single->Scale(1./double(h_eff->GetBinContent(1)));
      h_effgain_double->Scale(1./double(h_eff->GetBinContent(1)));
      h_eff->Scale(1./double(h_eff->GetBinContent(1))); //order important

      h_eff->SetMarkerSize(1.2);
      h_effgain_single->SetMarkerSize(1.2);
      h_effgain_double->SetMarkerSize(1.2);

      cerr << "WARNING" << endl;
      h_eff->GetXaxis()->SetBinLabel(7,"no. of tracks #geq 1");
      h_effgain_single->GetXaxis()->SetBinLabel(6,"no. of tracks #geq 1");
      h_effgain_double->GetXaxis()->SetBinLabel(5,"no. of tracks #geq 1");

      h_eff->GetXaxis()->LabelsOption("u");
      h_effgain_single->GetXaxis()->LabelsOption("u");
      h_effgain_double->GetXaxis()->LabelsOption("u");

      h_eff->SetMarkerStyle(21);
      h_effgain_single->SetMarkerStyle(21);
      h_effgain_double->SetMarkerStyle(21);

      h_eff->GetYaxis()->SetRangeUser(0.9,1.01);
      h_effgain_single->GetYaxis()->SetRangeUser(0.0,0.10);
      h_effgain_double->GetYaxis()->SetRangeUser(0.0,0.10);

      h_eff->SetTitle("EPOS-LHC;;selection efficiency");
      h_effgain_single->SetTitle("EPOS-LHC;;selection efficiency gain");
      h_effgain_double->SetTitle("EPOS-LHC;;selection efficiency gain");

      TPaveText* leg = new TPaveText(0.6,0.73,0.75,0.78,"NDC b t l");
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddText("EPOS-LHC");
      leg->SetTextSize(0.033);

      TCanvas* c1 = new TCanvas;
      h_eff->Draw("HIST P");
      leg->Draw();
      MCText(0,1);
      c1->SaveAs((string("plots/mc_eff_")+list[i]+string(".pdf")).c_str());

      TCanvas* c2 = new TCanvas;
      h_effgain_single->Draw("HIST P");
      leg->Draw();
      MCText(0,1);
      c2->SaveAs((string("plots/mc_effgain_single_")+list[i]+string(".pdf")).c_str());

      TCanvas* c3 = new TCanvas;
      h_effgain_double->Draw("HIST P");
      leg->Draw();
      MCText(0,1);
      c3->SaveAs((string("plots/mc_effgain_double_")+list[i]+string(".pdf")).c_str());
    }
}
