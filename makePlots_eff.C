void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void SetAttributes(TH1D* theGraph, int colour, int marker)
{
  theGraph->SetMarkerSize(1.8);
  theGraph->SetLineWidth(1.8);
  theGraph->SetFillColor(kWhite);
  theGraph->SetLineColor(colour);
  theGraph->SetMarkerColor(colour);
  theGraph->SetMarkerStyle(marker);
}

void makePlots_eff()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  name.push_back(string("HIJING 1.383"));
  name.push_back(string("QGSJetII-04"));

  TCanvas* c1 = new TCanvas;
  TCanvas* c2 = new TCanvas;
  TCanvas* c3 = new TCanvas;
  TCanvas* c4 = new TCanvas;

  TLegend* legan = new TLegend(0.23,0.23,0.43,0.45);
  TLegend* legpas = new TLegend(0.23,0.23,0.43,0.38);
  SetLegAtt(legan,1.5);
  SetLegAtt(legpas,1.5);
  for(int i=0; i<int(list.size()); i++)
    {
      cout << i+1 << "/" << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");

      TH1D* h_eff=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_eff")).c_str());
      TH1D* h_effgain_single=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_effgain_single")).c_str());
      TH1D* h_effgain_double=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_effgain_double")).c_str());

      int colour = kBlack;
      int marker = 21;
      if(list[i] == string("Epos"))
        {
         colour = kGreen-1;
         marker = 22;
        }
      if(list[i] == string("QGSJetII"))
        {
         colour = kBlue;
         marker = 34;
        }

      SetAttributes(h_eff,colour,marker);
      SetAttributes(h_effgain_single,colour,marker);
      SetAttributes(h_effgain_double,colour,marker);

      h_effgain_single->Scale(1./double(h_eff->GetBinContent(1)));
      h_effgain_double->Scale(1./double(h_eff->GetBinContent(1)));
      h_eff->Scale(1./double(h_eff->GetBinContent(1))); //order important

      h_eff->GetXaxis()->LabelsOption("u");
      h_effgain_single->GetXaxis()->LabelsOption("u");
      h_effgain_double->GetXaxis()->LabelsOption("u");

      h_eff->GetYaxis()->SetRangeUser(0.8,1.01);
      h_effgain_single->GetYaxis()->SetRangeUser(0.0,0.15);
      h_effgain_double->GetYaxis()->SetRangeUser(0.0,0.15);

      h_eff->SetTitle(name[i].c_str());
      h_effgain_single->SetTitle(name[i].c_str());
      h_effgain_double->SetTitle(name[i].c_str());

      h_eff->SetTitle(";;selection efficiency");
      h_effgain_single->SetTitle(";;selection efficiency gain");
      h_effgain_double->SetTitle(";;selection efficiency gain");

      legan->AddEntry(h_eff,name[i].c_str(),"p");
      if(list[i] != string("Hijing"))
        legpas->AddEntry(h_eff,name[i].c_str(),"p");

      c1->cd();
      h_eff->Draw(i<1?"HIST P":"HIST P SAME");

      c2->cd();
      h_effgain_single->Draw(i<1?"HIST P":"HIST P SAME");

      c3->cd();
      h_effgain_double->Draw(i<1?"HIST P":"HIST P SAME");

      c4->cd();
      h_eff->GetXaxis()->SetRange(1,h_eff->GetNbinsX()-1);
      if(list[i] != string("Hijing"))
        h_eff->Draw(i<1?"HIST P":"HIST P SAME");
    }
  c1->cd();
  CMSText(0,0,1);
  legan->Draw();
  c1->SaveAs((string("plots/mc_eff")+string(".pdf")).c_str());

  c2->cd();
  CMSText(0,0,1);
  legan->Draw();
  c2->SaveAs((string("plots/mc_effgain_single")+string(".pdf")).c_str());

  c3->cd();
  CMSText(0,0,1);
  legan->Draw();
  c3->SaveAs((string("plots/mc_effgain_double")+string(".pdf")).c_str());

  c4->cd();
  CMSText(0,0,1);
  legpas->Draw();
  c4->SaveAs((string("plots/mc_eff")+string("_PAS.pdf")).c_str());
}
