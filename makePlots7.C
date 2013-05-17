void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void makePlots7()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> list;
  //list.push_back(string("Epos"));
  //list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));

  
  for(int i=0; i<int(list.size()); i++)
    {
      cout << i << " " << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");

      cout << string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_single")) << endl;
      TH1D* a=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_single")).c_str());
      TH1D* b=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_double")).c_str());
      TH1D* c1=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_SD1")).c_str());
      TH1D* c2=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_SD2")).c_str());
      TH1D* d=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_DD")).c_str());
      TH1D* e=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_CD")).c_str());
      TH1D* f=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_ND")).c_str());
      TH1D* g=file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_all")).c_str());


      a->Scale(1./double(g->GetEntries()));
      b->Scale(1./double(g->GetEntries()));
      c1->Scale(1./double(g->GetEntries()));
      c2->Scale(1./double(g->GetEntries()));
      d->Scale(1./double(g->GetEntries()));
      e->Scale(1./double(g->GetEntries()));
      f->Scale(1./double(g->GetEntries()));
      g->Scale(1./double(g->GetEntries()));

      a->SetMarkerSize(1.2);
      b->SetMarkerSize(1.2);
      a->SetLineWidth(2.5);
      b->SetLineWidth(2.5);
      c1->SetLineWidth(2.5);
      c2->SetLineWidth(2.5);
      d->SetLineWidth(2.5);
      e->SetLineWidth(2.5);
      f->SetLineWidth(2.5);
      g->SetLineWidth(2.5);


      a->SetMarkerColor(kBlue+2);
      b->SetMarkerColor(kRed+2);
      c1->SetMarkerColor(kBlue+3);
      c2->SetMarkerColor(kBlue+1);
      d->SetMarkerColor(kRed+2);
      e->SetMarkerColor(kGreen+3);
      f->SetMarkerColor(kBlack);

      a->SetLineColor(kBlue+2);
      b->SetLineColor(kRed+2);
      c1->SetLineColor(kBlue+3);
      c2->SetLineColor(kBlue+1);
      d->SetLineColor(kRed+2);
      e->SetLineColor(kGreen+3);
      f->SetLineColor(kBlack);

      a->SetFillColor(kBlue+2);
      b->SetFillColor(kRed+2);
      c1->SetFillColor(kBlue+3);
      c2->SetFillColor(kBlue+1);
      d->SetFillColor(kRed+2);
      e->SetFillColor(kGreen+3);
      f->SetFillColor(kBlack);

      a->SetTitle("single arm > 3 GeV");
      b->SetTitle("double arm > 1.5 GeV");
      c1->SetTitle("SD1");
      c2->SetTitle("SD2");
      d->SetTitle("DD");
      e->SetTitle("CD");
      f->SetTitle("ND");
      g->SetTitle("all");

      g->GetYaxis()->SetRangeUser(0.0001,0.03);

      TCanvas* can1 = new TCanvas;
      THStack* hs = new THStack("hs","EPOS");
      hs->Add(c1);
      hs->Add(c2);
      hs->Add(d);
      hs->Add(e);
      hs->Add(f);
      hs->Draw("HIST");
      can1->SetLogy();
      hs->GetYaxis()->SetRangeUser(0.0001,0.03);
      hs->SetTitle(";lg(#xi);events (normalised)");
      TLegend* leg = can1->BuildLegend(0.25,0.65,0.45,0.85);
      leg->SetFillColor(kWhite);
      leg->Draw();
      can1->SaveAs((string("plots/diff_1_")+list[i]+string(".eps")).c_str());
      can1->SaveAs((string("plots/diff_1_")+list[i]+string(".pdf")).c_str());
      can1->SaveAs((string("plots/diff_1_")+list[i]+string(".png")).c_str());
  

      TCanvas* can2 = new TCanvas;
      THStack* hs2 = new THStack("hs","EPOS");
      hs2->Add(c1);
      hs2->Add(c2);
      hs2->Add(d);
      hs2->Add(e);
      hs2->Add(f);
      hs2->Draw("HIST");
      hs2->SetTitle(";lg(#xi);events (normalised)");
      leg->Draw();
      can2->SaveAs((string("plots/diff_2_")+list[i]+string(".eps")).c_str());
      can2->SaveAs((string("plots/diff_2_")+list[i]+string(".pdf")).c_str());
      can2->SaveAs((string("plots/diff_2_")+list[i]+string(".png")).c_str());


      TCanvas* c3 = new TCanvas;
      g->Draw("HIST LF2");
      a->Draw("HIST P SAME");
      b->Draw("HIST P SAME");
      c3->SetLogy();
      g->GetYaxis()->SetRangeUser(0.0001,0.03);
      g->SetTitle("Inelastic events;lg(#xi);events (normalised)");
      TLegend* leg = c3->BuildLegend(0.25,0.65,0.65,0.85);
      leg->SetFillColor(kWhite);
      leg->Draw();
      c3->SaveAs((string("plots/diff_3_")+list[i]+string(".eps")).c_str());
      c3->SaveAs((string("plots/diff_3_")+list[i]+string(".pdf")).c_str());
      c3->SaveAs((string("plots/diff_3_")+list[i]+string(".png")).c_str());


      TCanvas* c4 = new TCanvas;
      TH1D* single2 = a->Clone("single2");
      TH1D* double2 = b->Clone("double2");

      single2->SetTitle("single-arm > 3 GeV;lg(#xi);efficiency #epsilon");
      double2->SetTitle("double-arm > 3 GeV;lg(#xi);efficiency #epsilon");

      for(int n=0; n<=single2->GetNbinsX(); n++)
        {
          double ys = single2->GetBinContent(n);
          double yd = double2->GetBinContent(n);
          double y = g->GetBinContent(n);
          if(y)
            {
              single2->SetBinContent(n,ys/y);
              double2->SetBinContent(n,yd/y);
            }
          else
            {
              single2->SetBinContent(n,0);
              double2->SetBinContent(n,0);
            }
        }
      single2->Draw();
      double2->Draw("SAME");
      TLegend* leg = c4->BuildLegend(0.22,0.65,0.5,0.85);
      leg->SetFillColor(kWhite);
      leg->Draw();
      c4->SaveAs((string("plots/diff_4_")+list[i]+string(".eps")).c_str());
      c4->SaveAs((string("plots/diff_4_")+list[i]+string(".pdf")).c_str());
      c4->SaveAs((string("plots/diff_4_")+list[i]+string(".png")).c_str());
    }
}
