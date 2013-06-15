void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void makePlots7()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> list;
  list.push_back(string("Epos"));
  list.push_back(string("Hijing"));
  list.push_back(string("QGSJetII"));


  for(int i=0; i<int(list.size()); i++)
    {
      cout << i << " " << int(list.size()) << endl;
      TFile* file = TFile::Open("histos_old.root");
      file->cd();

      cout << string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_single")) << endl;
      TH1D* a=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_single")).c_str());
      TH1D* b=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_double")).c_str());
      TH1D* c1=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_SD1")).c_str());
      TH1D* c2=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_SD2")).c_str());
      TH1D* d=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_DD")).c_str());
      TH1D* e=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_CD")).c_str());
      TH1D* f=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_ND")).c_str());
      TH1D* g=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diffraction_all")).c_str());


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

      b->SetMarkerStyle(24);

      a->SetMarkerColor(kBlue-2);
      b->SetMarkerColor(kRed-2);
      c1->SetMarkerColor(kBlue-3);
      c2->SetMarkerColor(kBlue-1);
      d->SetMarkerColor(kRed-2);
      e->SetMarkerColor(kGreen-3);
      f->SetMarkerColor(kGray);

      a->SetLineColor(a->GetMarkerColor());
      b->SetLineColor(b->GetMarkerColor());
      c1->SetLineColor(c1->GetMarkerColor());
      c2->SetLineColor(c2->GetMarkerColor());
      d->SetLineColor(d->GetMarkerColor());
      e->SetLineColor(e->GetMarkerColor());
      f->SetLineColor(f->GetMarkerColor());

      a->SetFillColor(a->GetMarkerColor());
      b->SetFillColor(b->GetMarkerColor());
      c1->SetFillColor(c1->GetMarkerColor());
      c2->SetFillColor(c2->GetMarkerColor());
      d->SetFillColor(d->GetMarkerColor());
      e->SetFillColor(e->GetMarkerColor());
      f->SetFillColor(f->GetMarkerColor());

      a->SetTitle("single-arm E>6 GeV");
      b->SetTitle("double-arm E>2 GeV");
      c1->SetTitle("SD1");
      c2->SetTitle("SD2");
      d->SetTitle("DD");
      e->SetTitle("CD");
      f->SetTitle("ND");
      g->SetTitle("all");

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
      hs->SetTitle(";log_{10}(#xi);events (normalised)");
      TLegend* leg = can1->BuildLegend(0.25,0.65,0.45,0.85);
      SetLegAtt(leg);
      leg->Draw();
      MCText(0);
      can1->SaveAs((string("plots/diff_1_")+list[i]+string(".pdf")).c_str());


      TCanvas* can2 = new TCanvas;
      THStack* hs2 = new THStack("hs","EPOS");
      hs2->Add(c1);
      hs2->Add(c2);
      hs2->Add(d);
      hs2->Add(e);
      hs2->Add(f);
      hs2->Draw("HIST");
      hs2->SetTitle(";log_{10}(#xi);events (normalised)");
      leg->Draw();
      MCText(0);
      can2->SaveAs((string("plots/diff_2_")+list[i]+string(".pdf")).c_str());


      TCanvas* c3 = new TCanvas;
      g->Draw("HIST LF2");
      a->Draw("HIST P SAME");
      b->Draw("HIST P SAME");
      c3->SetLogy();
      //g->GetYaxis()->SetRangeUser(0.0001,0.03);
      g->SetTitle("Inelastic events;log_{10}(#xi);events (normalised)");
      TLegend* leg = new TLegend(0.22,0.70,0.62,0.87);
      SetLegAtt(leg);
      leg->AddEntry(g,"","l");
      leg->AddEntry(a,"","p");
      leg->AddEntry(b,"","p");
      leg->Draw();
      MCText(0,0);
      c3->SaveAs((string("plots/diff_3_")+list[i]+string(".pdf")).c_str());


      TCanvas* c4 = new TCanvas;
      TH1D* single2 = (TH1D*)a->Clone("single2");
      TH1D* double2 = (TH1D*)b->Clone("double2");

      single2->SetTitle("single-arm E>6 GeV;log_{10}(#xi);efficiency #epsilon");
      double2->SetTitle("double-arm E>2 GeV;log_{10}(#xi);efficiency #epsilon");

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
      single2->GetYaxis()->SetRangeUser(0.,1.3);
      single2->Draw();
      double2->Draw("SAME");
      TLegend* leg = new TLegend(0.22,0.77,0.62,0.87);
      SetLegAtt(leg);
      leg->AddEntry(a,"","p");
      leg->AddEntry(b,"","p");
      leg->Draw();
      MCText(0,0);
      c4->SaveAs((string("plots/diff_4_")+list[i]+string(".pdf")).c_str());
    }
}
