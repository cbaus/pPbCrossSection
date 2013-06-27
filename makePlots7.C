#include <TColor.h>

void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,TH1D* g);

void makePlots7()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  vector<string> list;
  list.push_back(string("Epos"));
  //list.push_back(string("Hijing"));
  //list.push_back(string("QGSJetII"));
  vector<string> name;
  name.push_back(string("EPOS-LHC"));
  //name.push_back(string("HIJING 1.383"));
  //name.push_back(string("QGSJetII-04"));


  for(int i=0; i<int(list.size()); i++)
    {
      cout << i << " " << int(list.size()) << endl;
      TFile* file = TFile::Open("histos.root");
      file->cd();

      cout << string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_single")) << endl;
      TH1D* a=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_single")).c_str());
      TH1D* b=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_double")).c_str());
      TH1D* c1=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_SD1")).c_str());
      TH1D* c2=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_SD2")).c_str());
      TH1D* d=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_DD")).c_str());
      TH1D* e=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_CD")).c_str());
      TH1D* f=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_ND")).c_str());
      TH1D* g=(TH1D*)file->Get(string(list[i]+string("/")+list[i]+string("_h_mc_diff_energy_all")).c_str());


      a->Scale(1./double(g->GetEntries()));
      b->Scale(1./double(g->GetEntries()));
      c1->Scale(1./double(g->GetEntries()));
      c2->Scale(1./double(g->GetEntries()));
      d->Scale(1./double(g->GetEntries()));
      e->Scale(1./double(g->GetEntries()));
      f->Scale(1./double(g->GetEntries()));
      g->Scale(1./double(g->GetEntries()));

      //a->SetBinLabel(1,"E=0 in acceptance");
      //b->SetBinLabel(1,"E=0 in acceptance");
      //c1->SetBinLabel(1,"E=0 in acceptance");
      //c2->SetBinLabel(1,"E=0 in acceptance");
      //d->SetBinLabel(1,"E=0 in acceptance");
      //e->SetBinLabel(1,"E=0 in acceptance");
      //f->SetBinLabel(1,"E=0 in acceptance");
      //g->SetBinLabel(1,"E=0 in acceptance");


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
      e->SetMarkerColor(TColor::GetColor(120,255,120));
      d->SetMarkerColor(TColor::GetColor(255,90,90));
      c2->SetMarkerColor(TColor::GetColor(65,65,255));
      c1->SetMarkerColor(TColor::GetColor(10,10,255));
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
      hs->GetXaxis()->SetBinLabel(1,"E=0 in acceptance");
      hs->GetXaxis()->SetBinLabel(hs->GetXaxis()->FindBin(0),"0");
      hs->GetXaxis()->SetBinLabel(hs->GetXaxis()->FindBin(1),"1");
      hs->GetXaxis()->SetBinLabel(hs->GetXaxis()->FindBin(2),"2");
      hs->GetXaxis()->SetBinLabel(hs->GetXaxis()->FindBin(3),"3");
      can1->SetLogy();
      hs->GetYaxis()->SetRangeUser(0.000001,0.03);
      hs->SetTitle(";log_{10}(#xi);events (normalised)");
      TLegend* leg = can1->BuildLegend(0.25,0.63,0.45,0.93);
      SetLegAtt(leg);
      leg->Draw();
      CMSText(0,0,1,);
      can1->SaveAs((string("plots/diff_1_")+list[i]+string(".pdf")).c_str());


      TCanvas* can2 = new TCanvas;
      THStack* hs2 = new THStack("hs","EPOS");
      hs2->Add(c1);
      hs2->Add(c2);
      hs2->Add(d);
      hs2->Add(e);
      hs2->Add(f);
      hs2->Draw("HIST");
      hs2->GetYaxis()->SetRangeUser(0,5e-3);
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
      TLegend* leg = new TLegend(0.22,0.75,0.62,0.93);
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

      single2->SetTitle("single-arm E>8 GeV;log_{10}(#xi);efficiency #epsilon");
      double2->SetTitle("double-arm E>2.5 GeV;log_{10}(#xi);efficiency #epsilon");

      for(int n=0; n<=single2->GetNbinsX(); n++)
        {
          double ys = single2->GetBinContent(n);
          double yd = double2->GetBinContent(n);
          double y = g->GetBinContent(n);
          double ys_e = a->GetBinError(n);
          double yd_e = b->GetBinError(n);
          double y_e = g->GetBinError(n);
          if(y)
            {
              double eff_single = ys/y;
              double eff_double = yd/y;
              double error_single=sqrt(pow(ys_e/y,2) + pow(y_e*ys/y/y,2));
              double error_double=sqrt(pow(yd_e/y,2) + pow(y_e*yd/y/y,2));
              single2->SetBinContent(n,ys/y);
              double2->SetBinContent(n,yd/y);
              if(eff_single < 0.99)
                single2->SetBinError(n,error_single);
              if(eff_single < 0.99)
                double2->SetBinError(n,error_double);
            }
          else
            {
              single2->SetBinContent(n,0);
              double2->SetBinContent(n,0);
            }
        }

      double thres_single=-10;
      double thres_double=-10;
      for (double x=0; x<single2->GetBinLowEdge(single2->GetNbinsX()); x+=0.001)
        {
          if(single2->Interpolate(x) > 0.95 && thres_single == -10)
            thres_single = x;
          if(double2->Interpolate(x) > 0.95 && thres_double == -10)
            thres_double = x;
        }

      single2->GetYaxis()->SetRangeUser(0.,1.3);
      single2->Draw();
      double2->Draw("SAME");
      TLegend* leg = new TLegend(0.22,0.78,0.62,0.93);
      SetLegAtt(leg);
      leg->AddEntry(a,"","p");
      leg->AddEntry(b,"","p");
      leg->Draw();
      CMSText(0,0,0,name[i]);
      c4->SaveAs((string("plots/diff_4_")+list[i]+string(".pdf")).c_str());
    }
}
