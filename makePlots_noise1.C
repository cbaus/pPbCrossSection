void Show(TH1D* a,TH1D* b,TH1D* c,TH1D* d, TH1D* e, TH1D* f);

void makePlots_noise1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  TFile* file = TFile::Open("histos_noise.root");

  {
  TCanvas* c1 = new TCanvas;
  TCanvas* c0;
  c1->Divide(4,3);

  vector<double> x;
  vector<double> xerr;
  vector<double> y;
  vector<double> yerr;
  vector<double> yeff;
  
  string add = "noise";
  int count=0;
  for (int i=-41; i<=41;)
    {
      cout << " Processing ... " << i << " count: " << count << endl;
      count++;
      ostringstream hf_str, hf_cut_str;
      hf_str  << add << "h_hf_rings_single_" << i;
      hf_cut_str  << add << "h_hf_rings_cut_single_" << i;
      TH1D* h_hf_rings_single = NULL;
      TH1D* h_hf_rings_cut_single = NULL;
      h_hf_rings_single = (TH1D*)(file->Get((add+string("/")+hf_str.str()).c_str()));
      h_hf_rings_cut_single = (TH1D*)(file->Get((add+string("/")+hf_cut_str.str()).c_str()));
      
      if(h_hf_rings_single != NULL && h_hf_rings_cut_single != NULL)
        {

          c1->cd(count)->SetLogy();
          h_hf_rings_single->Scale(1./double(h_hf_rings_single->GetEntries()));
          h_hf_rings_single->SetLineWidth(2);
          h_hf_rings_single->SetLineColor(kBlack-2);
          h_hf_rings_single->SetTitle((hf_str.str()+string(";HF tower E / GeV; events (normalised)")).c_str());
          h_hf_rings_single->GetXaxis()->SetRangeUser(0.,5.);
          h_hf_rings_single->GetYaxis()->SetRangeUser(0.001,1.001);
          h_hf_rings_single->Draw();
          
          x.push_back(i);
          xerr.push_back(0.);
          y.push_back(h_hf_rings_single->GetMean());
          yerr.push_back(h_hf_rings_single->GetRMS());
          yeff.push_back(h_hf_rings_cut_single->Interpolate(3)/h_hf_rings_cut_single->GetBinContent(1));
        }
      else
        {
          cout << "skipping ... " << i << endl;
        }
      

      if(i==-30)
        {
          c0 = c1;
          c1 = new TCanvas();
          c1->Divide(4,3);
          count = 0;
          i=30;
        }
      else
        i++;
    }
  }
  c0->SaveAs((string("plots/hf_m_noise_")+string(".eps")).c_str());
  c0->SaveAs((string("plots/hf_m_noise_")+string(".pdf")).c_str());
  c0->SaveAs((string("plots/hf_m_noise_")+string(".png")).c_str());
  c1->SaveAs((string("plots/hf_p_noise_")+string(".eps")).c_str());
  c1->SaveAs((string("plots/hf_p_noise_")+string(".pdf")).c_str());
  c1->SaveAs((string("plots/hf_p_noise_")+string(".png")).c_str());


  TCanvas* c2 = new TCanvas;
  c2->cd();
  TGraphErrors* graph = new TGraphErrors(x.size(),&x.front(),&y.front(),&xerr.front(),&yerr.front());
  graph->SetName("noisepertower");
  graph->SetTitle(";ieta;<HF tower E>_{Noise} / GeV");
  graph->Draw("AP");
  c2->SaveAs((string("plots/hf_meannoise_")+string(".eps")).c_str());
  c2->SaveAs((string("plots/hf_meannoise_")+string(".png")).c_str());
  c2->SaveAs((string("plots/hf_meannoise_")+string(".pdf")).c_str());


  
  TCanvas* c3 = new TCanvas;
  TGraph* graph2 = new TGraph(x.size(),&x.front(),&yeff.front());
  graph2->SetName("noiseeffpertower");
  graph2->SetTitle(";ieta;event fraction (3 GeV cut)");
  graph2->Draw("AP");
  c3->SaveAs((string("plots/hf_noisecut_")+string(".eps")).c_str());
  c3->SaveAs((string("plots/hf_noisecut_")+string(".pdf")).c_str());
  c3->SaveAs((string("plots/hf_noisecut_")+string(".png")).c_str());


}
