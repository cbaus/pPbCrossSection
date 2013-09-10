///Plot the stacked hit distribution of (noise,em + data/mc) for HF


int RingToIeta(int ring);
int RingToCanvas(int ring);
void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,int ring);

TCanvas* c1 = NULL;
void makePlots_perf_hits_rings()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  for(int i=0; i<24; i++)
    {
      if(i==0 || i==12)
        {
          ostringstream title; title << "c" << i;
          c1 = new TCanvas(title.str().c_str(),title.str().c_str(),1200,1600);
          c1->Divide(3,4);
        }
      ostringstream ss_ring; ss_ring << i;
      TFile* file2 = TFile::Open("histos.root");
      TFile* file = TFile::Open("histos.root");
      TH1D* a=(TH1D*)file2->Get((string("data210885/data210885_h_hf_hits_rings_")+ss_ring.str()).c_str());
      TH1D* a2=(TH1D*)file2->Get((string("data210885/data210885_h_hf_hits_rings_noise_")+ss_ring.str()).c_str());
      TH1D* b=(TH1D*)file->Get((string("Hijing/Hijing_h_hf_hits_rings_")+ss_ring.str()).c_str());
      TH1D* c=(TH1D*)file->Get((string("Epos/Epos_h_hf_hits_rings_")+ss_ring.str()).c_str());
      //TH1D* c=(TH1D*)file->Get((string("Epos_SL/Epos_SL_h_hf_hits_rings_")+ss_ring.str()).c_str());
      TH1D* d=(TH1D*)file->Get((string("QGSJetII/QGSJetII_h_hf_hits_rings_")+ss_ring.str()).c_str());
      TH1D* e=(TH1D*)file->Get((string("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_rings_")+ss_ring.str()).c_str());
      TH1D* f=(TH1D*)file->Get((string("Starlight_Pythia/Starlight_Pythia_h_hf_hits_rings_")+ss_ring.str()).c_str());

      ShowStack(a,a2,b,c,d,e,f,i);
      if(i==0 || i==11)
        {
          c1->SaveAs((string("plots/hf_ring_signal.pdf")).c_str());
        }
    }

}

void ShowStack(TH1D* data,TH1D* noise,TH1D* b,TH1D* c,TH1D* d,TH1D* sl1,TH1D* sl2,int ring)
{
  //RESCALING-------------------------------
  // TH1D** toBeRescaled = &b;
  // double fac = 1./1.3;
  // double* newbinning = new double[(*toBeRescaled)->GetNbinsX()+1];
  // for (int i=0; i<=(*toBeRescaled)->GetNbinsX();i++)
  //   newbinning[i]=(*toBeRescaled)->GetBinLowEdge(i+1)*fac;
  // newbinning[0]=(*toBeRescaled)->GetBinLowEdge(1);
  // newbinning[(*toBeRescaled)->GetNbinsX()]=(*toBeRescaled)->GetBinLowEdge((*toBeRescaled)->GetNbinsX()+1);
  // TH1D* bs = new TH1D("rescaled","rescaled",(*toBeRescaled)->GetNbinsX(),newbinning);
  // for (int i=1; i<=bs->GetNbinsX();i++)
  //   {
  //     bs->SetBinContent(i,(*toBeRescaled)->GetBinContent(i));
  //     bs->SetBinError(i,(*toBeRescaled)->GetBinError(i));
  //   }
  //(*toBeRescaled)=bs;
  //RESCALING-------------------------------

  const int normbin = data->FindBin(20);
  const int normendbin = data->GetNbinsX();
  double eposscale=double(c->Integral())/double(sl1->Integral());
  double eposscale2=double(c->Integral())/double(sl2->Integral());
  noise->Scale(1./double(data->Integral()) * 48.216 / (9.9845*(1.-(84.+296.)/3568.))); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet
  data->Scale(1./double(data->Integral()));
  b->Scale(data->Integral(normbin,normendbin)/b->Integral(normbin,normendbin));
  double eposnorm=data->Integral(normbin,normendbin)/c->Integral(normbin,normendbin);
  c->Scale(eposnorm);
  d->Scale(data->Integral(normbin,normendbin)/d->Integral(normbin,normendbin));
  sl1->Scale(eposscale*195./2130.*eposnorm);
  sl2->Scale(eposscale2*122./2130.*eposnorm);
  sl1->SetBit(TH1::kIsAverage);
  sl2->SetBit(TH1::kIsAverage);
  if(sl1->Add(sl2) == kFALSE) {new TCanvas; sl1->Draw(); new TCanvas; sl2->Draw(); cerr << "add failed" << endl; return;} //if error
  TH1D* sl = sl1;

  //data->Add(noise,-1);
  data->SetMarkerSize(1);

  // data->SetLineWidth(2);
  // noise->SetLineWidth(1.5);
  // b->SetLineWidth(2);
  // c->SetLineWidth(2);
  // d->SetLineWidth(2);
  // sl->SetLineWidth(2);


  data->SetMarkerColor(kBlack);
  noise->SetMarkerColor(kGreen-10);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);
  sl->SetMarkerColor(kOrange-9);

  noise->SetMarkerStyle(34);
  b->SetMarkerStyle(4);
  c->SetMarkerStyle(25);
  d->SetMarkerStyle(28);
  sl->SetMarkerStyle(22);

  data->SetLineColor(data->GetMarkerColor());
  noise->SetLineColor(noise->GetMarkerColor());
  b->SetLineColor(b->GetMarkerColor());
  c->SetLineColor(c->GetMarkerColor());
  d->SetLineColor(d->GetMarkerColor());
  sl->SetLineColor(sl->GetMarkerColor());

  data->SetFillColor(data->GetMarkerColor());
  noise->SetFillColor(noise->GetMarkerColor());
  b->SetFillColor(b->GetMarkerColor());
  c->SetFillColor(c->GetMarkerColor());
  d->SetFillColor(d->GetMarkerColor());
  sl->SetFillColor(sl->GetMarkerColor());


  data->SetTitle("Data");
  noise->SetTitle("Noise");
  b->SetTitle("HIJING 1.383");
  c->SetTitle("EPOS-LHC");
  d->SetTitle("QGSJetII-04");
  sl->SetTitle("#gamma-p (STARLIGHT+DPMJET/PYTHIA)");

  data->GetXaxis()->SetLimits(1,200);
  data->GetYaxis()->SetRangeUser(1e-6,5e1);
  data->GetXaxis()->SetTitle("E_{HF} [GeV]");
  data->GetYaxis()->SetTitle("events (normalised)");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  THStack* h_s_bg = new THStack("h_s_gb","events");
  h_s_bg->Add(noise,"HIST F");
  h_s_bg->Add(sl,"HIST F");
  // TList *list = new TList;
  // list->Add(noise);
  // list->Add(sl);
  // TH1D* h_s_bg = (TH1D*)b->Clone("h_s_bg");
  // h_s_bg->Reset();
  // h_s_bg->Merge(list);

  //TCanvas* c1 = new TCanvas;
  cout << "Ring: " << ring << "  Canvas: " << RingToCanvas(ring) << "  Ieta " << RingToIeta(ring) << endl;
  c1->cd(RingToCanvas(ring));
  data->Draw("P");
  h_s_bg->Draw("SAME");
  b->Draw("SAME");
  c->Draw("SAME");
  d->Draw("SAME");
  data->Draw("SAME P");
  data->Draw("SAME AXIS");

  TLegend* leg = new TLegend(0.23,0.72,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(data,"","p");
  leg->AddEntry(b,"","p");
  leg->AddEntry(c,"","p");
  leg->AddEntry(d,"","p");
  leg->AddEntry(sl,"","f");
  leg->AddEntry(noise,"","f");
  leg->Draw();
  c1->cd(RingToCanvas(ring))->SetLogy();
  c1->cd(RingToCanvas(ring))->SetLogx();
  ostringstream txt; txt << "Ring (ieta=" << RingToIeta(ring) << ")";
  CMSText(1,0,1,txt.str());

}

int RingToIeta(int ring)
///converts internal ring numbering to original ieta
{
  if(ring < 12)
    return ring - 41;
  else
    return ring + 18;
}

int RingToCanvas(int ring)
///converts internal ring numbering to original ieta
{
  if(ring < 12)
    return ring + 1;
  else
    return ring - 11;
}
