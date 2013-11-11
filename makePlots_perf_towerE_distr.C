void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);

void makePlots_crosscheck2()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  {
    TFile* file2 = TFile::Open("histos_old.root");
    TFile* file = TFile::Open("histos_mc.root");
  TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hf_hitdistr_coll");
  TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_hitdistr_noise");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_hitdistr_coll");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hitdistr_coll");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_hitdistr_coll");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hitdistr_coll");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hitdistr_coll");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hitdistr_coll");
  ShowStack(a,a2,b,c,d,e,f,"both");
  }
}

void ShowStack(TH1D* data,TH1D* noise,TH1D* b,TH1D* c,TH1D* d,TH1D* sl1,TH1D* sl2, string type)
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

    for(int k=0; k<=data->GetNbinsX();k++)
      {
        data->SetBinContent(k,data->GetBinContent(k)/data->GetBinWidth(k));
        noise->SetBinContent(k,noise->GetBinContent(k)/noise->GetBinWidth(k));
        b->SetBinContent(k,b->GetBinContent(k)/b->GetBinWidth(k));
        c->SetBinContent(k,c->GetBinContent(k)/c->GetBinWidth(k));
        d->SetBinContent(k,d->GetBinContent(k)/d->GetBinWidth(k));
        sl1->SetBinContent(k,sl1->GetBinContent(k)/sl1->GetBinWidth(k));
        sl2->SetBinContent(k,sl2->GetBinContent(k)/sl2->GetBinWidth(k));
      }

  const int normbin = data->FindBin(3);
  const int normendbin = data->GetNbinsX();
  double eposscale=double(c->Integral())/double(sl1->Integral()); //scale number of events in the two samples
  double eposscale2=double(c->Integral())/double(sl2->Integral());
  noise->Scale(1./double(data->Integral()) * 48.2 / 9.98); //HLT PAAccept rate is 48.2HZ compared to 9.98Hz of noise
  data->Scale(1./double(data->Integral()));
  b->Scale(data->Integral(normbin,normendbin)/b->Integral(normbin,normendbin));
  double eposnorm=data->Integral(normbin,normendbin)/c->Integral(normbin,normendbin);
  c->Scale(eposnorm);
  d->Scale(data->Integral(normbin,normendbin)/d->Integral(normbin,normendbin));
  sl1->Scale(eposscale*195./2130.*eposnorm);
  sl2->Scale(eposscale2*122./2130.*eposnorm);
  sl1->SetBit(TH1::kIsAverage);
  sl2->SetBit(TH1::kIsAverage);
  if(sl1->Add(sl2) == kFALSE) {new TCanvas; sl1->Draw(); new TCanvas; sl2->Draw(); cerr << "add failed" << endl; return;}
  TH1D* sl = sl1;

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
  data->GetYaxis()->SetRangeUser(1e-9,5e1);
  data->GetXaxis()->SetTitle("E_{HF} [GeV]");
  data->GetYaxis()->SetTitle("#frac{1}{N}#frac{dN}{dE_{HF}}");
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

  TCanvas* c1 = new TCanvas;
  data->Draw("P");
  //h_s_bg->Draw("SAME");
  b->Draw("SAME");
  c->Draw("SAME");
  d->Draw("SAME");
  data->Draw("SAME P");
  //data->Draw("SAME AXIS");

  TLegend* leg = new TLegend(0.23,0.72,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(data,"","p");
  leg->AddEntry(b,"","p");
  leg->AddEntry(c,"","p");
  leg->AddEntry(d,"","p");
  //leg->AddEntry(sl,"","f");
  //leg->AddEntry(noise,"","f");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  CMSText(1,0,1,"HF+ and HF-");

  c1->SaveAs((string("plots/hf_hitdistr_") + type + string(".pdf")).c_str());

}

