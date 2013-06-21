void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);
void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);

void makePlots3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  // {
  // TFile* file2 = TFile::Open("histos_old.root");
  // TFile* file = TFile::Open("histos_old.root");
  // TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfp_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hfp_hits_noise");
  // TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfp_hits_coll");
  // TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfp_hits_coll");
  // TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfp_hits_coll");
  // TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfp_hits_coll");
  // TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfp_hits_coll");
  // ShowStack(a,a2,b,c,d,e,f,"p");
  // }
  // {
  // TFile* file2 = TFile::Open("histos_old.root");
  // TFile* file = TFile::Open("histos_old.root");
  // TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfm_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hfm_hits_noise");
  // TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfm_hits_coll");
  // TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfm_hits_coll");
  // TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfm_hits_coll");
  // TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfm_hits_coll");
  // TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfm_hits_coll");
  // ShowStack(a,a2,b,c,d,e,f,"m");
  // }
  {
    TFile* file = TFile::Open("histos_skiprings.root");
    TFile* file2 = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file->Get("data210885/data210885_h_hf_hits_coll_single");
  TH1D* a2=(TH1D*)file->Get("data210885/data210885_h_hf_hits_noise_single");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_hits_coll_single");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_single");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_hits_coll_single");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_single");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_single");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_single");
  ShowStack(a,a2,b,c,d,e,f,"single");
  }
  {
  TFile* file = TFile::Open("histos_skiprings.root");
  TFile* file2 = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file->Get("data210885/data210885_h_hf_hits_coll_double");
  TH1D* a2=(TH1D*)file->Get("data210885/data210885_h_hf_hits_noise_double");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_hits_coll_double");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_double");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_hits_coll_double");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_double");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_double");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_double");
  ShowStack(a,a2,b,c,d,e,f,"double");
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

  const int normbin = data->FindBin(20);
  const int normendbin = data->GetNbinsX();
  double eposscale=double(c->Integral())/double(sl1->Integral());
  double eposscale2=double(c->Integral())/double(sl2->Integral());
  data->Scale(1./data->GetBinWidth(1));
  noise->Scale(1./double(data->Integral()) / noise->GetBinWidth(1) * 48.2 / 9.98); //HLT PAAccept rate is 48.2HZ compared to 9.98Hz of noise
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
  data->GetYaxis()->SetRangeUser(1e-6,5e1);
  data->GetXaxis()->SetTitle("E_{HF} [GeV]");
  data->GetYaxis()->SetTitle("events (normalised)");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  // THStack* h_s_b = new THStack("h_s","events");
  // h_s_b->Add(noise,"HIST F");
  // h_s_b->Add(sl,"HIST F");
  // h_s_b->Add(b,"p");
  TList *list = new TList;
  list->Add(noise);
  list->Add(sl);
  list->Add(b);
  TH1D* h_s_b = (TH1D*)b->Clone("h_s_b");
  h_s_b->Reset();
  h_s_b->Merge(list);
  THStack* h_s_c = new THStack("h_s_c","events");
  h_s_c->Add(noise,"HIST F");
  h_s_c->Add(sl,"HIST F");
  h_s_c->Add(c,"p");
  THStack* h_s_d = new THStack("h_s_d","events");
  h_s_d->Add(noise,"HIST F");
  h_s_d->Add(sl,"HIST F");
  h_s_d->Add(d,"p");

  TCanvas* c1 = new TCanvas;
  data->Draw("P");
  h_s_b->Draw("SAME");
  h_s_c->Draw("SAME");
  h_s_d->Draw("SAME");
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
  c1->SetLogy();
  c1->SetLogx();
  DataText(0,1);

  TLine* line = new TLine(type=="single"?6:2,1e-6,type=="single"?6:2,0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("SAME");

  c1->SaveAs((string("plots/hf_") + type + string("_signal")+string(".pdf")).c_str());

}

