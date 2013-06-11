void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string prefix);
void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string prefix);

void makePlots3()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  // {
  // TFile* file2 = TFile::Open("histos_old.root");
  // TFile* file = TFile::Open("histos_old.root");
  // TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfp_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_hits_noise");
  // TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfp_hits_coll");
  // TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfp_hits_coll");
  // TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfp_hits_coll");
  // TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfp_hits_coll");
  // TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfp_hits_coll");
  // Show(a,a2,b,c,d,e,f,"p");
  // }
  // {
  // TFile* file2 = TFile::Open("histos_old.root");
  // TFile* file = TFile::Open("histos_old.root");
  // TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hfm_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_hits_noise");
  // TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfm_hits_coll");
  // TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfm_hits_coll");
  // TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfm_hits_coll");
  // TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfm_hits_coll");
  // TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfm_hits_coll");
  // Show(a,a2,b,c,d,e,f,"m");
  // }
  {
    TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file->Get("data210885/data210885_h_hf_hits_coll_single");
  TH1D* a2=(TH1D*)file->Get("data210885/data210885_h_hf_hits_noise_single");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_hits_coll_single");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_single");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_single");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_single");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_single");
  ShowStack(a,a2,b,c,d,e,f,"single");
  }
  {
  TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file->Get("data210885/data210885_h_hf_hits_coll_double");
  TH1D* a2=(TH1D*)file->Get("data210885/data210885_h_hf_hits_noise_double");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_hits_coll_double");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_double");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_double");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_double");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_double");
  ShowStack(a,a2,b,c,d,e,f,"double");
  }
}

void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f, string prefix)
  {
  const int normbin = a->FindBin(20);
  const int normendbin = a->GetNbinsX();
  double eposscale=double(c->Integral())/double(e->Integral());
  double eposscale2=double(c->Integral())/double(f->Integral());
  a->Scale(1./a->GetBinWidth(1)/3329./20.);
  a2->Scale(1./double(a->Integral()) / a2->GetBinWidth(1) / 60 * 48.2 / 9.98); //HLT PAAccept rate is 48.2HZ compared to 9.98Hz of noise
  a->Scale(1./double(a->Integral()));
  b->Scale(a->Integral(normbin,normendbin)/b->Integral(normbin,normendbin));
  double eposnorm=a->Integral(normbin,normendbin)/c->Integral(normbin,normendbin);
  c->Scale(eposnorm);
  d->Scale(a->Integral(normbin,normendbin)/d->Integral(normbin,normendbin));
  e->Scale(eposscale*195./2130.*eposnorm);
  f->Scale(eposscale2*122./2130.*eposnorm);

  a->SetMarkerSize(1.25);
  e->SetMarkerSize(1.25);
  f->SetMarkerSize(1.25);

  // a->SetLineWidth(2);
  // a2->SetLineWidth(1.5);
  // b->SetLineWidth(2);
  // c->SetLineWidth(2);
  // d->SetLineWidth(2);
  // e->SetLineWidth(2);
  // f->SetLineWidth(2);


  a->SetMarkerColor(kBlack);
  a2->SetMarkerColor(kCyan+2);
  b->SetMarkerColor(kRed);
  c->SetMarkerColor(kBlue);
  d->SetMarkerColor(kGreen+2);
  e->SetMarkerColor(kMagenta+3);
  f->SetMarkerColor(kMagenta-5);

  a2->SetMarkerStyle(34);
  b->SetMarkerStyle(4);
  c->SetMarkerStyle(4);
  d->SetMarkerStyle(4);
  e->SetMarkerStyle(22);
  f->SetMarkerStyle(23);

  a->SetLineColor(kBlack);
  a2->SetLineColor(kCyan+2);
  b->SetLineColor(kRed);
  c->SetLineColor(kBlue);
  d->SetLineColor(kGreen+2);
  e->SetLineColor(kMagenta+3);
  f->SetLineColor(kMagenta-5);

  a->SetTitle("zero bias");
  a2->SetTitle("noise");
  b->SetTitle("HIJING");
  c->SetTitle("EPOS");
  d->SetTitle("QGSJetII");
  e->SetTitle("SL+DPMJet");
  f->SetTitle("SL+Pythia");

  a->GetXaxis()->SetLimits(0.5,200);
  a->GetYaxis()->SetRangeUser(1e-7,1.01);
  a->GetXaxis()->SetTitle("hottest HF tower E / GeV");
  a->GetYaxis()->SetTitle("events (normalised)");

//   double fac = 1.25;
//   TH1D* bs = new TH1D ("asd","asd",a->GetNbinsX(),0,200*fac);
//   for (int i=0; i<a->GetNbinsX()+1;i++)
//     {
//       bs->SetBinContent(i,a->GetBinContent(i));
//     }
//   a=bs;

  TCanvas* c1 = new TCanvas;
  a->Draw("P");
  a2->Draw("L SAME");
  b->Draw("L SAME");
  c->Draw("L SAME");
  d->Draw("L SAME");
  e->Draw("L SAME");
  f->Draw("L SAME");
  TLegend* leg = c1->BuildLegend(0.6,0.67,0.78,0.88);
  SetLegAtt(leg);
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  DataText();
  c1->SaveAs((string("plots/hf_") + prefix + string("_signal")+string(".pdf")).c_str());

}

void ShowStack(TH1D* data,TH1D* noise,TH1D* b,TH1D* c,TH1D* d,TH1D* sl1,TH1D* sl2, string prefix)
  {
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
  sl->SetTitle("#gammaP");

  data->GetXaxis()->SetLimits(1,200);
  data->GetYaxis()->SetRangeUser(1e-6,1.01);
  data->GetXaxis()->SetTitle("E_{HF} [GeV]");
  data->GetYaxis()->SetTitle("events (normalised)");

//   double fac = 1.25;
//   TH1D* bs = new TH1D ("asd","asd",data->GetNbinsX(),0,200*fac);
//   for (int i=0; i<data->GetNbinsX()+1;i++)
//     {
//       bs->SetBinContent(i,data->GetBinContent(i));
//     }
//   a=bs;

  THStack* h_s_b = new THStack("h_s","events");
  h_s_b->Add(noise,"HIST F");
  h_s_b->Add(sl,"HIST F");
  h_s_b->Add(b,"p");
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

  TLegend* leg = new TLegend(0.6,0.67,0.78,0.88);
  SetLegAtt(leg);
  leg->AddEntry(data,"","p");
  leg->AddEntry(noise,"","f");
  leg->AddEntry(sl,"","f");
  leg->AddEntry(b,"","p");
  leg->AddEntry(c,"","p");
  leg->AddEntry(d,"","p");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  DataText();
  c1->SaveAs((string("plots/hf_") + prefix + string("_signal")+string(".pdf")).c_str());

}

