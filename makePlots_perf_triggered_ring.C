///Plot the stacked hit distribution of (noise,em + data/mc) for HF


void Show(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);
void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,string type);
TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal);

void makePlots_perf_triggered_ring()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();

  {
    TFile* file2 = TFile::Open("histos.root");
    TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hf_triggered_ring_single");
  TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_triggered_ring_noise_single");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_triggered_ring_single");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_triggered_ring_single");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_triggered_ring_single");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_triggered_ring_single");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_triggered_ring_single");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_triggered_ring_single");
  ShowStack(a,a2,b,c,d,e,f,"single");
  }
  {
  TFile* file2 = TFile::Open("histos.root");
  TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file2->Get("data210885/data210885_h_hf_triggered_ring_double");
  TH1D* a2=(TH1D*)file2->Get("data210885/data210885_h_hf_triggered_ring_noise_double");
  TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hf_triggered_ring_double");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_triggered_ring_double");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_triggered_ring_double");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_triggered_ring_double");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_triggered_ring_double");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_triggered_ring_double");
  ShowStack(a,a2,b,c,d,e,f,"double");
  }
}

void ShowStack(TH1D* data,TH1D* noise,TH1D* b,TH1D* c,TH1D* d,TH1D* sl1,TH1D* sl2, string type)
  {
    cout << data << " " << noise << " " << b << " " << c << " " << d << " " << sl1 << " " << sl2 << endl;
  double eposscale=double(c->Integral())/double(sl1->Integral());
  double eposscale2=double(c->Integral())/double(sl2->Integral());
  noise->Scale(1./double(data->Integral()) * 48.216 / (9.9845*(1.-(84.+296.)/3568.))); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet
  data->Scale(1./double(data->Integral()));
  b->Scale(data->Integral()/b->Integral());
  double eposnorm=data->Integral()/c->Integral();
  c->Scale(eposnorm);
  d->Scale(data->Integral()/d->Integral());
  sl1->Scale(eposscale*195./2130.*eposnorm);
  sl2->Scale(eposscale2*122./2130.*eposnorm);
  sl1->SetBit(TH1::kIsAverage);
  sl2->SetBit(TH1::kIsAverage);
  if(sl1->Add(sl2) == kFALSE) {new TCanvas; sl1->Draw(); new TCanvas; sl2->Draw(); cerr << "add failed" << endl; return;}
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

  //data->GetXaxis()->SetLimits(1,);
  data->GetYaxis()->SetRangeUser(0,0.2);
  data->GetXaxis()->SetTitle("HF Ring (ieta)");
  data->GetYaxis()->SetTitle("events (normalised)");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  THStack* h_s_bg = new THStack("h_s_gb","events");
  h_s_bg->Add(noise,"HIST F");
  h_s_bg->Add(sl,"HIST F");

  if(0)
    {
      b = merge(noise,sl,b);
      c = merge(noise,sl,c);
      d = merge(noise,sl,d);
    }

  TCanvas* c1 = new TCanvas;
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
  CMSText(1,0,1,type=="single"?"single-arm selection":"double-arm selection");

  data->GetXaxis()->SetRangeUser(-41,-29);
  c1->SaveAs((string("plots/hf_") + type + string("_triggered_rings_m")+string(".pdf")).c_str());
  data->GetXaxis()->SetRangeUser(29,41);
  c1->SaveAs((string("plots/hf_") + type + string("_triggered_rings_p")+string(".pdf")).c_str());

}

TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal)
{
  TList *list = new TList;
  list->Add(bg1);
  list->Add(bg2);
  list->Add(signal);
  ostringstream newname; newname << signal->GetName() << "_merged";
  TH1D* newsignal = (TH1D*)signal->Clone(newname.str().c_str());
  newsignal->Reset();
  newsignal->Merge(list);
  return newsignal;
}
