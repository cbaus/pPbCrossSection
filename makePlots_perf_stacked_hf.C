///Plot the stacked hit distribution of (noise,em + data/mc) for HF

void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* mc4, TH1D* e,TH1D* f,string type);
TH1D* merge(TH1D* bg1, TH1D* bg2, TH1D* signal);

void makePlots_perf_stacked_hf()
{
  gROOT->ProcessLine(" .L style.cc");
  style();

  // {
  // TFile* file2 = TFile::Open("histos_old.root");
  // TFile* file = TFile::Open("histos_old.root");
  // TH1D* a=(TH1D*)file2->Get("data247324/data247324_h_hfp_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data247324/data247324_h_hfp_hits_noise");
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
  // TH1D* a=(TH1D*)file2->Get("data247324/data247324_h_hfm_hits_coll");
  // TH1D* a2=(TH1D*)file2->Get("data247324/data247324_h_hfm_hits_noise");
  // TH1D* b=(TH1D*)file->Get("Hijing/Hijing_h_hfm_hits_coll");
  // TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hfm_hits_coll");
  // TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hfm_hits_coll");
  // TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hfm_hits_coll");
  // TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hfm_hits_coll");
  // ShowStack(a,a2,b,c,d,e,f,"m");
  // }
  {
    TFile* file2 = TFile::Open("histos.root");
    TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file2->Get("data247324/data247324_h_hf_hits_coll_single");
  TH1D* a2=(TH1D*)file2->Get("data247324/data247324_h_hf_hits_noise_single");
  TH1D* b=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_hf_hits_coll_single");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_single");
  TH1D* dpm=(TH1D*)file->Get("DPMJet/DPMJet_h_hf_hits_coll_single");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_hits_coll_single");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_single");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_single");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_single");
  ShowStack(a,a2,b,0,0,0,0,0,"single");
  }
  {
  TFile* file2 = TFile::Open("histos.root");
  TFile* file = TFile::Open("histos.root");
  TH1D* a=(TH1D*)file2->Get("data247324/data247324_h_hf_hits_coll_double");
  TH1D* a2=(TH1D*)file2->Get("data247324/data247324_h_hf_hits_noise_double");
  TH1D* b=(TH1D*)file->Get("PythiaZ2Star/PythiaZ2Star_h_hf_hits_coll_double");
  TH1D* c=(TH1D*)file->Get("Epos/Epos_h_hf_hits_coll_double");
  //TH1D* c=(TH1D*)file->Get("Epos_SL/Epos_SL_h_hf_hits_coll_double");
  TH1D* d=(TH1D*)file->Get("QGSJetII/QGSJetII_h_hf_hits_coll_double");
  TH1D* dpm=(TH1D*)file->Get("DPMJet/DPMJet_h_hf_hits_coll_double");
  TH1D* e=(TH1D*)file->Get("Starlight_DPMJet/Starlight_DPMJet_h_hf_hits_coll_double");
  TH1D* f=(TH1D*)file->Get("Starlight_Pythia/Starlight_Pythia_h_hf_hits_coll_double");
  ShowStack(a,a2,b,0,0,0,0,0,"double");
  }
}

void ShowStack(TH1D* data,TH1D* noise,TH1D* mc1,TH1D* mc2,TH1D* mc3,TH1D* mc4, TH1D* sl1,TH1D* sl2, string type)
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
  double mc1scale;
  double mc1scale2;
  if(sl1)
    {
      mc1scale=double(mc1->Integral())/double(sl1->Integral());
      mc1scale2=double(mc1->Integral())/double(sl2->Integral());
    }
  noise->Scale(1./double(data->Integral()) * 260.96 / (9.59*(1.-18./3568.)) / 240 * 58.97); //HLT PAAccept rate is 48.2216HZ compared to 9.9845Hz of noise //(84+296)/3568 get skipped because pf BPTX_quiet

  data->Scale(1./double(data->Integral()));
  double mc1norm=data->Integral(normbin,normendbin)/mc1->Integral(normbin,normendbin);
  mc1->Scale(mc1norm);
  if(mc2) mc2->Scale(data->Integral(normbin,normendbin)/mc2->Integral(normbin,normendbin));
  if(mc2) mc3->Scale(data->Integral(normbin,normendbin)/mc3->Integral(normbin,normendbin));
  if(mc2) mc4->Scale(data->Integral(normbin,normendbin)/mc4->Integral(normbin,normendbin));
  if(mc2) sl1->Scale(mc1scale*195./2085.*mc1norm);
  if(mc2) sl2->Scale(mc1scale2*122./2085.*mc1norm);
  if(sl1) sl1->SetBit(TH1::kIsAverage);
  if(sl2) sl2->SetBit(TH1::kIsAverage);
  if(mc2) if(sl1->Add(sl2) == kFALSE) {new TCanvas; sl1->Draw(); new TCanvas; sl2->Draw(); cerr << "add failed" << endl; return;}
  TH1D* sl = 0;
  if(sl1) sl = sl1;

  //data->Add(noise,-1);
  data->SetMarkerSize(1);

  // data->SetLineWidth(2);
  // noise->SetLineWidth(1.5);
  // mc2->SetLineWidth(2);
  // mc1->SetLineWidth(2);
  // mc3->SetLineWidth(2);
  // sl->SetLineWidth(2);


  data->SetMarkerColor(kBlack);
  noise->SetMarkerColor(kGreen-9);
  mc1->SetMarkerColor(kRed);
  if(mc2) mc2->SetMarkerColor(kBlue);
  if(mc2) mc3->SetMarkerColor(kGreen+2);
  if(mc2) mc4->SetMarkerColor(kMagenta);
  if(mc2) sl->SetMarkerColor(kOrange+1);

  noise->SetMarkerStyle(34);
  mc1->SetMarkerStyle(4);
  if(mc2) mc2->SetMarkerStyle(25);
  if(mc2) mc3->SetMarkerStyle(28);
  if(mc2) mc4->SetMarkerStyle(26);
  if(mc2) sl->SetMarkerStyle(22);

  data->SetLineColor(data->GetMarkerColor());
  noise->SetLineColor(noise->GetMarkerColor());
  mc1->SetLineColor(mc1->GetMarkerColor());
  if(mc2) mc2->SetLineColor(mc2->GetMarkerColor());
  if(mc2) mc3->SetLineColor(mc3->GetMarkerColor());
  if(mc2) mc4->SetLineColor(mc4->GetMarkerColor());
  if(mc2) sl->SetLineColor(sl->GetMarkerColor());

  data->SetFillColor(data->GetMarkerColor());
  noise->SetFillColor(noise->GetMarkerColor());
  mc1->SetFillColor(mc1->GetMarkerColor());
  if(mc2) mc2->SetFillColor(mc2->GetMarkerColor());
  if(mc2) mc3->SetFillColor(mc3->GetMarkerColor());
  if(mc2) mc4->SetFillColor(mc4->GetMarkerColor());
  if(mc2) sl->SetFillColor(sl->GetMarkerColor());

  data->SetFillStyle(1001);
  noise->SetFillStyle(1001);
  mc1->SetFillStyle(1001);
  if(mc2) mc2->SetFillStyle(1001);
  if(mc2) mc3->SetFillStyle(1001);
  if(mc2) mc4->SetFillStyle(1001);
  if(mc2) sl->SetFillStyle(1001);

  data->SetTitle("Data");
  noise->SetTitle("Noise");
  mc1->SetTitle("Pythia6 Z2*");
  if(mc2) mc2->SetTitle("EPOS-LHC");
  if(mc2) mc3->SetTitle("QGSJETII-04");
  if(mc2) mc4->SetTitle("DPMJET3.06");
  if(mc2) sl->SetTitle("#gammap (STARLIGHT+DPMJET/PYTHIA)");

  data->GetXaxis()->SetLimits(1,200);
  data->GetYaxis()->SetRangeUser(1e-6,2e2);
  data->GetXaxis()->SetTitle("#it{E}_{HF} [GeV]");
  data->GetYaxis()->SetTitle("Events (normalised)");
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset()*1.1);


  THStack* h_s_bg = new THStack("h_s_gb","events");
  h_s_bg->Add(noise,"HIST F");
  if(sl1) h_s_bg->Add(sl,"HIST F");

  if(0)
    {
      if(sl1) mc1 = merge(noise,sl,mc1);
      if(mc2) mc2 = merge(noise,sl,mc2);
      if(mc2) mc3 = merge(noise,sl,mc3);
    }

  TCanvas* c1 = new TCanvas;
  data->Draw("P");
  h_s_bg->Draw("SAME HIST F");
  mc1->Draw("SAME");
  if(mc2) mc2->Draw("SAME");
  if(mc2) mc3->Draw("SAME");
  if(mc2) mc4->Draw("SAME");
  data->Draw("SAME P");
  data->Draw("SAME AXIS");

  TLegend* leg = new TLegend(0.23,0.70,0.43,0.93);
  SetLegAtt(leg);
  leg->AddEntry(data,"","p");
  leg->AddEntry(mc1,"","p");
  if(mc2) leg->AddEntry(mc2,"","p");
  if(mc2) leg->AddEntry(mc3,"","p");
  if(mc2) leg->AddEntry(mc4,"","p");
  if(mc2) leg->AddEntry(sl,"","f");
  leg->AddEntry(noise,"","f");
  leg->Draw();
  c1->SetLogy();
  c1->SetLogx();
  CMSText(3,0,1,type=="single"?"Single-arm selection":"Double-arm selection");

  TLine* line = new TLine(type=="single"?8:4,1e-6,type=="single"?8:4,0.1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("SAME");

  c1->SaveAs((string("plots/hf_") + type + string("_signal_paper")+string(".pdf")).c_str());

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
