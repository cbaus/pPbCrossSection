///Plot the stacked hit distribution of (noise,em + data/mc) for HF


int RingToIeta(int ring);
int RingToCanvas(int ring);
void ShowStack(TH1D* a,TH1D* a2,TH1D* b,TH1D* c,TH1D* d,TH1D* e,TH1D* f,int ring);
void DividegPad(Int_t nx, Int_t ny, Float_t l, Float_t r, Float_t t, Float_t b);

TCanvas* c1 = NULL;
void makePlots_perf_hits_rings()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();


  gStyle->SetPadTopMargin(0.);
  gStyle->SetPadBottomMargin(0.);
  gStyle->SetPadLeftMargin(0.);
  gStyle->SetPadRightMargin(0.); 
  gStyle->SetFrameBorderMode(0);
  for(int i=0; i<24; i++)
    {
      if(i==0 || i==24)
        {
          ostringstream title; title << "c" << i;
          c1 = new TCanvas(title.str().c_str(),title.str().c_str(),1200,1600);
          Double_t mlb = 0.3;
          Double_t mrt = 0.1;
          Double_t nx  = 3;
          Double_t ny  = 4;
          DividegPad(nx,ny,mlb,mrt,mrt,mlb);
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
      if(i==11) c1->SaveAs((string("plots/hf_ring_signal_m.pdf")).c_str());
      if(i==23) c1->SaveAs((string("plots/hf_ring_signal_p.pdf")).c_str());
    }

}

void ShowStack(TH1D* data,TH1D* noise,TH1D* b,TH1D* c,TH1D* d,TH1D* sl1,TH1D* sl2,int ring)
{

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

  TPad* refPad = c1->cd(3);
  TPad* curPad = c1->cd(RingToCanvas(ring));
  double heightRatio = curPad->GetHNDC()/refPad->GetHNDC();
  double widthRatio = curPad->GetWNDC()/refPad->GetWNDC();
  double smallestRatio = TMath::Min(heightRatio,widthRatio);
  data->GetXaxis()->SetTitleSize(data->GetXaxis()->GetTitleSize() / widthRatio);
  data->GetYaxis()->SetTitleSize(data->GetYaxis()->GetTitleSize() / heightRatio);
  data->GetXaxis()->SetLabelSize(data->GetXaxis()->GetLabelSize() / widthRatio);
  data->GetYaxis()->SetLabelSize(data->GetYaxis()->GetLabelSize() / heightRatio);
  data->GetYaxis()->SetTitleOffset(data->GetYaxis()->GetTitleOffset() * heightRatio);
  data->GetXaxis()->SetTitleOffset(data->GetXaxis()->GetTitleOffset() * widthRatio);
  data->GetYaxis()->SetLabelOffset(data->GetYaxis()->GetLabelOffset() * heightRatio);
  data->GetXaxis()->SetLabelOffset(data->GetXaxis()->GetLabelOffset() * widthRatio);
  

  cout << "Ring: " << ring << "  Canvas: " << RingToCanvas(ring) << "  Ieta " << RingToIeta(ring) << endl;
  data->Draw("P");
  h_s_bg->Draw("SAME");
  b->Draw("SAME");
  c->Draw("SAME");
  d->Draw("SAME");
  data->Draw("SAME P");
  data->Draw("SAME AXIS");

  if(ring == 1 || ring == 13)
    {
      TLegend* leg = new TLegend(0.23,0.62,0.63,0.87);
      SetLegAtt(leg);
      leg->AddEntry(data,"","p");
      leg->AddEntry(b,"","p");
      leg->AddEntry(c,"","p");
      leg->AddEntry(d,"","p");
      leg->AddEntry(sl,"","f");
      leg->AddEntry(noise,"","f");
      leg->SetTextFont(62);
      leg->Draw();
    }
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
  if(ring >= 12) ring -=12;

  if(ring == 0) return 10;
  else if(ring == 1) return 11;
  else if(ring == 2) return 12;
  else if(ring == 3) return 7;
  else if(ring == 4) return 8;
  else if(ring == 5) return 9;
  else if(ring == 6) return 4;
  else if(ring == 7) return 5;
  else if(ring == 8) return 6;
  else if(ring == 9) return 1;
  else if(ring ==10) return 2;
  else if(ring ==11) return 3;
  else cerr << "asd" << endl;
}

//   if(ring < 12)
//     return ring + 1;
//   else
//     return ring - 11;
// }

void DividegPad(Int_t nx, Int_t ny, Float_t l, Float_t r, Float_t t, Float_t b)
{
   Int_t ix, iy, n=0;
   Double_t x1, x2, y1, y2;
   Double_t dx = ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l);
   Double_t dl = dx/(1-l);
   Double_t dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
   Double_t db = dy/(1-b);
   char *name  = new char [strlen(gPad->GetName())+6];

   y1 = 0;
   y2 = db;
   for (iy=0; iy<ny; iy++) {
      x1 = 0;
      x2 = dl;
      for (ix=0;ix<nx;ix++) {
         if (x1 > x2) continue;
         n++;
         sprintf(name,"%s_%d",gPad->GetName(),n);
         pad = new TPad(name,name,x1,y1,x2,y2,0);
         if (ix==0)    pad->SetLeftMargin(l);
         if (ix==nx-1) pad->SetRightMargin(r);
         if (iy==ny-1) pad->SetTopMargin(t);
         if (iy==0)    pad->SetBottomMargin(b);
         x1 = x2;
         if (ix==nx-2) x2 = 1;
         else          x2 = x1+dx;
         pad->SetNumber(n);
         pad->Draw();
      }
      y1 = y2;
      if (iy==ny-2) y2 = 1;
      else          y2 = y1+dy;
   }
}
