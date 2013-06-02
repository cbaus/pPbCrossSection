void makePlots_concl1()
{
  gROOT->ProcessLine(" .L style.cc+");
  style();
  
  TCanvas* can1 = new TCanvas;
  can1->Divide(2);

  TGraphErrors* pp_epos  = new TGraphErrors(1);
  pp_epos->SetName("pp_epos");
  pp_epos->SetTitle("EPOS");
  pp_epos->SetPoint(0,0.9,64.783);
  pp_epos->SetMarkerColor(kBlue);
  pp_epos->SetMarkerStyle(22);
  TGraphErrors* pb_epos  = new TGraphErrors(1);
  pb_epos->SetName("pb_epos");
  pb_epos->SetTitle("EPOS");
  pb_epos->SetPoint(0,0.9,2096.168);
  pb_epos->SetMarkerColor(kBlue);
  pb_epos->SetMarkerStyle(22);

  TGraphErrors* pp_eposlhc  = new TGraphErrors(1);
  pp_eposlhc->SetName("pp_eposlhc");
  pp_eposlhc->SetTitle("EPOSLHC");
  pp_eposlhc->SetPoint(0,0.95,67.448);
  pp_eposlhc->SetMarkerColor(kGreen-2);
  pp_eposlhc->SetMarkerStyle(23);
  TGraphErrors* pb_eposlhc  = new TGraphErrors(1);
  pb_eposlhc->SetName("pb_eposlhc");
  pb_eposlhc->SetTitle("EPOSLHC");
  pb_eposlhc->SetPoint(0,0.95,2085.057);
  pb_eposlhc->SetMarkerColor(kGreen-2);
  pb_eposlhc->SetMarkerStyle(23);


  TGraphErrors* pp_qgsjetii  = new TGraphErrors(1);
  pp_qgsjetii->SetName("pp_qgsjetii");
  pp_qgsjetii->SetTitle("QGSJETII-04");
  pp_qgsjetii->SetPoint(0,1,69.422);
  pp_qgsjetii->SetMarkerColor(kMagenta);
  pp_qgsjetii->SetMarkerStyle(24);
  TGraphErrors* pb_qgsjetii  = new TGraphErrors(1);
  pb_qgsjetii->SetName("pb_qgsjetii");
  pb_qgsjetii->SetTitle("QGSJETII-04");
  pb_qgsjetii->SetPoint(0,1,2176.422);
  pb_qgsjetii->SetMarkerColor(kMagenta);
  pb_qgsjetii->SetMarkerStyle(24);


  TGraphErrors* pp_fit  = new TGraphErrors(1);
  pp_fit->SetName("pp_fit");
  pp_fit->SetTitle("COMPETE and Glauber");
  pp_fit->SetPoint(0,1.05,70);
  pp_fit->SetMarkerColor(kBlack);
  pp_fit->SetMarkerStyle(20);
  TGraphErrors* pb_fit  = new TGraphErrors(1);
  pb_fit->SetName("pb_fit");
  pb_fit->SetTitle("COMPETE and Glauber");
  pb_fit->SetPoint(0,1,2110);
  pb_fit->SetMarkerColor(kBlack);
  pb_fit->SetMarkerStyle(20);


  TGraphErrors* pp_this  = new TGraphErrors(1);
  pp_this->SetName("pp_this");
  pp_this->SetTitle("THIS WORK");
  pp_this->SetPoint(0,1.05,0);
  pp_this->SetMarkerColor(kRed);
  pp_this->SetMarkerStyle(21);
  TGraphErrors* pb_this  = new TGraphErrors(1);
  pb_this->SetName("pb_this");
  pb_this->SetTitle("THIS WORK");
  pb_this->SetPoint(0,1.05,2173);
  pb_this->SetPointError(0,0,2173*0.081);
  pb_this->SetMarkerColor(kRed);
  pb_this->SetMarkerStyle(21);


  ///RANGE
  pp_epos->GetYaxis()->SetTitle("sigma_{inel} / b");
  pp_epos->GetXaxis()->SetNdivisions(0);
  pp_epos->GetXaxis()->SetLimits(0.5,1.5);
  pp_epos->GetHistogram()->SetMaximum(80);
  pp_epos->GetHistogram()->SetMinimum(65);

  pb_epos->GetXaxis()->SetLimits(0.5,1.5);
  pb_epos->GetXaxis()->SetNdivisions(0);
  pb_epos->GetHistogram()->SetMaximum(2400);
  pb_epos->GetHistogram()->SetMinimum(1800);
  ///RANGE

  can1->cd(1);
  pp_epos->Draw("AP");
  pp_eposlhc->Draw("P");
  pp_qgsjetii->Draw("P");
  pp_fit->Draw("P");
  pp_this->Draw("P");

  can1->cd(2);
  pb_epos->Draw("AP");
  pb_eposlhc->Draw("P");
  pb_qgsjetii->Draw("P");
  pb_fit->Draw("P");
  pb_this->Draw("P");

  TLegend* leg1 = ((TCanvas*)(can1->GetPad(1)))->BuildLegend(0.55,0.7,0.85,0.9);
  leg1->SetFillColor(kWhite);
  
  can1->cd(1);
  DataText(true,true,"pp, #sqrt{s_{_{NN}}}=5.02 TeV");
  can1->cd(2);
  CMSPreliminary();
  DataText(true,true);
  can1->SaveAs((string("plots/concl1.pdf")+string(".pdf")).c_str());
}
