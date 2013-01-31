#include <style.h>

#include <TGaxis.h>
#include <TStyle.h>
#include <TROOT.h>

void style() {

  gROOT->ProcessLine(" .L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.15); //0.40
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.20);
  //gStyle->SetStatH(0.4);
  gStyle->SetStatY(0.82);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.3);
  //gStyle->SetStatH(0.7);
  gStyle->SetOptStat(0000); //0000 //112210
  gStyle->SetOptFit(0000);
  //
  //Can you please also increase the number of labels in the x-axis ?
  //(SetMoreLogLabels, SetNoExponent)
  //Brighter colors ? (e.g. 51 (violet), 95
  //(orange), 33-38-40 (greyish) ...)
  //(Barely readable) x,y-axes titles and labels ...
  //Can you please increase those too ? (SetLabelSize(0.05) and
  //SetTitleSize(0.05) at least)
  //
  gStyle->SetTitleFontSize(0.15);
  gStyle->SetOptTitle(0);

  gROOT->ForceStyle();

  TGaxis::SetMaxDigits(3);
}

