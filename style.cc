#include <style.h>

#include <TGaxis.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPaveText.h>

void style() {

  gROOT->ProcessLine(" .L tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.15); //0.40
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.20);
  gStyle->SetPadRightMargin(0.05); //0.40
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
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

// #################################################################################
// helper function to specify that proton/lead data is used
void CMSText(const bool data, const bool left, const bool top, const std::string str3, const std::string str2, const std::string str1) {
  const double marg = 0.005, margY = 0.03;
  std::string text;
  if(data)
    text=std::string("CMS Preliminary");
  else
    text=std::string("CMS Simulation");
  double size = 0.05;
  if(str1.length() != 0)
    size += 0.05;
  if(str2.length() != 0)
    size += 0.05;
  if(str3.length() != 0)
    size += 0.05;

  TPaveText* txt
    = new TPaveText(left ? gStyle->GetPadLeftMargin()+marg : 1-gStyle->GetPadRightMargin()-marg-0.30, //xlow
                    top ? 1-gStyle->GetPadTopMargin()-margY-size : gStyle->GetPadBottomMargin() + margY, //ylow
                    left ? gStyle->GetPadLeftMargin()+marg+0.30 : 1-gStyle->GetPadRightMargin()-marg, //xhigh
                    top ? 1-gStyle->GetPadTopMargin()-margY : gStyle->GetPadBottomMargin() + margY + size, //yhigh
                    "NDC b t l");
    txt->SetTextAlign(10*(left ? 1 : 3) + (top ? 1 : 3));
    txt->SetTextFont(42);
    txt->SetFillStyle(0);
    txt->SetTextColor(kBlack);
    txt->SetTextSize(0.033);
    txt->SetBorderSize(0);

    txt->AddText(text.c_str());
    if(str1.length()!=0)
      txt->AddText(str1.c_str());
    if(str3.length()!=0)
      txt->AddText(str3.c_str());
    if(str2.length()!=0)
      txt->AddText(str2.c_str());

    txt->Draw();
}

// #################################################################################
// helper function for legend
void SetLegAtt(TLegend* leg,const double txtsize)
{
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.033*txtsize);
}
