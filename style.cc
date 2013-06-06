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

// ########################################################
void CMSPreliminary() {
  TLatex* txt = new TLatex();
  txt->SetNDC();
  txt->SetTextFont(62);
  txt->SetTextColor(kGray+1);
  txt->SetTextSize(0.05);
  txt->SetTextAlign(32);
  txt->DrawLatex(1-gStyle->GetPadRightMargin(), 0.95, "CMS PRELIMINARY");
}
// #################################################################################
// helper function to specify that proton/lead data is used
void DataText(const bool left, const bool top, const std::string str) {
  const double marg = 0.005, margY = 0.03;
  TPaveText* txt
    = new TPaveText(left ? gStyle->GetPadLeftMargin()+marg : 1-gStyle->GetPadRightMargin()-marg-0.30,
                    top ? 1-gStyle->GetPadTopMargin()-margY-0.04 : gStyle->GetPadBottomMargin() + margY,
                    left ? gStyle->GetPadLeftMargin()+marg+0.30 : 1-gStyle->GetPadRightMargin()-marg,
                    top ? 1-gStyle->GetPadTopMargin()-margY : gStyle->GetPadBottomMargin() + margY + 0.04,
                    "NDC b t l");
    txt->SetTextAlign(10*(left ? 1 : 3) + (top ? 1 : 3));
    txt->SetTextFont(42);
    txt->SetFillStyle(0);
    txt->SetTextColor(kBlack);
    txt->SetTextSize(0.033);
    txt->SetBorderSize(0);
    txt->AddText(str.c_str());
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
  leg->SetTextSize(txtsize);
}
