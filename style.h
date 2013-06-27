#ifndef _include_style_h_
#define _include_style_h_

#include <TLegend.h>

void style();

// #################################################################################
// helper function to specify that lead lead data is used
void CMSText(const bool data=1, const bool left=true, const bool top=true, const std::string str3="", const std::string str2="", const std::string str1="pPb, #sqrt{s_{_{NN}}}=5.02 TeV");

// #################################################################################
// helper function for legend
void SetLegAtt(TLegend* leg, const double txtsize = 1.);

#endif
