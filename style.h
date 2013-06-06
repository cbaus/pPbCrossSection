#ifndef _include_style_h_
#define _include_style_h_

#include <TLegend.h>

void style();


// #################################################################################
// helper function to mark plot as "CMS PRELIMINARY"
void CMSPreliminary();

// #################################################################################
// helper function to specify that lead lead data is used
void DataText(const bool left=true, const bool top=true, const std::string txt="pPb, #sqrt{s_{_{NN}}}=5.02 TeV");

// #################################################################################
// helper function for legend
void SetLegAtt(TLegend* leg, const double txtsize = 0.033);

#endif
