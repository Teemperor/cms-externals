#ifndef __TAUSPINNER_EXAMPLE_COMMON_DEF__
#define __TAUSPINNER_EXAMPLE_COMMON_DEF__

// ROOT headers
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"

void draw3(TCanvas *c, TH1D *h1, TH1D *h2, TH1D *hdiff, bool isW = 0, TString fit = "", Double_t x0 = 0, Double_t x1 = 1);

#endif //__TAUSPINNER_EXAMPLE_COMMON_DEF__
