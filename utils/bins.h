#ifndef BINS_H
#define BINS_H

const int n_xf_bins = 4;
Float_t xf_bins[] = {0., 0.04, 0.07, 0.10, 1.0};
//const int n_cost_bins = 10;
//Float_t cost_bins[] = {-1.0, -.8, -.6, -.4, -.2, 0., 0.2, 0.4, 0.6, 0.8, 1.0};
const int n_cost_bins = 8;
Float_t cost_bins[] = {-1.0, -.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.0};
const int n_cost_ss_bins = n_cost_bins/2;
Float_t cost_ss_bins[] = {-1.0, -0.75, -0.5, -0.25, 0.0};
const int n_m_bins = 8;
Double_t m_bins[] = {150, 171, 200,  250, 320, 510, 700, 1000, 14000};

Double_t amc_alpha[n_m_bins] =        {0.056, 0.056, 0.047, 0.055, 0.042, 0.030, 0.018, 0.012};
Double_t amc_alpha_unc[n_m_bins] =    {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.007, 0.007};



#endif