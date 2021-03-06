#define STAND_ALONE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFitter.h"
#include "TSystem.h"
#include "fit_gen_cost.C"
#include "fit_amc_gen_cost.C"
#include "../../utils/PlotUtils.C"




void draw_generator_cmp(){
    gStyle->SetOptStat(0);
    int year = 2017;
    TFile *f_mad = TFile::Open("../analyze/output_files/DY17_gen_level_april6.root");
    TTree *t_mad_mu = (TTree *) f_mad->Get("T_gen_mu");
    TTree *t_mad_el = (TTree *) f_mad->Get("T_gen_el");

    //TFile *f_pwg = TFile::Open("../generator_stuff/root_files/powheg_m150_april30.root");
    TFile *f_pwg = TFile::Open("../generator_stuff/root_files/powheg_m150_evts.root");
    TTree *t_pwg = (TTree *)f_pwg->Get("T_lhe");

    float pt_low = 0.;
    float pt_high = 10000.;
    float rap_low = -1.;
    float rap_high = 100.;
    float m_low = 150;
    float m_high = 200.;
    bool phot_ind = false;

    char title[80];
    char plt_title1[80], plt_title2[80], plt_title3[80], plt_title4[80]; 

    sprintf(title, "POWHEG vs. aMC@NLO (2017) (%.0f < M < %.0f)", m_low, m_high);
    sprintf(plt_title1, "pwg_vs_mad%i_m%.0f_cost_st_cmp.pdf", year - 2000, m_low);
    sprintf(plt_title2, "pwg_vs_mad%i_m%.0f_pt_cmp.pdf", year - 2000, m_low);
    sprintf(plt_title3, "pwg_vs_mad%i_m%.0f_xf_cmp.pdf", year - 2000, m_low);
    sprintf(plt_title4, "pwg_vs_mad%i_m%.0f_cost_r_cmp.pdf", year - 2000, m_low);

    TH1F *h_mad_m = new TH1F("h_mad_m", "", 30, 150,1000);
    TH1F *h_pwg_m = new TH1F("h_pwg_m", "", 30, 150,1000);

    TH1F *h_pwg_cost_st = new TH1F("h_pwg_cost_st", title, 20, -1., 1.);
    TH1F *h_pwg_cost_r = new TH1F("h_pwg_cost_r", title, 20, -1., 1.);
    TH1F *h_pwg_pt = new TH1F("h_pwg_pt", title, 20, 0., 300.);
    TH1F *h_pwg_xf = new TH1F("h_pwg_xf", title, 20, 0., 0.5);

    TH1F *h_mad_cost_st = new TH1F("h_mad_cost_st", title, 20, -1., 1.);
    TH1F *h_mad_cost_r = new TH1F("h_mad_cost_r", title, 20, -1., 1.);
    TH1F *h_mad_pt = new TH1F("h_mad_pt", title, 20, 0., 300.);
    TH1F *h_mad_xf = new TH1F("h_mad_xf", title, 20, 0.,.5);


    make_amc_gen_cost(t_mad_mu,  h_mad_m, h_mad_cost_st, h_mad_cost_r, h_mad_pt, h_mad_xf, m_low, m_high, pt_low, pt_high, rap_low, rap_high);
    make_amc_gen_cost(t_mad_el,  h_mad_m, h_mad_cost_st, h_mad_cost_r, h_mad_pt, h_mad_xf, m_low, m_high, pt_low, pt_high, rap_low, rap_high);
    make_gen_cost(t_pwg,  h_pwg_cost_st, h_pwg_cost_r, h_pwg_pt, h_pwg_xf, m_low, m_high, phot_ind);

    h_mad_cost_st->Scale(1./h_mad_cost_st->Integral());
    h_pwg_cost_st->Scale(1./h_pwg_cost_st->Integral());

    h_mad_xf->Scale(1./h_mad_xf->Integral());
    h_pwg_xf->Scale(1./h_pwg_xf->Integral());

    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());


    h_mad_pt->Scale(1./h_mad_pt->Integral());
    h_pwg_pt->Scale(1./h_pwg_pt->Integral());


    h_mad_cost_r->Scale(1./h_mad_cost_r->Integral());
    h_pwg_cost_r->Scale(1./h_pwg_cost_r->Integral());


    make_ratio_plot(plt_title1, h_mad_cost_st, "amc@NLO",h_pwg_cost_st, "POWHEG", "aMC/POWHEG", "cos(#theta_{*})", false);
    make_ratio_plot(plt_title2, h_mad_pt, "amc@NLO",h_pwg_pt, "POWHEG", "aMC/POWHEG", "dilepton p_{T} (GeV)", false);
    make_ratio_plot(plt_title3, h_mad_xf, "amc@NLO",h_pwg_xf, "POWHEG", "aMC/POWHEG", "xF", false);
    make_ratio_plot(plt_title4, h_mad_cost_r, "amc@NLO",h_pwg_cost_r, "POWHEG", "aMC/POWHEG", "cos(#theta_{r})", false);

    return;
}

