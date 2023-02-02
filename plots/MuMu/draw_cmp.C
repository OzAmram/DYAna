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
#include "Math/Functor.h"
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/HistMaker.C"
#include "../../utils/PlotUtils.C"
#include "../../utils/root_files.h"
#include "../../utils/Colors.h"



const int type = FLAG_MUONS;
//const int year = 2016;

const bool write_out = true;
char *plot_dir = "AN_plots/prefit_kinematics_test_noptrw/";

//const bool write_out = false;
//char *plot_dir = "Paper_plots/prefit_kinematics/";

char *plot_label = "";




void draw_cmp(){


for(int year = 2016; year<=2018; year++){    
    setTDRStyle();
    init(year);
    init_indv_bkgs(year);
    setup_all_SFs(year);

    int n_pt_bins1 = 7;
    Float_t pt_bins1[] = {0., 10., 20., 30., 50., 70., 100., 300.};

    TH1F *dy_pt = new TH1F("dy_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *dy_tautau_pt = new TH1F("dy_tautau_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *data_pt = new TH1F("data_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *ttbar_pt = new TH1F("ttbar_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *diboson_pt = new TH1F("diboson_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *wt_pt = new TH1F("wt_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *QCD_pt = new TH1F("QCD_pt", "dy signal", n_pt_bins1, pt_bins1);
    TH1F *gg_pt = new TH1F("gg_pt", "dy signal", n_pt_bins1, pt_bins1);

    int n_xf_bins1 = 5;
    float xf_bins1[] = {0.,0.04, 0.07, 0.1, 0.2, 0.5};
    TH1F *dy_xf = new TH1F("dy_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *dy_tautau_xf = new TH1F("dy_tautau_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *data_xf = new TH1F("data_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *ttbar_xf = new TH1F("ttbar_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *diboson_xf = new TH1F("diboson_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *wt_xf = new TH1F("wt_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *QCD_xf = new TH1F("QCD_xf", "dy signal", n_xf_bins1,  xf_bins1);
    TH1F *gg_xf = new TH1F("gg_xf", "dy signal", n_xf_bins1,  xf_bins1);

    int n_m_bins = 10;
    float mbin_base = 10.;
    //Float_t mbins1[] = {170.,200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1800., 2400.};
    Float_t mbins1[] = {500., 550., 600., 700., 800., 900., 1000., 1200., 1400., 1800., 2400.};
    TH1F *data_m = new TH1F("data_m", "Data Dimuon Mass Distribution", n_m_bins, mbins1);
    TH1F *dy_m = new TH1F("dy_m", "dy Signal (qqbar, qglu, qbarglu)", n_m_bins, mbins1);
    TH1F *dy_tautau_m = new TH1F("dy_tautau_m", "dy no signal (qq, gluglu qbarqbar)", n_m_bins, mbins1);
    TH1F *ttbar_m = new TH1F("ttbar_m", "TTBar Background", n_m_bins, mbins1);
    TH1F *diboson_m = new TH1F("diboson_m", "DiBoson (WW, WZ, ZZ)", n_m_bins, mbins1);
    TH1F *QCD_m = new TH1F("QCD_m", "QCD", n_m_bins, mbins1);
    TH1F *gg_m = new TH1F("gg_m", "QCD", n_m_bins, mbins1);
    TH1F *WJets_m = new TH1F("WJets_m", "WJets", n_m_bins, mbins1);
    TH1F *wt_m = new TH1F("wt_m", "tw + #bar{t}w", n_m_bins, mbins1);

    int n_cost_bins = 10;
    float cost_bin_size = 2./n_cost_bins;
    TH1F *data_cost = new TH1F("data_cost", "Data", n_cost_bins, -1.,1.);
    TH1F *dy_cost = new TH1F("dy_cost", "dy Signal (qqbar, qglu, qbarglu)", n_cost_bins, -1,1);
    TH1F *dy_tautau_cost = new TH1F("dy_tautau_cost", "dy no signal (qq, gluglu qbarqbar)", n_cost_bins, -1.,1.);
    TH1F *ttbar_cost = new TH1F("ttbar_cost", "TTbar Background", n_cost_bins, -1.,1.);
    TH1F *diboson_cost = new TH1F("diboson_cost", "DiBoson (WW, WZ,ZZ)", n_cost_bins, -1,1);
    TH1F *QCD_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1,1);
    TH1F *gg_cost = new TH1F("gg_cost", "QCD", n_cost_bins, -1,1);
    TH1F *WJets_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1,1);
    TH1F *wt_cost = new TH1F("wt_cost", "tw + #bar{t}w", n_cost_bins, -1,1);

    int n_phi_bins = 20;
    TH1F *data_phi = new TH1F("data_phi", "Data", n_phi_bins, -4.,4.);
    TH1F *dy_phi = new TH1F("dy_phi", "dy Signal (qqbar, qglu, qbarglu)", n_phi_bins, -4,4);
    TH1F *dy_tautau_phi = new TH1F("dy_tautau_phi", "dy no signal (qq, gluglu qbarqbar)", n_phi_bins, -4.,4.);
    TH1F *ttbar_phi = new TH1F("ttbar_phi", "TTbar Background", n_phi_bins, -4.,4.);
    TH1F *diboson_phi = new TH1F("diboson_phi", "DiBoson (WW, WZ,ZZ)", n_phi_bins, -4,4);
    TH1F *QCD_phi = new TH1F("QCD_phi", "QCD", n_phi_bins, -4,4);
    TH1F *gg_phi = new TH1F("gg_phi", "QCD", n_phi_bins, -4,4);
    TH1F *WJets_phi = new TH1F("WJets_phi", "WJets", n_phi_bins, -4,4);
    TH1F *wt_phi = new TH1F("wt_phi", "tw + #bar{t}w", n_phi_bins, -4,4);

    int n_rap_bins = 20;
    float rap_bin_size = 5. / n_rap_bins;
    TH1F *data_rap = new TH1F("data_rap", "Data", n_rap_bins, -2.5,2.5);
    TH1F *dy_rap = new TH1F("dy_rap", "dy Signal (qqbar, qglu, qbarglu)", n_rap_bins, -2.5,2.5);
    TH1F *dy_tautau_rap = new TH1F("dy_tautau_rap", "dy no signal (qq, gluglu qbarqbar)", n_rap_bins, -2.5,2.5);
    TH1F *ttbar_rap = new TH1F("ttbar_rap", "TTbar Background", n_rap_bins, -2.5,2.5);
    TH1F *diboson_rap = new TH1F("diboson_rap", "DiBoson (WW, WZ,ZZ)", n_rap_bins, -2.5,2.5);
    TH1F *QCD_rap = new TH1F("QCD_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *gg_rap = new TH1F("gg_rap", "QCD", n_rap_bins, -2.5,2.5);
    TH1F *WJets_rap = new TH1F("WJets_rap", "WJets", n_rap_bins, -2.5,2.5);
    TH1F *wt_rap = new TH1F("wt_rap", "tw + #bar{t}w", n_rap_bins, -2.5,2.5);

    dy_tautau_cost->SetFillColor(tautau_c);
    dy_tautau_m->SetFillColor(tautau_c);
    dy_tautau_pt->SetFillColor(tautau_c);
    dy_tautau_xf->SetFillColor(tautau_c);
    dy_tautau_phi->SetFillColor(tautau_c);
    dy_tautau_rap->SetFillColor(tautau_c);

    dy_cost->SetFillColor(DY_c);
    dy_m->SetFillColor(DY_c);
    dy_pt->SetFillColor(DY_c);
    dy_xf->SetFillColor(DY_c);
    dy_phi->SetFillColor(DY_c);
    dy_rap->SetFillColor(DY_c);

    ttbar_cost->SetFillColor(ttbar_c);
    ttbar_m->SetFillColor(ttbar_c);
    ttbar_pt->SetFillColor(ttbar_c);
    ttbar_xf->SetFillColor(ttbar_c);
    ttbar_phi->SetFillColor(ttbar_c);
    ttbar_rap->SetFillColor(ttbar_c);


    wt_cost->SetFillColor(wt_c);
    wt_m->SetFillColor(wt_c);
    wt_pt->SetFillColor(wt_c);
    wt_xf->SetFillColor(wt_c);
    wt_phi->SetFillColor(wt_c);
    wt_rap->SetFillColor(wt_c);

    diboson_cost->SetFillColor(diboson_c);
    diboson_m->SetFillColor(diboson_c);
    diboson_pt->SetFillColor(diboson_c);
    diboson_xf->SetFillColor(diboson_c);
    diboson_phi->SetFillColor(diboson_c);
    diboson_rap->SetFillColor(diboson_c);

    QCD_cost->SetFillColor(qcd_c);
    QCD_m->SetFillColor(qcd_c);
    QCD_pt->SetFillColor(qcd_c);
    QCD_xf->SetFillColor(qcd_c);
    QCD_phi->SetFillColor(qcd_c);
    QCD_rap->SetFillColor(qcd_c);

    gg_cost->SetFillColor(gamgam_c);
    gg_m->SetFillColor(gamgam_c);
    gg_pt->SetFillColor(gamgam_c);
    gg_xf->SetFillColor(gamgam_c);
    gg_phi->SetFillColor(gamgam_c);
    gg_rap->SetFillColor(gamgam_c);






    float m_low = 500.;
    float m_high = 13000;
    bool ss = false;
    bool dy_ptrw = false;

    make_m_cost_pt_xf_hist(t_mumu_data, data_m, data_cost, data_pt, data_xf, data_phi, data_rap, true, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_mc, dy_m, dy_cost, dy_pt, dy_xf, dy_phi, dy_rap,              false, type,   year, m_low, m_high, false, dy_ptrw);
    make_m_cost_pt_xf_hist(t_mumu_tautau, dy_tautau_m, dy_tautau_cost, dy_tautau_pt, dy_tautau_xf, dy_tautau_phi, dy_tautau_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_ttbar, ttbar_m, ttbar_cost, ttbar_pt, ttbar_xf, ttbar_phi, ttbar_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_wt, wt_m, wt_cost, wt_pt, wt_xf, wt_phi, wt_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_gamgam, gg_m, gg_cost, gg_pt, gg_xf, gg_phi, gg_rap, false, type,  year, m_low, m_high);
    make_m_cost_pt_xf_hist(t_mumu_diboson, diboson_m, diboson_cost, diboson_pt, diboson_xf, diboson_phi, diboson_rap, false, type,   year, m_low, m_high);

    symmetrize1d(gg_cost);

    //gg_cost->Scale(0.);
    //gg_m->Scale(0.);

    bool ss_qcd = false;
    make_fakerate_est(t_mumu_WJets, t_mumu_QCD, t_mumu_WJets_contam, t_mumu_QCD_contam, QCD_m, QCD_cost, QCD_pt, QCD_xf, QCD_phi, QCD_rap, type, year, m_low, m_high, ss_qcd);
    //QCD_cost->Scale(0.);
    //



    printf("Data integral is %.2f \n", data_m->Integral());
    printf("DY integral is %.2f \n", dy_cost->Integral());
    printf("ttbar integral is %.2f \n", ttbar_cost->Integral());
    printf("diboson integral is %.2f \n", diboson_cost->Integral());
    printf("fakes integral is %.2f \n", QCD_cost->Integral());



    setHistError(QCD_m, qcd_sys_unc);
    setHistError(QCD_cost, qcd_sys_unc);
    setHistError(QCD_pt, qcd_sys_unc);
    setHistError(QCD_xf, qcd_sys_unc);
    setHistError(QCD_rap, qcd_sys_unc);
    setHistError(QCD_phi, qcd_sys_unc);

    setHistError(dy_m, dy_sys_unc);
    setHistError(dy_cost, dy_sys_unc);
    setHistError(dy_pt, dy_sys_unc);
    setHistError(dy_xf, dy_sys_unc);
    setHistError(dy_rap, dy_sys_unc);
    setHistError(dy_phi, dy_sys_unc);

    setHistError(diboson_m, diboson_sys_unc);
    setHistError(diboson_cost, diboson_sys_unc);
    setHistError(diboson_pt, diboson_sys_unc);
    setHistError(diboson_xf, diboson_sys_unc);
    setHistError(diboson_rap, diboson_sys_unc);
    setHistError(diboson_phi, diboson_sys_unc);

    setHistError(ttbar_m, top_sys_unc);
    setHistError(ttbar_cost, top_sys_unc);
    setHistError(ttbar_pt, top_sys_unc);
    setHistError(ttbar_xf, top_sys_unc);
    setHistError(ttbar_rap, top_sys_unc);
    setHistError(ttbar_phi, top_sys_unc);

    setHistError(wt_m, top_sys_unc);
    setHistError(wt_cost, top_sys_unc);
    setHistError(wt_pt, top_sys_unc);
    setHistError(wt_xf, top_sys_unc);
    setHistError(wt_rap, top_sys_unc);
    setHistError(wt_phi, top_sys_unc);


    setHistError(gg_m, gam_sys_unc);
    setHistError(gg_cost, gam_sys_unc);
    setHistError(gg_pt, gam_sys_unc);
    setHistError(gg_xf, gam_sys_unc);
    setHistError(gg_rap, gam_sys_unc);
    setHistError(gg_phi, gam_sys_unc);

    binwidth_normalize(data_m, mbin_base);
    binwidth_normalize(diboson_m, mbin_base);
    binwidth_normalize(QCD_m, mbin_base);
    binwidth_normalize(wt_m, mbin_base);
    binwidth_normalize(ttbar_m, mbin_base);
    binwidth_normalize(gg_m, mbin_base);
    binwidth_normalize(dy_tautau_m, mbin_base);
    binwidth_normalize(dy_m, mbin_base);


    THStack *m_stack = new THStack("m_stack", "MuMu Mass Distribution: Data vs MC ; m_{#mu^{+}#mu^{-}} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(wt_m);
    m_stack->Add(ttbar_m);
    m_stack->Add(gg_m);
    m_stack->Add(dy_tautau_m);
    m_stack->Add(dy_m);


    THStack *cost_stack = new THStack("cost_stack", "Cos(#theta) Distribution: Data vs MC; Cos(#theta)_{r}");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(wt_cost);
    cost_stack->Add(ttbar_cost);
    cost_stack->Add(gg_cost);
    cost_stack->Add(dy_tautau_cost);
    cost_stack->Add(dy_cost);

    THStack *pt_stack = new THStack("pt_stack", "ElEl Pt Distribution: Data vs dy; DiElectron Pt (GeV); Events / GeV");
    pt_stack->Add(diboson_pt);
    pt_stack->Add(QCD_pt);
    pt_stack->Add(wt_pt);
    pt_stack->Add(ttbar_pt);
    pt_stack->Add(gg_pt);
    pt_stack->Add(dy_tautau_pt);
    pt_stack->Add(dy_pt);

    binwidth_normalize(pt_stack);
    binwidth_normalize(data_pt);

    THStack *xf_stack = new THStack("xf_stack", "Di-electron x_F Distribution: Data vs MC; x_F");
    xf_stack->Add(diboson_xf);
    xf_stack->Add(QCD_xf);
    xf_stack->Add(wt_xf);
    xf_stack->Add(ttbar_xf);
    xf_stack->Add(gg_xf);
    xf_stack->Add(dy_tautau_xf);
    xf_stack->Add(dy_xf);

    THStack *phi_stack = new THStack("phi_stack", "DiElectron Phi Distribution: Data vs MC; #phi");
    phi_stack->Add(diboson_phi);
    phi_stack->Add(QCD_phi);
    phi_stack->Add(wt_phi);
    phi_stack->Add(ttbar_phi);
    phi_stack->Add(gg_phi);
    phi_stack->Add(dy_tautau_phi);
    phi_stack->Add(dy_phi);

    THStack *rap_stack = new THStack("rap_stack", "DiElectron Rapidity Distribution: Data vs MC; y");
    rap_stack->Add(diboson_rap);
    rap_stack->Add(QCD_rap);
    rap_stack->Add(wt_rap);
    rap_stack->Add(ttbar_rap);
    rap_stack->Add(gg_rap);
    rap_stack->Add(dy_tautau_rap);
    rap_stack->Add(dy_rap);



    gStyle->SetLegendBorderSize(0);
    float x_size = 0.4;
    float y_size = 0.3;
    float x_start = 0.3;
    float x_end = 0.7;
    float y_start = 0.7;
    float y_end = 1.0;

    //TLegend *leg1 = new TLegend(x_center - x_size/2, y_center - y_size/2, x_center + x_size/2, y_center + y_size/2);
    TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
    leg1->SetNColumns(2);
    leg1->SetHeader("Dimuon Signal Region");
    leg1->AddEntry(data_m, "data", "p");
    leg1->AddEntry(dy_m, "DY Signal", "f");
    leg1->AddEntry(dy_tautau_m, "DY #rightarrow #tau#tau", "f");
    leg1->AddEntry(gg_m, "#gamma#gamma #rightarrow #mu#mu", "f");
    leg1->AddEntry(ttbar_m, "t#bar{t}", "f");
    leg1->AddEntry(wt_m, "tW + #bar{t}W", "f");
    leg1->AddEntry(QCD_m, "QCD + WJets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");
    TLegend *leg3 = (TLegend *) leg1->Clone("leg3");
    TLegend *leg4 = (TLegend *) leg1->Clone("leg4");
    TLegend *leg5 = (TLegend *) leg1->Clone("leg5");
    TLegend *leg6 = (TLegend *) leg1->Clone("leg6");

    leg1->SetX1NDC(0.7);
    leg1->SetX2NDC(0.7);

 
    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    TCanvas *c_m, *c_cost, *c_pt, *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_pt, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    writeExtraText = true;
    char plt_file[100], y_ax_label[100];



    bool logy = true;

    bool logx = false;
    bool draw_sys_uncs = true;
    float ratio_range = 0.5;



    float x_start_m = 0.5;
    float y_start_m = 0.6;
    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);


    sprintf(y_ax_label, "Events/%.0f GeV", mbin_base);
    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack, leg1, "m", "M_{#mu#mu} (GeV)",y_ax_label, plot_label, -1., logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_m, year, 33 );
    sprintf(plt_file, "%sMuMu%i_m_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_m->Print(plt_file);

    
    float x_start_c = 0.57 - x_size/2;
    float y_start_c = 0.14;
    leg2->SetX1(x_start_c);
    leg2->SetX2(x_start_c+x_size);
    leg2->SetY1(y_start_c);
    leg2->SetY2(y_start_c+y_size);

    
    logy = false;
    sprintf(y_ax_label, "Events/%.1f", cost_bin_size);

    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost, cost_stack, leg2, "cost", "cos#theta_{r}",y_ax_label, plot_label,  -1., logy,logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_cost, year, 33);
    sprintf(plt_file, "%sMuMu%i_cost_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_cost->Print(plt_file);

    logy = true;
    //sprintf(y_ax_label, "Events/%.0f GeV", pt_bin_size);
    std::tie(c_pt, p_pt) = make_stack_ratio_plot(data_pt, pt_stack, leg3, "pt", "dimuon pt (GeV)","", plot_label, -1., logy, logx, false, ratio_range);
    CMS_lumi(p_pt, year, 33);
    sprintf(plt_file, "%sMuMu%i_pt_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_pt->Print(plt_file);

    std::tie(c_xf, p_xf) = make_stack_ratio_plot(data_xf, xf_stack, leg4, "xf", "x_F (GeV)","", plot_label, -1., logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_xf, year, 33);
    sprintf(plt_file, "%sMuMu%i_xf_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_xf->Print(plt_file);

    std::tie(c_phi, p_phi) = make_stack_ratio_plot(data_phi, phi_stack, leg5, "phi", "dimuon #phi","", plot_label, -1., logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_phi, year, 33);
    sprintf(plt_file, "%sMuMu%i_phi_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_phi->Print(plt_file);

    leg6->SetX1(x_start_c);
    leg6->SetX2(x_start_c+x_size);
    leg6->SetY1(y_start_c);
    leg6->SetY2(y_start_c+y_size);


    sprintf(y_ax_label, "Events/%.2f", rap_bin_size);
    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack, leg6, "rap", "dimuon rapidity",y_ax_label, plot_label, -1., logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_rap, year, 33);
    sprintf(plt_file, "%sMuMu%i_rap_cmp.pdf", plot_dir, year % 2000);
    if(write_out) c_rap->Print(plt_file);


}
}
    
    
