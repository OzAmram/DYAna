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

const bool write_out = true;
bool prelim = false;
char *plot_dir = "Paper_plots/EMu_plots/";
//char *plot_dir = "PAS_plots/EMu_plots/";
//char *plot_dir = "Misc_plots/mu_prefire_check_before/";
//char *plot_label = "e#mu Control Region";
//char *plot_label = "#bf{e#mu control region}";
char *plot_label = "";
char *fin_name = "EMu/saved_hists.root";



void draw_cmp_from_saved(){


    TH1F *data_m ,*diboson_m, *QCD_m, *top_m, *dy_m, *wt_m;
    TH1F *data_cost, *diboson_cost, *QCD_cost, *top_cost, *dy_cost, *wt_cost;
    TH1F *data_rap, *diboson_rap, *QCD_rap, *top_rap, *dy_rap, *wt_rap;



    TFile *fin = new TFile(fin_name, "READ");


    int year_start = 2016;
    int year_stop = 2018;

    for(int year = year_start; year <= year_stop; year++){

        char year_str[80];
        sprintf(year_str, "y%i", year);


        fin->cd(year_str);


        TH1F *data_m_ = (TH1F *) gDirectory->Get("data_m");
        TH1F *diboson_m_ = (TH1F *) gDirectory->Get("diboson_m");
        TH1F *QCD_m_ = (TH1F *) gDirectory->Get("QCD_m");
        TH1F *ttbar_m_ = (TH1F *) gDirectory->Get("ttbar_m");
        TH1F *dy_m_ = (TH1F *) gDirectory->Get("dy_m");
        TH1F *wt_m_ = (TH1F *) gDirectory->Get("wt_m");

        TH1F *data_cost_ = (TH1F *) gDirectory->Get("data_cost");
        TH1F *diboson_cost_ = (TH1F *) gDirectory->Get("diboson_cost");
        TH1F *QCD_cost_ = (TH1F *) gDirectory->Get("QCD_cost");
        TH1F *ttbar_cost_ = (TH1F *) gDirectory->Get("ttbar_cost");
        TH1F *dy_cost_ = (TH1F *) gDirectory->Get("dy_cost");
        TH1F *wt_cost_ = (TH1F *) gDirectory->Get("wt_cost");


        TH1F *data_rap_ = (TH1F *) gDirectory->Get("data_rap");
        TH1F *diboson_rap_ = (TH1F *) gDirectory->Get("diboson_rap");
        TH1F *QCD_rap_ = (TH1F *) gDirectory->Get("QCD_rap");
        TH1F *ttbar_rap_ = (TH1F *) gDirectory->Get("ttbar_rap");
        TH1F *dy_rap_ = (TH1F *) gDirectory->Get("dy_rap");
        TH1F *wt_rap_ = (TH1F *) gDirectory->Get("wt_rap");



        TH1F *top_cost_ = (TH1F *) ttbar_cost_->Clone("top_cost");
        TH1F *top_m_ = (TH1F *) ttbar_m_->Clone("top_m");
        TH1F *top_rap_ = (TH1F *) ttbar_rap_->Clone("top_rap");

        top_cost_->Add(wt_cost_);
        top_m_->Add(wt_m_);
        top_rap_->Add(wt_rap_);

        if(year == year_start){
            printf("%i \n", year);
            data_m = data_m_; diboson_m = diboson_m_; QCD_m = QCD_m_; top_m = top_m_; dy_m = dy_m_; wt_m = wt_m_;
            data_cost = data_cost_; diboson_cost = diboson_cost_; QCD_cost = QCD_cost_; top_cost = top_cost_; dy_cost = dy_cost_; wt_cost = wt_cost_;
            data_rap = data_rap_; diboson_rap = diboson_rap_; QCD_rap = QCD_rap_; top_rap = top_rap_; dy_rap = dy_rap_; wt_rap = wt_rap_;
        }
        else{
            printf("%i \n", year);
            data_m->Add(data_m_);
            diboson_m->Add(diboson_m_);
            QCD_m->Add(QCD_m_);
            top_m->Add(top_m_);
            dy_m->Add(dy_m_);
            wt_m->Add(wt_m_);


            data_cost->Add(data_cost_);
            diboson_cost->Add(diboson_cost_);
            QCD_cost->Add(QCD_cost_);
            top_cost->Add(top_cost_);
            dy_cost->Add(dy_cost_);
            wt_cost->Add(wt_cost_);

            data_rap->Add(data_rap_);
            diboson_rap->Add(diboson_rap_);
            QCD_rap->Add(QCD_rap_);
            top_rap->Add(top_rap_);
            dy_rap->Add(dy_rap_);
            wt_rap->Add(wt_rap_);
        }
    }

    setHistError(QCD_m, qcd_sys_unc);
    setHistError(QCD_cost, qcd_sys_unc);
    setHistError(QCD_rap, qcd_sys_unc);

    setHistError(diboson_m, diboson_sys_unc);
    setHistError(diboson_cost, diboson_sys_unc);
    setHistError(diboson_rap, diboson_sys_unc);

    setHistError(dy_m, dy_sys_unc);
    setHistError(dy_cost, dy_sys_unc);
    setHistError(dy_rap, dy_sys_unc);

    setHistError(top_m, top_sys_unc);
    setHistError(top_cost, top_sys_unc);
    setHistError(top_rap, top_sys_unc);

    setHistError(wt_m, top_sys_unc);
    setHistError(wt_cost, top_sys_unc);
    setHistError(wt_rap, top_sys_unc);



    TH1F *total_mc = (TH1F *) diboson_cost->Clone("total_mc");
    total_mc->Add(QCD_cost);
    total_mc->Add(dy_cost);
    total_mc->Add(wt_cost);
    total_mc->Add(top_cost);


    Double_t data_B = data_cost->Integral(1,n_cost_bins/2);
    Double_t data_F = data_cost->Integral(n_cost_bins/2 + 1, n_cost_bins);
    Double_t data_AFB = (data_F - data_B)/(data_F+data_B);
   
    Double_t data_dAFB = AFB_counting_unc(data_F, data_B, sqrt(data_F), sqrt(data_B));
    printf("F %.0f, B %.0f \n", data_F, data_B);
    printf("AFB %.3f +/- %.3f \n", data_AFB, data_dAFB);



    Double_t total_mc_dF, total_mc_dB;
    Double_t total_mc_B = total_mc->IntegralAndError(1,n_cost_bins/2, total_mc_dB);
    Double_t total_mc_F = total_mc->IntegralAndError(n_cost_bins/2 + 1, n_cost_bins, total_mc_dF);
    Double_t total_mc_AFB = (total_mc_F - total_mc_B)/(total_mc_F+total_mc_B);
   
    Double_t total_mc_dAFB = AFB_counting_unc(total_mc_F, total_mc_B, total_mc_dF, total_mc_dB);
    printf("total_mc F %.0f, B %.0f \n", total_mc_F, total_mc_B);
    printf("total_mc AFB %.3f +/- %.3f \n", total_mc_AFB, total_mc_dAFB);



    float mbin_base = 10.;
    binwidth_normalize(data_m, mbin_base);
    binwidth_normalize(diboson_m, mbin_base);
    binwidth_normalize(QCD_m, mbin_base);
    binwidth_normalize(wt_m, mbin_base);
    binwidth_normalize(top_m, mbin_base);
    binwidth_normalize(dy_m, mbin_base);





    dy_cost->SetFillColor(DY_c);
    dy_m->SetFillColor(DY_c);
    dy_rap->SetFillColor(DY_c);

    top_cost->SetFillColor(ttbar_c);
    top_m->SetFillColor(ttbar_c);
    top_rap->SetFillColor(ttbar_c);


    wt_cost->SetFillColor(wt_c);
    wt_m->SetFillColor(wt_c);
    wt_rap->SetFillColor(wt_c);

    diboson_cost->SetFillColor(diboson_c);
    diboson_m->SetFillColor(diboson_c);
    diboson_rap->SetFillColor(diboson_c);

    QCD_cost->SetFillColor(qcd_c);
    QCD_m->SetFillColor(qcd_c);
    QCD_rap->SetFillColor(qcd_c);



    dy_cost->SetLineColor(DY_c);
    dy_m->SetLineColor(DY_c);
    dy_rap->SetLineColor(DY_c);

    top_cost->SetLineColor(ttbar_c);
    top_m->SetLineColor(ttbar_c);
    top_rap->SetLineColor(ttbar_c);


    wt_cost->SetLineColor(wt_c);
    wt_m->SetLineColor(wt_c);
    wt_rap->SetLineColor(wt_c);

    diboson_cost->SetLineColor(diboson_c);
    diboson_m->SetLineColor(diboson_c);
    diboson_rap->SetLineColor(diboson_c);

    QCD_cost->SetLineColor(qcd_c);
    QCD_m->SetLineColor(qcd_c);
    QCD_rap->SetLineColor(qcd_c);






    dy_cost->SetFillStyle(DY_style);
    dy_m->SetFillStyle(DY_style);
    dy_rap->SetFillStyle(DY_style);

    top_cost->SetFillStyle(ttbar_style);
    top_m->SetFillStyle(ttbar_style);
    top_rap->SetFillStyle(ttbar_style);



    diboson_cost->SetFillStyle(diboson_style);
    diboson_m->SetFillStyle(diboson_style);
    diboson_rap->SetFillStyle(diboson_style);

    QCD_cost->SetFillStyle(qcd_style);
    QCD_m->SetFillStyle(qcd_style);
    QCD_rap->SetFillStyle(qcd_style);



    setTDRStyle();
    gStyle->SetLegendBorderSize(0);
    //gStyle->SetErrorX(0);



    float cost_bin_size = 2./8;
    int n_rap_bins = 20;
	float rap_bin_low = -2.4;
    float rap_bin_high = 2.4;
    float rap_bin_size = (rap_bin_high - rap_bin_low) / n_rap_bins;



    THStack *m_stack = new THStack("m_stack", "EMu Mass Distribution: Data vs MC ; m_{e#mu} (GeV)");
    m_stack->Add(diboson_m);
    m_stack->Add(QCD_m);
    m_stack->Add(top_m);
    //m_stack->Add(dy_m);
    //


    THStack *cost_stack = new THStack("cost_stack", "EMu Cos(theta) Distribution: Data vs MC ; cos(#theta)");
    cost_stack->Add(diboson_cost);
    cost_stack->Add(QCD_cost);
    cost_stack->Add(top_cost);

    //cost_stack->Add(dy_cost);





    THStack *rap_stack = new THStack("rap_stack", "EMu Cos(theta) Distribution: Data vs MC ; cos(#theta)");
    rap_stack->Add(diboson_rap);
    rap_stack->Add(QCD_rap);
    rap_stack->Add(top_rap);
    //rap_stack->Add(dy_rap);

    gStyle->SetLegendBorderSize(0);
    float x_size = 0.6;
    float y_size = 0.4;
    TLegend *leg1 = new TLegend(x_size, y_size);

    leg1->SetNColumns(2);
    //leg1->SetHeader("#bf{e#mu} control region");
    leg1->SetHeader("e#mu control region");
    //different data error bars
    TLegend *leg2 = (TLegend *) leg1->Clone("leg2");

    data_m->SetLineWidth(2);
    leg1->AddEntry(data_m, "Data", "lpe");
    data_m->SetLineWidth(2);
    leg2->AddEntry(data_m, "Data", "pe");

    //leg1->AddEntry(dy_m, "DY #rightarrow #tau#tau", "f");
    leg1->AddEntry(top_m, "t#bar{t} + single t", "f");
    leg1->AddEntry(QCD_m, "QCD and W+jets", "f");
    leg1->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg1->SetTextSize(0.058);


    leg2->AddEntry(top_m, "t#bar{t} + single t", "f");
    leg2->AddEntry(QCD_m, "QCD and W+jets", "f");
    leg2->AddEntry(diboson_m, "WW + WZ + ZZ", "f");
    leg2->SetTextSize(0.058);

    TLegend *leg3 = (TLegend *) leg2->Clone("leg3");

    TCanvas *c_m, *c_cost,  *c_xf, *c_phi, *c_rap;
    TPad *p_m, *p_cost, *p_xf, *p_phi, *p_rap;
    int iPeriod = 4; 
    
    if(prelim) writeExtraText = true;
    else writeExtraText = false;
    char plt_file[100], y_ax_label[100];
    
    bool logy = true;
    bool logx = false;
    bool draw_sys_uncs = true;
    float ratio_range = 0.4;


    char file_label[20];
    int year = -1;
    if(year_start == 2016 && year_stop == 2018){
        sprintf(file_label, "%s", "Comb");
    }
    else if(year_stop == 2016){
        sprintf(file_label, "%s", "16");
        year = 2016;
    }
    else if(year_stop == 2017){
        sprintf(file_label, "%s", "17");
        year = 2017;
    }
    else if(year_stop == 2018){
        sprintf(file_label, "%s", "18");
        year = 2018;
    }




    float x_start_m = 0.34;
    float y_start_m = 0.5;
    leg1->SetX1(x_start_m);
    leg1->SetX2(x_start_m+x_size);
    leg1->SetY1(y_start_m);
    leg1->SetY2(y_start_m+y_size);



    float hmax = 400000;
    /*
    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;
        */

    sprintf(y_ax_label, "Events / %.0f GeV", mbin_base);
    std::tie(c_m, p_m) = make_stack_ratio_plot(data_m, m_stack,  leg1, "m", "m_{e#mu} (GeV)", y_ax_label, plot_label, hmax, logy,logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_m, year, 11 );
    sprintf(plt_file, "%sEMu%s_m_cmp.pdf", plot_dir, file_label);
    if(write_out) c_m->Print(plt_file);
    sprintf(plt_file, "%sEMu%s_m_cmp.png", plot_dir, file_label);
    if(write_out) c_m->Print(plt_file);


    float x_start_c = 0.34;
    float y_start_c = 0.5;
    leg2->SetX1(x_start_c);
    leg2->SetX2(x_start_c+x_size);
    leg2->SetY1(y_start_c);
    leg2->SetY2(y_start_c+y_size);

    hmax = 50000;
    /*
    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;
        */


    logy = false;
    sprintf(y_ax_label, "Events / %.1f", cost_bin_size);
    std::tie(c_cost, p_cost) = make_stack_ratio_plot(data_cost,  cost_stack, leg2, "cost", "cos #theta_{R}", y_ax_label, plot_label, hmax, logy,logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_cost, year, 11);
    sprintf(plt_file, "%sEMu%s_cost_cmp.pdf", plot_dir, file_label);
    if(write_out) c_cost->Print(plt_file);
    sprintf(plt_file, "%sEMu%s_cost_cmp.png", plot_dir, file_label);
    if(write_out) c_cost->Print(plt_file);

    hmax = 35000;
    /*
    if(year == 2016)
        hmax *= 0.625;
    if(year == 2017)
        hmax *= 0.75;
        */

    leg3->SetX1(x_start_c);
    leg3->SetX2(x_start_c+x_size);
    leg3->SetY1(y_start_c);
    leg3->SetY2(y_start_c+y_size);

    logy = false;
    sprintf(y_ax_label, "Events / %.1f", rap_bin_size);
    std::tie(c_rap, p_rap) = make_stack_ratio_plot(data_rap, rap_stack,  leg3, "rap", "dilepton rapidity", y_ax_label,  plot_label, hmax, logy, logx, draw_sys_uncs, ratio_range);
    CMS_lumi(p_rap, year, 11);
    sprintf(plt_file, "%sEMu%s_rap_cmp.pdf", plot_dir, file_label);
    if(write_out) c_rap->Print(plt_file);
    sprintf(plt_file, "%sEMu%s_rap_cmp.png", plot_dir, file_label);
    if(write_out) c_rap->Print(plt_file);


}
