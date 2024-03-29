#ifndef PLOT_UTILS
#define PLOT_UTILS
#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TObject.h"
#include "TRatioPlot.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"



float dy_sys_unc = 0.12;
float qcd_sys_unc = 0.55;
float diboson_sys_unc = 0.15;
float top_sys_unc = 0.13;
float gam_sys_unc = 0.13;

float computeChi2(TH1 *h){
    // only use on ratio plots, with expected value of 1
    float sum = 0.;
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        float err = h->GetBinError(i);
        if (val > 0. && err > 0.){
            sum += std::pow((val-1)/err,2);
        }
    }
    return sum;
}

void unzero_bins(TH1 *h){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = std::max(h->GetBinContent(i), 1e-8);
        h->SetBinContent(i,val);
    }
}

void binwidth_normalize(TH1 *h, float base = 1.){
    for (int i=1; i <= h->GetNbinsX(); i++){
        float content = h->GetBinContent(i);
        float error = h->GetBinError(i);
        float width = h->GetBinWidth(i);
        h->SetBinContent(i, base*content/width);
        h->SetBinError(i, base*error/width);
    }
}


void binwidth_normalize(THStack *h_stack){
    for(auto h: *h_stack->GetHists()){
        binwidth_normalize( (TH1 *) h);
    }
}

TCanvas *draw_single_plot(TString plt_fname, TH1F *h, TString xlabel, TString ylabel, TString plot_label, bool logy = true, bool write_out = true){ 
    TCanvas *c = new TCanvas(plt_fname, "Histograms", 800,800);
    float hmax = 1.2 * h->GetMaximum();
    if(logy) hmax *=2;
    h->SetMaximum(hmax);
    h->SetMinimum(1e-5);
    if(logy && hmax > 10) h->SetMinimum(1e-1);
    h->SetLineColor(kBlue);
    h->SetLineWidth(3);
    h->Draw("histe");
    h->GetXaxis()->SetTitle(xlabel);
    h->GetYaxis()->SetTitle(ylabel);

    //h->GetYaxis()->SetNdivisions(505);
    //h->GetYaxis()->SetTitleSize(20);
    //h->GetYaxis()->SetTitleFont(43);
    //h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(25);
    // X axis h plot settings
    //h->GetXaxis()->SetTitleSize(20);
    //h->GetXaxis()->SetTitleFont(43);
    //h->GetXaxis()->SetTitleOffset(3.);
    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(25);



    TLatex latext; 
    latext.SetNDC();
    latext.SetTextColor(kBlack);
    latext.SetTextAlign(22); //centered
    latext.SetTextFont(42);
    latext.SetTextSize(0.04);    


    float H = c->GetWh();
    float W = c->GetWw();
    float l = c->GetLeftMargin();
    float t = c->GetTopMargin();
    float r = c->GetRightMargin();
    float b = c->GetBottomMargin();

    float w = 1-l-r;

    //printf("l %.3f,r %.3f,t %.3f,b %.3f,W %.3f,H %.3f \n",l,r,t,b,W,H);


    latext.DrawLatex(l + w/2, 1-2*t ,plot_label);

    if(write_out) c->Print(plt_fname);


    return c;
}








TCanvas *draw_ratio_plot(std::string title, TH1F *h, TH1F *ratio, char axis_label[80], char ratio_label[80], float ratio_min = 0.01, float ratio_max = 2.){
    TCanvas *c = new TCanvas(title.c_str(), "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad((title+"p1").c_str(), "pad1", 0.,0.3,0.98,1.);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    h->SetMinimum(1e-5);
    h->Draw();


    c->cd();
    TPad *pad2 = new TPad((title+"p2").c_str(), "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    ratio->SetMinimum(ratio_min);
    ratio->SetMaximum(ratio_max);
    ratio->Sumw2();
    ratio->SetStats(0);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();

    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(25);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(25);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(3.);
    ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(20);

    return c;
}



TCanvas* make_ratio_plot(std::string title, TH1* h1, char h1_label[80], TH1* h2, char h2_label[80], char ratio_label[80], 
        char axis_label[80], bool logy=false, bool write_out = true, float ratio_min = 0.5, float ratio_max = 1.5, char plot_label[100] = "",
        TH1F *custom_ratio = nullptr){
    //ratio is done as h1/h2

    unzero_bins(h1);
    unzero_bins(h2);

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);

    h1->SetLineWidth(3);
    h2->SetLineWidth(3);




    TCanvas *c = new TCanvas(title.c_str(), "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad((title+"p1").c_str(), "pad1", 0.,0.3,0.98,1.);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0);
    pad1->SetRightMargin(0.03);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    float hmax = 1.2 * std::max(h1->GetMaximum(), h2->GetMaximum());
    if(logy) hmax *=2;
    h1->SetMaximum(hmax);
    h1->SetMinimum(1e-5);
    if(logy && hmax > 10) h1->SetMinimum(1e-1);
    h1->Draw("hist E");
    gStyle->SetEndErrorSize(0);
    h2->Draw("hist E same");
    pad1->cd();
    TLatex latext; 
    latext.SetNDC();
    latext.SetTextColor(kBlack);
    latext.SetTextAlign(22); //centered
    latext.SetTextFont(42);
    latext.SetTextSize(0.07);    


    float H = pad1->GetWh();
    float W = pad1->GetWw();
    float l = pad1->GetLeftMargin();
    float t = pad1->GetTopMargin();
    float r = pad1->GetRightMargin();
    float b = pad1->GetBottomMargin();

    float w = 1-l-r;

    //printf("l %.3f,r %.3f,t %.3f,b %.3f,W %.3f,H %.3f \n",l,r,t,b,W,H);


    latext.DrawLatex(l + w/2, 1-2*t ,plot_label);
    //latext.DrawLatex( l + W/2, 0.9*top, plot_label); 
    pad1->Update();



    gStyle->SetLegendBorderSize(0);
    TLegend *leg1 = new TLegend(0.4, 0.15);
    //TLegend *leg1 = new TLegend(0.4, 0.55, 0.75, 0.8);
    leg1->AddEntry(h1, h1_label, "l");
    leg1->AddEntry(h2, h2_label, "l");
    leg1->Draw();

    c->cd();
    //gPad->BuildLegend();
    c->cd();
    TPad *pad2 = new TPad((title+"p2").c_str(), "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();


    TH1F *ratio;
    if(custom_ratio == nullptr){
        ratio = (TH1F *) h1->Clone("h_ratio");
        ratio->Sumw2();
        ratio->SetStats(0);
        ratio->Divide(h2);
    }
    else{
        ratio = custom_ratio;
        ratio->SetStats(0);
    }
    ratio->SetMinimum(ratio_min);
    ratio->SetMaximum(ratio_max);
    ratio->SetMarkerStyle(21);
    ratio->SetLineColor(kBlack);
    ratio->Draw("ep");
    c->cd();


    float rTS = 0.2;
    float rLS = 0.12;
    float rTOffset = 0.2;



    /*
    ratio->SetTitle("");
    // Y axis ratio plot settings
    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(TS);
    //ratio->GetYaxis()->SetTitleFont(43);
    //ratio->GetYaxis()->SetTitleOffset(1.2);
    //ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetYaxis()->SetLabelSize(15);
    // X axis ratio plot settings
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(TS);
    //ratio->GetXaxis()->SetTitleFont(43);
    //ratio->GetXaxis()->SetTitleOffset(3.);
    //ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio->GetXaxis()->SetLabelSize(LS);
    */

    ratio->GetYaxis()->SetTitle(ratio_label);
    ratio->GetYaxis()->SetNdivisions(205);
    ratio->GetYaxis()->SetTitleSize(rTS);
    ratio->GetYaxis()->SetLabelSize(rLS);
    ratio->GetYaxis()->SetTitleOffset(rTOffset);

    ratio->GetXaxis()->SetNdivisions(808);
    ratio->GetXaxis()->SetTickLength(0.06);
    ratio->GetXaxis()->SetTitle(axis_label);
    ratio->GetXaxis()->SetTitleSize(rTS);
    ratio->GetXaxis()->SetTitleOffset(rTOffset + 0.45);
    ratio->GetXaxis()->SetLabelSize(rLS);

    //lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    //int iPeriod = 4; 
    //CMS_lumi(pad1, iPeriod, 33 );
    if(write_out) c->Print(title.c_str());
    return c;
}


std::tuple<TCanvas*, TPad*> make_stack_ratio_plot(TH1F *h_data,  THStack *h_stack, TLegend *leg, TString label, TString xlabel,  TString ylabel, TString plot_label,
        float hmax =-1., bool logy = true, bool logx= false, bool draw_sys_unc = false, float ratio_range = 1.0, bool draw_chi2 = false, float hmin = -1.){



    TList *stackHists = h_stack->GetHists();
    TH1* sum = (TH1*)stackHists->At(0)->Clone();
    sum->Reset();


    for (int i=0;i<stackHists->GetSize();++i) {
        sum->Add((TH1*)stackHists->At(i));
    }
    //make separate tgraph for uncertainty on 'stacked' hists
    int nbins = sum->GetNbinsX();
    TGraphErrors *ratio_unc = new TGraphErrors(nbins+2);
    for(int i= 0; i<= nbins+1; i++){
        float x,y,ex,ey,content;
        int idx = i;
        float eps = 0.001;

        if(i ==0){
            idx = i+1; 
            x = sum->GetXaxis()->GetBinLowEdge(idx) + eps;
            ex = 0;
        }
        else if(i==nbins+1){
            idx = i-1;
            x = sum->GetXaxis()->GetBinUpEdge(idx) - eps;
            ex = 0;
        }
        else{
            x = sum->GetXaxis()->GetBinCenter(idx);
            ex = sum->GetXaxis()->GetBinWidth(idx);
        }

        content = sum->GetBinContent(idx);
        y = 1.;
        if(content > 0) ey = sum->GetBinError(idx) / sum->GetBinContent(idx);
        else ey = 0.;

        ratio_unc->SetPoint(i, x,y);
        ratio_unc->SetPointError(i, ex,ey);

    }

    TH1* h_stack_err = (TH1F *) sum->Clone();


    //for(int i= 1; i<= nbins; i++) sum->SetBinError(i, 0.);

    unzero_bins(sum);

    auto h_ratio = (TH1F *) h_data->Clone("h_ratio" + label);
    h_ratio->Sumw2();
    h_ratio->SetStats(0);
    bool do_diff = false;

    float center = 1.0;
    if(do_diff){
        center = 0.0;
        h_ratio->Add(sum, -1.);
    }
    h_ratio->Divide(sum);


    TCanvas *c = new TCanvas("c_" + label, "Histograms", 200, 10, 900, 700);
    TPad *pad1 = new TPad("pad1" + label, "pad1", 0.,0.3,0.98,1.);
    pad1->SetTopMargin(0.07);
    pad1->SetBottomMargin(0);
    pad1->SetRightMargin(0.03);
    pad1->Draw();
    pad1->cd();
    if(logy) pad1->SetLogy();
    if(logx) pad1->SetLogx();
    h_stack->Draw("hist");



    if(hmin < 0){
        hmin = 0.1;
        if(logy && h_stack->GetMinimum() > 0.1) hmin = 10.;
    }
    if(hmax <= 0. ) hmax = 1.2 * std::max(h_stack->GetMaximum(), h_data->GetMaximum());
    if(logy) hmax *=2;
    h_stack->SetMaximum(hmax);
    h_stack->SetMinimum(hmin);

    if(draw_sys_unc){
        gStyle->SetHatchesLineWidth(2);
        //gStyle->SetHatchesSpacing(1);
        h_stack_err->SetLineColor(kWhite);
        h_stack_err->SetFillColor(kBlack);
        h_stack_err->SetFillStyle(3353);
        h_stack_err->SetMarkerStyle(0);
        h_stack_err->SetMarkerSize(0.001);
        h_stack_err->Draw("e2 same");

         leg->AddEntry(h_stack_err, "Sys. unc.", "f");
    }


    bool const_size = h_data->GetXaxis()->GetBinWidth(1) == h_data->GetXaxis()->GetBinWidth(h_data->GetXaxis()->GetNbins());

    gStyle->SetEndErrorSize(0);
    h_data->SetMarkerStyle(kFullCircle);
    h_data->SetMarkerColor(1);
    h_data->SetMarkerSize(1.3);
    h_data->SetLineWidth(2.5);
    if(const_size){
        printf("const size\n");
        //gStyle->SetErrorX(0);
        //
        //
        char cname[100];
        sprintf(cname, "clone_%s", h_data->GetName());
        TH1F *h_data_pois = (TH1F *) h_data->Clone(cname);
        h_data_pois->Reset();
        for (int b =0;  b <= h_data->GetXaxis()->GetNbins() + 1; b++){
            h_data_pois->SetBinContent(b, h_data->GetBinContent(b));
        }
        h_data_pois->SetBinErrorOption(::TH1::kPoisson);
        h_data_pois->Draw("PEX0 same");
        h_data = h_data_pois;
        //gStyle->SetErrorX(1);

    }
    else{
        printf("not const size\n");
        h_data->Draw("P E0 same");
    }

    leg->Draw();


    TLatex latext; 
    latext.SetNDC();
    latext.SetTextColor(kBlack);
    latext.SetTextAlign(22); //center
    latext.SetTextFont(42);
    latext.SetTextSize(0.05);    


    float H = pad1->GetWh();
    float W = pad1->GetWw();
    float l = pad1->GetLeftMargin();
    float t = pad1->GetTopMargin();
    float r = pad1->GetRightMargin();
    float b = pad1->GetBottomMargin();

    float w = 1-l-r;


    latext.DrawLatex(l + w/2, 1-2*t ,plot_label);

    float y_title_offset = 1.5;


    float TS = 0.1;
    float TOffset = 0.8;
    float LS = 0.06;

    float rTS = TS * 0.7/0.3;
    float rLS = LS * 0.7/0.3;
    float rTOffset = TOffset * 0.3 / 0.7 - 0.05;

    h_stack->GetYaxis()->SetTitle(ylabel);
    h_stack->GetYaxis()->SetNdivisions(305);
    h_stack->GetYaxis()->SetTitleSize(TS);
    h_stack->GetYaxis()->SetTitleOffset(TOffset);
    h_stack->GetYaxis()->SetLabelSize(LS);

    h_stack->GetXaxis()->SetNdivisions(808);


    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0,.98,0.3);
    //pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.03);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    if(logx) pad2->SetLogx();



    ratio_unc->SetFillColor(kGray);
    //ratio_unc->SetFillStyle(3010);
    //ratio_unc->Print("all");




    float ratio_min = center - ratio_range;
    float ratio_max = center + ratio_range;
    h_ratio->GetYaxis()->SetRangeUser(ratio_min, ratio_max);

    h_ratio->SetTitle("");
    h_ratio->GetYaxis()->SetTitle("Obs/exp");
    h_ratio->GetYaxis()->SetNdivisions(205);
    h_ratio->GetYaxis()->SetTitleSize(rTS);
    h_ratio->GetYaxis()->SetLabelSize(rLS);
    h_ratio->GetYaxis()->SetTitleOffset(rTOffset);

    h_ratio->GetXaxis()->SetNdivisions(808);
    h_ratio->GetXaxis()->SetTickLength(0.06);
    h_ratio->GetXaxis()->SetTitle(xlabel);
    h_ratio->GetXaxis()->SetTitleSize(rTS);
    h_ratio->GetXaxis()->SetTitleOffset(rTOffset + 0.45);
    h_ratio->GetXaxis()->SetLabelSize(rLS);




    //h_ratio->SetMarkerStyle(21);
    if(const_size) h_ratio->Draw("PE0X0 same");
    else h_ratio->Draw("PE0 same");
    if(draw_sys_unc){
        ratio_unc->Draw("3 same");
        float line_start = h_ratio->GetXaxis()->GetXmin();
        float line_stop = h_ratio->GetXaxis()->GetXmax();
        TLine *l1 = new TLine(line_start,1,line_stop,1);
        l1->SetLineStyle(7);
        l1->SetLineWidth(2);
        l1->Draw();
    }
    if(const_size) h_ratio->Draw("EPX0 same");
    else h_ratio->Draw("EP same");
    gPad->RedrawAxis();

    float chi2 = computeChi2(h_ratio);
    int n_bins = h_ratio->GetNbinsX();


    if(draw_chi2){
        char chi2_label[100];
        sprintf(chi2_label, "#chi^{2} = %.1f", chi2);

        TLatex latext; 
        latext.SetNDC();
        latext.SetTextColor(kBlack);
        latext.SetTextAlign(32); //right
        latext.SetTextFont(42);
        latext.SetTextSize(0.1);    


        float H = pad2->GetWh();
        float W = pad2->GetWw();
        float l = pad2->GetLeftMargin();
        float t = pad2->GetTopMargin();
        float r = pad2->GetRightMargin();
        float b = pad2->GetBottomMargin();

        pad2->cd();
        latext.DrawLatex(1-1.2*r, 1-2.2*t ,chi2_label);
    }

    printf("Made ratio plot for label %s chi2/dof = %.1f/%i \n", label.Data(), chi2, n_bins);
    return std::make_pair(c, pad1);
}

#endif
