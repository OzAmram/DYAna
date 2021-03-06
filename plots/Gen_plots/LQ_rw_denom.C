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
#include "../../utils/root_files.h"
#include "../../utils/HistUtils.C"
#include "../../utils/PlotUtils.C"

int make_amc_gen_cost(TTree *t_gen, TH3D *h_3d, TH3D *h_3d_up, TH3D *h_3d_down, float m_low, float m_high){
    TLorentzVector *gen_lep_p(0), *gen_lep_m(0), cm;
    float gen_weight, m, cost, cost_st;
    int inc_id1, inc_id2;
    Bool_t sig_event(1);
    t_gen->SetBranchAddress("gen_p", &gen_lep_p);
    t_gen->SetBranchAddress("gen_m", &gen_lep_m);
    t_gen->SetBranchAddress("inc_id1", &inc_id1);
    t_gen->SetBranchAddress("inc_id2", &inc_id2);
    //t_gen->SetBranchAddress("gen_mu_p", &gen_lep_p);
    //t_gen->SetBranchAddress("gen_mu_m", &gen_lep_m);
    t_gen->SetBranchAddress("m", &m);
    t_gen->SetBranchAddress("cost", &cost);
    t_gen->SetBranchAddress("cost_st", &cost_st);
    t_gen->SetBranchAddress("gen_weight", &gen_weight);
    t_gen->SetBranchAddress("sig_event", &sig_event);





    int nEvents=0;

    for (int i=0; i<t_gen->GetEntries(); i++) {
        t_gen->GetEntry(i);
        //bool pass = abs(gen_lep_p->Eta()) < 2.4 && abs(gen_lep_m->Eta()) < 2.4 && max(gen_lep_m->Pt(), gen_lep_p->Pt()) > 30.;
        if(m >= m_low && m <= m_high && sig_event){
            nEvents++;
            cm = *gen_lep_p + *gen_lep_m;
            float pt = cm.Pt();
            float rap = abs(cm.Rapidity());
            /*
            float my_cost = get_cost(*gen_lep_p, *gen_lep_m);
            if(cost_st > 0) my_cost = abs(my_cost);
            else my_cost = -abs(my_cost);
            */
            float my_cost = cost_st;
            
            //fill in units of pb
            h_3d->Fill(m, rap,  my_cost, 0.5*gen_weight);
            h_3d->Fill(m, rap, -my_cost, 0.5*gen_weight);

            if(abs(inc_id1) == 1 && abs(inc_id2) == 1 && inc_id1 * inc_id2 < 0){ //d dbar
                h_3d_down->Fill(m, rap,  my_cost, 0.5*gen_weight);
                h_3d_down->Fill(m, rap, -my_cost, 0.5*gen_weight);
            }
            if(abs(inc_id1) == 2 && abs(inc_id2) == 2 && inc_id1 * inc_id2 < 0){ //u ubar
                h_3d_up->Fill(m, rap,  my_cost, 0.5* gen_weight);
                h_3d_up->Fill(m, rap, -my_cost, 0.5* gen_weight);
            }


        }
    }
    printf("selected %i events \n", nEvents);

    return nEvents;

}

void normalize(TH3D *h){
    for(int i=1; i<=h->GetNbinsX(); i++){
        for(int j=1; j<=h->GetNbinsY(); j++){
            for(int k=1; k<=h->GetNbinsZ(); k++){
                float xw = h->GetXaxis()->GetBinWidth(i);
                float yw = h->GetYaxis()->GetBinWidth(j);
                float zw = h->GetYaxis()->GetBinWidth(k);
                float content = h->GetBinContent(i,j,k);
                float new_content = content/(xw*yw *zw);

                float err = h->GetBinError(i,j,k);
                float new_err = err/(xw*yw*zw);
                //printf("i,j xw, yw %i %i %.2f %.2f \n", i,j, xw,yw);

                h->SetBinContent(i,j,k, new_content);
                h->SetBinError(i,j,k, new_err);
                if(content < 1e-12 || new_content < 0. || content * 1.5 < err){
                    printf("WARNING bin %i %i (m = %.0f) has content %.3e +/- %.3e \n", i,j, h->GetXaxis()->GetBinCenter(i), content, err);
                }

            }
        }
    }
}




void LQ_rw_denom(){

    bool write_out = true;
    char *out_file = "../analyze/SFs/2016/LQ_rw.root";
    TFile *f_gen = TFile::Open("../analyze/output_files/DY16_gen_level_aug4.root");

    TFile * f_out;
    if(write_out)
        f_out = TFile::Open(out_file, "RECREATE");
    

    //TFile *f_gen = TFile::Open("../MuMu17_dy_gen.root");
    TTree *t_gen_mu = (TTree *) f_gen->Get("T_gen_mu");
    TTree *t_gen_el = (TTree *) f_gen->Get("T_gen_el");


    char plot_dir[] = "Misc_plots/A0_fits/";



    int n_rap_bins = 4;
    float rap_bins[] = {0., 0.6, 1.0, 1.5, 2.4};
    int n_cost_bins = 20;
    //float cost_bins[] = {-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.0};
    float cost_bins[] = {-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};
    //float cost_bins[] = {0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0};

    int n_LQ_pt_bins = 1;
    //float LQ_pt_bins[] = {0., 20., 60., 100., 10000};
    float LQ_pt_bins[] = {0., 10000};

    int n_LQ_m_bins = 28;
    float m_LQ_bins[] = {350., 375., 400., 425., 450., 475., 500., 525., 550., 575., 600., 625., 650., 700., 750., 800., 850., 900., 950., 1000.,  1100.,  
        1200.,  1300., 1400.,  1600.,  1800.,  2000.,  2500.,   3000., };

    TH3D *h = new TH3D("h", "", n_LQ_m_bins, m_LQ_bins, n_rap_bins, rap_bins,    n_cost_bins, cost_bins);

    TH3D *h_up = new TH3D("h_up", "", n_LQ_m_bins, m_LQ_bins, n_rap_bins, rap_bins,   n_cost_bins, cost_bins);
    TH3D *h_down = new TH3D("h_down", "", n_LQ_m_bins, m_LQ_bins,   n_rap_bins, rap_bins, n_cost_bins, cost_bins);

    TH1D *h_1d = new TH1D("h1", "", n_cost_bins, cost_bins);


    int nEvents = 0;
    float m_low = m_LQ_bins[0];
    float m_high = 14000.;

    make_amc_gen_cost(t_gen_mu,  h, h_up, h_down, m_low,m_high);
    make_amc_gen_cost(t_gen_el, h, h_up, h_down,  m_low, m_high);
    
    /*
    TH1D *h_el_cost_pre = h_el->ProjectionY("h_el_cost_pre", 10,10);
    TH1D *h_mu_cost_pre = h_mu->ProjectionY("h_mu_cost_pre", 10,10);
    make_ratio_plot("LQ_cos_theta_pre.png", h_el_cost_pre, "El",h_mu_cost_pre, "Mu", "El/Mu", "cos(#theta_{*})", false, true);
    */

    normalize(h);

    normalize(h_up);
    normalize(h_down);
    
    //h_el->Print("range");


    /*
    TH1D *h_el_cost = h_el->ProjectionY("h_el_cost", 10,10);
    TH1D *h_mu_cost = h_mu->ProjectionY("h_mu_cost", 10,10);
    make_ratio_plot("LQ_cos_theta.png", h_el_cost, "El",h_mu_cost, "Mu", "El/Mu", "cos(#theta_{*})", false, true);

    TH1D *h_el_m = h_el->ProjectionX("h_el_m");
    TH1D *h_mu_m = h_mu->ProjectionX("h_mu_m");
    make_ratio_plot("LQ_m.png", h_el_m, "El",h_mu_m, "Mu", "El/Mu", "M (GeV)", false, true);


    h_el_cost->SetLineColor(kBlue);
    h_mu_cost->SetLineColor(kRed);
    h_el_cost->Draw();
    h_mu_cost->Draw("same");
    */

    if(write_out){
        f_out->cd();
        h->Write();
        h_up->Write();
        h_down->Write();
        f_out->Print();
        f_out->Close();
    }

}
