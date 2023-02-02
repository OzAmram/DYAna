#ifndef DYAFB_SFS
#define DYAFB_SFS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "TFile.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "HistUtils.C"
#include "bins.h"





typedef struct {
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
    TH2D *ISO_SF;
    TH2D *ID_SF;
    TH2D *ISO_SF_SYS;
    TH2D *ID_SF_SYS;
} mu_SFs;


typedef struct {
    TH1D *data_pileup;
    TH1D *pileup_ratio;
} pileup_SFs;

typedef struct {
    TH1D *data_pileup_up;
    TH1D *data_pileup_down;
    TH1D *data_pileup_nom;

    TH1D *ratio_pileup_up;
    TH1D *ratio_pileup_down;
    TH1D *ratio_pileup_nom;
} pileup_systematics;

typedef struct{
    TH2D *ID_SF;
    TH2D *RECO_SF;
    TH2D *HLT_SF;
    TH2D *HLT_MC_EFF;
} el_SFs;


typedef struct{
  TFile *f_;
  TF1* parametrization0p0To0p2_;
  TF1* parametrization0p2To0p3_;
  TF1* parametrization0p3To0p55_;
  TF1* parametrization0p55To0p83_;
  TF1* parametrization0p83To1p24_;
  TF1* parametrization1p24To1p4_;
  TF1* parametrization1p4To1p6_;
  TF1* parametrization1p6To1p8_;
  TF1* parametrization1p8To2p1_;
  TF1* parametrization2p1To2p25_;
  TF1* parametrization2p25To2p4_;
  TF1* parametrizationHotSpot_;
} MuPrefireHelper;


typedef struct{
    TH2F *jet_rate;
    TH2F *el_rate;
    MuPrefireHelper mu_era1_helper;
} prefire_SFs;

typedef struct{
    TH1F *ratios[n_m_bins][n_pt_bins][n_rap_bins];
    RooRealVar *A0_fits[n_m_bins][n_pt_bins][n_rap_bins];
} A0_helpers;

typedef struct{
    TH1F *el_rw[n_m_bins];
    TH1F *mu_rw[n_m_bins];
    TH1F *comb_rw[n_m_bins];
    TH1F *el_data_sub[n_m_bins];
    TH1F *mu_data_sub[n_m_bins];
    TH1F *comb_data_sub[n_m_bins];
} ptrw_helper;

typedef struct{
    TH1F *el_rw;
    TH1F *mu_rw;
} fakes_costrw_helper;

typedef struct{
    TH1F *rw[n_emu_rw_m_bins];
} emu_costrw_helper;

typedef struct{
    TH3D *h;
    TH3D *h_up;
    TH3D *h_down;
} LQ_rw_helper;

typedef struct{
    TProfile *h_F_up;
    TProfile *h_R_up;
    TProfile *h_RF_up;
    TProfile *h_F_down;
    TProfile *h_R_down;
    TProfile *h_RF_down;
    TProfile *h_pdfs[60];
} RF_pdf_norm_helper;



    


double get_var(Float_t vals[100]){
    float mean(0.), var(0.);
    int n_vars = 100;
    int n_entries = n_vars;
    for(int i =0; i< n_vars; i++){
        //printf("%.2f \n", vals[i]);
        if(std::isnan((float)vals[i])) n_entries--;
        else{
            //printf("val %.2f \n", vals[i]);
            mean += vals[i];
        }
        //printf("%.3e \n", vals[i]);
    }
    mean = mean / n_entries;
    //printf("mean %.3f n_entries %i\n", mean, n_entries);

    for(int  i=0; i< n_vars; i++){
        if(std::isnan((float)vals[i])) continue;
        else var += pow(vals[i] - mean, 2);
    }
    var = var/(n_entries -1);
    //printf("std %.3f \n\n\n", sqrt(var));
    return var;
}

Float_t get_pileup_SF(Int_t n_int, TH1D *h){
    TAxis* x_ax =  h->GetXaxis();
    Double_t max = x_ax->GetXmax();
    if(n_int > max) return 1;

    int xbin = x_ax->FindBin(n_int);

    Float_t result = h->GetBinContent(xbin);
    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}


float get_LQ_reweighting_denom(LQ_rw_helper h_LQ, int FLAG1, int FLAG2, float m, float rap, float cost){
    TH3D *h_rw;
    rap = abs(rap);
    if(FLAG2 == 0){ // everything
        h_rw = h_LQ.h;
    }
    else if(FLAG2 == 1){ // down quarks
        h_rw = h_LQ.h_down;
    }
    else if(FLAG2 == 2){ // up quarks
        h_rw = h_LQ.h_up;
    }
    else{
        printf("Invalid LQ reweighting flag %i! Options are 0, 1 or 2 \n", FLAG2);
        exit(1);
    }

    int xbin = h_rw->GetXaxis()->FindBin(m);
    int ybin = h_rw->GetYaxis()->FindBin(abs(rap));
    int zbin = h_rw->GetZaxis()->FindBin(cost);
    float weight = h_rw->GetBinContent(xbin, ybin, zbin);
   // if(weight < 1e-8){
      //  printf("m %.2f rap %.2f cost %.2f, xbin %i ybin %i zbin %i,  weight %f \n", m, rap, cost, xbin, ybin, zbin, weight);
        //weight = 1e-6;
   // }
    return weight;

}


float get_reweighting_denom(A0_helpers h, float cost, float m, float pt, float rap, int systematic = 0){
    if(m <= m_bins[0]) m = m_bins[0] + 0.1;
    int m_bin = find_bin(m_bins, m);
    int pt_bin = find_bin(pt_bins, pt);
    int rap_bin = find_bin(rap_bins, rap);
    float A0_ = h.A0_fits[m_bin][pt_bin][rap_bin]->getValV();
    if(systematic !=0){
        float err = h.A0_fits[m_bin][pt_bin][rap_bin]->getError();
        A0_ += systematic * err;
    }
    TH1F *h_correction = h.ratios[m_bin][pt_bin][rap_bin];
    TAxis* x_ax =  h_correction->GetXaxis();
    int bin = x_ax->FindBin(cost);
    float correction = h_correction->GetBinContent(bin);
    
    float denom = (3./8.*(1.+cost*cost + 0.5 * A0_ * (1. - 3. *cost*cost)));
    if(denom < 0  || isnan(denom)){
        printf("Denom %.3f, cost,m,pt,rap: %.2f, %.0f, %.0f %.1f \n", denom, cost, m, pt, rap);
        A0_ = 0.05;
        denom = (3./8.*(1.+cost*cost + 0.5 * A0_ * (1. - 3. *cost*cost)));
    }
    return denom;
}


float get_emu_costrw_SF(emu_costrw_helper h, float cost, float m, int systematic = 0){
    if(m < emu_rw_m_bins[0]) m = emu_rw_m_bins[0] + 1.;
    if(m > emu_rw_m_bins[n_emu_rw_m_bins - 1]) m = emu_rw_m_bins[n_emu_rw_m_bins] - 1.;
    int m_bin = find_bin(emu_rw_m_bins, m);
    TH1F *h_rw;
    if(m_bin < n_emu_rw_m_bins && m_bin >= 0) h_rw = h.rw[m_bin];
    else{
        printf("EMu mbin lookup error for rw. mbin %i, m %.1f \n", m_bin, m);
        exit(1);
    }

    int bin = h_rw->FindBin(cost);
    float correction = h_rw->GetBinContent(bin);
    //one systematic for every |cos(theta)| bin (nbins / 2)
    if(systematic != 0){
        float stat_err = h_rw->GetBinError(bin);
        float error = stat_err;

        int sys_bin = abs(systematic);
        int opp_bin = (h_rw->GetNbinsX() + 1) -sys_bin;
        if(bin == sys_bin || bin == opp_bin){
            //shift the reweighting in this bin by the error
            if(systematic >0) correction += error;
            if(systematic <0) correction -= error;
        }

    }


    
    return correction;
}

float get_LQ_ptrw_SF(ptrw_helper h, float m, float pt, int systematic = 0, bool muons = true){
    TH1F *h_rw, *h_N;
    //only use single mass bin, use OPPOSITE channel for correction
    if(muons){
        h_rw = h.el_rw[5];
        h_N = h.el_data_sub[5];
    }

    else{
        h_rw = h.mu_rw[5];
        h_N = h.el_data_sub[5];
    }

    TAxis* x_ax =  h_rw->GetXaxis();
    int bin = h_rw->FindBin(pt);
    float correction = h_rw->GetBinContent(bin);
    //one systematic for each pt bin
    if(systematic != 0){
        int sys_bin = abs(systematic);
        float stat_err = h_rw->GetBinError(sys_bin);
        //Low stat bins should not go crazy
        stat_err = max(stat_err, 0.1f);
        //float sys_correction = h_rw->GetBinContent(sys_bin);
        //float sys_err = 0.2 * std::fabs( sys_correction - 1.);
        //float error = pow(stat_err * stat_err + sys_err * sys_err, 0.5);
        float error = stat_err;

        if(bin == abs(systematic)){
            //shift the reweighting in this bin by the error
            if(systematic >0) correction += error;
            if(systematic <0) correction -= error;
        }
        else{
            //shift other bins to account for changing normalization
            float delta_N;
            if(systematic > 0) delta_N = error * h_N->GetBinContent(sys_bin);
            if(systematic < 0) delta_N = -error * h_N->GetBinContent(sys_bin);
            float integral = h_N->Integral();
            float ratio = integral/(integral + delta_N);
            correction *= ratio;
        }


    }


    
    return correction;
}


float get_ptrw_SF(ptrw_helper h, float m, float pt, int systematic = 0){
    TH1F *h_rw, *h_N;
    if(m > m_bins[n_m_bins-2]) m = m_bins[n_m_bins-2] - 0.1;
    int m_bin = find_bin(m_bins, m);
    h_rw = h.comb_rw[m_bin];
    h_N = h.comb_data_sub[m_bin];
    /*
    if(flag == FLAG_MUONS){
        h_rw = h.mu_rw[m_bin];
        h_N = h.mu_data_sub[m_bin];
    }
    else{
        h_rw = h.el_rw[m_bin];
        h_N = h.el_data_sub[m_bin];
    }
    */
    TAxis* x_ax =  h_rw->GetXaxis();
    int bin = h_rw->FindBin(pt);
    float correction = h_rw->GetBinContent(bin);
    //one systematic for each pt bin
    if(systematic != 0){
        int sys_bin = abs(systematic);
        float stat_err = h_rw->GetBinError(sys_bin);
        //Low stat bins should not go crazy
        stat_err = max(stat_err, 0.1f);
        //float sys_correction = h_rw->GetBinContent(sys_bin);
        //float sys_err = 0.2 * std::fabs( sys_correction - 1.);
        //float error = pow(stat_err * stat_err + sys_err * sys_err, 0.5);
        float error = stat_err;

        if(bin == abs(systematic)){
            //shift the reweighting in this bin by the error
            if(systematic >0) correction += error;
            if(systematic <0) correction -= error;
        }
        else{
            //shift other bins to account for changing normalization
            float delta_N;
            if(systematic > 0) delta_N = error * h_N->GetBinContent(sys_bin);
            if(systematic < 0) delta_N = -error * h_N->GetBinContent(sys_bin);
            float integral = h_N->Integral();
            float ratio = integral/(integral + delta_N);
            correction *= ratio;
        }


    }


    
    return correction;
}

Float_t get_Mu_trk_SF(Float_t eta, TGraphAsymmErrors *h, int systematic = 0){
    eta = abs(eta);
    Float_t result = h->Eval(eta);
    
    if(systematic !=0){
        TAxis* x_ax =  h->GetXaxis();
        int xbin = x_ax->FindBin(eta);
        Float_t err;
        if(systematic >0) err = h->GetErrorYhigh(xbin);
        else if(systematic <0) err = h->GetErrorYlow(xbin);
        result += err*systematic;
    }

    //if(result < 0.0001) printf("0 pileup SF for %i vertices\n", n_int);
    return result;
}

void get_pdf_avg_std_dev(Float_t pdf_Weights[100], Float_t *pdf_avg, Float_t *pdf_std_dev){
    Float_t sum = 0;
    for (int i=0; i<100; i++){
        sum+= pdf_Weights[i];
    }
    *pdf_avg = sum/100.;
    Float_t var = 0;
    for (int i=0; i<100; i++){
        var += pow(*pdf_avg - pdf_Weights[i], 2);
    }
    *pdf_std_dev = sqrt(var/99.);
    return;
}

Float_t get_prefire_rate(float pt, float eta, TH2F *map, int systematic = 0){
    //based on
    //https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/PrefireCorr.py
    float min_pt = 20.;
    float max_pt = 499.;
    float min_eta = 2.0;
    float max_eta = 3.0;
    if((pt < min_pt) || (std::abs(eta) < min_eta) || (std::abs(eta) > max_eta)) return 0.;

    int bin = map->FindBin(eta, std::min(pt, max_pt));
    float prefire_prob = map->GetBinContent(bin);

    float stat_err = map->GetBinError(bin);
    float syst_err = 0.2 * prefire_prob;
    if(systematic == -1){
        prefire_prob = std::max(0., prefire_prob - std::pow(stat_err*stat_err + syst_err*syst_err, 0.5));
    }
    else if(systematic == 1){
        prefire_prob = std::min(1., prefire_prob + std::pow(stat_err*stat_err + syst_err*syst_err, 0.5) );
    }
    //printf("pt %.0f eta %.2f, prob %.3f \n", pt, eta, prefire_prob);
    return prefire_prob;
}

Float_t get_mu_prefire_rate(float pt, float eta, float phi, MuPrefireHelper helper, int year, int systematic = 0){
    //https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1PrefiringWeightProducer.cc#L373-L440
  float prefrate;
  float statuncty;
  if (year == 2016 && (eta > 1.24 && eta < 1.6) &&
      (phi > 2.44346 && phi < 2.79253)) {
    prefrate = helper.parametrizationHotSpot_->Eval(pt);
    statuncty = helper.parametrizationHotSpot_->GetParError(2);
  } else if (std::abs(eta) < 0.2) {
    prefrate = helper.parametrization0p0To0p2_->Eval(pt);
    statuncty = helper.parametrization0p0To0p2_->GetParError(2);
  } else if (std::abs(eta) < 0.3) {
    prefrate = helper.parametrization0p2To0p3_->Eval(pt);
    statuncty = helper.parametrization0p2To0p3_->GetParError(2);
  } else if (std::abs(eta) < 0.55) {
    prefrate = helper.parametrization0p3To0p55_->Eval(pt);
    statuncty = helper.parametrization0p3To0p55_->GetParError(2);
  } else if (std::abs(eta) < 0.83) {
    prefrate = helper.parametrization0p55To0p83_->Eval(pt);
    statuncty = helper.parametrization0p55To0p83_->GetParError(2);
  } else if (std::abs(eta) < 1.24) {
    prefrate = helper.parametrization0p83To1p24_->Eval(pt);
    statuncty = helper.parametrization0p83To1p24_->GetParError(2);
  } else if (std::abs(eta) < 1.4) {
    prefrate = helper.parametrization1p24To1p4_->Eval(pt);
    statuncty = helper.parametrization1p24To1p4_->GetParError(2);
  } else if (std::abs(eta) < 1.6) {
    prefrate = helper.parametrization1p4To1p6_->Eval(pt);
    statuncty = helper.parametrization1p4To1p6_->GetParError(2);
  } else if (std::abs(eta) < 1.8) {
    prefrate = helper.parametrization1p6To1p8_->Eval(pt);
    statuncty = helper.parametrization1p6To1p8_->GetParError(2);
  } else if (std::abs(eta) < 2.1) {
    prefrate = helper.parametrization1p8To2p1_->Eval(pt);
    statuncty = helper.parametrization1p8To2p1_->GetParError(2);
  } else if (std::abs(eta) < 2.25) {
    prefrate = helper.parametrization2p1To2p25_->Eval(pt);
    statuncty = helper.parametrization2p1To2p25_->GetParError(2);
  } else if (std::abs(eta) < 2.4) {
    prefrate = helper.parametrization2p25To2p4_->Eval(pt);
    statuncty = helper.parametrization2p25To2p4_->GetParError(2);
  } else {
      std::cout << "Muon outside of |eta| <= 2.4. Prefiring weight set to 0." << std::endl;
    return 0.;
  }
  
  float prefiringRateSysUnc = 0.11;
  float systuncty = prefiringRateSysUnc * prefrate;

  if (systematic == 1)
    prefrate = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
  else if (systematic == -1)
    prefrate = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));

  if (prefrate > 1.) {
      std::cout << "Found a prefiring probability > 1. Setting to 1." << std::endl;
    return 1.;
  }
  return prefrate;
}







Float_t get_mu_SF(Float_t pt, Float_t eta, int year, TH2D *h, int systematic_barrel = 0, int systematic_endcap = 0){
    if(year <2016 || year > 2018) printf("Year is not from 2016-2018. This is bad!! \n");
    //stay in range of histogram
    float sys_unc_mult = 1.0;

    //not used
    int pt_systematic_sep = 0;

    TAxis *a_pt;
    if(year == 2016){
        a_pt=h->GetYaxis();
    }
    else{
        a_pt=h->GetXaxis();
    }


    float min_pt = a_pt->GetBinLowEdge(1);
    float max_pt = a_pt->GetBinUpEdge(a_pt->GetNbins());
    float mid_pt = 60.;


    if( pt <= min_pt) pt = min_pt + 1.;
    if (pt >= max_pt){
        pt = max_pt - 1.;
        sys_unc_mult = 1.5;
    }

    TAxis* x_ax =  h->GetXaxis();
    TAxis *y_ax =  h->GetYaxis();

    int xbin(0),ybin(0);
    if(year == 2016){
        //2016 hists are eta and pt
        xbin = x_ax->FindBin(eta);
        ybin = y_ax->FindBin(pt);
    }
    else{
        //2017 & 18 hists are pt and abs(eta)
        xbin = x_ax->FindBin(pt);
        ybin = y_ax->FindBin(abs(eta));
    }
    bool correct_pt_range = true;
    if(pt_systematic_sep != 0){
        correct_pt_range = (pt >= mid_pt && pt_systematic_sep > 0) || (pt < mid_pt && pt_systematic_sep < 0);
        //correct_pt_range = (ybin > mid_bin && pt_systematic_sep > 0) || (ybin <= mid_bin && pt_systematic_sep < 0);
    }


    Float_t result = h->GetBinContent(xbin, ybin);
    int systematic = systematic_barrel;
    if(abs(eta) > 1.4) systematic = systematic_endcap;
    if(systematic != 0 && correct_pt_range){
        Float_t err = sys_unc_mult * h->GetBinError(xbin, ybin);
        err =abs(err);
        result += (systematic * err);
    }
    if(result < 0.001){ 
        printf("0 muon SF for Pt %.1f, Eta %1.2f \n", pt, eta);
        result = 1;
    }
    return result;
}

Float_t get_el_SF(Float_t pt, Float_t eta, TH2D *h, int systematic_barrel = 0, int systematic_endcap = 0, int pt_systematic_sep = 0){
    float sys_unc_mult = 1.0;
    //get SF for eta's in overlap region (already vetoed superclusters in the
    //region)
    if(!goodElEta(eta)) eta = 1.6;

    int n_bins_pt = h->GetNbinsY();

    float min_pt = h->GetYaxis()->GetBinLowEdge(1);
    float max_pt = h->GetYaxis()->GetBinUpEdge(n_bins_pt);

    float mid_pt = 100.;
    int mid_bin = h->GetYaxis()->FindBin(mid_pt);

    if( pt <= min_pt) pt = min_pt + 1.;
    if (pt >= max_pt){
        pt = max_pt - 1.;
        sys_unc_mult = 1.5;
    }
    TAxis* x_ax =  h->GetXaxis();

    TAxis *y_ax =  h->GetYaxis();
    int xbin = x_ax->FindBin(eta);
    int ybin = y_ax->FindBin(pt);
    Float_t result = h->GetBinContent(xbin, ybin);
    bool correct_pt_range = true;
    if(pt_systematic_sep != 0){
        correct_pt_range = (pt >= mid_pt && pt_systematic_sep > 0) || (pt < mid_pt && pt_systematic_sep < 0);
        //correct_pt_range = (ybin > mid_bin && pt_systematic_sep > 0) || (ybin <= mid_bin && pt_systematic_sep < 0);
    }
    int systematic = systematic_barrel;
    if(abs(eta) > 1.4) systematic = systematic_endcap;
    if(systematic != 0 && correct_pt_range){
        Float_t err = sys_unc_mult * h->GetBinError(xbin, ybin);
        err =abs(err);
        if(err > 0.9 * result) err = 0.;
        result += (systematic * err);
    }
    //printf("eta %.2f sys_b %i sys_e %i, sys %i \n", eta, systematic_barrel, systematic_endcap, systematic);
    

    if(result < 0.01 || result > 10.){
        printf("%.2f el SF for Pt %.1f, Eta %1.2f \n", result, pt, eta);
        result = 1;
    }
    return result;
}

Float_t get_nnpdf_normfix(float gen_m){

    const float m_bins_mc[] = {100, 200,  400, 500, 700, 800, 1000, 1500, 2000, 3000}; //match mbins of MC samples
    const float norms[] = { 0.950, 0.954, 0.977, 0.984, 0.993, 1.000, 1.014, 1.039, 1.074};
    //float norms[] = { 0.949, 0.951, 0.953, 0.956, 0.966, 0.984, 0.996, 1.019, };

    if(gen_m<=m_bins_mc[0]) return norms[0];
    if(gen_m>=m_bins_mc[9]) return norms[8];

    int bin = find_bin(m_bins_mc, gen_m);
    float res = norms[bin];
    if(isnan(res) || res > 2. || res < 0.5 ){
        printf("m, i, norm: %.1f %i %.3f \n", gen_m, bin, res);
        res = 1.;
    }

    return res;
}


Float_t get_HLT_SF_1mu(Float_t mu1_pt, Float_t mu1_eta, TH2D *h_SF){
    //get HLT SF for event with just 1 muon
    //stay in range of histogram
    if (mu1_pt >= 350.) mu1_pt = 350.;
    mu1_eta = abs(mu1_eta);
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();
    int xbin1_SF = x_ax_SF->FindBin(std::abs(mu1_eta));
    int ybin1_SF = y_ax_SF->FindBin(mu1_pt);


    Float_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);

    Float_t result = SF1;
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan HLT SF 1 mu for Pt %.1f, Eta %1.2f \n", mu1_pt, mu1_eta);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}

Float_t get_HLT_SF(Float_t lep1_pt, Float_t lep1_eta, Float_t lep2_pt, Float_t lep2_eta, TH2D *h_SF, TH2D *h_MC_EFF, 
       int systematic_barrel = 0, int systematic_endcap = 0, int pt_systematic_sep = 0){
    float sys1_unc_mult = 1.0;
    float sys2_unc_mult = 1.0;
    //printf("Getting HLT for %.2f %.2f %.2f %.2f \n", lep1_pt, lep1_eta, lep2_pt, lep2_eta);
    //Get HLT SF for event with 2 leptons
    //stay in range of histogram
    TAxis *x_ax_SF =  h_SF->GetXaxis();
    TAxis *y_ax_SF =  h_SF->GetYaxis();

    int pt_max_bin = y_ax_SF->GetLast();
    double pt_max = y_ax_SF->GetBinUpEdge(pt_max_bin);
    double pt_overflow = y_ax_SF->GetBinCenter(pt_max_bin);

    float mid_pt = 75.;
    int mid_bin = y_ax_SF->FindBin(mid_pt);

    if (lep1_pt >= pt_max){
        sys1_unc_mult = 1.5;
        lep1_pt = pt_overflow;
    }
    if (lep2_pt >= pt_max){
        sys2_unc_mult = 1.5;
        lep2_pt = pt_overflow;
    }
    lep1_eta = abs(lep1_eta);
    lep2_eta = abs(lep2_eta);
    int xbin1_SF = x_ax_SF->FindBin(std::fabs(lep1_eta));
    int ybin1_SF = y_ax_SF->FindBin(lep1_pt);

    int xbin2_SF = x_ax_SF->FindBin(std::fabs(lep2_eta));
    int ybin2_SF = y_ax_SF->FindBin(lep2_pt);

    Float_t SF1 = h_SF->GetBinContent(xbin1_SF, ybin1_SF);
    Float_t SF2 = h_SF->GetBinContent(xbin2_SF, ybin2_SF);

    bool correct_pt_range1  = true;
    bool correct_pt_range2  = true;
    if(pt_systematic_sep != 0){
        correct_pt_range1 = (lep1_pt >= mid_pt && pt_systematic_sep > 0) || (lep1_pt < mid_pt && pt_systematic_sep < 0);
        correct_pt_range1 = (lep2_pt >= mid_pt && pt_systematic_sep > 0) || (lep2_pt < mid_pt && pt_systematic_sep < 0);
    }

    int systematic1 = systematic_barrel;
    int systematic2 = systematic_barrel;
    if(abs(lep1_eta) > 1.4) systematic1 = systematic_endcap;
    if(abs(lep2_eta) > 1.4) systematic2 = systematic_endcap;

    if(systematic1 != 0 && correct_pt_range1){
        Float_t SF1_err = h_SF->GetBinError(xbin1_SF, ybin1_SF);
        SF1 += sys1_unc_mult * SF1_err * systematic1;
    }
    if(systematic2 != 0 && correct_pt_range2){
        Float_t SF2_err = h_SF->GetBinError(xbin2_SF, ybin2_SF);
        SF2 += sys2_unc_mult * SF2_err * systematic2;
    }
    //printf("eta1 %.2f eta2 %.2f sys_b %i sys_e %i sys1 %i sys2 %i \n", lep1_eta, lep2_eta, systematic_barrel, systematic_endcap, systematic1, systematic2);



    TAxis *x_ax_MC_EFF =  h_MC_EFF->GetXaxis();
    TAxis *y_ax_MC_EFF =  h_MC_EFF->GetYaxis();
    int xbin1_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(lep1_eta));
    int ybin1_MC_EFF = y_ax_MC_EFF->FindBin(lep1_pt);

    int xbin2_MC_EFF = x_ax_MC_EFF->FindBin(std::fabs(lep2_eta));
    int ybin2_MC_EFF = y_ax_MC_EFF->FindBin(lep2_pt);

    Float_t MC_EFF1 = h_MC_EFF->GetBinContent(xbin1_MC_EFF, ybin1_MC_EFF);
    Float_t MC_EFF2 = h_MC_EFF->GetBinContent(xbin2_MC_EFF, ybin2_MC_EFF);


    Float_t result = (1 - (1-MC_EFF1*SF1)*(1-MC_EFF2*SF2))/
                      (1 - (1-MC_EFF1)*(1-MC_EFF2));
    //printf("%.1f %.1f SF: %.3f %.3f, MC EFF: %.3f %.3f  comb %.3f \n", lep1_pt, lep2_pt, SF1, SF2, MC_EFF1, MC_EFF2, result);
    if(result < 0.01) printf("0 HLT SF for Pt %.1f, Eta %1.2f \n", lep1_pt, lep1_eta);
    if(TMath::IsNaN(result)){ 
        printf("Nan lep HLT SF for Pt1 %.2f, Eta1 %.2f Pt2 %.2f Eta2 %.2f  %.2f %.2f \n", lep1_pt, lep1_eta, lep2_pt, lep2_eta, SF1, SF2);
        result = 1;
    }
    //printf("Result, SF1 = (%0.3f, %0.3f) \n", result, SF1);
    return result;
}



void setup_A0_helper(A0_helpers *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/a0_fits.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/a0_fits.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/a0_fits.root");
    for (int i=0; i< n_m_bins; i++){
        for (int j=0; j< n_pt_bins; j++){
            for (int k=0; k< n_pt_bins; k++){

                char title[100];
                sprintf(title, "amc_fit_ratio_y%i_m%i_pt%i_rap%i", year -2000, i, j, k);
                h->ratios[i][j][k] = (TH1F *) f->Get(title)->Clone();
                h->ratios[i][j][k]->SetDirectory(0);
                sprintf(title, "a0_y%i_m%i_pt%i_rap%i", year -2000, i, j, k);
                h->A0_fits[i][j][k] = (RooRealVar *) f->Get(title)->Clone();
            }
        }
    }
}



void setup_RF_pdf_norm_helper(RF_pdf_norm_helper *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/RF_pdf_weights.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/RF_pdf_weights.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/RF_pdf_weights.root");
    for (int i=0; i< n_m_bins; i++){
        for (int j=0; j< n_pt_bins; j++){

            h->h_R_up = (TProfile *) f->Get("h_R_up")->Clone();
            h->h_R_up->SetDirectory(0);
            h->h_F_up = (TProfile *) f->Get("h_F_up")->Clone();
            h->h_F_up->SetDirectory(0);
            h->h_RF_up = (TProfile *) f->Get("h_RF_up")->Clone();
            h->h_RF_up->SetDirectory(0);

            h->h_R_down = (TProfile *) f->Get("h_R_down")->Clone();
            h->h_R_down->SetDirectory(0);
            h->h_F_down = (TProfile *) f->Get("h_F_down")->Clone();
            h->h_F_down->SetDirectory(0);
            h->h_RF_down = (TProfile *) f->Get("h_RF_down")->Clone();
            h->h_RF_down->SetDirectory(0);

            char title[100];
            for(int k = 0; k<60; k++){
                sprintf(title, "h_pdf%i", k);
                h->h_pdfs[k] = (TProfile *) f->Get(title);
                h->h_pdfs[k]->SetDirectory(0);
            }


        }
    }
}


void setup_LQ_rw_helper(LQ_rw_helper *h_lq, int year){
    TFile *f;

    if(year == 2016)      f = TFile::Open("../analyze/SFs/2016/LQ_rw.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/LQ_rw.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/LQ_rw.root");



    h_lq->h = (TH3D *) f->Get("h")->Clone();
    h_lq->h->SetDirectory(0);


    h_lq->h_up = (TH3D *) f->Get("h_up")->Clone();
    h_lq->h_up->SetDirectory(0);


    h_lq->h_down = (TH3D *) f->Get("h_down")->Clone();
    h_lq->h_down->SetDirectory(0);

    f->Close();
}

void setup_LQ_ptrw_helper(ptrw_helper *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/pt_rw_LQ.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/pt_rw_LQ.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/pt_rw_LQ.root");
    int i = 5;
    char h_name[100];
    sprintf(h_name, "elel%i_m%i_pt_ratio", year % 2000, i);
    h->el_rw[i] = (TH1F *) f->Get(h_name)->Clone();
    h->el_rw[i]->SetDirectory(0);
    sprintf(h_name, "mumu%i_m%i_pt_ratio", year % 2000, i);
    h->mu_rw[i] = (TH1F *) f->Get(h_name)->Clone();
    h->mu_rw[i]->SetDirectory(0);
    sprintf(h_name, "comb%i_m%i_pt_ratio", year % 2000, i);
    h->comb_rw[i] = (TH1F *) f->Get(h_name)->Clone();
    h->comb_rw[i]->SetDirectory(0);

    sprintf(h_name, "elel%i_m%i_pt_data_sub", year % 2000, i);
    h->el_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
    h->el_data_sub[i]->SetDirectory(0);
    sprintf(h_name, "mumu%i_m%i_pt_data_sub", year % 2000, i);
    h->mu_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
    h->mu_data_sub[i]->SetDirectory(0);
    sprintf(h_name, "comb%i_m%i_pt_data_sub", year % 2000, i);
    h->comb_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
    h->comb_data_sub[i]->SetDirectory(0);
}




void setup_ptrw_helper(ptrw_helper *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/pt_rw.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/pt_rw.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/pt_rw.root");
    for (int i=0; i< n_m_bins -1; i++){
        char h_name[100];
        sprintf(h_name, "elel%i_m%i_pt_ratio", year % 2000, i);
        h->el_rw[i] = (TH1F *) f->Get(h_name)->Clone();
        h->el_rw[i]->SetDirectory(0);
        sprintf(h_name, "mumu%i_m%i_pt_ratio", year % 2000, i);
        h->mu_rw[i] = (TH1F *) f->Get(h_name)->Clone();
        h->mu_rw[i]->SetDirectory(0);
        sprintf(h_name, "comb%i_m%i_pt_ratio", year % 2000, i);
        h->comb_rw[i] = (TH1F *) f->Get(h_name)->Clone();
        h->comb_rw[i]->SetDirectory(0);

        sprintf(h_name, "elel%i_m%i_pt_data_sub", year % 2000, i);
        h->el_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
        h->el_data_sub[i]->SetDirectory(0);
        sprintf(h_name, "mumu%i_m%i_pt_data_sub", year % 2000, i);
        h->mu_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
        h->mu_data_sub[i]->SetDirectory(0);
        sprintf(h_name, "comb%i_m%i_pt_data_sub", year % 2000, i);
        h->comb_data_sub[i] = (TH1F *) f->Get(h_name)->Clone();
        h->comb_data_sub[i]->SetDirectory(0);
    }
}

void setup_emu_costrw_helper(emu_costrw_helper *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/emu_cost_rw.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/emu_cost_rw.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/emu_cost_rw.root");
    else{
        printf("Year is %i ?? ", year);
        exit(1);
    }
    char name[100];
    for(int i=0; i<n_emu_rw_m_bins; i++){
        sprintf(name, "emu%i_mbin0_cost_ratio", year % 2000);
        h->rw[i] = (TH1F *) f->Get(name)->Clone();
        h->rw[i]->SetDirectory(0);
    }

}


void setup_fakes_costrw_helper(fakes_costrw_helper *h, int year){
    TFile *f;

    if(year == 2016) f = TFile::Open("../analyze/SFs/2016/fakes_cost_rw.root");
    else if(year == 2017) f = TFile::Open("../analyze/SFs/2017/fakes_cost_rw.root");
    else if(year == 2018) f = TFile::Open("../analyze/SFs/2018/fakes_cost_rw.root");
    char el_name[100], mu_name[100];
    sprintf(el_name, "elel%i_ss_cost_ratio", year % 2000);
    sprintf(mu_name, "mumu%i_ss_cost_ratio", year % 2000);
    h->el_rw = (TH1F *) f->Get(el_name)->Clone();
    h->el_rw->SetDirectory(0);

    h->mu_rw = (TH1F *) f->Get(mu_name)->Clone();
    h->mu_rw->SetDirectory(0);
}

void setup_mu_prefire_helper(TFile *f_, TString era_name_, MuPrefireHelper &era_helper){
    //https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1PrefiringWeightProducer.cc#L373-L440
    era_helper.f_ = f_;
    TString paramName = "L1prefiring_muonparam_HotSpot_" + era_name_;
    era_helper.parametrizationHotSpot_ = (TF1 * ) f_->Get(paramName);

    paramName = "L1prefiring_muonparam_0.0To0.2_" + era_name_;
    era_helper.parametrization0p0To0p2_ = (TF1 *) f_->Get(paramName);
    paramName = "L1prefiring_muonparam_0.2To0.3_" + era_name_;
    era_helper.parametrization0p2To0p3_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_0.3To0.55_" + era_name_;
    era_helper.parametrization0p3To0p55_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_0.55To0.83_" + era_name_;
    era_helper.parametrization0p55To0p83_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_0.83To1.24_" + era_name_;
    era_helper.parametrization0p83To1p24_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_1.24To1.4_" + era_name_;
    era_helper.parametrization1p24To1p4_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_1.4To1.6_" + era_name_;
    era_helper.parametrization1p4To1p6_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_1.6To1.8_" + era_name_;
    era_helper.parametrization1p6To1p8_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_1.8To2.1_" + era_name_;
    era_helper.parametrization1p8To2p1_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_2.1To2.25_" + era_name_;
    era_helper.parametrization2p1To2p25_ = (TF1 *)f_->Get(paramName);
    paramName = "L1prefiring_muonparam_2.25To2.4_" + era_name_;
    era_helper.parametrization2p25To2p4_ = (TF1 *)f_->Get(paramName);
 

    if(era_helper.parametrizationHotSpot_ == nullptr)
        era_helper.parametrizationHotSpot_ = era_helper.parametrization1p4To1p6_ ;


    //printf("Finish setting up era %s \n", era_name_.Data());



}



void setup_prefire_SFs(prefire_SFs *pre_SF,  int year){


    TFile *f_mu_prefire = TFile::Open("../analyze/SFs/2016/L1MuonPrefiringParametriations.root");


    if(year == 2016){

        TString era1_name = "2016";

        setup_mu_prefire_helper(f_mu_prefire, era1_name, pre_SF->mu_era1_helper);
        pre_SF->mu_era1_helper.f_->Print();

        TFile *f_el = TFile::Open("../analyze/SFs/2016/L1prefiring_photonpt_2016BtoH.root");
        TH2F *h1 = (TH2F *) f_el->Get("L1prefiring_photonpt_2016BtoH")->Clone();
        h1->SetDirectory(0);
        pre_SF->el_rate= h1;
        f_el->Close();

        TFile *f_jet = TFile::Open("../analyze/SFs/2016/L1prefiring_jetpt_2016BtoH.root");
        TH2F *h2 = (TH2F *) f_jet->Get("L1prefiring_jetpt_2016BtoH")->Clone();
        h2->SetDirectory(0);
        pre_SF->jet_rate= h2;
        f_jet->Close();
    }
    else if(year == 2017 || year == 2018){

        TString era1_name = "20172018";

        setup_mu_prefire_helper(f_mu_prefire, era1_name, pre_SF->mu_era1_helper);

        TFile * f_el = TFile::Open("../analyze/SFs/2017/L1prefiring_photonpt_2017BtoF.root");
        TH2F *h1 = (TH2F *) f_el->Get("L1prefiring_photonpt_2017BtoF")->Clone();
        h1->SetDirectory(0);
        pre_SF->el_rate= h1;
        f_el->Close();

        TFile *f_jet = TFile::Open("../analyze/SFs/2017/L1prefiring_jetpt_2017BtoF.root");
        TH2F *h2 = (TH2F *) f_jet->Get("L1prefiring_jetpt_2017BtoF")->Clone();
        h2->SetDirectory(0);
        pre_SF->jet_rate= h2;
        f_jet->Close();
    }
    //f_mu_prefire->Close();
    f_mu_prefire->Print();
    pre_SF->mu_era1_helper.parametrization0p0To0p2_->Print();

}





void setup_mu_SFs(mu_SFs *era1, mu_SFs *era2, int year){
    TH1::AddDirectory(kFALSE);

    TFile *f1, *f2, *f3, *f4, *f5, *f6;


    if(year == 2016){
        f1 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_HLT.root");
        f1->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("../analyze/SFs/2016/Mu_GH_HLT.root");
        f4->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF ->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF ->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_1->SetDirectory(0);
        era1->ID_SF = ID_1;
        TH2D *ID_1_SYS = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_eta_pt_syst")->Clone();
        ID_1_SYS->SetDirectory(0);
        era1->ID_SF_SYS = ID_1_SYS;
        f2->Close();

        f3 = TFile::Open("../analyze/SFs/2016/Mu_BCDEF_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_1->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        TH2D *ISO_1_SYS = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_syst")->Clone();
        ISO_1_SYS->SetDirectory(0);
        era1->ISO_SF_SYS = ISO_1_SYS;
        f3->Close();


        f5 = TFile::Open("../analyze/SFs/2016/Mu_GH_ID.root");
        TH2D *ID_2 = (TH2D *) f5->Get("NUM_TightID_DEN_genTracks_eta_pt")->Clone();
        ID_2->SetDirectory(0);
        era2->ID_SF = ID_2;
        TH2D *ID_2_SYS = (TH2D *) f5->Get("NUM_TightID_DEN_genTracks_eta_pt_syst")->Clone();
        ID_2_SYS->SetDirectory(0);
        era2->ID_SF_SYS = ID_2_SYS;
        f5->Close();

        f6 = TFile::Open("../analyze/SFs/2016/Mu_GH_ISO.root");
        TH2D *ISO_2 = (TH2D *) f6->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")->Clone();
        ISO_2->SetDirectory(0);
        era2->ISO_SF = ISO_2;
        TH2D *ISO_2_SYS = (TH2D *) f6->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_syst")->Clone();
        ISO_2_SYS->SetDirectory(0);
        era2->ISO_SF_SYS = ISO_2_SYS;
        f6->Close();

    }
    else if(year == 2017){
        printf("hlt\n");
        f1 = TFile::Open("../analyze/SFs/2017/Mu_HLT.root");
        f1->cd("IsoMu27_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();
        era2->HLT_SF = era1->HLT_SF;
        era2->HLT_MC_EFF = era1->HLT_MC_EFF;

        printf("id\n");
        f2 = TFile::Open("../analyze/SFs/2017/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        TH2D *ID_1_SYS = (TH2D *) f2->Get("NUM_TightID_DEN_genTracks_pt_abseta_syst")->Clone();
        ID_1_SYS->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        era1->ID_SF_SYS = ID_1_SYS;
        era2->ID_SF_SYS = ID_1_SYS;
        f2->Close();

        printf("iso\n");
        f3 = TFile::Open("../analyze/SFs/2017/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        TH2D *ISO_1_SYS = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst")->Clone();
        ISO_1_SYS->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        era1->ISO_SF_SYS = ISO_1_SYS;
        era2->ISO_SF_SYS = ISO_1_SYS;
        f3->Close();


    }
    else if(year == 2018){
        f1 = TFile::Open("../analyze/SFs/2018/Mu_per1_HLT.root");
        f1->cd("IsoMu24_PtEtaBins");
        era1->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era1->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era1->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era1->HLT_MC_EFF->SetDirectory(0);
        f1->Close();

        f4 = TFile::Open("../analyze/SFs/2018/Mu_per2_HLT.root");
        f4->cd("IsoMu24_PtEtaBins");
        era2->HLT_SF = (TH2D *) gDirectory->Get("abseta_pt_ratio")->Clone();
        era2->HLT_SF->SetDirectory(0);
        gDirectory->cd("efficienciesMC");
        era2->HLT_MC_EFF = (TH2D *) gDirectory->Get("abseta_pt_MC")->Clone();
        era2->HLT_MC_EFF->SetDirectory(0);
        f4->Close();

        f2 = TFile::Open("../analyze/SFs/2018/Mu_ID.root");
        TH2D *ID_1 = (TH2D *) f2->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta")->Clone();
        ID_1->SetDirectory(0);
        TH2D *ID_1_SYS = (TH2D *) f2->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta_syst")->Clone();
        ID_1_SYS->SetDirectory(0);
        era1->ID_SF = ID_1;
        era2->ID_SF = ID_1;
        era1->ID_SF_SYS = ID_1_SYS;
        era2->ID_SF_SYS = ID_1_SYS;
        f2->Close();

        f3 = TFile::Open("../analyze/SFs/2018/Mu_ISO.root");
        TH2D *ISO_1 = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone();
        ISO_1->SetDirectory(0);
        TH2D *ISO_1_SYS = (TH2D *) f3->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst")->Clone();
        ISO_1_SYS->SetDirectory(0);
        era1->ISO_SF = ISO_1;
        era2->ISO_SF = ISO_1;
        era1->ISO_SF_SYS = ISO_1_SYS;
        era2->ISO_SF_SYS = ISO_1_SYS;
        f3->Close();


    }

    if(era1->HLT_SF == NULL || era1->ID_SF == NULL || era1->ISO_SF == NULL  || era1->ID_SF_SYS == NULL || era1->ISO_SF_SYS == NULL ||
       era2->HLT_SF == NULL || era2->ID_SF == NULL || era2->ISO_SF == NULL  || era2->ID_SF_SYS == NULL || era2->ISO_SF_SYS == NULL)
        printf("Something wrong setup muon SF'S !!!! \n\n\n");
}

void setup_el_SF(el_SFs *sf, int year){
    //Setup electron SF's
    TFile *f_id, *f_reco, *f_hlt, *f_hlt2;
    if(year == 2016){
        f_hlt = TFile::Open("../analyze/SFs/2016/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2016/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2016/El_RECO.root");
    }
    else if(year == 2017){
        f_hlt = TFile::Open("../analyze/SFs/2017/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2017/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2017/El_RECO.root");
    }
    else if(year == 2018){
        f_hlt = TFile::Open("../analyze/SFs/2018/El_HLT.root");
        f_id = TFile::Open("../analyze/SFs/2018/El_ID.root");
        f_reco = TFile::Open("../analyze/SFs/2018/El_RECO.root");
    }

    TH2D *h1 = (TH2D *) f_id->Get("EGamma_SF2D")->Clone();
    h1->SetDirectory(0);
    sf->ID_SF = h1;
    f_id->Close();

    TH2D *h2 = (TH2D *) f_reco->Get("EGamma_SF2D")->Clone();
    h2->SetDirectory(0);
    sf->RECO_SF = h2;
    f_reco->Close();

    TH2D *h_hltsf1 = (TH2D *) f_hlt->Get("EGamma_SF2D")->Clone();
    h_hltsf1->SetDirectory(0);
    sf->HLT_SF = h_hltsf1;

    TH2D *h_hlt_mceff = (TH2D *) f_hlt->Get("EGamma_MCEff")->Clone();
    h_hlt_mceff->SetDirectory(0);
    sf->HLT_MC_EFF = h_hlt_mceff;

    
    if(sf->HLT_SF == NULL || sf->HLT_MC_EFF == NULL || sf->ID_SF == NULL || sf->RECO_SF == NULL ) printf("Something wrong setup electron SF'S !!!! \n\n\n");
}


void setup_pileup_systematic(pileup_systematics *pu_sys, int year){

    TFile *fin;
    if(year == 2016){
        fin = TFile::Open("../analyze/SFs/2016/pu_SF.root");
    }
    if(year == 2017){
        fin = TFile::Open("../analyze/SFs/2017/pu_SF.root");
    }
    if(year == 2018){
        fin = TFile::Open("../analyze/SFs/2018/pu_SF.root");
    }


    pu_sys->ratio_pileup_nom = (TH1D *) fin->Get("pu_ratio")->Clone();
    pu_sys->ratio_pileup_nom->SetDirectory(0);

    pu_sys->ratio_pileup_up = (TH1D *) fin->Get("pu_ratio_up")->Clone();
    pu_sys->ratio_pileup_up->SetDirectory(0);

    pu_sys->ratio_pileup_down = (TH1D *) fin->Get("pu_ratio_down")->Clone();
    pu_sys->ratio_pileup_down->SetDirectory(0);
}





void setup_pileup_systematic_old(pileup_systematics *pu_sys, int year){


    TFile *f7, *f8, *f9;
    if(year == 2016){
        f7 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2016/Data16PileupHistogram_72383.root");
    }
    if(year == 2017){
        f7 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2017/Data17PileupHistogram_72383.root");
    }
    if(year == 2018){
        f7 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_69200.root");
        f8 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_66017.root");
        f9 = TFile::Open("../analyze/SFs/2018/Data18PileupHistogram_72383.root");
    }


    TH1D * pileup_data_nom = (TH1D *) f7->Get("pileup")->Clone();
    pileup_data_nom->Scale(1./pileup_data_nom->Integral());
    pileup_data_nom->SetDirectory(0);

    TH1D * pileup_data_up = (TH1D *) f8->Get("pileup")->Clone();
    pileup_data_up->Scale(1./pileup_data_up->Integral());
    pileup_data_up->SetDirectory(0);

    TH1D * pileup_data_down = (TH1D *) f9->Get("pileup")->Clone();
    pileup_data_down->Scale(1./pileup_data_down->Integral());
    pileup_data_down->SetDirectory(0);

    //Printf(" Ints are %.4e %.4e %.4e \n", pileup_data_nom->Integral(), pileup_data_up->Integral(), pileup_data_down->Integral());
    pu_sys->data_pileup_up = (TH1D *) pileup_data_up->Clone("data_pu_up");
    pu_sys->data_pileup_down = (TH1D *) pileup_data_down->Clone("data_pu_down");
    pu_sys->data_pileup_nom = (TH1D *) pileup_data_nom->Clone("data_pu_nom");
    pu_sys->data_pileup_up->SetDirectory(0);
    pu_sys->data_pileup_down->SetDirectory(0);
    pu_sys->data_pileup_nom->SetDirectory(0);

    pu_sys->ratio_pileup_up = (TH1D *) pileup_data_up->Clone("ratio_pu_up");
    pu_sys->ratio_pileup_down = (TH1D *) pileup_data_down->Clone("ratio_pu_down");
    pu_sys->ratio_pileup_nom = (TH1D *) pileup_data_nom->Clone("ratio_pu_nom");
    pu_sys->ratio_pileup_up->SetDirectory(0);
    pu_sys->ratio_pileup_down->SetDirectory(0);
    pu_sys->ratio_pileup_nom->SetDirectory(0);
    f7->Close();
    f8->Close();
    f9->Close();
}
#endif
