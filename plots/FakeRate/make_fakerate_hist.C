#define STAND_ALONE
//perform fits to Reconstructed MuMu data to extract Asym

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
#include "../../utils/HistMaker.C"
//

void SetErrors(TH2D *h_rate, TH2D *h_total){
    int nBins_x = h_rate->GetXaxis()->GetNbins();
    int nBins_y = h_rate->GetYaxis()->GetNbins();
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){
            //binomial distribution divided by n
            Double_t n = h_total->GetBinContent(i,j);
            Double_t r = h_rate->GetBinContent(i,j);
            Double_t r_prev = h_rate->GetBinContent(i, j-1);
            Double_t err = sqrt(r*(1-r)/n);
            h_rate->SetBinError(i,j, err);
        }
    }
    return;
}

void SetErrors(TH1D *h_rate, TH1D *h_total){
    int nBins = h_rate->GetXaxis()->GetNbins();
    for (int i=0; i <= nBins; i++){
        //binomial distribution divided by n
        Double_t n = h_total->GetBinContent(i);
        Double_t p = h_rate->GetBinContent(i);
        Double_t err = sqrt(p*(1-p)/n);
        //printf("Err is %.2f \n", err);
        if(p == 0) err = 0.2;
        h_rate->SetBinError(i, err);
    }
    return;
}

void cleanup(TH2D *h_rate){
    int nBins_x = h_rate->GetXaxis()->GetNbins();
    int nBins_y = h_rate->GetYaxis()->GetNbins();
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){
            Double_t r = h_rate->GetBinContent(i,j);
            Double_t r_prev = h_rate->GetBinContent(i, j-1);
            if(r <=0.005 && r_prev > 2. * r){
                if (j != 0){
                    printf("Substituting 0 rate with 1 lower pt bin \n");
                    r = h_rate->GetBinContent(i, j-1);
                    Double_t prev_err = h_rate->GetBinError(i, j-1);
                    Double_t err = 1.4 * prev_err;
                    h_rate->SetBinContent(i,j,r);
                    h_rate->SetBinError(i,j,err);
                    printf(" i j prev error %i %i %.2f \n", i,j, prev_err );
                }

            }

        }
    }
}

void SetCorrectedRate(TH2D *h_data_rate, TH2D *h_data_total, TH2D *h_contam_rate, TH2D* h_contam_total){
    int nBins_x = h_data_rate->GetXaxis()->GetNbins();
    int nBins_y = h_data_rate->GetYaxis()->GetNbins();
    printf("nx,ny = %i %i \n", nBins_x, nBins_y);
    //printf("Get size %i \n", nBins);
    for (int i=1; i <= nBins_x; i++){
        for (int j=1; j <= nBins_y; j++){
            Double_t r_old = h_data_rate->GetBinContent(i,j);
            Double_t n_old = h_data_total->GetBinContent(i,j);
            Double_t p_old = n_old * r_old;
            Double_t p_old_err = sqrt(p_old);
            if(r_old == 0 || n_old == 0) continue;

            Double_t r_contam = h_contam_rate->GetBinContent(i,j);
            Double_t r_contam_err = h_contam_rate->GetBinError(i,j);
            Double_t n_contam = h_contam_total->GetBinContent(i,j);
            Double_t p_contam = r_contam * n_contam;
            Double_t p_contam_err = r_contam_err * n_contam;

            Double_t eta_center = h_data_total->GetXaxis()->GetBinCenter(i);
            Double_t pt_center = h_data_total->GetYaxis()->GetBinCenter(j);

            Double_t p_new = p_old - p_contam;
            Double_t p_new_err = sqrt(p_old_err*p_old_err + p_contam_err*p_contam_err);
            if(p_new <= 0) p_new = 0;
            Double_t n_new = n_old - n_contam;
            Double_t r_new = p_new/n_new;
            Double_t r_new_err = p_new_err / n_new;
            
            printf("Bin (%.3f, %.0f): (r, n,p) Old (%.2f, %.0f, %.0f) New (%.2f, %.0f, %.0f) \n", eta_center, pt_center, r_old, n_old, p_old, r_new, n_new, p_new);
            printf("p_old, p_old_err, p_contam, p_contam_err, p_new, p_new_err: %.2f %.2f %.2f  %.2f %.2f %.2f  \n", p_old, p_old_err, p_contam, p_contam_err, p_new, p_new_err);

            h_data_rate->SetBinContent(i,j, r_new);
            h_data_rate->SetBinError(i,j,r_new_err);
            h_data_total->SetBinContent(i,j, n_new);
        }
    }
        

}

void construct_fakerate_template(TH2D *h_rate, TH2D *h_total, TTree *t, int flag = FLAG_MUONS, bool isData = true, int year = 2016, bool mt_cut = false, float z_mass_window = -1., float met_cut = 100.){
    Double_t mt, pt, eta, gen_weight, z_mass, met, m01, m02, m12;
    Bool_t pass;
    Double_t lumi = 0;
    t->SetBranchAddress("mt", &mt);
    t->SetBranchAddress("met", &met);
    t->SetBranchAddress("z_mass", &z_mass);
    t->SetBranchAddress("m01", &m01);
    t->SetBranchAddress("m02", &m02);
    t->SetBranchAddress("m12", &m12);
    z_mass = 0.;
    mt = 0.;
    if(flag == FLAG_MUONS){
        t->SetBranchAddress("mu_pt", &pt);
        t->SetBranchAddress("mu_eta", &eta);
        if(year == 2016) lumi = mu_lumi16;
        else if(year == 2017) lumi = mu_lumi17;
        else if(year == 2018) lumi = mu_lumi18;
        else printf("Year is %i \n", year);
    }
    else{
        t->SetBranchAddress("el_pt", &pt);
        t->SetBranchAddress("el_eta", &eta);
        if(year == 2016) lumi = el_lumi16;
        else if(year == 2017) lumi = el_lumi17;
        else if(year == 2018) lumi = el_lumi18;
        else printf("Year is %i \n", year);
    }
    t->SetBranchAddress("pass", &pass);
    if(!isData) t->SetBranchAddress("gen_weight", &gen_weight);
    Long64_t size  =  t->GetEntries();
    for (int i=0; i<size; i++) {
        t->GetEntry(i);
        Double_t fill_weight = 1;
        if(!isData) fill_weight = gen_weight;
        bool cut = (z_mass_window < 0. || fabs(91.2 - z_mass) < z_mass_window) && (!mt_cut || (mt < 30.)) && (met < met_cut);
        cut = cut && m01 < 500. && m02 < 500. && m12 < 500.;
        if(cut){
            h_total->Fill(eta,pt, fill_weight);
            if(pass) h_rate->Fill(eta,pt, fill_weight);
            }
    }
    h_rate->Divide(h_total);
    if(!isData) h_total->Scale(1000*lumi);
    return;
}
    
        
    


void make_fakerate_hist(){

    bool write_out = true;
    int year = 2016;
    //int FLAG = FLAG_ELECTRONS;
    int FLAG = FLAG_MUONS;
    int n_pt_bins = 4;
    TFile *f, *f_mc, *f_new; 
    Float_t *pt_bins;
    Float_t eta_bins[] = {0,1.479,2.4};
    int n_eta_bins = 2;
    bool mt_cut = false;
    float met_cut = 100.;

    float z_mass_window = 15.;

    char f_plot[100], f_ewk_plot[100];

    if(FLAG == FLAG_ELECTRONS){
        Float_t pt_bins_[] = {10,20,30,50, 100, 1000};
        pt_bins = pt_bins_;
        n_pt_bins = 5;
        sprintf(f_plot, "Paper_plots/El%i_fakerate.pdf", year % 2000);
        sprintf(f_ewk_plot, "Paper_plots/El%i_ewk_contam_fakerate.pdf", year % 2000);
        if(year == 2018){
            f = TFile::Open("../analyze/output_files/FakeRate18_el_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate18_el_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2018/SingleElectron18_data_fakerate.root", "RECREATE");
        }
        else if(year == 2017){
            f = TFile::Open("../analyze/output_files/FakeRate17_el_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate17_el_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2017/SingleElectron17_data_fakerate.root", "RECREATE");
        }
        else if(year == 2016){
            f = TFile::Open("../analyze/output_files/FakeRate16_el_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate16_el_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2016/SingleElectron16_data_fakerate.root", "RECREATE");
        }
        else{
            printf("Somethign went wrong.. \n");
            exit(0);
        }
    }

    else{
        Float_t pt_bins_[] = {10,20,30,50, 100, 1000};
        pt_bins = pt_bins_;
        n_pt_bins = 5;
        sprintf(f_plot, "Paper_plots/Mu%i_fakerate.pdf", year % 2000);
        sprintf(f_ewk_plot, "Paper_plots/Mu%i_ewk_contam_fakerate.pdf", year % 2000);
        if(year == 2018){
            f = TFile::Open("../analyze/output_files/FakeRate18_mu_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate18_mu_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2018/SingleMuon18_data_fakerate.root", "RECREATE");
        }
        else if(year == 2017){
            f = TFile::Open("../analyze/output_files/FakeRate17_mu_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate17_mu_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2017/SingleMuon17_data_fakerate.root", "RECREATE");
        }
        else if(year == 2016){
            f = TFile::Open("../analyze/output_files/FakeRate16_mu_data_jan23_2022.root");
            f_mc = TFile::Open("../analyze/output_files/FakeRate16_mu_contam_jan23_2022.root");
            if(write_out) f_new = TFile::Open("../analyze/FakeRate/root_files/2016/SingleMuon16_data_fakerate.root", "RECREATE");
        }
        else{
            printf("Somethign went wrong.. \n");
            exit(0);
        }
    }



    TH2D *h_rate = new TH2D("h_rate", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_total = new TH2D("h_total", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TH2D *h_contam_rate = new TH2D("h_contam_rate", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);
    TH2D *h_contam_total = new TH2D("h_contam_total", "passing iso cuts", n_eta_bins, eta_bins, n_pt_bins, pt_bins);

    TTree *t_data = (TTree *) f->Get("T_data");
    TTree *t_mc = (TTree *) f_mc->Get("T_data");



    construct_fakerate_template(h_rate, h_total, t_data, FLAG, true, year, mt_cut, z_mass_window, met_cut);
    construct_fakerate_template(h_contam_rate, h_contam_total, t_mc, FLAG,  false, year, mt_cut, z_mass_window, met_cut);
    h_contam_rate->Print("range");
    h_contam_total->Print("range");
    //h_contam_total->Scale(0.5);

    //TH2D* h_rate = (TH2D *)f->Get("h_rate"); 
    //TH2D* h_total = (TH2D *)f->Get("h_total"); 

    TH2D* h_rate_new = (TH2D *)h_rate->Clone("h_rate_new");
    TH2D* h_total_new = (TH2D *)h_total->Clone("h_total_new");

    TH2D* h_ewk_frac = (TH2D *) h_contam_total->Clone("h_ewk_frac");
    h_ewk_frac->Divide(h_total);


    /*
    TH2D* h_contam_pass = (TH2D *)f_mc->Get("h_pass"); 
    TH2D* h_contam_total = (TH2D *)f_mc->Get("h_total"); 
    h_contam_pass->Scale(1000*tot_lumi);
    h_contam_total->Scale(1000*tot_lumi);
    TH2D* h_contam_rate = h_contam_pass->Clone("h_conam_rate");
    h_contam_rate->Divide(h_contam_total);
    */

    SetErrors(h_rate_new, h_total_new);

    SetCorrectedRate(h_rate_new, h_total_new, h_contam_rate, h_contam_total);
    printf("Before cleanup\n");
    h_rate_new->Print("range");
    cleanup(h_rate_new);
    printf("After\n");
    h_rate_new->Print("range");
    if(write_out){
        f_new->cd();
        h_rate_new->Write();
        h_total_new->Write();
        printf("Wrote out corrected rate \n");
    }

    TH1D *rate_barrel = h_rate_new->ProjectionY("rate_barrel", 1,1, "e");
    TH1D *rate_endcap = h_rate_new->ProjectionY("rate_endcap", 2,2, "e");

    TH1D *total_barrel = h_total_new->ProjectionY("total_bar", 1,1, "e");
    TH1D *total_endcap = h_total_new->ProjectionY("total_endcap", 2,2, "e");


    

    //SetErrors(rate_barrel, total_barrel);
    //SetErrors(rate_endcap, total_endcap);

    char title[100], title2[100]; 
    if(FLAG == FLAG_MUONS) sprintf(title, "Muons Fakerate: %i ", year);
    else sprintf(title, "Electrons Fakerate: %i ", year);


    if(FLAG == FLAG_MUONS) sprintf(title2, "Muons EWK Contam: %i ", year);
    else sprintf(title2, "Electrons EWK Contam: %i ", year);

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 900, 700);
    c1->SetLogx();




    rate_barrel->SetTitle(title);
    rate_barrel->SetStats(0);
    rate_barrel->SetLineWidth(3);
    rate_barrel ->Draw("E1");
    rate_barrel->GetXaxis()->SetTitle("p_t (Gev)");
    rate_barrel->GetXaxis()->SetRangeUser(10,500);
    if(FLAG == FLAG_MUONS) rate_barrel->GetYaxis()->SetRangeUser(0, 0.3);
    else rate_barrel->GetYaxis()->SetRangeUser(0, 0.9);
    c1->Update();

    //TCanvas *c2 = new TCanvas("c2", "Histograms", 200, 10, 900, 700);
    //rate_endcap->SetTitle("Fake-Rate Muons with Trigger: End Cap");
    rate_endcap->SetStats(0);
    rate_endcap->SetLineWidth(3);
    rate_endcap->Draw("E1 same");
    rate_endcap->SetLineColor(kRed);
    //rate_endcap->GetXaxis()->SetRangeUser(10,100);
    rate_endcap->GetYaxis()->SetRangeUser(0, 0.5);


    TLegend *leg1 = new TLegend(0.25,0.15);
    leg1->AddEntry(rate_barrel, "Barrel", "f");
    leg1->AddEntry(rate_endcap, "Endcap",  "f");
    leg1->Draw();

    if(write_out){
        c1->Print(f_plot);
    }

    TH1D *h_ewk_frac_barrel = h_ewk_frac->ProjectionY("ewk_frac_barrel", 1,1, "e");
    TH1D *h_ewk_frac_endcap = h_ewk_frac->ProjectionY("ewk_frac_endcap", 2,2, "e");


    float ymax = std::min(1.0, h_ewk_frac_barrel->GetMaximum() * 2.0);
    TCanvas *c2 = new TCanvas("c2", "", 900, 700);
    c2->SetLogx();
    h_ewk_frac_barrel->SetTitle(title2);
    h_ewk_frac_barrel->SetStats(0);
    h_ewk_frac_barrel->SetLineWidth(3);
    h_ewk_frac_barrel->SetLineColor(kGreen + 3);
    h_ewk_frac_barrel->Draw("E1");
    h_ewk_frac_barrel->GetYaxis()->SetRangeUser(0., 1.0);
    h_ewk_frac_barrel->GetXaxis()->SetRangeUser(10., 500.);
    h_ewk_frac_barrel->GetYaxis()->SetTitle("Frac. of Ewk. Contam.");
    h_ewk_frac_barrel->GetXaxis()->SetTitle("p_t (Gev)");


    h_ewk_frac_endcap->SetLineColor(kMagenta+ 3);
    h_ewk_frac_endcap->Draw("E1 same");
    h_ewk_frac_endcap->SetLineWidth(3);

    TLegend *leg2 = new TLegend(0.25,0.15);
    leg2->AddEntry(h_ewk_frac_barrel, "Barrel", "f");
    leg2->AddEntry(h_ewk_frac_endcap, "Endcap",  "f");
    leg2->Draw();

    if(write_out){ c2->Print(f_ewk_plot);}

        
    



    /*
    TCanvas *c3 = new TCanvas("c3", "Histograms", 200, 10, 900, 700);
    total_barrel->SetTitle("Number of events");
    total_barrel->SetStats(0);
    total_barrel->SetLineWidth(3);
    total_barrel ->Draw("hist");
    total_barrel->GetXaxis()->SetTitle("p_t (Gev)");
    total_barrel->GetXaxis()->SetRangeUser(10,100);
    c3->SetLogy();
    c3->SetLogx();

    total_endcap->SetTitle("Number of events");
    total_endcap->SetStats(0);
    total_endcap->SetLineWidth(3);
    total_endcap->SetLineColor(kRed);
    total_endcap->Draw("hist same");
    total_endcap->GetXaxis()->SetTitle("p_t (Gev)");
    total_endcap->GetXaxis()->SetRangeUser(10,100);
    printf("Totals: Bar %.0f End %.0f \n", total_barrel->Integral(), total_endcap->Integral());
    TLegend *leg2 = new TLegend(0.5, 0.65, 0.75, 0.8);
    leg2->AddEntry(total_barrel, "Barrel",  "f");
    leg2->AddEntry(total_endcap, "Endcap" , "f");
    leg2->Draw();
    */
}

    
    
