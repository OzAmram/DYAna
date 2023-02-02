
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TH2D.h"
#include "Math/Functor.h"
#include "TempMaker.C"

using namespace std;



void cleanup_hist(TH1 *h){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        if(val<0.) h->SetBinContent(i,0.);
    }
}

void setHistError(TH1 *h, float e, bool add_err = true){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        float err = h->GetBinError(i);
        float new_err = val*e;
        if(add_err) new_err = pow(err*err + new_err*new_err, 0.5);
        h->SetBinError(i, new_err);
    }
}


void setHistMassDepError(TH1 *h, bool add_err = true){
    int nBins = h->GetNbinsX();
    for(int i=1; i<= nBins; i++){
        float val = h->GetBinContent(i);
        float m = h->GetXaxis()->GetBinCenter(i);

        float e = 0.;

        if(m<200) e = 0.1;
        else if(m > 200 && m < 320) e = 0.13; 
        else if(m > 320 && m < 510) e = 0.15; 
        else if(m > 510 && m < 700) e = 0.20; 
        else if(m > 700 && m < 1000) e = 0.20; 
        else if( m > 1000) e = 0.26; 

        float err = h->GetBinError(i);
        float new_err = val*e;
        if(add_err) new_err = pow(err*err + new_err*new_err, 0.5);
        h->SetBinError(i, new_err);
    }
}



void symmetrize1d(TH1F *h){
    int n_cost_bins = h->GetNbinsX();

        for(int j=1; j<= n_cost_bins/2; j++){
            float content = h->GetBinContent(j);
            float error = h->GetBinError(j);

            int opp_j = (n_cost_bins + 1) -j;
            float content2 = h->GetBinContent(opp_j);
            float error2 = h->GetBinError(opp_j);

            float new_content = (content + content2)/2.0;
            float new_error = pow((error*error + error2*error2), 0.5)/2.0;
            h->SetBinContent(j, new_content);
            h->SetBinContent(opp_j, new_content);

            h->SetBinError(j, new_error);
            h->SetBinError(opp_j, new_error);


        }
    
}




void make_emu_m_cost_pt_rap_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt,   TH1F *h_rap, bool is_data = false, 
        int year=2016, float m_low = 150., float m_high = 999999., bool ss = false, bool costrw = false){
    Long64_t size  =  t1->GetEntries();
    TempMaker tm(t1, is_data, year);
    tm.do_emu = true;
    if(!is_data) tm.do_emu_costrw = costrw;

    tm.setup();
    int nEvents=0;

    for (int i=0; i<tm.nEntries; i++) {
        tm.getEvent(i);
        bool pass  = tm.m >= m_low && tm.m <= m_high; //&& tm.met_pt < met_cut && tm.has_no_bjets;
        if(pass){

            nEvents++;
            tm.doCorrections();
            tm.getEvtWeight();

            h_m->Fill(tm.m, tm.evt_weight);
            if(ss) h_cost->Fill(abs(tm.cost), tm.evt_weight);
            else h_cost->Fill(tm.cost, tm.evt_weight);
            h_pt->Fill(tm.cm.Pt(), tm.evt_weight);
            h_rap->Fill(tm.cm.Rapidity(), tm.evt_weight);


        }
    }
    printf("Selected %i events \n", nEvents);
    tm.finish();
}





void make_m_cost_pt_xf_hist(TTree *t1, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf, TH1F *h_phi, TH1F *h_rap,
        bool is_data=false, int flag1 = FLAG_MUONS,
        int year = 2016, Double_t m_low = 150., Double_t m_high = 9999999., bool ss = false, bool ptrw = false){
        h_m->Sumw2();
        h_cost->Sumw2();
        h_pt->Sumw2();
        h_xf->Sumw2();
        h_phi->Sumw2();
        h_rap->Sumw2();
    //read event data
        TempMaker tm(t1, is_data, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;
        tm.do_ptrw = ptrw;


        tm.setup();
        int nEvents=0;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.not_cosmic;//&& tm.met_pt < met_cut  && tm.has_no_bjets ;
            //bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.has_no_bjets && tm.not_cosmic;
            //bool good_eta = (fabs(tm.el1_eta) > 1.5 && fabs(tm.el2_eta) < 1.4) || (fabs(tm.el2_eta) > 1.5 && fabs(tm.el1_eta) < 1.4);
            //pass = pass && good_eta;


            if(pass){
                nEvents++;
                tm.doCorrections();
                tm.getEvtWeight();
                //if(flag1 == FLAG_ELECTRONS && (tm.el2_pt < 20 && tm.el2_eta > 1.8)) tm.evt_weight = 0;

                h_m->Fill(tm.m,tm.evt_weight);
                h_pt->Fill(tm.pt,tm.evt_weight);
                h_xf->Fill(tm.xF, tm.evt_weight);
                if(h_phi != nullptr) h_phi->Fill(tm.cm.Phi(), tm.evt_weight);
                if(h_rap != nullptr) h_rap->Fill(tm.cm.Rapidity(), tm.evt_weight);
                //if(h_phi != nullptr) h_phi->Fill(tm.lep_p->Phi(), tm.evt_weight);

                if(ss){
                    h_cost->Fill(-tm.cost, 0.5* tm.evt_weight);
                    h_cost->Fill(tm.cost, 0.5* tm.evt_weight);
                }
                else h_cost->Fill(tm.cost, tm.evt_weight);
            }


        
        }
        printf("Selected %i events \n", nEvents);


    tm.finish();
}

void make_fakerate_est(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_contam, TTree *t_QCD_contam, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_xf, TH1F *h_phi, TH1F *h_rap,
        int flag1 = FLAG_MUONS, int year = 2016, float m_low=150., float m_high = 99999., bool ss = false, bool reweight = false, bool sys_errors = false){
    TLorentzVector *lep_p=0;
    TLorentzVector *lep_m=0;
    Double_t pt;
    FakeRate FR;
    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    h_cost->Sumw2();
    if(flag1 == FLAG_MUONS) setup_new_mu_fakerate(&FR, year);
    else setup_new_el_fakerate(&FR, year);
    TH2D *h_err = (TH2D *) FR.h->Clone("h_err");
    h_err->Reset();
    for (int l=0; l<=3; l++){
        TTree *t;
        bool is_data, is_one_iso;
        if (l==0){
            t = t_WJets;
            is_data =true;
            is_one_iso = true;
        }
        if (l==1){
            t = t_QCD;
            is_data =true;
            is_one_iso = false;
        }
        if (l==2){
            t = t_WJets_contam;
            is_data =false;
            is_one_iso = true;
        }
        if (l==3){
            t = t_QCD_contam;
            is_data =false;
            is_one_iso = false;
        }
        TempMaker tm(t, is_data, year);
        tm.is_one_iso = is_one_iso;
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;
        tm.do_RC = true;
        tm.setup();
        int n= 0;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);

            bool pass = (tm.m >= m_low && tm.m <= m_high); //&& tm.met_pt < met_cut  && tm.has_no_bjets;

            //bool good_eta = (fabs(tm.el1_eta) > 1.5 && fabs(tm.el2_eta) < 1.4) || (fabs(tm.el2_eta) > 1.5 && fabs(tm.el1_eta) < 1.4);
            //pass = pass && good_eta;

            bool opp_sign = false;
            if(flag1 == FLAG_MUONS) opp_sign = ((abs(tm.mu1_charge - tm.mu2_charge)) > 0.01);
            else opp_sign = ((abs(tm.el1_charge - tm.el2_charge)) > 0.01);
            if(!ss) pass = pass && opp_sign;
            //else if(flag1 == FLAG_ELECTRONS) pass = pass && !opp_sign;
            //else pass = pass && !opp_sign;
            if(pass){
                double evt_reweight = 0.;

                if(flag1 == FLAG_MUONS){
                    Double_t mu1_fakerate, mu2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);

                        evt_reweight = mu1_fakerate/(1-mu1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 1);
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), 1);
                    }

                    if(l==1){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));

                        h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f),  -0.5*(mu2_fakerate/(1-mu2_fakerate)));
                        h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f),  -0.5*(mu1_fakerate/(1-mu1_fakerate)));
                    }
                    if(l==2){
                        if(tm.iso_lep ==1) mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        if(tm.iso_lep ==0) mu1_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = -(mu1_fakerate )/(1-mu1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), -tm.getEvtWeight());
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), -tm.getEvtWeight());
                    }
                    if(l==3){
                        mu1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, FR.h);
                        mu2_fakerate = get_new_fakerate_prob(tm.mu2_pt, tm.mu2_eta, FR.h);
                        evt_reweight = (mu1_fakerate/(1-mu1_fakerate)) * (mu2_fakerate/(1-mu2_fakerate));
                        h_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 0.5*tm.getEvtWeight() * (mu2_fakerate/(1-mu2_fakerate)));
                        h_err->Fill(min(abs(tm.mu2_eta), 2.3f), min(tm.mu2_pt, 150.f), 0.5*tm.getEvtWeight() * (mu1_fakerate/(1-mu1_fakerate)));
                    }
                }
                else{
                    Double_t el1_fakerate, el2_fakerate; 
                    if(l==0){
                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = el1_fakerate/(1-el1_fakerate);

                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), 1);
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), 1);
                    }
                    if(l==1){
                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));

                        h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), -0.5*(el2_fakerate/(1-el2_fakerate)));
                        h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), -0.5*(el1_fakerate/(1-el1_fakerate)));

                    }
                    if(l==2){

                        if(tm.iso_lep ==1) el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        if(tm.iso_lep ==0) el1_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = -(el1_fakerate )/(1-el1_fakerate);


                        if(tm.iso_lep ==1) h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), -tm.getEvtWeight());
                        if(tm.iso_lep ==0) h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f), -tm.getEvtWeight());
                    }
                    if(l==3){

                        el1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, FR.h);
                        el2_fakerate = get_new_fakerate_prob(tm.el2_pt, tm.el2_eta, FR.h);
                        evt_reweight = (el1_fakerate/(1-el1_fakerate)) * (el2_fakerate/(1-el2_fakerate));

                        h_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f),  0.5*tm.getEvtWeight() * (el2_fakerate/(1-el2_fakerate)));
                        h_err->Fill(min(abs(tm.el2_eta), 2.3f), min(tm.el2_pt, 150.f),  0.5*tm.getEvtWeight() * (el1_fakerate/(1-el1_fakerate)));
                    }
                }
                double tot_evt_weight = 0.;
                if(is_data) tot_evt_weight = evt_reweight; 
                else{
                    tot_evt_weight = evt_reweight * tm.getEvtWeight();
                    //printf("%.3e %.3e \n", tm.getEvtWeight(), tot_evt_weight);
                }


                if(!ss){
                    h_cost->Fill(tm.cost, tot_evt_weight);
                    h_rap->Fill(tm.cm.Rapidity(), tot_evt_weight);
                }
                else{
                    h_cost->Fill(tm.cost, 0.5*tot_evt_weight);
                    h_cost->Fill(-tm.cost, 0.5*tot_evt_weight);
                    h_rap->Fill(tm.cm.Rapidity(), 0.5 * tot_evt_weight);
                    h_rap->Fill(-tm.cm.Rapidity(), 0.5 * tot_evt_weight);
                }
                
                h_m->Fill(tm.m, tot_evt_weight);
                h_pt->Fill(tm.cm.Pt(), tot_evt_weight);
                h_xf->Fill(tm.xF, tot_evt_weight);
                h_phi->Fill(tm.cm.Phi(), tot_evt_weight);
                //if(h_phi != nullptr) h_phi->Fill(tm.lep_p->Phi(), tot_evt_weight);
                if(opp_sign) tot_weight_os += tot_evt_weight;
                else tot_weight_ss += tot_evt_weight;


            }
        }
        printf("After iter %i current fakerate est is %.0f \n", l, h_cost->Integral());
        //printf("N is %i \n", n);

    }
    cleanup_hist(h_m);
    cleanup_hist(h_pt);
    cleanup_hist(h_xf);
    cleanup_hist(h_cost);
    cleanup_hist(h_rap);

    if(ss){
        float scaling = 1.0;
        scaling = tot_weight_ss / (tot_weight_ss + tot_weight_os);
        printf("Scaling is %.2f \n", scaling);
        h_m->Scale(scaling);
        h_pt->Scale(scaling);
        h_xf->Scale(scaling);
        h_cost->Scale(scaling);
        h_phi->Scale(scaling);
        h_rap->Scale(scaling);
        h_err->Scale(scaling);
    }


    if(sys_errors) set_fakerate_errors(h_err, FR.h, h_cost);
    if(reweight){
        fakes_costrw_helper h_rw;
        setup_fakes_costrw_helper(&h_rw, year);
        //h_cost->Print("range");
        if(flag1 == FLAG_MUONS) fakes_cost_reweight(h_cost, h_rw.mu_rw);
        else fakes_cost_reweight(h_cost, h_rw.el_rw);
        //h_cost->Print("range");
    }


    //if(ss && flag1 != FLAG_MUONS){
    Double_t err;
    printf("Total fakerate est is %.0f +/- %.0f \n", h_cost->IntegralAndError(1, h_cost->GetNbinsX(), err), err);
    return;
}


void Fakerate_est_emu(TTree *t_WJets, TTree *t_QCD, TTree *t_WJets_MC, TTree *t_QCD_MC, TH1F *h_m, TH1F *h_cost, TH1F *h_pt, TH1F *h_rap, int flag1 = FLAG_MUONS, 
        int year=2016, float m_low = 150., float m_high = 10000., bool reweight = false, bool sys_errors = false){
    FakeRate el_FR, mu_FR;
    //TH2D *FR;
    setup_new_el_fakerate(&el_FR, year);
    setup_new_mu_fakerate(&mu_FR, year);
    TH2D *h_mu_err = (TH2D *) mu_FR.h->Clone("h_err");
    h_mu_err->Reset();
    TH2D *h_el_err = (TH2D *) el_FR.h->Clone("h_err");
    h_el_err->Reset();
    //FR.h->Print();
    bool is_data = false;
    bool is_one_iso = false;

    float tot_weight_os = 0.;
    float tot_weight_ss = 0.;
    for (int l=0; l<=3; l++){
        TTree *t;
        if (l==0){
            t = t_WJets;
            is_data = true;
            is_one_iso = true;
        }
        if (l==1){
            t = t_QCD;
            is_data = true;
            is_one_iso = false;
        }
        if (l==2){
            t = t_WJets_MC;
            is_data = false;
            is_one_iso = true;
        }
        if (l==3){
            t = t_QCD_MC;
            is_data = false;
            is_one_iso = false;
        }
        TempMaker tm(t, is_data, year);
        tm.do_emu = true;
        tm.is_one_iso = is_one_iso;

        tm.setup();



        Long64_t size  =  t->GetEntries();
        int n=0;
        float avg_mc_weight = 0.;
        for (int i=0; i<size; i++) {
            tm.getEvent(i);

            bool opp_sign =  ((abs(tm.mu1_charge - tm.el1_charge)) > 0.01);
            bool pass = tm.m>= m_low && tm.m <= m_high && opp_sign &&//tm.met_pt < met_cut  && tm.has_no_bjets && opp_sign &&

                ((flag1 == FLAG_MUONS && tm.mu1_pt > 27.) || (flag1 == FLAG_ELECTRONS && tm.el1_pt > 29.));
            if(pass){
                double evt_reweight = 0.;
                float lep1_fakerate, lep2_fakerate;
                if(l==0){
                    //iso_lep: 0=muon passed iso, 1=electron passed iso
                    if(tm.iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, mu_FR.h);
                    else  lep1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, el_FR.h);

                    evt_reweight = lep1_fakerate/(1-lep1_fakerate);

                    if(tm.iso_lep ==1) h_mu_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 1);
                    if(tm.iso_lep ==0) h_el_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), 1);
                }
                if(l==1){
                    lep1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, mu_FR.h);
                    lep2_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, el_FR.h);
                    evt_reweight = -(lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));


                        h_mu_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f),  -0.5*(lep2_fakerate/(1-lep2_fakerate)));
                        h_el_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f),  -0.5*(lep1_fakerate/(1-lep1_fakerate)));
                }
                if(l==2){
                    if(tm.iso_lep ==1) lep1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, mu_FR.h);
                    else            lep1_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, el_FR.h);
                    evt_reweight = -(lep1_fakerate)/(1-lep1_fakerate);

                        if(tm.iso_lep ==1) h_mu_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), -tm.getEvtWeight());
                        if(tm.iso_lep ==0) h_el_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), -tm.getEvtWeight());

                }
                if(l==3){
                    lep1_fakerate = get_new_fakerate_prob(tm.mu1_pt, tm.mu1_eta, mu_FR.h);
                    lep2_fakerate = get_new_fakerate_prob(tm.el1_pt, tm.el1_eta, el_FR.h);
                    evt_reweight = (lep1_fakerate/(1-lep1_fakerate)) * (lep2_fakerate/(1-lep2_fakerate));

                        h_mu_err->Fill(min(abs(tm.mu1_eta), 2.3f), min(tm.mu1_pt, 150.f), 0.5*tm.getEvtWeight() * (lep2_fakerate/(1-lep2_fakerate)));
                        h_el_err->Fill(min(abs(tm.el1_eta), 2.3f), min(tm.el1_pt, 150.f), 0.5*tm.getEvtWeight() * (lep1_fakerate/(1-lep1_fakerate)));

                    }

                    double tot_evt_weight = 0.;
                    if(is_data) tot_evt_weight = evt_reweight; 
                    else{
                        tot_evt_weight = evt_reweight * tm.getEvtWeight();
                        //printf("%.3e %.3e \n", tm.getEvtWeight(), tot_evt_weight);
                    }

                    h_m->Fill(tm.m, tot_evt_weight);
                    h_pt->Fill(tm.cm.Pt(), tot_evt_weight);
                    h_cost->Fill(tm.cost, tot_evt_weight);
                    h_rap->Fill(tm.cm.Rapidity(),tot_evt_weight);
                    if(opp_sign) tot_weight_os += tot_evt_weight;
                    else tot_weight_ss += tot_evt_weight;
                }
            }

        printf("After iter %i current fakerate est is %.0f \n", l, h_m->Integral());
    }
    cleanup_hist(h_m);
    cleanup_hist(h_cost);
    cleanup_hist(h_pt);
    cleanup_hist(h_rap);


    if(sys_errors){
        set_fakerate_errors(h_mu_err, mu_FR.h, h_cost);
        set_fakerate_errors(h_mu_err, mu_FR.h, h_m);

        set_fakerate_errors(h_el_err, el_FR.h, h_cost);
        set_fakerate_errors(h_el_err, el_FR.h, h_m);
    }

    if(reweight){
        fakes_costrw_helper h_rw;
        setup_fakes_costrw_helper(&h_rw, year);

        float prev_weight = h_cost->Integral();
        fakes_cost_reweight(h_cost, h_rw.el_rw);
        float new_weight = h_cost->Integral();

        //scale other histograms to account for this reweight

        h_m->Scale(new_weight/prev_weight);
        h_pt->Scale(new_weight/prev_weight);
        h_rap->Scale(new_weight/prev_weight);
    }


    printf("Total fakerate est is %.0f \n", h_m->Integral());

}



void make_pileup_hist(TTree *t1, TH1F *h_before, TH1F *h_after, bool is_data=false, int year = 2016, int flag1 = FLAG_MUONS){
    //read event data
        TempMaker tm(t1, is_data, year);
        if(flag1 == FLAG_MUONS) tm.do_muons = true;
        else tm.do_electrons = true;

        tm.setup();
        int nEvents=0;
        double m_low = 150.;
        double m_high = 100000.;

        for (int i=0; i<tm.nEntries; i++) {
            tm.getEvent(i);
            bool pass = (tm.m >= m_low && tm.m <= m_high) && tm.met_pt < met_cut  && tm.has_no_bjets && tm.not_cosmic;

            if(pass){
                nEvents++;
                    Double_t evt_weight = tm.getEvtWeight();
                    if(tm.pu_SF < 1e-4) tm.pu_SF = 1.0;
                    h_before->Fill(tm.pu_NtrueInt, evt_weight/tm.pu_SF);
                    h_after->Fill(tm.pu_NtrueInt, evt_weight);

                }

        
        }
        printf("Selected %i events \n", nEvents);


    tm.finish();
}


