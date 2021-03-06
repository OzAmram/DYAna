#define STAND_ALONE
#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "../../utils/root_files.h"
#include "../../analyze/combine/TemplateUtils.h"



void draw_templates(){
        gStyle->SetOptStat(0);
        gROOT->SetBatch(1);
    
        int year = 2018;
        init(year);
        setup_all_SFs(year);
        //char *plot_dir = "Paper_plots/template_plots";
        char *plot_dir = "Misc_plots/template_plots";
        //setup_all_SFs(year);
        string sys_label = "";

        int num_bins = n_m_bins;
        const int n_rap_bins = 4;
        float rap_bins[] = {0., 0.6, 1., 1.5,  2.4};

        for(int i=0; i<num_bins; i++){
            double m_low = m_bins[i];
            double m_high = m_bins[i+1];

            char title[100];
            auto h_mumu_sym = new TH2F(title, "Symmetric template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_sym->SetDirectory(0);
            sprintf(title, "mumu%i_alpha%s", year %2000, sys_label.c_str());
            auto h_mumu_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_alpha->SetDirectory(0);
            sprintf(title, "mumu%i_asym%s", year %2000, sys_label.c_str());
            auto h_mumu_asym = new TH2F(title, "Asymmetric template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_mumu_asym->SetDirectory(0);

            auto h_elel_sym = new TH2F(title, "Symmetric template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_sym->SetDirectory(0);
            sprintf(title, "elel%i_alpha%s", year %2000, sys_label.c_str());
            auto h_elel_alpha = new TH2F(title, "Gauge boson polarization template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_alpha->SetDirectory(0);
            sprintf(title, "elel%i_asym%s", year %2000, sys_label.c_str());
            auto h_elel_asym = new TH2F(title, "Asymmetric template of mc",
                    n_rap_bins, rap_bins, n_cost_bins, cost_bins);
            h_elel_asym->SetDirectory(0);



            bool ss = false;
            bool use_xF = false;



            gen_mc_template(t_mumu_mc, h_mumu_sym, h_mumu_asym, h_mumu_alpha, year, m_low, m_high, FLAG_MUONS, use_xF, "");

            gen_mc_template(t_elel_mc, h_elel_sym, h_elel_asym, h_elel_alpha, year, m_low, m_high, FLAG_ELECTRONS, use_xF, "");


            char mu_title[100], el_title[100];
            sprintf(mu_title, "Muons Mass %.0f to %.0f", m_low, m_high);
            sprintf(el_title, "Electrons Mass %.0f to %.0f", m_low, m_high);

            char mu_fname1[100], mu_fname2[100], el_fname1[100], el_fname2[100];

            sprintf(mu_fname1, "%s/MuMu%i_M%.0f_sym_temps.png", plot_dir, year, m_low);
            sprintf(mu_fname2, "%s/MuMu%i_M%.0f_fit_temps.png", plot_dir, year, m_low);
            sprintf(el_fname1, "%s/ElEl%i_M%.0f_sym_temps.png", plot_dir, year, m_low);
            sprintf(el_fname2, "%s/ElEl%i_M%.0f_fit_temps.png", plot_dir, year, m_low);


            auto h_mumu_pl = *h_mumu_sym + *h_mumu_asym;
            auto h_mumu_mn = *h_mumu_sym - *h_mumu_asym;
            h_mumu_pl.Scale(0.5);
            h_mumu_mn.Scale(0.5);
            Double_t alpha= 0.05;
            double norm = 3./4./(2.+alpha);
            h_mumu_alpha->Scale(norm);

            //h_mumu_pl.Print("range");
            auto h1_mumu_pl = convert2d(&h_mumu_pl);
            //h1_mumu_pl->Print("range");
            auto h1_mumu_mn = convert2d(&h_mumu_mn);
            auto h1_mumu_alpha = convert2d(h_mumu_alpha);
            auto h1_mumu_sym = convert2d(h_mumu_sym);
            auto h1_mumu_asym = convert2d(h_mumu_asym);
            


            h1_mumu_alpha->SetLineColor(kGreen +3);
            h1_mumu_sym->SetLineColor(kBlue);
            h1_mumu_asym->SetLineColor(kRed+1);
            h1_mumu_pl->SetLineColor(kOrange +7);
            h1_mumu_mn->SetLineColor(kRed + 1);

            h1_mumu_alpha->SetLineWidth(2);
            h1_mumu_sym->SetLineWidth(2);
            h1_mumu_asym->SetLineWidth(2);
            h1_mumu_pl->SetLineWidth(2);
            h1_mumu_mn->SetLineWidth(2);

            h1_mumu_asym->SetMaximum(h1_mumu_sym->GetMaximum()*1.2);

            TCanvas *c_mumu1 = new TCanvas("c_mumu", "Histograms", 200, 10, 900, 700);
            h1_mumu_asym->SetTitle(mu_title); 
            h1_mumu_asym->Draw("hist");
            h1_mumu_sym->Draw("hist same ");

            float x_start = 0.75;
            float x_end = 0.9;

            float y_start = 0.75;
            float y_end = 0.9;

            TLegend *leg1 = new TLegend(x_start, y_start, x_end, y_end);
            leg1->AddEntry(h1_mumu_asym, "Asym Template", "l");
            leg1->AddEntry(h1_mumu_sym, "Sym Template", "l");
            leg1->Draw();
            c_mumu1->Print(mu_fname1);

            TCanvas *c_mumu2 = new TCanvas("c_mumu2", "Histograms", 200, 10, 900, 700);
            h1_mumu_pl->SetTitle(mu_title);
            h1_mumu_pl->Draw("hist");
            h1_mumu_mn->Draw("hist same");
            h1_mumu_alpha->Draw("hist same");


            c_mumu2->cd();
            TLegend *leg2 = new TLegend(x_start, y_start, x_end, y_end);
            leg2->AddEntry(h1_mumu_pl, "Plus Template", "l");
            leg2->AddEntry(h1_mumu_mn, "Minus Template", "l");
            leg2->AddEntry(h1_mumu_alpha, "#alpha Template", "l");
            leg2->Draw();

            c_mumu2->Print(mu_fname2);


            auto h_elel_pl = *h_elel_sym + *h_elel_asym;
            auto h_elel_mn = *h_elel_sym - *h_elel_asym;
            h_elel_pl.Scale(0.5);
            h_elel_mn.Scale(0.5);
            h_elel_alpha->Scale(norm);

            auto h1_elel_pl = convert2d(&h_elel_pl);
            auto h1_elel_mn = convert2d(&h_elel_mn);
            auto h1_elel_alpha = convert2d(h_elel_alpha);
            auto h1_elel_sym = convert2d(h_elel_sym);
            auto h1_elel_asym = convert2d(h_elel_asym);
            


            h1_elel_alpha->SetLineColor(kGreen +3);
            h1_elel_sym->SetLineColor(kBlue);
            h1_elel_asym->SetLineColor(kRed+1);
            h1_elel_pl->SetLineColor(kOrange +7);
            h1_elel_mn->SetLineColor(kRed + 1);

            h1_elel_alpha->SetLineWidth(2);
            h1_elel_sym->SetLineWidth(2);
            h1_elel_asym->SetLineWidth(2);
            h1_elel_pl->SetLineWidth(2);
            h1_elel_mn->SetLineWidth(2);

            h1_elel_asym->SetMaximum(h1_elel_sym->GetMaximum()*1.2);


            TCanvas *c_elel1 = new TCanvas("c_elel", "Histograms", 200, 10, 900, 700);
            h1_elel_asym->SetTitle(el_title);
            h1_elel_asym->Draw("hist");
            h1_elel_sym->Draw("hist same ");

            leg1->Draw();

            c_elel1->Print(el_fname1);

            TCanvas *c_elel2 = new TCanvas("c_elel2", "Histograms", 200, 10, 900, 700);
            h1_elel_pl->SetTitle(el_title);
            h1_elel_pl->Draw("hist");
            h1_elel_mn->Draw("hist same");
            h1_elel_alpha->Draw("hist same");


            leg2->Draw();

            c_elel2->Print(el_fname2);
        }
}




