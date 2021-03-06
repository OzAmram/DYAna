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

const int type = FLAG_ELECTRONS;
int year = 2018;
char *plot_dir = "Misc_plots/";



void draw_fakes_ss_os_ratio(){
    setTDRStyle();
    init(year);

    TH1F * dummy = new TH1F("dummy", "", 100, 0, 100.);


    int n_cost_bins = 8;
    TH1F *ss_fakes_cost = new TH1F("QCD_cost", "QCD", n_cost_bins, -1.,1);
    TH1F *os_fakes_cost = new TH1F("WJets_cost", "WJets", n_cost_bins, -1.,1);

    bool cost_reweight = false;
    
    bool ss = false;

    float m_low = 150.;
    float m_high = 10000.;

    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, dummy, os_fakes_cost, dummy,dummy,dummy,dummy, type, year, m_low, m_high, ss,  cost_reweight);

    dummy->Reset();

    ss = true;
    make_fakerate_est(t_elel_WJets, t_elel_QCD, t_elel_WJets_contam, t_elel_QCD_contam, dummy, ss_fakes_cost, dummy,dummy,dummy,dummy, type, year, m_low, m_high, ss,  cost_reweight);

    symmetrize1d(os_fakes_cost);


    Double_t n_os = os_fakes_cost->Integral();
    Double_t n_ss = ss_fakes_cost->Integral();
    bool normalize = true;
    
    if(normalize){
        os_fakes_cost->Scale(1./n_os);
        ss_fakes_cost->Scale(1./n_ss);
    }
    ss_fakes_cost->Print();
    os_fakes_cost->Print();



    bool scale_error=false;
    float qcd_err = 0.5;
    float diboson_err = 0.05;
    bool add_err = true;

    


    char plot_file[100];
    TCanvas *c_ratio_plot = make_ratio_plot(string("elel_ss_cost_comparison"), os_fakes_cost, "Opposite Sign Fakes Est.", ss_fakes_cost, "Same Sign Fakes Est.", "ratio", "ee cos(#theta)", false, false, 0.01, 2.0);
    sprintf(plot_file, "%sy%i_elel_fakes_ss_os_cmp.png", plot_dir, year - 2000);
    c_ratio_plot->Print(plot_file);


}

    
    
