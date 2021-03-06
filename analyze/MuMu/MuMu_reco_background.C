#include "../../utils/NTupleReader.C"




void MuMu_reco_background(int nJobs =1, int iJob = 0, string fin ="", int year = -1)
{

    if(fin == "") fin = string("EOS_files/2016/WT_files.txt");
    NTupleReader nt(fin.c_str(),"output_files/test.root", false);
    if (year == -1) year = 2016;
    nt.year = year;
    printf("Year is %i \n", nt.year);

    nt.nJobs = nJobs;
    nt.iJob = iJob;
    nt.do_muons = true;
    nt.do_SFs = true;
    nt.do_RC = true;
    nt.RC_from_gen = false;
    if(fin.find("TTbar") != string::npos){
        printf("Turning on top pt reweighting \n");
        nt.do_top_ptrw = true;
    }

    nt.setupSFs();

    if(fin.find("LPair") != string::npos){
        nt.btag_mc_eff_idx = 1;
    }
    nt.setupRC();

    nt.setupOutputTree("T_sig");
    nt.setupOutputTree("T_WJets");
    nt.setupOutputTree("T_QCD");
    nt.setupOutputTree("T_ss");

    int iso_mu;
    nt.outTrees[1]->Branch("iso_mu", &iso_mu); 


    while(nt.getNextFile()){


        for (int i=0; i<nt.tin_nEntries; i++) {
            nt.getEvent(i);
            if(nt.good_trigger && nt.loose_dimuon_id && nt.cm_m > 50. ){
                //if(nt.nEvents > 10) exit(1);
                nt.fillEvent();

                nt.fillEventSFs();
                nt.fillEventRC();
                bool one_tight = nt.mu_tight_id0 ^ nt.mu_tight_id1;

                //pick the category
                if(nt.opp_sign && nt.tight_dimuon_id){ //signal region
                    nt.outTrees[0]->Fill();
                }
                else if(!nt.opp_sign && nt.tight_dimuon_id){ //samesign region
                    nt.outTrees[3]->Fill();
                }
                else if(one_tight){ //wjets control region
                    if(nt.mu_tight_id0) iso_mu = 0;
                    else           iso_mu = 1;
                    nt.outTrees[1]->Fill();
                }
                else if(!nt.mu_tight_id0 && !nt.mu_tight_id1){ //qcd control region
                    nt.outTrees[2]->Fill();
                }

            }
        } 

        printf("moving on to next file, currently %i events \n\n", nt.nEvents);


    }
    nt.finish();

    return;
}

