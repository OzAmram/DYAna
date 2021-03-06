
import ROOT
from ROOT import *
from utils import *
import numpy as np

gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--Afb",  default=0.6, type='float', help="Afb value to inject")
parser.add_option("--A0",  default=0.05, type='float', help="A0 value to inject")
parser.add_option("--nToys",  default=100, type='int', help="How many toys to run")
parser.add_option("--mbin",  default=0, type='int', help="Which mass bin to run on ")
parser.add_option("-o", "--odir", default="", help = "output directory")
parser.add_option("--reuse_fit", default=False, action="store_true", help="Reuse initial fit from previous run to save time")
parser.add_option("--prefit", default=False, action="store_true", help="Sample toys from prefit uncs")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")

(options, args) = parser.parse_args()

gStyle.SetOptFit(1) 
chan = "combined"


workspace = "workspaces/%s_fit_bias_tests_%i.root" % (chan, options.mbin)
make_workspace(workspace, options.mbin, no_sys = options.no_sys)

extra_params = "--X-rtd MINIMIZER_no_analytic"


print("Will inject AFB %.2f A0 %.2f for all toys " %(options.Afb, options.A0))

res_afb = []
res_a0 = []
pull_afb = []
pull_a0 = []

if(not options.prefit):
    print("Sampling toys based on postfit")
    if(not options.reuse_fit):
        print_and_do("combine -M MultiDimFit -d %s --saveFit --saveWorkspace --robustFit 1 %s" % (workspace, extra_params))

for i in range(options.nToys):

    if(not options.prefit):
        fitted_afb, fitted_a0 = setSnapshot(Afb_val = -1., mdf = True)
        if(options.no_sys):
            print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i --snapshotName initialFit --toysNoSystematics --bypassFrequentistFit --saveToys -t 1  --setParameters Afb=%.2f,A0=%.2f" 
                % (i, options.Afb, options.A0))

        else:
            print_and_do("combine -M GenerateOnly -d initialFitWorkspace.root -s %i --snapshotName initialFit --toysFrequentist --bypassFrequentistFit --saveToys -t 1  --setParameters Afb=%.2f,A0=%.2f" 
            % (i, options.Afb, options.A0))
    else:

        print_and_do("combine -M GenerateOnly -d %s -s %i --saveToys -t 1 --toysFrequentist --setParameters Afb=%.2f,A0=%.2f" 
            % (workspace, i, options.Afb, options.A0))

    print_and_do("combine -M MultiDimFit -d %s --saveWorkspace --saveFitResult --toysFile higgsCombineTest.GenerateOnly.mH120.%i.root --toysFrequentist  -t 1 --robustFit 1 --forceRecreateNLL %s" 
            %(workspace, i, extra_params))
    f_fit = TFile.Open("multidimfitTest.root")
    if f_fit:
        fr = f_fit.Get('fit_mdf')
        myargs = RooArgSet(fr.floatParsFinal())

        Afb = myargs.find("Afb").getVal()
        Afb_err = myargs.find("Afb").getError()
        A0 = myargs.find("A0").getVal()
        A0_err = myargs.find("A0").getError()
        print("Afb %.3f err %.3f " % (Afb, Afb_err))
        res_afb.append(Afb - options.Afb)
        res_a0.append(A0 - options.A0)
        if(Afb_err > 0.): pull_afb.append((Afb- options.Afb)/ Afb_err)
        if(A0_err > 0.):  pull_a0.append((A0- options.A0)/ A0_err)
        

        #h_res_afb.Fill(Afb - options.Afb)
        #h_res_a0.Fill(A0 - options.A0)
        #if(Afb_err > 0.): h_pull_afb.Fill((Afb- options.Afb)/ Afb_err)
        #if(A0_err > 0.):  h_pull_a0.Fill((A0- options.A0)/ A0_err)
        f_fit.Close()
        #print_and_do("rm higgsCombineTest.GenerateOnly.mH120.%i.root" % i)
    else:
        print("Can't open multidimfit. Looks like the fit failed.")



#f_ws = TFile.Open("higgsCombineTest.FitDiagnostics.mH120.123456.root")
#f_toys = TFile.Open("higgsCombineTest.GenerateOnly.mH120.123456.root")
#ws = f_ws.Get('w')
#toy = f_toys.Get('toys/toy_1')
#getattr(ws, 'import')(toy)
#ws.writeToFile("toy_ws.root")
#print_and_do("PostFitShapesFromWorkspace -w toy_ws.root --dataset model_sData  -f fitDiagnostics.root:fit_s -o toy_shapes.root --sampling --samples 100")
##print_and_do("python scripts/plot_postfit.py -i toy_ws.root -o test/ -m %i" % (mbin))


n_bins = 20
h_pull_afb = TH1F("h_pull_Afb", "", n_bins, -3.5, 3.5)
h_pull_a0 = TH1F("h_pull_A0", "", n_bins, -3.5, 3.5)

res_afb_range = max(3.5*np.std(res_afb), 0.15)
res_a0_range = max(3.5*np.std(res_a0), 0.15)

h_res_afb = TH1F("h_res_Afb", "", n_bins, -res_afb_range, res_afb_range)
h_res_a0 = TH1F("h_res_A0", "", n_bins, -res_a0_range, res_a0_range)


def fill_h(arr, h):
    #print(arr)
    for arg in arr:
        h.Fill(arg)

fill_h(pull_afb, h_pull_afb)
fill_h(pull_a0, h_pull_a0)
fill_h(res_afb, h_res_afb)
fill_h(res_a0, h_res_a0)



c1 = TCanvas("c1", "", 900, 900)
h_pull_afb.Fit("gaus")
fit_afb= h_pull_afb.GetFunction("gaus")
if(fit_afb): fit_afb.SetLineColor(kBlue)
h_pull_afb.Draw()
h_pull_afb.SetTitle("Signal Inject Test Mass bin %i, Inject AFB %.2f" % (options.mbin, options.Afb))
h_pull_afb.GetXaxis().SetTitle("Pull Afb")
c1.Print("%sbias_test_mbin%i_afb%.0f.png" %(options.odir, options.mbin, 100.* options.Afb))


c2 = TCanvas("c1", "", 900, 900)
h_pull_a0.Fit("gaus")
fit_a0= h_pull_a0.GetFunction("gaus")
if(fit_a0): fit_a0.SetLineColor(kBlue)
h_pull_a0.Draw()
h_pull_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
h_pull_a0.GetXaxis().SetTitle("Pull A0")
c2.Print("%sbias_test_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))


c3 = TCanvas("c1", "", 900, 900)
h_res_afb.Fit("gaus")
fit_afb= h_res_afb.GetFunction("gaus")
if(fit_afb): fit_afb.SetLineColor(kBlue)
h_res_afb.Draw()
h_res_afb.SetTitle("Signal Inject Test Mass bin %i, Inject AFB %.2f" % (options.mbin, options.Afb))
h_res_afb.GetXaxis().SetTitle("#Delta Afb")
c3.Print("%sbias_test_res_mbin%i_afb%.0f.png" %(options.odir, options.mbin, 100.* options.Afb))


c4 = TCanvas("c1", "", 900, 900)
h_res_a0.Fit("gaus")
fit_a0= h_res_a0.GetFunction("gaus")
if(fit_a0): fit_a0.SetLineColor(kBlue)
h_res_a0.Draw()
h_res_a0.SetTitle("Signal Inject Test Mass bin %i, Inject A0 %.2f" % (options.mbin, options.A0))
h_res_a0.GetXaxis().SetTitle("#Delta A0")
c4.Print("%sbias_test_res_mbin%i_Az%.0f.png" %(options.odir, options.mbin, 100.* options.A0))

