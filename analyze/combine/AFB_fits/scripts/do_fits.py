from utils import *
import ROOT
from ROOT import *

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--chan",  default="combined", type="string", help="What channels to run the fit over (combined, ee, or mumu)")
parser.add_option("--no_plot",  default=False, action="store_true", help="Don't make postfit plots")
parser.add_option("--no_sys",  default=False, action="store_true", help="Use fit template without any shape systematics")
parser.add_option("--fake_data",  default=False, action="store_true", help="Use fit template without any shape systematics and no fakes")
parser.add_option("--no_cleanup",  default=False, action="store_true", help="Don't remove root files created by fit")
parser.add_option("--mbin", default = -1, type='int', help="Only do fits for this single mass bin, default is all bins")
parser.add_option("-y", "--year", default = -1, type='int', help="Only do fits for this single year (2016,2017, or 2018), default is all years")
parser.add_option("-v", "--verbose", default = 0, type='int', help="Turn up verbosity on fits")
parser.add_option("--noSymMCStats", default = False, action="store_true",  help="Don't add constraints to mcStat nuisances")
parser.add_option("--diff", default = False, action="store_true",  help="Measure difference between electron and muon AFB's")
parser.add_option("--sep", default = False, action="store_true",  help="Measure separate electron and muon AFB's")
parser.add_option("--ratio", default = False, action="store_true",  help="Measure ratio between electron and muon AFB's")
parser.add_option("--all_mbins", default = False, action="store_true",  help="All da mbins")

(options, args) = parser.parse_args()


extra_params = ""

if(options.diff and not (options.chan == "combined")):
    print("invalid options. Selected diff but channel is %s " % options.chan)
    sys.exit(1)

if(options.chan == "ee"):
    print("Chan is ee, will mask mumu channels")
    if(options.year < 0): 
        extra_params += " --setParameters mask_Y16_mumu16=1,mask_Y17_mumu17=1,mask_Y18_mumu18=1" 
    else:
        extra_params += " --setParameters mask_Y%i_mumu%i=1" % (options.year % 2000, options.year % 2000)
elif(options.chan == "mumu"):
    print("Chan is mumu, will mask ee channels")
    if(options.year < 0):
        extra_params += " --setParameters mask_Y16_ee16=1,mask_Y17_ee17=1,mask_Y18_ee18=1" 
    else:
        extra_params += " --setParameters mask_Y%i_ee%i=1" % (options.year % 2000, options.year % 2000, options.year % 2000, options.year % 2000)

if(options.verbose > 0):
    extra_params +=" --verbose %i" % options.verbose

#No analytic minimization of MC stats nuisances
extra_params += " --X-rtd MINIMIZER_no_analytic"



bin_start = 1
bin_stop = 8

fit_name = options.chan
if(options.diff): fit_name += "_diff"
if(options.sep): fit_name += "_sep"
if(options.ratio): fit_name += "_ratio"
if(options.no_sys): 
    fit_name +="_nosys"
    extra_params += " --freezeParameters allConstrainedNuisances"
if(options.noSymMCStats):
    fit_name += "_noSymMC"
if(options.fake_data): fit_name +="_fake_data"
if(options.year > 0): fit_name +="_y%i" % (options.year % 2000)


if(options.mbin >= 0):
    print("Will only do fit for bin %i " % options.mbin)
    bin_start = options.mbin
    bin_stop = bin_start + 1

mbins = [i for i in range(bin_start, bin_stop)]

if(options.all_mbins):
    extra_params += " --cminApproxPreFitTolerance 1.0 --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerStrategy 0 "
    mbins = [9]




for mbin in mbins:
    print(" \n \n Starting fit for bin %i \n\n" % mbin)

    workspace="workspaces/%s_fit_%i.root" % (fit_name, mbin)

    if(options.all_mbins):
        make_all_mbins_workspace(workspace, year = options.year, symMCStats = not (options.noSymMCStats), diff = options.diff , ratio = options.ratio, sep = options.sep)
    else:
        make_workspace(workspace, mbin, fake_data = options.fake_data, year = options.year, symMCStats = not (options.noSymMCStats), 
                diff = options.diff , sep = options.sep, ratio = options.ratio)

    plotdir="postfit_plots/%s_mbin%i" % (fit_name, mbin)
    print_and_do("[ -e %s ] && rm -r %s" % (plotdir, plotdir))
    print_and_do("mkdir %s" % (plotdir))

    plotdir_swap="postfit_plots/%s_swap_axis_mbin%i" % (fit_name, mbin)
    print_and_do("[ -e %s ] && rm -r %s" % (plotdir_swap, plotdir_swap))
    print_and_do("mkdir %s" % (plotdir_swap))
    print_and_do("combine %s -M MultiDimFit  --saveWorkspace --saveFitResult --robustFit 1 %s" %(workspace, extra_params))

    if(not options.no_plot ):
        print_and_do("PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfitTest.root:fit_mdf --postfit -o %s_fit_shapes_mbin%i.root --samples 100"
                % (fit_name, mbin))
        print_and_do("python scripts/sample_comb_postfit.py -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfitTest.root:fit_mdf -o %s_fit_shapes_mbin%i.root -n 100"
                % (fit_name, mbin))
        extra_args = ""
        if(options.year > 0): extra_args = " -y %i " % options.year
        print_and_do("python scripts/plot_postfit.py -i %s_fit_shapes_mbin%i.root -o %s -m %i %s" % (fit_name, mbin, plotdir, mbin, extra_args))
        print_and_do("python scripts/plot_comb_postfit.py -i %s_fit_shapes_mbin%i.root -o %s -m %i" % (fit_name, mbin, plotdir, mbin))
        #print_and_do("python scripts/plot_swaped_axis_postfit.py -i %s_fit_shapes_mbin%i.root -o %s -m %i %s" % (fit_name, mbin, plotdir_swap, mbin, extra_args))
        print_and_do("combine %s -M FitDiagnostics --skipBOnlyFit %s" % (workspace, extra_params)) #only to get prefit, probably a better way
        if(not options.diff and not options.ratio and not options.sep): 
            print_and_do("python scripts/my_diffNuisances.py multidimfitTest.root --multidim --prefit fitDiagnosticsTest.root -p Afb --skipFitB -g %s" % (plotdir))
        print_and_do("mv %s_fit_shapes_mbin%i.root %s" %(fit_name, mbin, plotdir))
        if(not options.no_cleanup): print_and_do("rm fitDiagnosticsTest.root higgsCombineTest.FitDiagnostics.mH120.root")


    print_and_do("""echo "fit_mdf->Print();" > cmd.txt""")
    print_and_do("""echo ".q" >> cmd.txt """)
    print_and_do("root -l -b multidimfitTest.root < cmd.txt > fit_results/%s_fit_results_mbin%i.txt" % (fit_name, mbin))
    print_and_do("rm -f cards/sed*")
    if(not options.no_cleanup): print_and_do("rm cmd.txt combine_logger.out higgsCombineTest.MultiDimFit.mH120.root multidimfitTest.root")

