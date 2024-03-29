


import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup

years = [2016,2017,2018]
#prefixes = ["MuMu"]
prefixes = ["EMu"]
#years = [2018]
#prefixes = ["ElEl"]
labels = ["wt", "dy", "ttbar", "diboson"]
#labels = ["phot_ind"]

#labels = ["wt", "ttbar", "diboson", "dy", "phot_ind"]
#labels = ["dy_mlow", "data_mlow"]
#labels = ["dy"]
#labels = ["diboson", "phot_ind"]
#labels = ["phot_ind"]
#labels = ["data"]
#labels = ["data"]
ending = "_aug31"
redo_fakes = True

for year in years:
    for prefix in prefixes:
        f_names = []
        for label in labels:
            #os.system("rm output_files/" + str(year) +"/" + prefix + str(year % 2000) + "_" + label + "*.root") 
            #os.system("python doCondor.py -e -o output_files/" + str(year) +"/ -n " + prefix + str(year % 2000) + "_" + label + ending) 
            f_names.append("output_files/" + str(year) +"/" + prefix + str(year % 2000) + "_" + label + ending + ".root") 

        if(redo_fakes):
            rm_name = "output_files/" + str(year)+"/" + prefix +str(year %2000) + "_fakes_contam" + "*.root"
            os.system("rm %s" % rm_name)
            fakes_name = "output_files/" + str(year)+"/" + prefix +str(year %2000) + "_fakes_contam" + ending + ".root"
            cmd = "hadd " + fakes_name
            for name in f_names:
                cmd += " " + name
            print(cmd)
            os.system(cmd)
            os.system("./scripts/prune_root_file.sh %s" % fakes_name)

