"""
Abraham Tishelman-Charny
17 November 2021

The purpose of this notebook is to plot quantities from ETTAnalyzer outputs.

Example commands:
python3 ETTPlot.py --doSymLog

"""

import numpy as np
import os 
import argparse 
import uproot 
import pickle

from python.ETTPlot_Tools import * 

parser = argparse.ArgumentParser()
parser.add_argument("--variable", type = str, default = "RealVsEmu", help = "Variable to plot")
parser.add_argument("--selections", type = str, default = "clean", help = "Comma separated list of selections to apply")
parser.add_argument("--plotIndividuals", action="store_true", default = False, help = "Produce a plot for each selection")
parser.add_argument("--plotRatio", action="store_true", default = False, help = "Produce a ratio plot of tagged / all TPs")
parser.add_argument("--doSymLog", action="store_true", default = False, help = "Plot with symmetric log scale")
parser.add_argument("--addPlotText", type = int, default = 1, help = "Add text on plot")
parser.add_argument("--ymax", type = float, default = 256, help = "ymax")
parser.add_argument("--xmax", type = float, default = 256, help = "xmax")
parser.add_argument("--zmax", type = float, default = -1, help = "zmax")
parser.add_argument("--maxFiles", type = int, default = 1, help = "Max number of files to process")
args = parser.parse_args()

doSymLog = args.doSymLog 
ymax = args.ymax
xmax = args.xmax 
variable = args.variable
maxFiles = args.maxFiles 
selections = args.selections.split(',')
plotIndividuals = args.plotIndividuals
plotRatio = args.plotRatio 
addPlotText = args.addPlotText
zmax = args.zmax 

if(__name__ == '__main__'):

    varLabelDict = {
        "realVsEmu": "realVsEmu",
        "EnergyVsTimeOccupancy" : "EnergyVsTimeOccupancy"
    }

    varLabel = varLabelDict[variable]

    ##-- EB Occupancy
    direc = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346446_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211115_170649/0000_output/"
    n_files = 384 # number of input files 

    for selection in selections:
        values_0 = pickle.load(open("%s/%s_%s_values_0.p"%(direc, varLabel, selection), "rb"))
        values_1 = pickle.load(open("%s/%s_%s_values_1.p"%(direc, varLabel, selection), "rb"))

        total_values = np.add(values_0, values_1)

        for i in range(n_files):
            if(i > maxFiles):
                print("Max files reached: ",i)
                break 
            if(i == 0 or i == 1): continue 
            if(i%10 == 0): print("on file %s / %s"%(i, n_files))

            exec('thisValues = pickle.load(open("%s/%s_%s_values_%s.p"%(direc, varLabel, selection, i), "rb"))')  
            total_values = np.add(total_values, thisValues)

        exec("%s_%s_values = np.copy(total_values)"%(varLabel, selection))

    
    # add a second directory 
    #n_files_2 = 1000 
    direc_2 = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211116_115908/0000_output/"
    files_2 = ["%s/%s"%(direc_2, f) for f in os.listdir(direc_2) if varLabel in f]

    for file_i,file in enumerate(files_2):
        if(i > maxFiles):
            print("Max files reached: ",i)
            break         
        if(file_i%10 == 0): print("on file %s / %s"%(file_i, len(files_2)))

        thisValues = pickle.load(open(file, "rb"))
        # exec('thisValues = pickle.load(open("%s"%(file), "rb"))')  
        total_values = np.add(total_values, thisValues)     

    # Plot 
    vmaxAll = -1
    if(plotIndividuals):
        for selection in selections:
            print("On selection:",selection)
            exec("Values_array = %s_%s_values"%(varLabel, selection))
            isRatio = 0 
            plotText = selection
            if("all" in selection):
                vmaxAll = np.max(Values_array)
            plotText_params = [plotText, addPlotText]
            MakeETTPlot(Values_array, varLabel, selection, doSymLog, 0, plotText_params, vmaxAll, ymax, zmax)

    # take ratio 
    if(plotRatio):
        if(("clean_all" in selections) and ("clean_tagged" in selections)):

            exec("total_values = np.copy(%s_clean_all_values)"%(varLabel))
            exec("tagged_values = np.copy(%s_clean_tagged_values)"%(varLabel))

            total_values = np.array(total_values, dtype=float) 
            tagged_values = np.array(tagged_values, dtype=float)   

            # See which region of phase space is tagged 
            fraction = np.divide(tagged_values, total_values, out=np.zeros_like(tagged_values), where=total_values!=0)
            isRatio = 1 
            # doSymLog = 0 
            plotText = "Tagged / Total"
            plotText_params = [plotText, addPlotText]
            MakeETTPlot(fraction, "%s_ratio"%(varLabel), "clean", doSymLog, isRatio, plotText_params, vmaxAll, ymax, zmax)

    print("DONE")    