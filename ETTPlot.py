"""
Abraham Tishelman-Charny
17 November 2021

The purpose of this notebook is to plot quantities from ETTAnalyzer outputs.

Example commands:
python3 ETTPlot.py --directory /eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220209_125921/0000_output/ --variables oneMinusEmuOverRealvstwrADCCourseBinning --maxFiles 100
python3 ETTPlot.py --doSymLog

"""

import numpy as np
import os 
import argparse 
import uproot 
import pickle
from matplotlib import pyplot as plt 

from python.ETTPlot_Tools import * 

parser = argparse.ArgumentParser()
parser.add_argument("--directory", type = str, default = "NODIRECTORY", required=True, help = "Directory with input files")
parser.add_argument("--variables", type = str, default = "RealVsEmu", help = "Comma separated list of variables to plot") 
parser.add_argument("--times", type = str, default = "all,Early,inTime,Late,VeryLate", help = "Comma separated list of times to plot") 
parser.add_argument("--severities", type = str, default = "zero,three,four", help = "Comma separated list of severities to plot") 

# parser.add_argument("--selections", type = str, default = "clean", help = "Comma separated list of selections to apply")
parser.add_argument("--plotIndividuals", action="store_true", default = False, help = "Produce a plot for each selection")
parser.add_argument("--plotRatio", action="store_true", default = False, help = "Produce a ratio plot of tagged / all TPs")
parser.add_argument("--doSymLog", action="store_true", default = False, help = "Plot with symmetric log scale")
parser.add_argument("--addPlotText", type = int, default = 1, help = "Add text on plot")
# parser.add_argument("--ymax", type = float, default = 256, help = "ymax")
# parser.add_argument("--xmax", type = float, default = 256, help = "xmax")
# parser.add_argument("--zmax", type = float, default = -1, help = "zmax")
parser.add_argument("--maxFiles", type = int, default = 1, help = "Max number of files to process")
args = parser.parse_args()

directory = args.directory 
doSymLog = args.doSymLog 
# ymax = args.ymax
# xmax = args.xmax 
variables = args.variables.split(',')
times = args.times.split(',')
severities = args.severities.split(',')
maxFiles = args.maxFiles 
# selections = args.selections.split(',')
plotIndividuals = args.plotIndividuals
plotRatio = args.plotRatio 
addPlotText = args.addPlotText
# zmax = args.zmax 

if(__name__ == '__main__'):

    varLabelDict = {
        "realVsEmu": "realVsEmu",
        "EnergyVsTimeOccupancy" : "EnergyVsTimeOccupancy",
        "oneMinusEmuOverRealvstwrADCCourseBinning" : "oneMinusEmuOverRealvstwrADCCourseBinning"
    }

    for variable in variables:
        varLabel = varLabelDict[variable]
        for severity in severities:
            for time in times:
                if(time == "inTime" and severity == "three"):
                    print("Skipping in time severity three as they shouldn't exist by definition")
                    continue 
                print("On Var: %s, Sev: %s, time: %s"%(variable, severity, time))
                thisDirec = "%s/%s/%s/%s/"%(directory, variable, severity, time) # directory for a given variable, severity, time 
                files = [f for f in os.listdir(thisDirec)]
                n_files = len(files)
                values_0 = pickle.load(open("%s/%s_sev%s_%s_values_0.p"%(thisDirec, variable, severity, time), "rb"))
                values_1 = pickle.load(open("%s/%s_sev%s_%s_values_1.p"%(thisDirec, variable, severity, time), "rb"))
                total_values = np.add(values_0, values_1)
                for i in range(n_files):
                    if(i > maxFiles):
                        print("Max files reached: ",i)
                        break 
                    if(i == 0 or i == 1): continue 
                    if(i%100 == 0): print("on file %s / %s"%(i, n_files))

                    exec('thisValues = pickle.load(open("%s/%s_sev%s_%s_values_%s.p"%(thisDirec, variable, severity, time, i), "rb"))')  
                    total_values = np.add(total_values, thisValues)

                exec("%s_%s_%s_values = np.copy(total_values)"%(variable, severity, time)) # copy with unique name to plot together on same plot later 
                exec("These_Values = np.copy(total_values)") 

                averages, stdevs = MakeETTPlot(These_Values, variable, severity, time) # make plots and return averages 

                exec("%s_%s_averages = np.copy(averages)"%(severity, time))
                exec("%s_%s_stdevs = np.copy(stdevs)"%(severity, time))

        colorDict = {
            "all" : "blue",
            "Early" : "cyan",
            "inTime" : "green",
            "Late" : "orange",
            "VeryLate" : "red"
        }

        shapeDict = {
            "zero" : "o",
            "three" : "^",
            "four" : "s"
        }   

        ol = "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta2p5prime_WithOddPF_MultiFitReco/"
        error = 1 
        log = 1

        # plot average lines on same plots
        if(variable == "oneMinusEmuOverRealvstwrADCCourseBinning"):
            xbins, ybins = GetBins(variable)
            ##-- A plot for each timing 
            for time in times:
                fig, ax = plt.subplots()
                for sev_i, severity in enumerate(severities):
                    if(time == "inTime" and severity == "three"): 
                        print("Skipping in time severity 3")
                        continue     
                    shape = shapeDict[severity]
                    color = colorDict[time]

                    exec("averages_ = np.copy(%s_%s_averages)"%(severity, time))
                    exec("stdevs_ = np.copy(%s_%s_stdevs)"%(severity, time))

                    energy_bins = xbins
                    energies_ = xbins
                    xmin_, xmax_ = xbins[0], xbins[-1]
                    centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
                    xerrors_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) for i in range(len(energy_bins) - 1) ]   
                    
                    averages = np.array(averages_) 
                    stdevs = np.array(stdevs_) 
                    centered_energy_bins = np.array(centered_energy_bins_)
                    xerrors = np.array(xerrors_)

                    MASK = tuple([averages != -1])

                    centered_energy_bins = centered_energy_bins[MASK]
                    averages = averages[MASK]
                    stdevs = stdevs[MASK]  
                    xerrors = xerrors[MASK]

                    zero_errors = [0. for i in range(0, len(averages))]            
        
                    if(error):
                        plt.scatter(x = centered_energy_bins, y = averages, label = "Severity = %s, %s"%(severity, time), s = 15)
                        plt.errorbar(x = centered_energy_bins, y = averages, xerr = xerrors, yerr = zero_errors, fmt = " ")  
                    else:
                        plt.scatter(x = energies, y = averages, label = "Severity = %s, %s"%(severity, time), s = 10, marker = shape, color = color)
                plt.legend(loc = 'best', fontsize = 10)
                plt.ylim(0, 1.01)
                plt.xlim(xmin_, xmax_)
                plt.xlabel("Real data TP Et (ADC)", fontsize=15)
                plt.ylabel("Average 1 - (Emulated / Real)", fontsize=15)
                plt.grid()
                if(log):
                    plt.ylim(0.0001, 1)
                    plt.yscale('log')  
                plt.savefig("%s/Sev_all_Average_%s_%s.png"%(ol, varLabel, time), dpi = 300)
                plt.savefig("%s/Sev_all_Average_%s_%s.pdf"%(ol, varLabel, time), dpi = 300)    
                plt.close()
                
            ##-- A plot for each severity 
            for sev_i, severity in enumerate(severities):
                fig, ax = plt.subplots()

                for time in times:
                    if(time == "inTime" and severity == "three"): 
                        print("Skipping in time severity 3")
                        continue     

                    shape = shapeDict[severity]
                    color = colorDict[time]

                    exec("averages_ = np.copy(%s_%s_averages)"%(severity, time))
                    exec("stdevs_ = np.copy(%s_%s_stdevs)"%(severity, time))
                    
                    energy_bins = xbins
                    energies_ = xbins
                    xmin_, xmax_ = xbins[0], xbins[-1]
                    centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
                    xerrors_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) for i in range(len(energy_bins) - 1) ]   
                    
                    averages = np.array(averages_) 
                    stdevs = np.array(stdevs_) 
                    centered_energy_bins = np.array(centered_energy_bins_)
                    xerrors = np.array(xerrors_)

                    MASK = tuple([averages != -1])

                    centered_energy_bins = centered_energy_bins[MASK]
                    averages = averages[MASK]
                    stdevs = stdevs[MASK]  
                    xerrors = xerrors[MASK]

                    zero_errors = [0. for i in range(0, len(averages))]            

                    if(error):
                        plt.scatter(x = centered_energy_bins, y = averages, label = "Severity = %s, %s"%(severity, time), s = 15)
                        plt.errorbar(x = centered_energy_bins, y = averages, xerr = xerrors, yerr = zero_errors, fmt = " ")  
                    else: #yer
                        plt.scatter(x = energies_, y = averages, label = "Severity = %s, %s"%(severity, time), s = 10, marker = shape, color = color)

                plt.legend(loc = 'best', fontsize = 10)
                plt.ylim(0, 1.01)
                plt.xlim(xmin_, xmax_)
                plt.xlabel("Real data TP Et (ADC)", fontsize=15)
                plt.ylabel("Average 1 - (Emulated / Real)", fontsize=15)
                plt.grid()
                if(log):
                    plt.ylim(0.0001, 1)
                    plt.yscale('log')  
                plt.savefig("%s/Sev_%s_Average_%s_allTimings.png"%(ol, severity, varLabel), dpi = 300)
                plt.savefig("%s/Sev_%s_Average_%s_allTimings.pdf"%(ol, severity, varLabel), dpi = 300)    
                plt.close()    

    print("DONE")

    """

        # ##-- EB Occupancy
        # #direc = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346446_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211115_170649/0000_output/"
        # n_files = 384 # number of input files 

        # for selection in selections:
        #     values_0 = pickle.load(open("%s/%s_%s_values_0.p"%(direc, varLabel, selection), "rb"))
        #     values_1 = pickle.load(open("%s/%s_%s_values_1.p"%(direc, varLabel, selection), "rb"))

        #     total_values = np.add(values_0, values_1)

        #     for i in range(n_files):
        #         if(i > maxFiles):
        #             print("Max files reached: ",i)
        #             break 
        #         if(i == 0 or i == 1): continue 
        #         if(i%10 == 0): print("on file %s / %s"%(i, n_files))

        #         exec('thisValues = pickle.load(open("%s/%s_%s_values_%s.p"%(direc, varLabel, selection, i), "rb"))')  
        #         total_values = np.add(total_values, thisValues)

        #     exec("%s_%s_values = np.copy(total_values)"%(varLabel, selection))

        
        # # add a second directory 
        # #n_files_2 = 1000 
        # direc_2 = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211116_115908/0000_output/"
        # files_2 = ["%s/%s"%(direc_2, f) for f in os.listdir(direc_2) if varLabel in f]

        # for file_i,file in enumerate(files_2):
        #     if(i > maxFiles):
        #         print("Max files reached: ",i)
        #         break         
        #     if(file_i%10 == 0): print("on file %s / %s"%(file_i, len(files_2)))

        #     thisValues = pickle.load(open(file, "rb"))
        #     # exec('thisValues = pickle.load(open("%s"%(file), "rb"))')  
        #     total_values = np.add(total_values, thisValues)     

        # # Plot 
        # vmaxAll = -1
        # if(plotIndividuals):
        #     for selection in selections:
        #         print("On selection:",selection)
        #         exec("Values_array = %s_%s_values"%(varLabel, selection))
        #         isRatio = 0 
        #         plotText = selection
        #         if("all" in selection):
        #             vmaxAll = np.max(Values_array)
        #         plotText_params = [plotText, addPlotText]
        #         MakeETTPlot(Values_array, varLabel, selection, doSymLog, 0, plotText_params, vmaxAll, ymax, zmax)

        # # take ratio 
        # if(plotRatio):
        #     if(("clean_all" in selections) and ("clean_tagged" in selections)):

        #         exec("total_values = np.copy(%s_clean_all_values)"%(varLabel))
        #         exec("tagged_values = np.copy(%s_clean_tagged_values)"%(varLabel))

        #         total_values = np.array(total_values, dtype=float) 
        #         tagged_values = np.array(tagged_values, dtype=float)   

        #         # See which region of phase space is tagged 
        #         fraction = np.divide(tagged_values, total_values, out=np.zeros_like(tagged_values), where=total_values!=0)
        #         isRatio = 1 
        #         # doSymLog = 0 
        #         plotText = "Tagged / Total"
        #         plotText_params = [plotText, addPlotText]
        #         MakeETTPlot(fraction, "%s_ratio"%(varLabel), "clean", doSymLog, isRatio, plotText_params, vmaxAll, ymax, zmax)

        #print("DONE")    

        """