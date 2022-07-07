"""
7 July 2022 
Abraham Tishelman-Charny

The purpose of this module is to compare data and emulated double weights tagging probabilities.

Example usage:

python TaggingProbability.py --times inTime --datasets Run352912,PilotBeam2021  --severities zero
python TaggingProbability.py --times VeryLate --datasets Run352912,PilotBeam2021  --severities four

"""

import argparse 
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np 
from python.ETTPlot_Tools import GetBins, Add_CMS_Header
from matplotlib import pyplot as plt 

parser = argparse.ArgumentParser()
parser.add_argument("--times", type = str, default = "all,Early,inTime,Late,VeryLate", help = "Comma separated list of times to plot") 
parser.add_argument("--datasets", type = str, default = "zero,three,four", help = "Comma separated list of datasets to plot") 
parser.add_argument("--severities", type = str, default = "zero,three,four", help = "Comma separated list of severities to plot") 

args = parser.parse_args()

times = args.times.split(',')
datasets = args.datasets.split(',')
severities = args.severities.split(',')

vars = ["TaggingProbability", "EmulatedTaggingProbability"]

# get values 
for v in vars:
    for dataset in datasets:
        for severity in severities: 
            for time in times: 
                values_path = "output/{dataset}/{v}/{v}_Sev_{severity}_{time}-times_values.parquet".format(v=v, severity = severity, time = time, dataset=dataset)
                unc_path = "output/{dataset}/{v}/{v}_Sev_{severity}_{time}-times_statuncer.parquet".format(v=v, severity = severity, time = time, dataset=dataset)

                values_table = pq.read_table(values_path)
                yUnc_table = pq.read_table(unc_path)

                values_table_p = values_table.to_pandas()
                yUnc_table_p = yUnc_table.to_pandas()

                values = values_table_p.to_numpy().ravel()
                uncertainties = yUnc_table_p.to_numpy().ravel()

                exec("{v}_{dataset}_{severity}_{time}_values = np.copy(values)".format(v=v, severity = severity, time = time, dataset=dataset))
                exec("{v}_{dataset}_{severity}_{time}_unc = np.copy(uncertainties)".format(v=v, severity = severity, time = time, dataset=dataset))

# make ratio plot for each case: {dataset, severity, time}
verbose = 1 
for dataset in datasets:
    for severity in severities:
        for time in times:

            # get the values 
            exec("Data_Values = TaggingProbability_{dataset}_{severity}_{time}_values".format(severity = severity, time = time, dataset=dataset))
            exec("Data_Unc = TaggingProbability_{dataset}_{severity}_{time}_unc".format(severity = severity, time = time, dataset=dataset))
            exec("Emulator_Values = EmulatedTaggingProbability_{dataset}_{severity}_{time}_values".format(severity = severity, time = time, dataset=dataset))
            exec("Emulator_Unc = EmulatedTaggingProbability_{dataset}_{severity}_{time}_unc".format(severity = severity, time = time, dataset=dataset)) 

            if(verbose):
                print("dataset:",dataset)
                print("severity:",severity)
                print("time:",time)
                print(" ")
                print("Data values:",Data_Values)
                print("Emulator values:",Emulator_Values)

            # make the plot 
            fig, axarr = plt.subplots(2, 
                                        sharex=True, 
                                        gridspec_kw={
                                            'hspace': 0.15,
                                            'height_ratios': (0.7,0.3)
                                            }
                                        )    
            fig.set_size_inches(10, 7.5)
            upper = axarr[0]                   
            lower = axarr[1]       

            xbins, ybins = GetBins("EmulEnergyVsTimeOccupancy", "2021_2022_900GeVCollisions") # assuming same bins 
            energy_bins = ybins

            # assuming same rebinning as when the parquet files were made 
            rebin = 1
            ADCToGeV = 1 

            if(severity == "four"):
                reBins = [[1,32],[33,255]]
                # reBins = [[1,21], [22,23], [24,32], [33,255]]
            else: 
                reBins = [[5,32],[33,255]]                            

            if(rebin):             
                energy_bins_rebinned = np.copy(energy_bins)

                binEdgesRemoved = 0 
                for newBin in reBins:
                    newBinMin_i, newBinMax_i = newBin[0], newBin[1] 

                    # for bin edges, need to remove middle edges not including min and max
                    for b_i in range(newBinMin_i, newBinMax_i):
                        bin_i_toRemove = b_i + 1 - binEdgesRemoved
                        energy_bins_rebinned = np.delete(energy_bins_rebinned, bin_i_toRemove) 
                        binEdgesRemoved += 1 
                        if(b_i == newBinMax_i): continue # don't delete upper edge of rebin range 

                energy_bins = np.copy(energy_bins_rebinned)                                       

            if(ADCToGeV):
                energy_bins = [float(i)/2. for i in energy_bins]      

            centerTheXbins = 1

            if(centerTheXbins): centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
            else: centered_energy_bins_ = energy_bins[:-1]

            xerrors_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) for i in range(len(energy_bins) - 1) ]   
            
            # averages = np.array(averages_) 
            # yUnc = np.array(yUnc_) 
            centered_energy_bins = np.array(centered_energy_bins_)
            xerrors = np.array(xerrors_)

            # going to assume same mask for data and emulator which is not necessarily true 
            MASK = tuple([Data_Values != -1])

            centered_energy_bins = centered_energy_bins[MASK]
            Data_Values = Data_Values[MASK]
            Data_Unc = Data_Unc[MASK]
            Emulator_Values = Emulator_Values[MASK]
            Emulator_Unc = Emulator_Unc[MASK]            

            xerrors = xerrors[MASK]

            zero_errors = [0. for i in range(0, len(Data_Values))]            
        
            labelReplaceDict = {
                "MinDelta": "",
                "0p5prime" : "0.5",
                "2p5prime" : "2.5",

                "Sevzero" : "",
                "three" : "3",
                "Sevfour" : "Spike",
                "inTime" : "",
                "VeryLate" : ""
            }

            color = next(upper._get_lines.prop_cycler)['color']

            addXErr = 1 
            addLine = 1 
            error = 1 
            
            if(addXErr):
                xerr = xerrors 
                addLine = 0 
            else:
                xerr = zero_errors  

            # add values to upper plot 
            upper.scatter(x = centered_energy_bins, y = Data_Values, label = "Data", s = 15, color = "C0") 
            upper.errorbar(x = centered_energy_bins, y = Data_Values, xerr = xerr, yerr = Data_Unc, fmt = " ", color = "C0")  

            upper.scatter(x = centered_energy_bins, y = Emulator_Values, label = "Emulator", s = 15, color = "C1") 
            upper.errorbar(x = centered_energy_bins, y = Emulator_Values, xerr = xerr, yerr = Emulator_Unc, fmt = " ", color = "C1")              

            # add ratio to lower plot 

            EB_LABEL_XMIN = 0.18
            if(severity == "zero" and time == "inTime"):
                EB_LABEL_XMIN = 0.16 
                EB_LABEL_YMAX = 0.85
                loc = 'upper right'
                PlotTextLabel = "EM signal, |t|<3ns"
                xLabel = "Signal $E_{T}$ (GeV)"
            elif(severity == "four" and time == "VeryLate"):
                loc = 'lower right'
                EB_LABEL_YMAX = 0.45
                PlotTextLabel = "Spike, t$\geq$10 ns"
                xLabel = "Spike $E_{T}$ (GeV)"
            else:
                PlotTextLabel = "Sev %s, Times %s"%(severity, time)
                xLabel = "TP $E_{T}$ (GeV)"  

            plt.rcParams['legend.title_fontsize'] = 20
            upper.legend(loc = loc, fontsize = 20, prop={'size' : 14})
            plt.rcParams['legend.title_fontsize'] = 20

            upper.set_ylim(0, 1.01)

            xmin_, xmax_ = energy_bins[0], energy_bins[-1]
            xmax_ = 17
            print("Setting xmax to:",xmax_)

            plt.xlim(xmin_, xmax_)
            plt.xlabel(xLabel, fontsize=15)
            plt.grid()
            plt.xticks(fontsize = 15)
            plt.yticks(fontsize = 15)

            # upper.set_xlim(xmin_, xmax_)
            yLabel = "Tagging Probability"
            upper.set_ylabel(yLabel, fontsize=15)
            upper.grid()
            # upper.set_xticks_size(fontsize = 15)
            # upper.set_yticks_size(fontsize = 15)       
            upper.tick_params(axis='y', labelsize = 15)


            lower.set_ylabel(r"  $\frac{Data}{Emulator}$", fontsize=15)           
            lower.set_ylim(0, 2)      
            lower.plot([xmin_, xmax_],[1,1],linestyle=':', color = 'black')

            ratio_vals = np.divide(Data_Values, Emulator_Values, out=np.zeros_like(Emulator_Values), where=Emulator_Values!=0)

            # compute statistical uncertainty in ratio 
            ratio_unc = [] 

            for i, ratio_val in enumerate(ratio_vals):
                if(ratio_val == 0):
                    ratio_unc.append(0)
                else:
                    N_0_unc = float(Data_Unc[i])
                    N_1_unc = float(Emulator_Unc[i])

                    N_0 = float(Data_Values[i])
                    N_1 = float(Emulator_Values[i])

                    if(N_0 == 0):
                        N_0_unc_term = 0
                    else:
                        N_0_unc_term = (N_0_unc/N_0)**2

                    if(N_1 == 0):
                        N_1_unc_term = 0
                    else:
                        N_1_unc_term = (N_1_unc/N_1)**2                                      

                    rel_error = np.sqrt( N_0_unc_term + N_1_unc_term)
                    abs_error = ratio_val * float(rel_error)
                    ratio_unc.append(abs_error)                        

            if(error):
                lower.scatter(x = centered_energy_bins, y = ratio_vals, s = 15) ### [:-7] = remove the final 7 points (hack)
                lower.errorbar(x = centered_energy_bins, y = ratio_vals, xerr = xerr, yerr = ratio_unc, fmt = " ", color = 'C0')  
                #if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-', color = color)   
            else:
                lower.scatter(x = centered_energy_bins, y = ratio_vals, s = 15)
                #if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-')                        

            if(verbose):
                print("Data / Emulator values:",ratio_vals)

            upperRightText = "upperRightText"
            lumi = "lumi"
            addLumi = 1
            fontsize = 16
            text_xmin = 0.1                          

            unit = "{\mu}b"

            lumiDict = {
                "PilotBeam2021" : "82",
                "Run352912" : "145",
            }

            lumi = lumiDict[dataset]

            sqrts = "900 GeV"

            Add_CMS_Header(plt, upper, upperRightText, text_xmin, addLumi, lumi, fontsize, unit, sqrts)
            plt.text(
                EB_LABEL_XMIN, EB_LABEL_YMAX, u"ECAL Barrel",
                fontsize=14, fontweight='bold',
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=upper.transAxes
            )   

            plt.text(
                EB_LABEL_XMIN, EB_LABEL_YMAX - 0.075, PlotTextLabel,
                fontsize=14, fontweight='bold',
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=upper.transAxes
            )                           
            
            olDict = {
                "PilotBeam2021" : "/2021Collisions/Runs_346446_346447/Reemul_by_globalTag/",
                "Run352912" : "/2022Collisions/Run_352912/Reemul_by_globalTag/",
            }

            olLabel = olDict[dataset]
            ol = "/eos/user/a/atishelm/www/EcalL1Optimization/{olLabel}/".format(olLabel=olLabel)

            plt.tight_layout()
            plt.savefig("%s/Sev_%s_%s_DataOverEmulatorTaggingProbability_linear.png"%(ol, severity, time), dpi = 300)
            plt.savefig("%s/Sev_%s_%s_DataOverEmulatorTaggingProbability_linear.pdf"%(ol, severity, time), dpi = 300)   

            upper.set_ylim(0.00001, 1)
            upper.set_yscale('log')  
            plt.tight_layout()
            plt.savefig("%s/Sev_%s_%s_DataOverEmulatorTaggingProbability_log.png"%(ol, severity, time), dpi = 300)
            plt.savefig("%s/Sev_%s_%s_DataOverEmulatorTaggingProbability_log.pdf"%(ol, severity, time), dpi = 300)    

            plt.close()            

print("DONE")
