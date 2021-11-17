"""
Abraham Tishelman-Charny
17 November 2021

The purpose of this notebook is to plot quantities from ETTAnalyzer outputs.
"""

import numpy as np
import os 
import argparse 
import uproot 
import pickle
import matplotlib
from matplotlib import pyplot as plt 
from matplotlib.colors import LogNorm

from python.ETTPlot_Tools import * 

# parser = argparse.ArgumentParser()
# parser.add_argument('--nodes',default = "", required=False, type=str, help = "Comma separated list of nodes to run")
# parser.add_argument('--reweightNodes',default = "cttHH3", required=False, type=str, help = "Comma separated list of nodes to reweight to")
# parser.add_argument('--years',default = "2017", required=False, type=str, help = "Comma separated list of years to run")
# parser.add_argument('--NominalOnly',action="store_true",help = "Only run on nominal tree")
# parser.add_argument('--categorize', action="store_true", required=False, help = "Split trees into categories based on DNN score")
# args = parser.parse_args()

if(__name__ == '__main__'):

    # parameters 
    isWide = 0
    varLabel = "realVsEmu"
    upperRightText = "Pilot Beam 2021"
    xmin = 0.115

    ##-- EB Occupancy
    direc = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346446_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211115_170649/0000_output/"
    n_files = 384 # number of input files 

    selection = "clean"

    values_0 = pickle.load(open("%s/%s_%s_values_0.p"%(direc, varLabel, selection), "rb"))
    values_1 = pickle.load(open("%s/%s_%s_values_1.p"%(direc, varLabel, selection), "rb"))

    total_values = np.add(values_0, values_1)

    for i in range(n_files):
        if(i == 0 or i == 1): continue 
        if(i%10 == 0): print("on file %s / %s"%(i, n_files))

        exec('thisValues = pickle.load(open("%s/%s_%s_values_%s.p"%(direc, varLabel, selection, i), "rb"))')  
        #print("thisValues:",thisValues)
        total_values = np.add(total_values, thisValues)

    print("total_values:",total_values)

    ##-- Save array for plotting in next cell
    EBOcc_values = np.copy(total_values)
    exec("%s_%s_values = np.copy(total_values)"%(varLabel, selection))

    # print("Values:",total_values)

    ##-- Real vs emulated TPs
    isWide = 0

    ##-- Prepare figure and axes 
    fig, ax = plt.subplots()
    fig.set_dpi(100)
    # fig.set_size_inches(8, 6)
    #     fig.set_size_inches(15, 7.5)
    fig.set_size_inches(10, 7.5)
    cmap = plt.cm.Blues
    cmap.set_under(color='white') 

    exec("Values_array = %s_%s_values"%(varLabel, selection))    

    xbins = np.linspace(0,256,134) # coffea and dqm binning
    ybins = np.linspace(0,256,134) # coffea and dqm binning

    ##-- Custom binning 
    ymin_large = int(ymin * 1000)
    ybins = range(ymin_large, 1200, 25)
    Values_array = SetYMin(Values_array, ymin)

    # plot with colormesh 
    pos = ax.pcolormesh(xbins, ybins, Values_array.transpose(1,0), cmap = cmap, vmin = 1)
    cb = fig.colorbar(pos, 
                    ax=ax,
                )    

    """
    # plot with imshow 
    pos = ax.imshow(Values_array.transpose(1,0),
                interpolation='none',
                aspect="auto",
                extent = [0, 256, 256, 0],
    #                vmin=0.00000001, vmax=1,
    #            vmin=0, vmax=1,
    #            vmin=0.000001, vmax=1,
    #                 vmin=0,
                cmap = cmap,
    #            norm = LogNorm(vmin = 1, vmax = eval("sev_%s_all_maxValue"%(severity))) ##-- Set vmax to max from total number of TPs
    #                norm = LogNorm(vmin = 1) ##-- Set vmax to max from total number of TPs
                norm = LogNorm(vmin = 1) ##-- Set vmax to max from total number of TPs
                )
    plt.gca().invert_yaxis()

    cb = fig.colorbar(pos, 
                    ax=ax,
                )

    """

    plt.xlabel("Emulated TP Et (ADC)", fontsize=20)
    plt.ylabel("Real data TP Et (ADC)", fontsize=20)
    #     plt.xscale('log')

    Add_CMS_Header(plt, isWide, ax, upperRightText, xmin)

    # plt.grid()
    # plt.show()
    # ol = "/eos/user/a/atishelm/www/EcalL1Optimization/ZeroBias/"
    # plt.savefig("%s/EnergyVsTimePercentTagged_sev%s_2d.png"%(ol, sevLabel))
    # plt.savefig("%s/EnergyVsTimePercentTagged_sev%s_2d.pdf"%(ol, sevLabel))

    plt.savefig("RealVsEmu.png")
    plt.close()

    print("DONE")    