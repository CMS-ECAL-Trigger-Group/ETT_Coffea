"""
17 November 2021
Abraham Tishelman-Charny 

The purpose of this module is to provide tools for ETTPlot. 
"""

import numpy as np 
import matplotlib
from matplotlib import pyplot as plt 
from matplotlib.colors import LogNorm, SymLogNorm

##-- CMS header 
def Add_CMS_Header(plt, ax, upperRightText, xmin):
    ##-- Upper left plot text
    ##-- CMS 
    plt.text(
        0., 1., u"CMS ",
        fontsize=20, fontweight='bold',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )

    prelim_x = xmin
    
    ##-- Preliminary 
    plt.text(
        prelim_x, 0.998, u"$\it{Preliminary}$",
        fontsize=19,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )    

    # upper right text 
    plt.text(
        1., 1., upperRightText,
        fontsize=16, horizontalalignment='right', 
        verticalalignment='bottom', 
        transform=ax.transAxes
    )  

# for changing the axes of plots made from pickled arrays 
def SetYMin(Values_array_, ymin_, binWidth_):
    Transformed_List = []
    ymin_bin_i = (1./binWidth_) * ymin_ + 10./binWidth_ 
    for i, xBinVals in enumerate(Values_array_):
        nVals = len(xBinVals)
        MASK = [i > ((ymin_bin_i) - 1) for i, val in enumerate(xBinVals)]
        xBinVals = xBinVals[MASK]
        Transformed_List.append(xBinVals)
    Transformed_List = np.array(Transformed_List)
    return Transformed_List

def SetYMax(Values_array_, ymax_, binWidth_):
    Transformed_List = [] 
    ymax_bin_i = ymax_ # for [0, 256] with 1 spaced binnings 
    for i, xBinVals in enumerate(Values_array_):
        nVals = len(xBinVals)
        MASK = [i < ((ymax_bin_i) + 1) for i, val in enumerate(xBinVals)]
        xBinVals = xBinVals[MASK]
        Transformed_List.append(xBinVals)        
    Transformed_List = np.array(Transformed_List)
    return Transformed_List

def SetXMax(Values_array_, xmax_, binWidth_): # as you loop the Values_array, each object is a y slice starting from the left 
    Transformed_List = [] 
    xmax_bin_i = xmax_ # for [0, 256] with 1 spaced binnings 
    for i, xBinVals in enumerate(Values_array_):
        if(i > xmax_bin_i): continue 
        else: Transformed_List.append(xBinVals)        
    Transformed_List = np.array(Transformed_List)
    return Transformed_List    

def GetBins(varLabel_):
    binDict = {
        "realVsEmu" : [[0, 256, 256], [0, 256, 256]],
        "EnergyVsTimeOccupancy" : [[-50, 50, 100],[0, 35, 35]],
        "EnergyVsTimeOccupancy_ratio" : [[-50, 50, 100],[0, 35, 35]]

    }

    xmin, xmax, xbins = binDict[varLabel_][0]
    ymin, ymax, ybins = binDict[varLabel_][1]

    xbinning = np.linspace(xmin, xmax, xbins + 1)
    ybinning = np.linspace(ymin, ymax, ybins + 1)

    return [xbinning, ybinning]

def GetPlotLabels(varLabel_):
    labelDict = {
        "realVsEmu" : ["Emulated TP Et (ADC)", "Real data TP Et (ADC)"],
        "EnergyVsTimeOccupancy" : ["time (ns)", "Real data TP Et (ADC)"],
        "EnergyVsTimeOccupancy_ratio" : ["time (ns)", "Real data TP Et (ADC)"]

    }

    return labelDict[varLabel_]     

def MakeETTPlot(Values_array, varLabel, selection, doSymLog, isRatio, plotText, vmaxAll):

    # parameters 
    upperRightText = "Pilot Beam 2021"
    text_xmin = 0.1

    ##-- Prepare figure and axes 
    fig, ax = plt.subplots()
    
    fig.set_dpi(100)
    fig.set_size_inches(10, 7.5)
    # cmap = plt.cm.Blues
    # cmap = plt.cm.RdBu_r # for sym log 
    cmap = plt.cm.jet
    cmap.set_under(color='white') 

    #exec("Values_array = %s_%s_values"%(varLabel, selection))    

    xbins, ybins = GetBins(varLabel)

    # ##-- Custom binning 
    # ymin_large = int(ymin * 1000)
    # ymin = 0 
    # ybins = range()
    # ybins = range(ymin_large, 1200, 25)

    # ymin = 0 
    # binWidth = 1
    # Values_array = SetYMin(Values_array, ymin, binWidth)
    # Values_array = SetYMin(Values_array, ymin)
    # xmax = 50
    # ymax = 50

    # custom binning for realVsEmu    
    # binWidth = 1 
    # Values_array = SetYMax(Values_array, ymax, binWidth)
    # Values_array = SetXMax(Values_array, xmax, binWidth)
    #xbins = np.linspace(0, xmax, int(xmax)+1)
    #ybins = np.linspace(0, ymax, int(ymax)+1)

    # if("Tagged" in selection):
    #     vmax = vmaxAll
    # else:
    #     vmax = None

    if("tagged" in selection):
        print("Setting clim")
        vmax = vmaxAll
    else:
        vmax = None
        # plt.clim(0, vmaxAll)

    # plot with colormesh 
    if(isRatio): 
        vmin = 0.00000001
    else: 
        vmin = 1 
    if(doSymLog): 
        norm = norm = SymLogNorm(linthresh=0.03, vmin = vmin)
    else: 
        norm = None 
    pos = ax.pcolormesh(xbins, 
                        ybins, 
                        Values_array.transpose(1,0), 
                        cmap = cmap, 
                        vmin = vmin,
                        vmax = vmax,
                        norm = norm
                        )
    cb = fig.colorbar(pos, 
                    ax=ax,
                )    

    # object_methods = [method_name for method_name in dir(cb)
                #   if callable(getattr(cb, method_name))]

    # print("object_methods:",object_methods)

    # cmin, cmax = cb.get_clim()
    # if("tagged" in selection):
    #     print("Setting clim")
    #     plt.clim(0, vmaxAll)

    #ax.matshow(Values_array, cmap='seismic')

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

    xLabel, yLabel = GetPlotLabels(varLabel)

    plt.xlabel(xLabel, fontsize=25)
    plt.ylabel(yLabel, fontsize=25)

    Add_CMS_Header(plt, ax, upperRightText, text_xmin)

    plt.grid()
    plotText = plotText.replace("clean_", "")
    plt.text(
        0.1, 0.75, plotText,
        fontsize=30, fontweight='bold',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )

    # if(isRatio): plt.set_zlim(thisZmin,zmax_)
    
    ol = "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/"
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    fig.tight_layout()
    plt.savefig("{ol}/{varLabel}_{selection}.png".format(ol=ol, varLabel=varLabel, selection=selection))
    plt.savefig("{ol}/{varLabel}_{selection}.pdf".format(ol=ol, varLabel=varLabel, selection=selection))
    plt.close()    