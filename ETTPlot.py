from python.ETTPlot_Tools import * 
import uproot 

"""
Abraham Tishelman-Charny
17 November 2021

The purpose of this notebook is to plot quantities from ETTAnalyzer outputs.
"""

import matplotlib
import pickle
from matplotlib import pyplot as plt 
import numpy as np
import os 
from matplotlib.colors import LogNorm

##-- Real vs. Emu 
#realVsEmu_sevall_all_values_43.p
isWide = 0
varLabel = "realVsEmu"

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
"""



plt.xlabel("Emulated TP Et (ADC)", fontsize=20)
plt.ylabel("Real data TP Et (ADC)", fontsize=20)
#     plt.xscale('log')

# cb = fig.colorbar(pos, 
#                 ax=ax,
#             )

upperRightText = "Pilot Beam 2021"
xmin = 0.115
Add_CMS_Header(plt, isWide, ax, upperRightText, xmin)

#     plt.text(
#         0.4, 6., r"Emu < Real",
#         fontsize=20, fontweight='bold',
#         horizontalalignment='left',
#         verticalalignment='bottom',
#         transform=ax.transAxes
#     )

# plt.grid()
# plt.show()
# ol = "/eos/user/a/atishelm/www/EcalL1Optimization/ZeroBias/"
# plt.savefig("%s/EnergyVsTimePercentTagged_sev%s_2d.png"%(ol, sevLabel))
# plt.savefig("%s/EnergyVsTimePercentTagged_sev%s_2d.pdf"%(ol, sevLabel))

plt.savefig("RealVsEmu.png")
plt.close()

print("DONE")    