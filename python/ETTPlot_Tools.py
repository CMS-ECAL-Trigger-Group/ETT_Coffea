"""
17 November 2021
Abraham Tishelman-Charny 

The purpose of this module is to provide tools for ETTPlot. 
"""

##-- CMS header 
def Add_CMS_Header(plt, isWide, ax, upperRightText, xmin):
    ##-- Upper left plot text
    ##-- CMS 
    plt.text(
        # 0.05, 0.9, u"CMS $\it{Preliminary}$",
        0., 1., u"CMS ",
        fontsize=20, fontweight='bold',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )

    prelim_x = xmin
    
#     if(isWide):
#         prelim_x = 0.08
#     else:
#         prelim_x = 0.115
#         prelim_x = 0.135
    
    ##-- Preliminary 
    plt.text(
#         prelim_x, 0.998, u"$\it{Simulation}$ $\it{Preliminary}$",
        prelim_x, 0.998, u"$\it{Preliminary}$",
        fontsize=19,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )    

    ##-- Lumi 
    plt.text(
#         1., 1., r"%s fb$^{-1}$ (13 TeV)"%(lumi),
#         1., 1., "(13 TeV)",
#         1., 1., "(14 TeV)",
        1., 1., upperRightText,
        fontsize=16, horizontalalignment='right', 
        verticalalignment='bottom', 
        transform=ax.transAxes
    )  

# for changing the axes of plots made from pickled arrays 
def SetYMin(Values_array_, ymin_):
    binWidth = 0.025 
    Transformed_List = []
    ymin_bin_i = (1./binWidth) * ymin_ + 10./binWidth 
    for i, xBinVals in enumerate(Values_array_):
        nVals = len(xBinVals)
        MASK = [i > ((ymin_bin_i) - 1) for i, val in enumerate(xBinVals)]
        xBinVals = xBinVals[MASK]
        Transformed_List.append(xBinVals)
    Transformed_List = np.array(Transformed_List)
    return Transformed_List


def SetYMax():
    binWidth = 0.    