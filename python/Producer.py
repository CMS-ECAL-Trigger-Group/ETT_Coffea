"""
HistProducer.py
Produce histograms using coffea.
"""

from coffea.hist import Hist, Bin, export1d, plot2d
from coffea.processor import ProcessorABC, LazyDataFrame, dict_accumulator
from uproot3 import recreate
import numpy as np
import awkward as ak 
import copy 

class HistProducer(ProcessorABC):
    """
    A coffea Processor which produces a histogram
    This applies selections 
    """

    histograms = NotImplemented
    selection = NotImplemented

    def __init__(self, dim = "", do_syst=False, syst_var='', severities=["all"], weight_syst=False, haddFileName=None, flag=False):
        self._flag = flag
        self.dim = int(dim)
        self.do_syst = do_syst
        self.syst_var, self.syst_suffix = (syst_var, f'_sys_{syst_var}') if do_syst and syst_var else ('', '')
        self.severities = severities
        self.weight_syst = weight_syst

        # 1d histograms
        if(self.dim == 1):

            self._accumulator = dict_accumulator({
                name: Hist('Events', Bin(name=name, **axis))
                for name, axis in ((self.naming_schema(hist['name'], region), hist['axis'])
                                for _, hist in list(self.histograms_1d.items())
                                for region in hist['region'])
            })

        # 2d histograms 
        elif(self.dim == 2):
            self._accumulator = dict_accumulator({
                name: Hist('Events', Bin(name = axes['xaxis']['label'], **axes['xaxis']), Bin(name = axes['yaxis']['label'], **axes['yaxis'])) ##-- Make it 2d by specifying two Binnings 
                for name, axes in ((self.naming_schema(hist['name'], region), hist['axes'])
                                for _, hist in list(self.histograms_2d.items())
                                for region in hist['region'])
            })        

        self.outfile = haddFileName

    def __repr__(self):
        return f'{self.__class__.__name__}(do_syst: {self.do_syst}, syst_var: {self.syst_var}, severities : {self.severities}, weight_syst: {self.weight_syst}, output: {self.outfile})'

    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, df, *args):
        output = self.accumulator.identity()

        # 1d histograms 
        if(self.dim == 1):
            for h, hist in list(self.histograms_1d.items()):
                variable_name = hist['name']
                for region in hist['region']:

                    name = self.naming_schema(hist['name'], region)
                    selec = self.passbut(df, hist['target'], region, variable_name)

                    selectedValues = np.hstack(ak.to_list(df[hist['target']][selec])).flatten()

                    output[name].fill(**{
                        name: selectedValues
                    })

                    del selectedValues               


        elif(self.dim == 2):
            for h, hist in list(self.histograms_2d.items()):
                variable_name = hist['name']
                for region in hist['region']:
                    print("region:",region)

                    name = self.naming_schema(hist['name'], region)
                    selec = self.passbut(df, hist['target_x'], region, variable_name) ##-- Should the selection depend on target?

                    xax_lab = hist['target_x']
                    yax_lab = hist['target_y']

                    xVals = np.hstack(ak.to_list(df[hist['target_x']][selec])).flatten()

                    # for this variable, seems I need to specify the computation
                    if((variable_name == "oneMinusEmuOverRealvstwrADCCourseBinning") or (variable_name == "oneMinusEmuOverRealvstwrADCCourseBinningZoomed")):
                        #print("min value:",min(ak.to_list(df["twrADC"]).flatten()))
                        yVals = np.hstack(ak.to_list((np.subtract(1., np.divide(df["twrEmul3ADC"], df["twrADC"])))[selec])).flatten() ## 'target_y' : '1 - (twrEmul3ADC/twrADC)',

                        # yVals = np.hstack(ak.to_list((np.subtract(1., np.divide(df["twrEmul3ADC"], df["twrADC"], out=np.zeros_like(df["twrEmul3ADC"]), where=df["twrADC"]!=0)))[selec])).flatten() ## 'target_y' : '1 - (twrEmul3ADC/twrADC)',
                        # np.divide(a, b, out=np.zeros_like(a), where=b!=0)
                        #yVals = np.hstack(ak.to_list((np.divide(df["twrEmul3ADC"], df["twrADC"]))[selec])).flatten() ## 'target_y' : '(twrEmul3ADC/twrADC)',

                    else:
                        yVals = np.hstack(ak.to_list(df[hist['target_y']][selec])).flatten() # read from target_y 

                    output[name].fill(**{
                        xax_lab : xVals,
                        yax_lab : yVals 
                    }
                    )

                    del xVals
                    del yVals 
                    del name 
                    del selec 
                    del xax_lab
                    del yax_lab            

        return output

    def postprocess(self, accumulator):
        return accumulator

    def passbut(self, event: LazyDataFrame, excut: str, cat: str, variable_name: str):
        """Backwards-compatible passbut."""
        if( (variable_name == "oneMinusEmuOverRealvstwrADCCourseBinning") or (variable_name == "oneMinusEmuOverRealvstwrADCCourseBinningZoomed")): # add selection to not include TPs with ET = 0 
            evalStr = '&'.join('(' + cut.format(sys=('' if self.weight_syst else self.syst_suffix)) + ')' for cut in self.selection[cat] )
            evalStr += '&(event.twrADC != 0)'
            return eval(evalStr)#if excut not in cut))
        else: 
            return eval('&'.join('(' + cut.format(sys=('' if self.weight_syst else self.syst_suffix)) + ')' for cut in self.selection[cat] ))#if excut not in cut))

class ETT_NTuple(HistProducer):

    addRunSelection = 0 
    SelectionsToRun = []

    sevDict = {
        "zero" : "0",
        "three" : "3",
        "four" : "4"
    }

    severities = ["all", "zero", "three", "four"]
    times = ["all", "inTime", "Early", "Late", "VeryLate"]
    # times = ["all"]
    FGSelections = ["all", "Tagged"] # all: all TPs. Tagged: FGbit=1

    for severity in severities:
        for time in times:
            for FGSel in FGSelections:
                if(time == "inTime" and severity == "three"): continue ##-- skip in time severity 3 
                selname = "sev%s_%s_%s"%(severity, time, FGSel)
                SelectionsToRun.append(selname)   

    # 1d histograms 
    histograms_1d = {

        ##-- 1d histograms 
        # 'time': {
        #     'target': 'time',
        #     'name': 'time', 
        #     'region' : ['sevzero_all_', 'sevzero_MostlyZeroed',
        #                 'sevthree_all', 'sevthree_MostlyZeroed',
        #                 'sevfour_all', 'sevfour_MostlyZeroed'
        #                ],

        #     # 'region' : [
        #                 # 'sevzero_all_twrADC1to2', 'sevzero_MostlyZeroed_twrADC1to2',
        #                 # 'sevzero_all_twrADC2to32', 'sevzero_MostlyZeroed_twrADC2to32',
        #                 # 'sevzero_all_twrADC32to256', 'sevzero_MostlyZeroed_twrADC32to256',

        #                 # 'sevthree_all_twrADC1to2', 'sevthree_MostlyZeroed_twrADC1to2',
        #                 # 'sevthree_all_twrADC2to32', 'sevthree_MostlyZeroed_twrADC2to32',
        #                 # 'sevthree_all_twrADC32to256', 'sevthree_MostlyZeroed_twrADC32to256',

        #                 # 'sevfour_all_twrADC1to2', 'sevfour_MostlyZeroed_twrADC1to2',
        #                 # 'sevfour_all_twrADC2to32', 'sevfour_MostlyZeroed_twrADC2to32',
        #                 # 'sevfour_all_twrADC32to256', 'sevfour_MostlyZeroed_twrADC32to256',                                                

        #             #    ],

        #     # 'region': ['sevall_all', 'sevall_MostlyZeroed', 'sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
        #     # 'region': ['sevzero_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevfour_MostlyZeroed'],
        #     # 'axis': {'label': 'time', 'n_or_arr': 120, 'lo': -225, 'hi': 125}
        #     # 'axis': {'label': 'time', 'n_or_arr': 120, 'lo': -225, 'hi': 125}
        #     'axis': {'label': 'time', 'n_or_arr': 100, 'lo': -50, 'hi': 50}
        # }, 

        # 'twrADC': {
        #     'target': 'twrADC',
        #     'name': 'twrADC', 
        #     'region' : ['sevzero_all', 'sevzero_MostlyZeroed',
        #                 'sevthree_all', 'sevthree_MostlyZeroed',
        #                 'sevfour_all', 'sevfour_MostlyZeroed'
        #                ],
        #     'axis': {'label': 'twrADC', 'n_or_arr': 255, 'lo': 1, 'hi': 256}
        # }, 

        # 'twrEmul3ADC': {
        #     'target': 'twrEmul3ADC',
        #     'name': 'twrEmul3ADC', 
        #     'region' : ['sevzero_all', 'sevzero_MostlyZeroed',
        #                 'sevthree_all', 'sevthree_MostlyZeroed',
        #                 'sevfour_all', 'sevfour_MostlyZeroed'
        #                ],
        #     'axis': {'label': 'twrEmul3ADC', 'n_or_arr': 255, 'lo': 1, 'hi': 256}
        # },                 

        'twrADC': {
            'target': 'twrADC',
            'name': 'twrADC', 
            # 'region': ['sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
            # 'region': ['sevall_all', 'sevall_MostlyZeroed', 'sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
            'region' : SelectionsToRun,
            'axis': {'label': 'twrADC', 'n_or_arr': 256, 'lo': 0, 'hi': 256}
        }, 

        # 'twrEmul3ADC': {
        #     'target': 'twrADC',
        #     'name': 'twrADC', 
        #     # 'region': ['sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
        #     'region': ['sevall_all', 'sevall_MostlyZeroed', 'sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
        #     'axis': {'label': 'twrADC', 'n_or_arr': 256, 'lo': 0, 'hi': 256}
        # },                 

        # 'time': {
        #     'target': 'time',
        #     'name': 'time', 
        #     'region': ['sevzero_all', 'sevthree_all', 'sevfour_all', 'sevzero_MostlyZeroed', 'sevthree_MostlyZeroed', 'sevfour_MostlyZeroed'],
        #     'axes' : {
        #         'xaxis': {'label': 'time', 'n_or_arr': 120, 'lo': -225, 'hi': 125},
        #         'yaxis': {'label': 'time', 'n_or_arr': 120, 'lo': -225, 'hi': 125}
        #     }

        # },

    }

    # 2d histograms 
    histograms_2d = {

        ##-- emu/real 
        # 'emuOverRealvstwrADC': {
        #     'target_x' : 'twrADC',
        #     # 'target_y' : '(twrEmul3ADC/twrADC)',
        #     'target_y' : 'emuOverReal',
        #     'name': 'emuOverRealvstwrADC', 
        #     'region' : ['sevzero_all',
        #                 'sevthree_all',
        #                 'sevfour_all'
        #                ],
        #     'axes' : {
        #         'xaxis': {'label': 'twrADC', 'n_or_arr': 255, 'lo': 1, 'hi': 256},
        #         'yaxis': {'label': 'emuOverReal', 'n_or_arr': 48, 'lo': 0, 'hi': 1.2}                
        #     }

        # },  


        # """
        # Useful for exploring the data
        # """

        # 'EnergyVsTimeOccupancy': {
        #     # 'target': { 'x': 'twrADC', 'y' : 'twrEmul3ADC'},
        #     'target_x' : 'time',
        #     'target_y' : 'twrADC',
        #     'name': 'EnergyVsTimeOccupancy', 
        #     # 'region' : ['clean_all', 'clean_tagged'],
        #     'region' : SelectionsToRun,
        #     'axes' : {
        #         'xaxis': {'label': 'time', 'n_or_arr': 100, 'lo': -50, 'hi': 50},
        #         'yaxis': {'label': 'twrADC', 'n_or_arr': 256, 'lo': 0, 'hi': 256} # Full ET range           
        #         # 'yaxis': {'label': 'twrADC', 'n_or_arr': 35, 'lo': 0, 'hi' : 35} # Low ET range                 
        #     }

        # },  

        # """
        # Useful for evaluating effect of emulation, e.g. double weights, on digis
        # """

        # ##-- 1 - emu/real 
        # 'oneMinusEmuOverRealvstwrADCCourseBinning': {
        #     'target_x' : 'twrADC',
        #     # 'target_y' : '(twrEmul3ADC/twrADC)',
        #     'target_y' : 'oneMinusEmuOverRealvstwrADCCourseBinning',
        #     'name': 'oneMinusEmuOverRealvstwrADCCourseBinning', 
        #     'region' : SelectionsToRun, 
        #     'axes' : {
        #         # 'xaxis': {'label': 'twrADC', 'n_or_arr': 255, 'lo': 1, 'hi': 256}, 
        #         'xaxis': {'label': 'twrADC', 'n_or_arr': [1.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0, 80.0, 88.0, 96.0, 104.0, 112.0, 150.0, 256.0], 'lo': 1, 'hi': 256}, 
        #         # 'xaxis': {'label': 'twrADC', 'n_or_arr': [1.0, 8.0, 16.0, 24.0, 32.0, 40.0], 'lo': 1, 'hi': 40},
        #         # 'xaxis': {'label': 'twrADC', 'n_or_arr': 40, 'lo': 1, 'hi': 41}, # Low ET range 
        #         'yaxis': {'label': 'oneMinusEmuOverRealvstwrADCCourseBinning', 'n_or_arr': 48, 'lo': 0, 'hi': 1.2}                
        #         # 'yaxis': {'label': 'oneMinusEmuOverRealvstwrADCCourseBinning', 'n_or_arr': 88, 'lo': -1, 'hi': 1.2}                
        #         # 'yaxis': {'label': 'oneMinusEmuOverRealvstwrADCCourseBinning', 'n_or_arr': 128, 'lo': -2, 'hi': 1.2}                
        #         # 'yaxis': {'label': 'oneMinusEmuOverRealvstwrADCCourseBinning', 'n_or_arr': 448, 'lo': -10, 'hi': 1.2}    ##-- 0.025 space bins in y axis             
        #     }

        # },  

        # ##-- 1 - emu/real zoomed
        # 'oneMinusEmuOverRealvstwrADCCourseBinningZoomed': {
        #     'target_x' : 'twrADC',
        #     # 'target_y' : '(twrEmul3ADC/twrADC)',
        #     'target_y' : 'oneMinusEmuOverRealvstwrADCCourseBinningZoomed',
        #     'name': 'oneMinusEmuOverRealvstwrADCCourseBinningZoomed', 
        #     'region' : SelectionsToRun, 
        #     'axes' : {
        #         'xaxis': {'label': 'twrADC', 'n_or_arr': 40, 'lo': 1, 'hi': 41}, # Low ET range 
        #         'yaxis': {'label': 'oneMinusEmuOverRealvstwrADCCourseBinningZoomed', 'n_or_arr': 48, 'lo': 0, 'hi': 1.2}                
        #     }

        # },          

        # 'EBOcc': {
        #     'target_x' : 'iphi',
        #     'target_y' : 'ieta',
        #     'name': 'EBOcc', 
        #     'region' : ['all'],
        #     # 'region' : ['sevzero_all', 'sevzero_MostlyZeroed',
        #     #             'sevthree_all', 'sevthree_MostlyZeroed',
        #     #             'sevfour_all', 'sevfour_MostlyZeroed'
        #     #            ],
        #     'axes' : {
        #         'xaxis': {'label': 'iphi', 'n_or_arr': 80, 'lo': 0, 'hi': 80},
        #         'yaxis': {'label': 'ieta', 'n_or_arr': 36, 'lo': -18, 'hi': 18}                
        #     }

        # },          

        'realVsEmu': {
            # 'target': { 'x': 'twrADC', 'y' : 'twrEmul3ADC'},
            'target_x' : 'twrEmul3ADC',
            'target_y' : 'twrADC',
            'name': 'realVsEmu', 
            'region' : SelectionsToRun,
            'axes' : {
                'xaxis': {'label': 'twrEmul3ADC', 'n_or_arr': 256, 'lo': 0, 'hi': 256},
                'yaxis': {'label': 'twrADC', 'n_or_arr': 256, 'lo': 0, 'hi': 256}
                # 'xaxis': {'label': 'twrEmul3ADC', 'n_or_arr': 133, 'lo': 0, 'hi': 256},
                # 'yaxis': {'label': 'twrADC', 'n_or_arr': 133, 'lo': 0, 'hi': 256}                
            }
        },  

    }

    time_regions = {
        "all" : ["event.time != -999"],
        "inTime" : ["event.time < 3", "event.time >= -3"],
        "Early" : ["event.time < -3"],
        "Late" : ["event.time >= 3", "event.time < 10"],
        "VeryLate" : ["event.time >= 10"]
    }

    severities = ["all", "zero", "three", "four"]

    sevDict = {
        "all" : "all",
        "zero" : "0",
        "three" : "3",
        "four" : "4"
    }
    times = ["all", "inTime", "Early", "Late", "VeryLate"]
    # times = ["all"]
    FGSelections = ["all", "Tagged"] # all: all TPs. Tagged: FGbit=1

    # Define selections based on event values 
    selection = {}

    # basic ETT TP cleaning 
    """
    selection = {
        "clean_all" : ["event.time != -999", # TP matched to (highest energy) rec hit in tower 
                   "event.ttFlag != 4", # TP does not have TTF4. TTF4 is usually a masked or problematic TP 
        ],
        "clean_tagged" : [
            "event.time != -999",
            "event.ttFlag != 4"
            #"event.FineGrainBit == 1" # tagged by double weights in tagging mode 
        ]
    }
    """

    # based on severity and time 
    #"""
    selection = {}

    for severity in severities:
        sevNum = sevDict[severity]
        for time in times:
            for FGSel in FGSelections: 
                selec_selections = []
                if(FGSel == "Tagged"):
                    selec_selections.append("event.FineGrainBit == 1")
                selec_selections.append("event.ttFlag != 4") # add TTF4 cleaning for all severities, times. 
                selec_selections.append("event.time != -999") # don't include TPs which were not matched to a recHit 
                #selec_selections.append("event.twrADC != 0") # for anything with dividing, don't include TP = 0 ET entries
                if(sevNum != "all"): # if there is a sevnum, add the selection. If there is not, add no selection (be inclusive of all severities)
                    selec_selections.append("event.sevlv == %s"%(sevNum))
                else:
                    selec_selections.append("event.sevlv != -999")
                # if(addRunSelection): 
                #     print("Adding run selection")
                #     selec_selections.append("event.runNb == 324725")
                timeSels = time_regions[time]
                for timeSel in timeSels:
                    selec_selections.append(timeSel)
                Selection_Name = "sev%s_%s_%s"%(severity, time, FGSel)
                selection[Selection_Name] = (selec_selections)
    #""" 

    def weighting(self, event: LazyDataFrame):
        weight = 1.0
        try:
            weight = event.xsecscale
        except:
            return "ERROR: weight branch doesn't exist"
        return weight

    def naming_schema(self, name, region):
        return f'{name}_{region}{self.syst_suffix}'
