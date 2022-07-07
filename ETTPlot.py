"""
Abraham Tishelman-Charny
17 November 2021

The purpose of this module is to plot quantities from ETT_Coffea outputs.

Example commands:

conda activate higgs-dna # to create parquet files 

# Comparison of 2021 and 2022 900 GeV collisions

python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EmulEnergyVsTimeOccupancy --maxFiles 999999 --times all,inTime,VeryLate --severities zero,four --plotIndividuals --FGSelections all,EmulTagged
python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EBOcc --maxFiles 999999 --times all --severities all --plotIndividuals --FGSelections all,DataNotEqEmul
python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EBOcc --maxFiles 10 --times all --severities all --plotIndividuals --FGSelections all,Tagged,DataNotEqEmul
python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EnergyVsTimeOccupancy --maxFiles 999999 --times all,Early,inTime,Late,VeryLate --severities all,zero,three,four --plotIndividuals
python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EnergyVsTimeOccupancy --maxFiles 999999 --times all,Early,inTime,Late,VeryLate --severities all,zero,three,four --plotTogether --fromParquet

python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EnergyVsTimeOccupancy --maxFiles 999999 --times all --severities zero,three,four --plotIndividuals
python3 ETTPlot.py --dataset 2021_2022_900GeVCollisions --variables EnergyVsTimeOccupancy --maxFiles 100 --times all --severities zero  --plotTogether --fromParquet

# 2021 Collisions
python3 ETTPlot.py --dataset PilotBeam2021 --variables EnergyVsTimeOccupancy --maxFiles 1412 --times all --severities all,zero,three,four --plotIndividuals

# 2022 Collisions 
python3 ETTPlot.py --dataset Run352912 --variables EnergyVsTimeOccupancy,realVsEmu --maxFiles 3077 --times all,Early,inTime,Late,VeryLate --severities all,zero,three,four --plotIndividuals

# Full readout data from 2017/2018 
python3 ETTPlot.py --dataset FullReadoutData_2018 --variables oneMinusEmuOverRealvstwrADCCourseBinningZoomed  --maxFiles 1000000  --plotTogether --fromParquet --times inTime --severities zero
python3 ETTPlot.py --dataset FullReadoutData_2018 --variables EnergyVsTimeOccupancy,oneMinusEmuOverRealvstwrADCCourseBinning  --maxFiles 1000000  --plotTogether --fromParquet
python3 ETTPlot.py --dataset FullReadoutData_2018 --variables EnergyVsTimeOccupancy,oneMinusEmuOverRealvstwrADCCourseBinning  --maxFiles 1000000  --plotIndividual
python3 ETTPlot.py --variables oneMinusEmuOverRealvstwrADCCourseBinning --severities zero --times inTime --maxFiles 10 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2017_2018/" --plotIndividuals
python3 ETTPlot.py --dataset FullReadoutData_2017_2018 --directory --variables oneMinusEmuOverRealvstwrADCCourseBinning,EnergyVsTimeOccupancy --severities zero --times inTime --maxFiles 100 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2017_2018/MinDelta2p5prime_WithOddPF_WeightsReco/" --plotIndividuals

# 2021 Pilot Beam
python3 ETTPlot.py --variables oneMinusEmuOverRealvstwrADCCourseBinning  --plotTogether --fromParquet

python3 ETTPlot.py --variables oneMinusEmuOverRealvstwrADCCourseBinning --severities zero --times inTime --maxFiles 10000 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_MultiFitReco/" --plotIndividuals
python3 ETTPlot.py --variables oneMinusEmuOverRealvstwrADCCourseBinning --severities zero --times inTime --maxFiles 100000 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_MultiFitReco/" --plotTogether --fromParquet

python3 ETTPlot.py --directory /eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_094402/all_output/ --variables oneMinusEmuOverRealvstwrADCCourseBinning,EnergyVsTimeOccupancy --severities zero --times inTime --maxFiles 100 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_MultiFitReco/" --plotIndividuals --plotTogether --fromParquet

python3 ETTPlot.py --directory /eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_104615/all_output/ --variables oneMinusEmuOverRealvstwrADCCourseBinning,EnergyVsTimeOccupancy --severities zero --times inTime --maxFiles 10000 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_WeightsReco/"
python3 ETTPlot.py --directory /eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_094402/all_output/ --variables oneMinusEmuOverRealvstwrADCCourseBinning,EnergyVsTimeOccupancy --severities zero --times inTime --maxFiles 10000 --outputLocation "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_MultiFitReco/"
"""

import numpy as np
import os 
import argparse 
import uproot 
import pickle
from matplotlib import pyplot as plt 

from python.ETTPlot_Tools import * 

import pyarrow as pa
import pyarrow.parquet as pq

parser = argparse.ArgumentParser()
parser.add_argument("--variables", type = str, default = "RealVsEmu", help = "Comma separated list of variables to plot") 
parser.add_argument("--times", type = str, default = "all,Early,inTime,Late,VeryLate", help = "Comma separated list of times to plot") 
parser.add_argument("--severities", type = str, default = "zero,three,four", help = "Comma separated list of severities to plot") 
parser.add_argument("--FGSelections", type = str, default = "all,Tagged", help = "Comma separated list of FGBit scenarios to plot") 
parser.add_argument("--maxFiles", type = int, default = 1, help = "Max number of files to process")
#parser.add_argument("--outputLocation", type = str, required = True, help = "Output directory for beautiful plots") 
parser.add_argument("--dataset", type = str, required = True, help = "Dataset that was run over. Options: FullReadoutData_2017_2018, PilotBeam2021") 
parser.add_argument("--plotIndividuals", action="store_true", default = False, help = "Produce a plot for each selection")
parser.add_argument("--plotTogether", action="store_true", default = False, help = "Produce a plot for each selection")
parser.add_argument("--fromParquet", action="store_true", default = False, help = "Plot using parquet files")

args = parser.parse_args()

variables = args.variables.split(',')
times = args.times.split(',')
severities = args.severities.split(',')
FGSelections = args.FGSelections.split(',')
maxFiles = args.maxFiles 
dataset = args.dataset
plotIndividuals = args.plotIndividuals 
plotTogether = args.plotTogether
fromParquet = args.fromParquet

N_FGSels = len(FGSelections)

ol = "/eos/user/a/atishelm/www/EcalL1Optimization/{dataset}/".format(dataset=dataset)

if(dataset == "FullReadoutData_2018"):
    directories = [
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220214_092624/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220214_122937/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220214_122635/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220214_122824/all_output/",
    ]

    direc_ol_dict = {
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220214_092624/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2018/MinDelta2p5prime_WithOddPF_WeightsReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220214_122937/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2018/MinDelta2p5prime_WithOddPF_MultiFitReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220214_122635/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2018/MinDelta0p5prime_WithOddPF_WeightsReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220214_122824/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/FullReadoutData_2018/MinDelta0p5prime_WithOddPF_MultiFitReco/",
    }

# My name is Abe and I like making hacky python code to serve my short-term needs   
elif(dataset == "PilotBeam2021"):
    directories = [
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights_ReemulateFromGlobalTag/220621_091715/allFiles_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_0p5PrimeODDweights/220210_104023/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_2p5PrimeODDweights/220210_103954/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_104615/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_094402/all_output/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220209_125921/all_output/"
    ]

    direc_ol_dict = {
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights_ReemulateFromGlobalTag/220621_091715/allFiles_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/2021Collisions/Runs_346446_346447/Reemul_by_globalTag/"
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_0p5PrimeODDweights/220210_104023/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithoutOddPF_MultiFitReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_2p5PrimeODDweights/220210_103954/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta2p5prime_WithoutOddPF_MultiFitReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_104615/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_WeightsReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/220210_094402/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta0p5prime_WithOddPF_MultiFitReco/",
        # "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/220209_125921/all_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/PilotBeam2021/MinDelta2p5prime_WithOddPF_MultiFitReco/",    
    }

elif(dataset == "Run352912"):
    directories = [
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_352912/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights/220615_220151/allFiles_output/"
    ]

    direc_ol_dict = {
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_352912/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights/220615_220151/allFiles_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/2022Collisions/Run_352912/Reemul_by_globalTag/"
    }

elif(dataset == "2021_2022_900GeVCollisions"):
    directories = [
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_352912/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights/220615_220151/allFiles_output/",
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights_ReemulateFromGlobalTag/220621_091715/allFiles_output/"
    ]

    direc_ol_dict = {
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_352912/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights/220615_220151/allFiles_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/2022Collisions/Run_352912/Reemul_by_globalTag/",
        "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_3_0_DoubleWeights_ReemulateFromGlobalTag/220621_091715/allFiles_output/" : "/eos/user/a/atishelm/www/EcalL1Optimization/2021Collisions/Runs_346446_346447/Reemul_by_globalTag/"        
    }    

upperRightTextDict = {
    "PilotBeam2021" : ["Runs 346446, 346447", "82", "{\mu}b", "900 GeV"], # label, lumi 
    "FullReadoutData_2018" : ["FR 2018", "0.014"], # label, lumi 
    "Run352912" : ["Run 352912", "145", "{\mu}b", "900 GeV"]
}

if(__name__ == '__main__'):

    varLabelDict = {
        "realVsEmu": "realVsEmu",
        "EnergyVsTimeOccupancy" : "EnergyVsTimeOccupancy",
        "oneMinusEmuOverRealvstwrADCCourseBinning" : "oneMinusEmuOverRealvstwrADCCourseBinning",
        "oneMinusEmuOverRealvstwrADCCourseBinningZoomed" : "oneMinusEmuOverRealvstwrADCCourseBinningZoomed",
        "EBOcc" : "EBOcc",
        "EmulEnergyVsTimeOccupancy" : "EmulEnergyVsTimeOccupancy"
    }

    # Make output directory if it doesn't exist 
    mkdirCmd = "mkdir -p {ol}".format(ol=ol)
    cpIndCmd = "cp {ol}/../index.php {ol}".format(ol=ol)
    print("$",mkdirCmd)
    os.system(mkdirCmd)
    print("$",cpIndCmd)
    os.system(cpIndCmd)

    # go through variables to plot 
    for variable in variables:
        varLabel = varLabelDict[variable]
        if(plotIndividuals):
            for direc in directories:
                print("direc:",direc)
                if("Run_352912" in direc or "PilotBeam_2021" in direc):
                    WP, PF, RECO = "GT", "GT", "Weights"
                else: 
                    WP, PF, RECO = GetWorkingPointLabels(direc) # if varying working points 

                direc_ol = direc_ol_dict[direc]

                if("Run_352912" in direc):
                    sub_dataset = "Run352912"
                elif("Runs_346446_346447_PilotBeam_2021" in direc):
                    sub_dataset = "PilotBeam2021"
                else: sub_dataset = dataset 

                if(sub_dataset == "PilotBeam2021"): 
                    print("Skipping ",sub_dataset)
                    continue 

                for severity in severities:
                    print("On severity:",severity)
                    for time in times:
                        print("On time:",time)
                        if(time == "inTime" and severity == "three"):
                            print("Skipping in time severity three as they shouldn't exist by definition")
                            continue 

                        for FGSelection in FGSelections:
                            print("FGSelection:",FGSelection)
                            print("On Var: %s, Sev: %s, time: %s"%(variable, severity, time))
                            thisDirec = "%s/%s/%s/%s/"%(direc, variable, severity, time) # directory for a given variable, severity, time 
                            files = [f for f in os.listdir(thisDirec)]
                            n_files = len(files)
                            values_0 = pickle.load(open("%s/%s_sev%s_%s_%s_values_0.p"%(thisDirec, variable, severity, time, FGSelection), "rb")) # EnergyVsTimeOccupancy_sevall_all_all_values_3040.p
                            values_1 = pickle.load(open("%s/%s_sev%s_%s_%s_values_1.p"%(thisDirec, variable, severity, time, FGSelection), "rb"))
                            total_values = np.add(values_0, values_1)

                            if(maxFiles == 999999):
                                maxFilesDict = {
                                    "Run352912" : 3077,
                                    "PilotBeam2021" : 1412                                    
                                }
                                maxFiles = maxFilesDict[sub_dataset]
                                print("maxFiles:",maxFiles)

                            for i in range(n_files):
                                if(i > maxFiles):
                                    print("Max files reached: ",i)
                                    break 
                                if(i == 0 or i == 1): continue 
                                if(i%100 == 0): print("on file %s / %s"%(i, n_files))

                                exec('thisValues = pickle.load(open("%s/%s_sev%s_%s_%s_values_%s.p"%(thisDirec, variable, severity, time, FGSelection, i), "rb"))')  
                                total_values = np.add(total_values, thisValues)

                            maxFiles = args.maxFiles # reset maxFiles value 

                            exec("%s_%s_%s_values = np.copy(total_values)"%(variable, severity, time)) # copy with unique name to plot together on same plot later 
                            exec("These_Values = np.copy(total_values)") 
                            os.system("mkdir -p %s"%(direc_ol))
                            os.system("cp %s/../index.php %s"%(direc_ol, direc_ol))

                            upperRightText, lumi, unit, sqrts = upperRightTextDict[sub_dataset]
                            averages, stdevs = MakeETTPlot(These_Values, variable, severity, time, direc_ol, upperRightText, sub_dataset, lumi, unit, sqrts, FGSelection) # make plots and return averages 

                            exec("%s_Values = np.copy(These_Values)"%(FGSelection))           

                            # averageComputingVars = ["oneMinusEmuOverRealvstwrADCCourseBinning", "oneMinusEmuOverRealvstwrADCCourseBinningZoomed", "EmulEnergyVsTimeOccupancy"]                     
                            averageComputingVars = ["oneMinusEmuOverRealvstwrADCCourseBinning", "oneMinusEmuOverRealvstwrADCCourseBinningZoomed"]                     

                            if(variable in averageComputingVars):
                                # save averages and stdevs as parquet files (choosing parquet for fun / experience) to combine everything later 
                                # save averages for each case: var, WP, PF, RECO, Sev, Time

                                parquetOutDir = "output/%s/%s_%s_%s_%s/"%(sub_dataset, variable, WP, PF, RECO)
                                os.system("mkdir -p %s"%(parquetOutDir))
                                outName_averages = "%s/%s_%s_%s_%s_%s_%s_averages.parquet"%(parquetOutDir, variable, WP, PF, RECO, severity, time)
                                outName_stdevs = "%s/%s_%s_%s_%s_%s_%s_stdevs.parquet"%(parquetOutDir, variable, WP, PF, RECO, severity, time)

                                averages_table = pa.table({"data": averages})
                                stdevs_table = pa.table({"data": stdevs})
                                print("Saving parquet file:",outName_averages)   
                                print("Saving parquet file:",outName_stdevs)                   
                                pa.parquet.write_table(averages_table, outName_averages)                
                                pa.parquet.write_table(stdevs_table, outName_stdevs)    

                                exec("%s_%s_averages = np.copy(averages)"%(severity, time))
                                exec("%s_%s_stdevs = np.copy(stdevs)"%(severity, time))

                        if(variable != "EBOcc" and variable != "realVsEmu"):
                            if(variable == "EmulEnergyVsTimeOccupancy"): # emulated 
                                exec("Tagged_Values = np.copy(EmulTagged_Values)")
                            # Make ratio of tagged / all to see which region of phase space is tagged by double weights 
                            fraction = np.divide(Tagged_Values, all_Values, out=np.zeros_like(Tagged_Values), where=all_Values!=0)
                            averages, stdevs = MakeETTPlot(fraction, variable, severity, time, direc_ol, upperRightText, sub_dataset, lumi, unit, sqrts, "ratio") # make plots and return averages 

                        # If variable is energy vs time, compute Tagging probability vs. ET 
                        if(variable == "EnergyVsTimeOccupancy" or variable == "EmulEnergyVsTimeOccupancy"):
                            # Compute ratio of tagged / all per ET bin (integrate over times)
                            # integrate over time 

                            Total_Tagged_PerEnergyBin = np.sum(Tagged_Values, axis=0)
                            Total_PerEnergyBin = np.sum(all_Values, axis=0)

                            rebin = 1 
                            debug = 0

                            # if rebinning, need to manipulate arrays 
                            xbins, ybins = GetBins(variable, sub_dataset)
                            energy_bins = ybins             
                            # combine bins (0 indexed) (1 -> 9), (10 -> 20), (21 -> 30), (31 -> 40) 

                            if(severity == "four"):
                                reBins = [[1,32],[33,255]]
                                # reBins = [[1,21], [22,23], [24,32], [33,255]]
                            else: 
                                reBins = [[5,32],[33,255]]

                            if(debug):
                                print("Before rebin tagged:",Total_Tagged_PerEnergyBin)
                                print("Before rebin total:",Total_PerEnergyBin)
                                print("Bin edges before rebin:",energy_bins)
                                print("length:",len(Total_Tagged_PerEnergyBin))
                                print("length:",len(Total_PerEnergyBin))
                                print("length:",len(energy_bins))

                            # reBins = [[1,9], [10,20]] # inclusive bin ranges to combine, needs to be in order (i.e. cannot have [10, 20] before [1,9])          
                            # reBins = [[1,50]] # inclusive bin ranges to combine, needs to be in order (i.e. cannot have [10, 20] before [1,9])          
                            if(rebin):

                                # begin with copies of original 
                                Total_Tagged_PerEnergyBin_rebinned = np.copy(Total_Tagged_PerEnergyBin)
                                Total_PerEnergyBin_rebinned = np.copy(Total_PerEnergyBin)
                                energy_bins_rebinned = np.copy(energy_bins)

                                # keep track of number of elements removed to shift indicies 
                                binsRemoved = 0
                                binEdgesRemoved = 0 

                                # for each new rebin to make 
                                for newBin in reBins:
                                    newBinMin_i, newBinMax_i = newBin[0], newBin[1] # indices of range to consolidate into one bin 

                                    # indicies to sum 
                                    sumIndicesArr = [True if (i >= (newBinMin_i) and i <= (newBinMax_i)) else False for i, val in enumerate(Total_Tagged_PerEnergyBin)]
                                    Tagged_newBinSum = np.sum(Total_Tagged_PerEnergyBin, where=sumIndicesArr) # get sum from original array so you don't need to worry about already removed elements 
                                    Total_newBinSum = np.sum(Total_PerEnergyBin, where=sumIndicesArr)

                                    # replace first value in summed vals with sum of those vals
                                    Total_Tagged_PerEnergyBin_rebinned[newBinMin_i - binsRemoved] = Tagged_newBinSum
                                    Total_PerEnergyBin_rebinned[newBinMin_i - binsRemoved] = Total_newBinSum

                                    # remove subsequent vals since they're included in sum
                                    for b_i in range(newBinMin_i, newBinMax_i):
                                        bin_i_toRemove = b_i + 1 - binsRemoved
                                        Total_Tagged_PerEnergyBin_rebinned = np.delete(Total_Tagged_PerEnergyBin_rebinned, bin_i_toRemove)
                                        Total_PerEnergyBin_rebinned = np.delete(Total_PerEnergyBin_rebinned, bin_i_toRemove)
                                        binsRemoved += 1 

                                    # for bin edges, need to remove middle edges not including min and max
                                    for b_i in range(newBinMin_i, newBinMax_i):
                                        bin_i_toRemove = b_i + 1 - binEdgesRemoved
                                        energy_bins_rebinned = np.delete(energy_bins_rebinned, bin_i_toRemove) 
                                        binEdgesRemoved += 1 
                                        if(b_i == newBinMax_i): continue # don't delete upper edge of rebin range 

                                # update original arrays to rebinned versions
                                Total_Tagged_PerEnergyBin = np.copy(Total_Tagged_PerEnergyBin_rebinned)
                                Total_PerEnergyBin = np.copy(Total_PerEnergyBin_rebinned)
                                energy_bins = np.copy(energy_bins_rebinned)

                            if(debug):
                                print("After rebin tagged:",Total_Tagged_PerEnergyBin)
                                print("After rebin total:",Total_PerEnergyBin)
                                print("Bin edges after rebin:",energy_bins)
                                print("length:",len(Total_Tagged_PerEnergyBin))
                                print("length:",len(Total_PerEnergyBin))
                                print("length:",len(energy_bins))

                            NOENTRIESMASK = tuple([Total_PerEnergyBin != 0])
                            TaggingProbPerEnergyBin = np.divide(Total_Tagged_PerEnergyBin, Total_PerEnergyBin, out=np.zeros_like(Total_Tagged_PerEnergyBin), where=Total_PerEnergyBin!=0)
                            # TaggingProbPerEnergyBin[NOENTRIESMASK] = -1 # If there are no TPs, set tagging probability to -1

                            # compute statistical uncertainty
                            TaggingProbPerEnergyBin_STATUNC = [] 

                            for i, avg_val in enumerate(TaggingProbPerEnergyBin):
                                if(avg_val == -1):
                                    TaggingProbPerEnergyBin_STATUNC.append(0)
                                else:
                                    NTagged = Total_Tagged_PerEnergyBin[i]
                                    NTotal = Total_PerEnergyBin[i]

                                    if(NTagged == 0):
                                        NTagged_unc = 0
                                    else:
                                        NTagged_unc = 1/NTagged

                                    if(NTotal == 0):
                                        abs_error = 0 
                                        # NTotal_unc = 0
                                    else:
                                        # NTotal_unc = 1/NTotal     
                                        abs_error = np.sqrt(avg_val*(1-float(avg_val))/NTotal) # uncertainty for ratios which have common entries in numerator and denominator (not statistically independent)                            

                                    # rel_error = np.sqrt(NTagged_unc + NTotal_unc)
                                    # abs_error = avg_val * float(rel_error)

                                    # abs_error = np.sqrt(avg_val*(1-float(avg_val))/NTotal)

                                    TaggingProbPerEnergyBin_STATUNC.append(abs_error)
                            
                            TaggingProbPerEnergyBin_STATUNC = np.array(TaggingProbPerEnergyBin_STATUNC)
                            
                            # plot 
                            # Prepare figure and axes 
                            fig, ax = plt.subplots()
                            fig.set_dpi(100)
                            fig.set_size_inches(10, 7.5)

                            centerTheXbins = 0

                            if(centerTheXbins): centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
                            else: centered_energy_bins_ = energy_bins

                            centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
                            xerrors_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) for i in range(len(energy_bins) - 1) ]  
                            centered_energy_bins = np.array(centered_energy_bins_)
                            xerrors = np.array(xerrors_)

                            yVals_before_mask = np.copy(TaggingProbPerEnergyBin)
                            yUnc_before_mask = np.copy(TaggingProbPerEnergyBin_STATUNC)

                            centered_energy_bins = centered_energy_bins[NOENTRIESMASK]
                            yVals = TaggingProbPerEnergyBin[NOENTRIESMASK]
                            yUnc = TaggingProbPerEnergyBin_STATUNC[NOENTRIESMASK]  
                            xerrors = xerrors[NOENTRIESMASK]

                            zero_errors = [0. for i in range(0, len(yVals))]
                            error = 1
                            log = 0
                            xmin_, xmax_ = ybins[0], ybins[-1]

                            if(error):
                                plt.scatter(x = centered_energy_bins, y = yVals, label = "Severity = %s, %s"%(severity, time), s = 15)
                                plt.errorbar(x = centered_energy_bins, y = yVals, xerr = zero_errors, yerr = yUnc, fmt = " ") 
                            else:
                                plt.scatter(x = centered_energy_bins, y = yVals, label = "Severity = %s, %s"%(severity, time), s = 10)
                                plt.plot(x = centered_energy_bins, y = yVals, label = "__nolegend__", linestyle = '-')

                            LabelDict = {
                                "oneMinusEmuOverRealvstwrADCCourseBinning" : ["Average 1 - (Emulated / Real)", r"Real data TP $E_{T}$ (ADC)"],
                                "oneMinusEmuOverRealvstwrADCCourseBinningZoomed": ["Average 1 - (Emulated / Real)", r"Real data TP $E_{T}$ (ADC)"],
                                "EnergyVsTimeOccupancy" : ["Tagging probability", r"Real data TP $E_{T}$ (ADC)"],
                                "EmulEnergyVsTimeOccupancy" : ["Emulated tagging probability", r"Emulated TP $E_{T}$ (ADC)"]
                            }

                            varFileNameDict = {
                                "EnergyVsTimeOccupancy" : "TaggingProbability",
                                "EmulEnergyVsTimeOccupancy" : "EmulatedTaggingProbability"                                
                            }

                            outFileVarName = varFileNameDict[variable]

                            yLabel, xLabel = LabelDict[variable]

                            plt.xlabel(xLabel, fontsize=15)
                            plt.ylabel(yLabel, fontsize=15)    
                            plt.legend(loc = 'best', fontsize = 15)
                            plt.ylim(0, 1.01)

                            xmax_ = 17
                            print("Setting xmax to:",xmax_)

                            plt.xlim(xmin_, xmax_)
                            plt.grid()
                            if(log):
                                plt.ylim(0.0001, 1)
                                plt.yscale('log')  

                            outLocation = "%s/Sev_%s_%s_%s-times.png"%(direc_ol, severity, outFileVarName, time)

                            plt.savefig("%s/Sev_%s_%s_%s-times.png"%(direc_ol, severity, outFileVarName, time), dpi = 100)
                            plt.savefig("%s/Sev_%s_%s_%s-times.pdf"%(direc_ol, severity, outFileVarName, time), dpi = 300)    

                            plt.ylim(0.0001, 1)
                            plt.yscale('log')  

                            plt.savefig("%s/Sev_%s_%s_%s-times_log.png"%(direc_ol, severity, outFileVarName, time), dpi = 100)
                            plt.savefig("%s/Sev_%s_%s_%s-times_log.pdf"%(direc_ol, severity, outFileVarName, time), dpi = 300)                                                           

                            plt.close()                     

                            # save values as parquet files to eventually plot together with other datasets / working points   

                            parquetOutDir = "output/%s/%s/"%(sub_dataset, outFileVarName)
                            os.system("mkdir -p %s"%(parquetOutDir))
                            outName_yVals = "%s/%s_Sev_%s_%s-times_values.parquet"%(parquetOutDir, outFileVarName, severity, time)
                            outName_yUnc = "%s/%s_Sev_%s_%s-times_statuncer.parquet"%(parquetOutDir, outFileVarName, severity, time)

                            yVals_table = pa.table({"data": yVals_before_mask})
                            yUnc_table = pa.table({"data": yUnc_before_mask})
                            print("Saving parquet file:",outName_yVals)   
                            print("Saving parquet file:",outName_yUnc)                   
                            pa.parquet.write_table(yVals_table, outName_yVals)                
                            pa.parquet.write_table(yUnc_table, outName_yUnc)    

                            # exec("%s_%s_averages = np.copy(averages)"%(severity, time))
                            # exec("%s_%s_stdevs = np.copy(stdevs)"%(severity, time))                              
                                
        if(plotTogether):

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

            error = 1
            log = 1
            ADCToGeV = 1 # if true, divide x axis values by 2 to get GeV 

            # plot average lines on same plots
            avgPlotVars = ["oneMinusEmuOverRealvstwrADCCourseBinning", "oneMinusEmuOverRealvstwrADCCourseBinningZoomed", "EnergyVsTimeOccupancy", "EmulEnergyVsTimeOccupancy"]
            if(variable in avgPlotVars):
            
                for time in times:
                    print("On time:",time)
                    for sev_i, severity in enumerate(severities):
                        print("On severity:",severity)
                        if(time == "inTime" and severity == "three"): 
                            print("Skipping in time severity 3")
                            continue     
                        shape = shapeDict[severity]
                        color = colorDict[time]

                        # ratio plot 
                        fig, axarr = plt.subplots(2, 
                                                    sharex=True, 
                                                    gridspec_kw={
                                                        'hspace': 0.15,
                    #                                     'height_ratios': (0.8,0.2)
                                                        'height_ratios': (0.7,0.3)
                                                        }
                                                    )    
                        fig.set_size_inches(10, 7.5)
                        upper = axarr[0]                   
                        lower = axarr[1]

                        for direc_i, direc in enumerate(directories):
                            if(variable == "EnergyVsTimeOccupancy" or variable == "EmulEnergyVsTimeOccupancy"):
                                pass 
                            else: 
                                WP, PF, RECO = GetWorkingPointLabels(direc)                        

                            if("Run_352912" in direc):
                                sub_dataset = "Run352912"
                            elif("Runs_346446_346447_PilotBeam_2021" in direc):
                                sub_dataset = "PilotBeam2021"
                            else: sub_dataset = dataset 

                            # if(sub_dataset == "Run352912"): 
                                # print("Skipping ",sub_dataset)
                                # continue 

                            varFileNameDict = {
                                "EnergyVsTimeOccupancy" : "TaggingProbability",
                                "EmulEnergyVsTimeOccupancy" : "EmulatedTaggingProbability"                                
                            }

                            outFileVarName = varFileNameDict[variable]

                            # plot from previously saved values so that you don't need to reproduce plots by running over all files again 
                            if(fromParquet): 
                                # currently different treatment for EnergyVsTimeOccupancy w.r.t. oneMinus... vars.
                                if(variable == "EnergyVsTimeOccupancy" or variable == "EmulEnergyVsTimeOccupancy"):
                                    parquetOutDir = "output/%s/%s/"%(sub_dataset, outFileVarName)
                                    averages_path = "%s/%s_Sev_%s_%s-times_values.parquet"%(parquetOutDir, outFileVarName, severity, time)
                                    yUnc_path = "%s/%s_Sev_%s_%s-times_statuncer.parquet"%(parquetOutDir, outFileVarName, severity, time)

                                else:
                                    parquetOutDir = "output/%s/%s_%s_%s_%s/"%(dataset, variable, WP, PF, RECO)
                                    averages_path = "%s/%s_%s_%s_%s_%s_%s_averages.parquet"%(parquetOutDir, variable, WP, PF, RECO, severity, time)
                                    yUnc_path = "%s/%s_%s_%s_%s_%s_%s_stdevs.parquet"%(parquetOutDir, variable, WP, PF, RECO, severity, time)     
                       
                                averages_table = pq.read_table(averages_path)
                                yUnc_table = pq.read_table(yUnc_path)
                                data_averages = averages_table.to_pandas()
                                yUnc = yUnc_table.to_pandas()

                                averages_ = data_averages.to_numpy().ravel()
                                yUnc_ = yUnc.to_numpy().ravel()

                                # save values 
                                averages_path_txt = averages_path.replace(".parquet", ".txt")
                                yUnc_path_txt = yUnc_path.replace(".parquet", ".txt") 

                                np.savetxt(averages_path_txt, averages_, delimiter =', ')   
                                np.savetxt(yUnc_path_txt, yUnc_, delimiter =', ')     


                            else:
                                exec("averages_ = np.copy(%s_%s_averages)"%(severity, time))
                                exec("yUnc_ = np.copy(%s_%s_stdevs)"%(severity, time))


                            xbins, ybins = GetBins(variable, dataset)

                            if(variable == "EnergyVsTimeOccupancy" or variable == "EmulEnergyVsTimeOccupancy"): 
                                energy_bins = ybins
                                ol = "/eos/user/a/atishelm/www/EcalL1Optimization/{dataset}/".format(dataset=dataset)                    
                            else: 
                                energy_bins = xbins
                                ol = "/eos/user/a/atishelm/www/EcalL1Optimization/Dataset_Comparisons/"

                            rebin = 1

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
                            
                            averages = np.array(averages_) 
                            yUnc = np.array(yUnc_) 
                            centered_energy_bins = np.array(centered_energy_bins_)
                            xerrors = np.array(xerrors_)

                            MASK = tuple([averages != -1])

                            centered_energy_bins = centered_energy_bins[MASK]
                            averages = averages[MASK]
                            yUnc = yUnc[MASK]  
                            xerrors = xerrors[MASK]

                            zero_errors = [0. for i in range(0, len(averages))]            
                        
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

                            plotLabelDict = {
                                "Run352912" : "2.5",
                                "PilotBeam2021" : "0.5"
                            }

                            plotLabel = plotLabelDict[sub_dataset]

                            for key in labelReplaceDict:
                                val = labelReplaceDict[key]
                                plotLabel = plotLabel.replace(key, val)

                            color = next(upper._get_lines.prop_cycler)['color']

                            addXErr = 1 
                            addLine = 1 
                            if(addXErr):
                                xerr = xerrors 
                                addLine = 0 
                            else:
                                xerr = zero_errors

                            if(error):
                                upper.scatter(x = centered_energy_bins, y = averages, label = plotLabel, s = 15) ### [:-7] = remove the final 7 points (hack)
                                upper.errorbar(x = centered_energy_bins, y = averages, xerr = xerr, yerr = yUnc, fmt = " ", color = color)  
                                if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-', color = color)   
                            else:
                                upper.scatter(x = centered_energy_bins, y = averages, label = plotLabel, s = 15)
                                if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-')   
                            
                            plt.yticks(fontsize = 15)     

                            exec("y_vals_%s = np.copy(averages)"%(direc_i)) # save values for ratio later 
                            exec("y_Unc_%s = np.copy(yUnc)"%(direc_i)) # save values for ratio later 

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
                        upper.legend(loc = loc, title = r"$\delta_{min}$ (GeV)", fontsize = 20, prop={'size' : 14})
                        plt.rcParams['legend.title_fontsize'] = 20

                        upper.set_ylim(0, 1.01)

                        xmin_, xmax_ = energy_bins[0], energy_bins[-1]
                        xmax_ = 17
                        print("Setting xmax to:",xmax_)

                        yLabelDict = {
                            "oneMinusEmuOverRealvstwrADCCourseBinning" : "Average $E_{T}$ fraction subtracted",
                            "oneMinusEmuOverRealvstwrADCCourseBinningZoomed": "Average $E_{T}$ fraction subtracted",
                            "EnergyVsTimeOccupancy" : "Tagging probability",
                            "EmulEnergyVsTimeOccupancy" : "Emulated tagging probability"
                        }                        

                        yLabel = yLabelDict[variable]

                        plt.xlim(xmin_, xmax_)
                        plt.xlabel(xLabel, fontsize=15)
                        # plt.ylabel(yLabel, fontsize=15)
                        plt.grid()
                        plt.xticks(fontsize = 15)
                        plt.yticks(fontsize = 15)

                        # upper.set_xlim(xmin_, xmax_)
                        upper.set_ylabel(yLabel, fontsize=15)
                        upper.grid()
                        # upper.set_xticks_size(fontsize = 15)
                        # upper.set_yticks_size(fontsize = 15)       
                        upper.tick_params(axis='y', labelsize = 15)

                        lower.set_ylabel(r"  $\frac{\delta_{min} = 2.5 GeV}{\delta_{min} = 0.5 GeV}$", fontsize=15)           
                        lower.set_ylim(0, 2)      
                        lower.plot([xmin_, xmax_],[1,1],linestyle=':', color = 'black')

                        ratio_vals = np.divide(y_vals_0, y_vals_1, out=np.zeros_like(y_vals_1), where=y_vals_1!=0)

                        # compute statistical uncertainty in ratio 
                        ratio_unc = [] 

                        for i, ratio_val in enumerate(ratio_vals):
                            if(ratio_val == 0):
                                ratio_unc.append(0)
                            else:
                                N_0_unc = float(y_Unc_0[i])
                                N_1_unc = float(y_Unc_1[i])

                                N_0 = float(y_vals_0[i])
                                N_1 = float(y_vals_1[i])

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
                            lower.scatter(x = centered_energy_bins, y = ratio_vals, label = plotLabel, s = 15) ### [:-7] = remove the final 7 points (hack)
                            lower.errorbar(x = centered_energy_bins, y = ratio_vals, xerr = xerr, yerr = ratio_unc, fmt = " ", color = 'C0')  
                            #if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-', color = color)   
                        else:
                            lower.scatter(x = centered_energy_bins, y = ratio_vals, label = plotLabel, s = 15)
                            #if(addLine): upper.plot(centered_energy_bins, averages, label = "__nolegend__", linestyle = '-')                        

                        upperRightText = "upperRightText"
                        lumi = "lumi"
                        addLumi = 1
                        fontsize = 16
                        text_xmin = 0.1                          

                        if(dataset == "2021_2022_900GeVCollisions"):
                            unit = "{\mu}b"
                            lumi = "82-145"
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
                        
                        # fig.tight_layout()
                        plt.tight_layout()
                        plt.savefig("%s/Sev_%s_%s_Average_%s_linear.png"%(ol, severity, time, outFileVarName), dpi = 300)
                        plt.savefig("%s/Sev_%s_%s_Average_%s_linear.pdf"%(ol, severity, time, outFileVarName), dpi = 300)   

                        if(log):
                            upper.set_ylim(0.00001, 1)
                            upper.set_yscale('log')  
                            # fig.tight_layout()
                            plt.tight_layout()
                            plt.savefig("%s/Sev_%s_%s_Average_%s_log.png"%(ol, severity, time, outFileVarName), dpi = 300)
                            plt.savefig("%s/Sev_%s_%s_Average_%s_log.pdf"%(ol, severity, time, outFileVarName), dpi = 300)    
                        plt.close()

    print("DONE")

    # Population fraction plots:
    """

    # The purpose of this cell is to produce plots of population of sev X timing divided by total 
    # Previous plots showed us, e.g., that we are good at zeroing very late spikes.
    # However, it would help to know what percentage of the total spikes are very late in order to 
    # begin to get an idea of how we can lower the overall spike rate, as a function of TP energy. 

    ##-- stack plots are used when the total is of important interest in addition to the per piece 
    ##-- in this case we already know the total is 1. Line plots may be more useful to show. 

    ##-- Parameters
    log = 0
    xmin, xmax = 1, 256
    lumi = "XXX" 
    isWide = 0

    ##-- One plot per severity 
    for severity in severities:
        print("severity:",severity)
        colors = ["C%s"%(i) for i in range(0,4)]
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 7.5)
        for timing_i, timing in enumerate(timings):  
            if(timing == "inTime" and severity == "three"): continue
            print("timing:",timing)
            color = colors[timing_i]
            exec("N_all_entries_per_slice = np.copy(%s_N_entries_per_slice_all_timings)"%(severity))
            energy_bins = [1.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0, 80.0, 88.0, 96.0, 104.0, 112.0, 150.0, 256.0]
            centered_energy_bins_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) + energy_bins[i] for i in range(len(energy_bins) - 1) ]
            centered_energy_bins = np.array(centered_energy_bins_)
            xerrors_ = [ ((energy_bins[i+1] - energy_bins[i]) / 2.) for i in range(len(energy_bins) - 1) ]  
            xerrors = np.array(xerrors_)        
            exec("N_entries = np.copy(%s_%s_N_entries_per_slice)"%(severity, timing))

            MASK = tuple([N_entries != 0])
            
            zero_errors_ = [0. for i in range(0, len(zero_Early_N_entries_per_slice))] 
            zero_errors = np.array(zero_errors_)
            
            centered_energy_bins = centered_energy_bins[MASK]
            xerrors = xerrors[MASK]
            zero_errors = zero_errors[MASK]
            
            N_all_entries_per_slice = N_all_entries_per_slice[MASK]
            N_entries = N_entries[MASK]
            N_entries = np.divide(N_entries, N_all_entries_per_slice)
            
            plt.scatter(x = centered_energy_bins, y = N_entries, label = "%s"%(timing), s = 13, color = color)
            plt.plot(centered_energy_bins, N_entries, color = color )
            plt.errorbar(x = centered_energy_bins, y = N_entries, xerr = xerrors, yerr = zero_errors, fmt = " ", color = color )  
            
        plt.xlabel("Real data TP Et (ADC)", fontsize = 20)
        plt.ylabel("Population fraction", fontsize = 20)    
        legend = plt.legend(loc = 'best', fontsize = 15, title = 'Severity %s'%(severity))
        plt.setp(legend.get_title(),fontsize='15')
        plt.ylim(0, 1.1)
        plt.xlim(xmin, xmax)
        plt.grid()
        if(log):
            plt.ylim(0.0001, 1)
            plt.yscale('log')    
            
        Add_CMS_Header(plt, lumi, isWide, ax)
            
        plt.show()
        plt.close()
            
    print("DONE")    

    """

    # 2d plots 

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