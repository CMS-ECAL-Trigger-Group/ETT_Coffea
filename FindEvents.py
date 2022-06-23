"""
25 November 2021
Abraham Tishelman-Charny 

The purpose of this python module is to find ECAL TPs which have disagreeing real data and emulated TP energies.

values to save: runNb:lumiBlock:evtNb:twrADC:twrEmul3ADC:ieta:iphi:TCCid:TowerInTCC
"""

import argparse 
import pickle 
import uproot 
import numpy as np 

parser = argparse.ArgumentParser()
parser.add_argument("--inFile", type = str, help = "Input file")
parser.add_argument("--treeName", type = str, default = "tuplizer/ETTAnalyzerTree", help = "TDirectory / tree name of input root file")
parser.add_argument("--condor", action = "store_true", help = "Run over HTCondor")
args = parser.parse_args()

inFile = args.inFile 
treeName = args.treeName
condor = args.condor 

# # change import location based on running over condor or not 
# if(condor): 
#     from DummyPythonModule_Tools import MultiplyParameter
# else:
#     from python.DummyPythonModule_Tools import MultiplyParameter # make sure __init__.py file exists in python directory

print("Running FindEvents.py")
print("Input file: ",inFile)

f_u = uproot.open(inFile)
t = f_u[treeName]
vars_to_access = ["iphi"]

# create selections based on twrADC, twrEmul3ADC
twrADC_values = t["twrADC"].array().flatten()
twrEmul3ADC_values = t["twrEmul3ADC"].array().flatten()
MASK = np.array([twrADC_values != twrEmul3ADC_values]).flatten()

for v in vars_to_access:
    exec("{v}_values = t[v].array().flatten()".format(v=v))
    exec("{v}_values = {v}_values[MASK]".format(v=v))

print("")

#newValue = MultiplyParameter(InputParameter)

# just to test transfer of output files, save newValue as a pickle file 
#pickle.dump( newValue, open( 'DummyValue.p', "wb" ))  # on condor, this then gets transferred from the condor area to your output location     

print("Done running FindEvents.py")
