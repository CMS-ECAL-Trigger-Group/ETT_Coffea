# ETT Coffea

This is a repository originating from [SUEPPhysics](https://github.com/SUEPPhysics/SUEPCoffea), whose purpose is to process CMS ECAL trigger team analyzer output files. 

## Example run steps  

```
  git clone git@github.com:CMS-ECAL-Trigger-Group/ETT_Coffea.git
  cd ETT_Coffea
  singularity shell -B ${PWD} -B /afs -B /eos /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask:latest ##-- Mount /afs and /eos space to run on files in those locations. https://hsf-training.github.io/hsf-training-docker/10-singularity/index.html
  python3 RunProducer.py --inDir="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Run_346446_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeightsTaggingMode/211115_170649/oneFile/" --treename="tuplizer/ETTAnalyzerTree" --outDir="/eos/user/a/atishelm/www/EcalL1Optimization/ETT_Coffea_singleFile/" --condor="0"
```
