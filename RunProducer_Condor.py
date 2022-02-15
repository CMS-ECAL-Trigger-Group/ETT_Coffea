"""
18 November 2021 

The purpose of this module is to submit condor jobs to run the ETT coffea producer. 

Example commands:

# Fullreadout data analysis:

# To submit:

# Submitted:
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/" --tag=220214_122937 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_324725_306425_FullReadoutData/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/" --tag=220214_092624 -s

# 2021 pilot beam analysis:
# Submitted:
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithoutOddPeakFinder_2p5PrimeODDweights/" --tag=220210_104645 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithoutOddPeakFinder_0p5PrimeODDweights/" --tag=220210_104904 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/" --tag=220210_104547 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_0p5PrimeODDweights/" --tag=220210_104023 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithoutOddPeakFinder_2p5PrimeODDweights/" --tag=220210_103954 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_weightsRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/" --tag=220210_104615 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_0p5PrimeODDweights/" --tag=220210_094402 -s
python RunProducer_Condor.py --direc="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/Trigger/DoubleWeights/Runs_346446_346447_PilotBeam_2021/ETTAnalyzer_CMSSW_12_1_0_pre3_DoubleWeights_MultifitRecoMethod_StripZeroingMode_WithOddPeakFinder_2p5PrimeODDweights/" --tag=220209_125921 -s

Misc:

Thank you: https://research.cs.wisc.edu/htcondor/manual/v8.5/condor_submit.html
https://github.com/htcondor/htcondor/blob/abbf76f596e935d5f2c2645e439cb3bee2eef9a7/src/condor_starter.V6.1/docker_proc.cpp ##-- Docker/HTCondor under the hood 

"""

import os, sys
import argparse
import logging
import pwd
import subprocess
import shutil
import time
import glob

#!/usr/bin/python

logging.basicConfig(level=logging.DEBUG)

script_TEMPLATE = """#!/bin/bash
##source /cvmfs/cms.cern.ch/cmsset_default.sh

##source /cvmfs/cms.cern.ch/cmsset_default.sh
#source /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask:latest
export SCRAM_ARCH=slc7_amd64_gcc820

echo
echo $_CONDOR_SCRATCH_DIR
cd   $_CONDOR_SCRATCH_DIR
echo
echo "... start job at" `date "+%Y-%m-%d %H:%M:%S"`
echo "----- directory before running:"
ls -lR .
#echo "----- CMSSW BASE, python path, pwd:"
#echo "+ CMSSW_BASE  = $CMSSW_BASE"
echo "+ PYTHON_PATH = $PYTHON_PATH"
echo "+ PWD         = $PWD"

python RunProducer.py --jobNum=$1 --infile=$2 --treename="tuplizer/ETTAnalyzerTree" # command 

echo "----- directory after running :"
ls -lR .
echo " ------ THE END (everyone dies !) ----- "
"""

condor_TEMPLATE = """
#request_disk          = 1024
request_disk          = 2048
request_memory = 8000
executable            = {jobdir}/script.sh
arguments             = $(ProcId) $(jobid) 
transfer_input_files = {transfer_files}, $(jobid)

# environment = "LD_LIBRARY_PATH_STORED=/eos" 

output                = $(ClusterId).$(ProcId).out
error                 = $(ClusterId).$(ProcId).err
log                   = $(ClusterId).$(ProcId).log
initialdir            = {jobdir}

transfer_output_remaps = "{transfer_output_remaps}"

#Requirements = HasSingularity
+JobFlavour           = "{queue}"
#+SingularityImage = "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask:latest"

HasDocker = true
universe = docker
docker_image = coffeateam/coffea-dask:latest

queue jobid from {jobdir}/inputfiles.dat
"""

def main():
    parser = argparse.ArgumentParser(description='Famous Submitter')
    parser.add_argument("-t"   , "--tag"   , type=str, default="Exorcism"  , help="production tag", required=True)
    parser.add_argument("-q"   , "--queue" , type=str, default="espresso", help="")
    parser.add_argument("-f"   , "--force" , action="store_true"          , help="recreate files and jobs")
    parser.add_argument("-s"   , "--submit", action="store_true"          , help="submit only")
    parser.add_argument("-dry" , "--dryrun", action="store_true"          , help="running without submission")
    parser.add_argument("--direc", type = str , help="Directory with input files")

    options = parser.parse_args()

    indir = "{}/{}/".format(options.direc, options.tag)

    for sample in os.listdir(indir):
        if "merged" in sample:
            continue
        if "WS" in sample:
            continue
        if "output" in sample:
            continue 

        jobs_dir = '_'.join(['jobs', options.tag, sample])
        logging.info("-- sample_name : " + sample)

        if os.path.isdir(jobs_dir):
            if not options.force:
                logging.error(" " + jobs_dir + " already exist !")
                continue
            else:
                logging.warning(" " + jobs_dir + " already exists, forcing its deletion!")
                shutil.rmtree(jobs_dir)
                os.mkdir(jobs_dir)
        else:
            os.mkdir(jobs_dir)
        with open(os.path.join(jobs_dir, "inputfiles.dat"), 'w') as infiles:
            in_files = glob.glob("{indir}/{sample}/*.root".format(sample=sample, indir=indir))
            for name in in_files:
                infiles.write(name+"\n")
            infiles.close()

        vars = ["oneMinusEmuOverRealvstwrADCCourseBinning", "EnergyVsTimeOccupancy"]
        # vars = ["oneMinusEmuOverRealvstwrADCCourseBinning"]
        times = ["all", "inTime", "Early", "Late", "VeryLate"]
        severities = ["zero", "three", "four"]

        outdir = indir + sample + "_output/"
        os.system("mkdir -p {}".format(outdir))
        print("Making directories if they don't exist...")
        for var in vars:
            for sev in severities:
                for time in times:
                    if(time == "inTime" and sev == "three"): continue ##-- skip in time severity 3 
                    outdir_perVarSevTime = indir + sample + "_output/{var}/{sev}/{time}/".format(var=var, sev=sev, time=time)
                    print(outdir_perVarSevTime)
                    os.system("mkdir -p {}".format(outdir_perVarSevTime))

        with open(os.path.join(jobs_dir, "script.sh"), "w") as scriptfile:
            script = script_TEMPLATE.format(
                outputdir=outdir
            )
            scriptfile.write(script)
            scriptfile.close()

        with open(os.path.join(jobs_dir, "condor.sub"), "w") as condorfile:
            allFiles = [
                "../RunProducer.py",
                "../python/Producer.py",
                "../python/SumWeights.py"
            ]

            # output files 
            # depends on variables being output
            vars = ["oneMinusEmuOverRealvstwrADCCourseBinning", "EnergyVsTimeOccupancy"]
            # vars = ["oneMinusEmuOverRealvstwrADCCourseBinning"]
            times = ["all", "inTime", "Early", "Late", "VeryLate"]
            severities = ["zero", "three", "four"]

            transfer_files = []

            for var in vars:
                for sev in severities:
                    for time in times:
                        transfer_file = "{var}_sev{sev}_{time}_values.p={output_dir}/{var}/{sev}/{time}/{var}_sev{sev}_{time}_values_$(ProcId).p".format(var=var, sev=sev, time=time, output_dir = outdir)
                        transfer_files.append(transfer_file)

            transfer_output_remaps = str(";".join(transfer_files))

            condor = condor_TEMPLATE.format(
                transfer_files = ",".join(allFiles),
                transfer_output_remaps=transfer_output_remaps,
                jobdir=jobs_dir,
                queue=options.queue
            )
            condorfile.write(condor)
            condorfile.close()
        if options.dryrun:
            continue

        htc = subprocess.Popen(
            "condor_submit " + os.path.join(jobs_dir, "condor.sub"),
            shell  = True,
            stdin  = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            close_fds=True
        )
        out, err = htc.communicate()
        exit_status = htc.returncode
        logging.info("condor submission status : {}".format(exit_status))

if __name__ == "__main__":
    main()
    print("DONE")
