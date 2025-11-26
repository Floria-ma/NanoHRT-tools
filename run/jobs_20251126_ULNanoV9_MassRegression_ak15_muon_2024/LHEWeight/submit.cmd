universe              = vanilla
requirements          = (Arch == "X86_64") && (OpSys == "LINUX")
MY.WantOS = "el7"
request_memory        = 2000
request_disk          = 10000000
executable            = /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/run_postproc_condor.sh
arguments             = $(jobid)
transfer_input_files  = /afs/cern.ch/user/z/zima/CMSSW.tar.gz,/afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/processor.py,/afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/metadata.json,/afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/heavyFlavSFTree_cfg.json,/afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/keep_and_drop_input.txt,/afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/keep_and_drop_output_LHEweights.txt
output                = /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/$(jobid).out
error                 = /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/$(jobid).err
log                   = /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/$(jobid).log
use_x509userproxy     = true
Should_Transfer_Files = YES
initialdir            = /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight
WhenToTransferOutput  = ON_EXIT
want_graceful_removal = true
on_exit_remove        = (ExitBySignal == False) && (ExitCode == 0)
on_exit_hold          = ( (ExitBySignal == True) || (ExitCode != 0) )
on_exit_hold_reason   = strcat("Job held by ON_EXIT_HOLD due to ", ifThenElse((ExitBySignal == True), "exit by signal", strcat("exit code ",ExitCode)), ".")
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 10*60)
transfer_output_files = ""

+MaxRuntime = 24*60*60
+AccountingGroup = "group_u_CMST3.all" 

queue 5000 jobid from /afs/cern.ch/user/z/zima/CMSSW_15_0_6/src/PhysicsTools/NanoHRTTools/run/jobs_20251126_ULNanoV9_MassRegression_ak15_muon_2024/LHEWeight/submit.txt
