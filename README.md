# Version of NanoHRT-tools for JMS/JMR scale factor measurement in 2024 scouting data

This is a fork from the very similar [earlier version for 2022-2023](https://github.com/beatrizlopes/NanoHRT-tools/tree/dev-UL-0201), which runs on 2024 scouting NanoAOD.
It produces ntuples in a ttbar single-muon enriched region for a tag-and-probe like JMS/JMR scale factor measurement.


### How to install

Set up CMSSW:
```bash
cmsrel CMSSW_15_0_6
cd CMSSW_15_0_6/src
cmsenv
```

Set up NanoAOD-tools
```bash
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
```

Get this repository
```bash
git clone https://github.com/Floria-ma/NanoHRT-tools.git PhysicsTools/NanoHRTTools
```

Note: depending on your system and environment, you might need to make the following changes:
- In `NanoAODTools/scripts/haddnano.py`, change the first line from `#!/bin/env python` to `#!/bin/env python3`.

Compile
```bash
scram b -j8
```

### How to run:

Prepare all jobs (for simulation):
```bash
python3 runHeavyFlavTrees.py --scouting -i /eos/cms/store/cmst3/group/vhcc/ScoutingNanoAOD/2024 -o /eos/user/z/zima/HadronicVH/20251204_ULNanoV9_MassRegression --jet-type ak15 --channel muon --year 2024 -n 20 --condor-extras '+AccountingGroup = "group_u_CMST3.all" ' --po jec False --po jes False --po jer False --sfbdt 0 --run-syst
```

Test locally (optional but recommended after every change) by going into `jobs_<name of output>/mc`
and run `python3 postprocessor 0`
(where the last argument is the number of the job to run, can be anything between 0 and the number of prepared jobs).

You can check the total number of jobs in the `jobs_<name of output>/mc/submit.txt`.
If it is not unreasonably large, you can submit them with `condor_submit jobs_<name of output>/mc/submit.cmd`.

Some command line options:
- `-n` is the number of input files per job.
- you can do faster testing with fewer jobs by omitting the `--run-syst` (in order to skip systematics).
- For more command line options, see the documentation of the [upstream repository](https://github.com/beatrizlopes/NanoHRT-tools/tree/dev-UL-0201).
