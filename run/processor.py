#!/usr/bin/env python3
from __future__ import print_function

import os
import sys
import time
import json
import argparse
import subprocess
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def xrd_prefix(filepaths):
    prefix = ''
    allow_prefetch = False
    if not isinstance(filepaths, (list, tuple)):
        filepaths = [filepaths]
    filepath = filepaths[0]
    if filepath.startswith('/eos/cms'):
        prefix = 'root://eoscms.cern.ch/'
    elif filepath.startswith('/eos/user'):
        prefix = 'root://eosuser.cern.ch/'
    elif filepath.startswith('/eos/uscms'):
        prefix = 'root://cmseos.fnal.gov/'
    elif filepath.startswith('/store/'):
        # remote file
        import socket
        host = socket.getfqdn()
        if 'cern.ch' in host:
            prefix = 'root://xrootd-cms.infn.it//'
        else:
            prefix = 'root://cmsxrootd.fnal.gov//'
        allow_prefetch = True
    expanded_paths = [(prefix + '/' + f if prefix else f) for f in filepaths]
    return expanded_paths, allow_prefetch


def outputName(md, jobid):
    info = md['jobs'][jobid]
    return '{samp}_{idx}_tree.root'.format(samp=info['samp'], idx=info['idx'])


def main(args):

    # load job metadata
    print('Loading job metadata...')
    with open(args.metadata) as fp:
        md = json.load(fp)

    # load modules
    modules = []
    print('Loading modules...')
    for mod, names in md['imports']:
        import_module(mod)
        obj = sys.modules[mod]
        selnames = names.split(",")
        for name in dir(obj):
            if name[0] == "_":
                continue
            if name in selnames:
                print("Loading %s from %s " % (name, mod))
                modules.append(getattr(obj, name)())

    # remove any existing root files
    for f in os.listdir('.'):
        if f.endswith('.root'):
            print('Removing file %s' % f)
            os.remove(f)

    # get input and output settings
    inputfiles = args.files if len(args.files) else md['jobs'][args.jobid]['inputfiles']
    filepaths, allow_prefetch = xrd_prefix(inputfiles)
    outputname = outputName(md, args.jobid)

    # do printouts for easier checking what's going on
    print('Found following input files:')
    for filepath in filepaths: print(f'  - {filepath}')
    print('Found following settings for the PostProcessor:')
    print(f'  - Preselection cut: {md.get("cut")}')
    print(f'  - Preselection json: {md.get("json")}')
    print(f'  - Modules: {modules}')
    # add more info as needed in debugging...

    # run the postprocessor
    print('Running PostProcessor...')
    sys.stdout.flush()
    sys.stderr.flush()
    p = PostProcessor(outputDir='.',
                      inputFiles=filepaths,
                      cut=md.get('cut'),
                      branchsel=os.path.basename(md['branchsel_in']),
                      modules=modules,
                      compression=md.get('compression', 'LZMA:9'),
                      friend=md.get('friend', False),
                      postfix=md.get('postfix'),
                      jsonInput=md.get('json'),
                      provenance=md.get('provenance', False),
                      haddFileName=None,
                      prefetch=(allow_prefetch and md.get('prefetch', False)),
                      longTermCache=md.get('longTermCache', False),
                      maxEntries=md.get('maxEntries', None),
                      firstEntry=md.get('firstEntry', 0),
                      outputbranchsel=os.path.basename(md['branchsel_out'])
                      )
    p.run()
    sys.stdout.flush()
    sys.stderr.flush()

    # sort files in descending order by number of entries
    # (hadd crashes if the first file to merge has 0 events)
    output_files = [f for f in os.listdir() if f.endswith('.root')]
    nentries = []
    for output_file in output_files:
        f = ROOT.TFile.Open(output_file)
        tree = f.Get('Events')
        nentries.append(tree.GetEntries())
        f.Close()
    sorted_ids = sorted(range(len(nentries)), key=lambda idx: nentries[idx], reverse=True)
    output_files = [output_files[idx] for idx in sorted_ids]
    nentries = [nentries[idx] for idx in sorted_ids]
    print('Found following number of entries in output files:')
    for f, n in zip(output_files, nentries): print(f'  - {f}: {n}')

    # hadd files
    cmd = f'haddnano.py {outputname} ' + ' '.join(output_files)
    print('Running haddnano as follows:')
    print(cmd)
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    if p.returncode != 0:
        raise RuntimeError('Hadd failed!')

    # keep only the hadd file
    for f in os.listdir('.'):
        if f.endswith('.root') and f != outputname:
            os.remove(f)

    # stage out
    if md['outputdir'].startswith('/eos'):
        cmd = 'xrdcp --silent -p -f {outputname} {outputdir}/{outputname}'.format(
            outputname=outputname, outputdir=xrd_prefix(md['joboutputdir'])[0][0])
        print(cmd)
        success = False
        for count in range(args.max_retry):
            p = subprocess.Popen(cmd, shell=True)
            p.communicate()
            if p.returncode == 0:
                success = True
                break
            else:
                time.sleep(args.sleep)
        if not success:
            raise RuntimeError("Stage out FAILED!")

        # clean up
        os.remove(outputname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-m', '--metadata',
                        default='metadata.json',
                        help='Path to the metadata file. Default:%(default)s')
    parser.add_argument('--max-retry',
                        type=int, default=3,
                        help='Max number of retry for stageout. Default: %(default)s'
                        )
    parser.add_argument('--sleep',
                        type=int, default=120,
                        help='Seconds to wait before retry stageout. Default: %(default)s'
                        )
    parser.add_argument('--files',
                        nargs='*', default=[],
                        help='Run over the specified input file. Default:%(default)s')
    parser.add_argument('jobid', type=int, help='Index of the output job.')

    args = parser.parse_args()
    main(args)
