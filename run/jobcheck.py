#!/usr/bin/env python3

################################################
# Check the error and output log files of jobs #
################################################
# Use case:
#   Use this tool for automated checking of job failures.
#   The check happens on the basis of the _err_ files
#   that are created by the condor job scheduler.
#   Two separate checks are being performed:
#     - Scan for known error strings, see list below.
#     - Scan for 'starting' and 'done' tags.
#   For the latter check to work, your executable must write
#   the correct tags to stderr at the beginning and at the end of the job.
#   See existing examples for correct usage.
# Usage:
#   Run 'python jobCheck.py -h' for a list of options.
#   You can run the script from this directory and specify the job directory in the args,
#   or alternatively you can run the script from the job directory
#   (using 'python [path]/jobCheck.py) and leave the directory arg at its default.
#   If you have sourced the 'source.sh' script in the project's main directory,
#   you can simply run 'jobCheck [+args]' from anywhere, 
#   without specifying "python" or the path to this script.
# Note:
#   Should work for both condor and qsub log files.
#   The latter has not been used in a long time however, so not sure.


import sys
import os
import argparse
import glob


def check_start_done( filename, 
                      starting_tag='###starting###',
                      done_tag='###done###',
                      ntarget=None,
                      verbose=True ):
    ### check starting and done tags in a file.
    # returns 0 in case of no errors, which is defined as:
    #   - at least one starting tag is present in the file.
    #   - the number of done tags equals the number of starting tags.
    #   - the number of done tags equals the target number (if provided).
    # returns 1 in all other cases.

    # read the file content
    f = open(filename)
    filetext = f.read()
    f.close()
    
    # count number of starting tags
    nstarted = filetext.count(starting_tag)
    if(nstarted==0):
        if verbose:
            msg = 'WARNING in jobCheck.py: file {}'.format(filename)
            msg += ' does not contain a valid starting tag "{}"'.format(starting_tag)
            print(msg)
        return 1

    # count number of done tags
    ndone = filetext.count(done_tag)
    if ntarget is None: ntarget = ndone

    # return 0 if all is ok
    if(nstarted==ndone and ndone==ntarget): return 0

    # return 1 otherwise
    if verbose:
        msg = 'WARNING in jobCheck.py: found issue in file {}:\n'.format(filename)
        msg += '   {} commands were initiated.\n'.format(nstarted)
        msg += '   {} seem to have finished normally.\n'.format(ndone)
        msg += '   {} were expected.\n'.format(ntarget)
        print(msg)
    return 1


def check_error_content(filename, contentlist='default', verbose=True):
    ### check for known error messages in a file.
    # returns 0 if none of the elements of contentlist is present in the file;
    # returns 1 otherwise.

    # read the file content
    f = open(filename)
    filetext = f.read()
    f.close()

    # hard-coded default error content
    if( isinstance(contentlist,str) and contentlist=='default' ):
        contentlist = ([    'SysError',
                           '/var/torque/mom_priv/jobs',
                           'R__unzip: error',
                           'hadd exiting due to error in',
                           'Bus error',
                           'Exception:',
                           'Traceback (most recent call last):' ])
        contentlist.append('###error###') # custom error tag for targeted flagging

    # check if the file content contains provided error tags
    contains = []
    for idx,content in enumerate(contentlist):
        if filetext.count(content)>0:
            contains.append(idx)
    if len(contains)==0: return 0
    if verbose:
        msg = 'WARNING in jobCheck.py: found issue in file {}:\n'.format(filename)
        for idx in contains:
            msg += '   found sequence {}\n'.format(contentlist[idx])
        print(msg)
    return 1


if __name__=='__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Check job logs for errors.')
    parser.add_argument('-d', '--dir', default=[os.getcwd()], nargs='+',
                        help='Directory to scan for files (default: cwd)')
    args = parser.parse_args()

    # print arguments
    print('running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # find files
    print('finding files...')
    files = []
    for checkdir in args.dir:
        condorpattern = os.path.join(checkdir, '*.err')
        files += glob.glob(condorpattern)
    nfiles = len(files)
    print('found {} error log files.'.format(nfiles))
    print('start scanning...')

    # loop over files
    nerror = 0
    for fname in files:
        # initialize
        error_content = 0
        # error checking
        error_content = check_error_content(fname)
        if error_content > 0: nerror += 1

    # print results
    print('number of files scanned: {}'.format(nfiles))
    print('number of files with error: {}'.format(nerror))
    print('number of files without apparent error: {}'.format(nfiles-nerror))
