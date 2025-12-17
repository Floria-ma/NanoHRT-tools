# Simple utility script to print the available branches in a root file


import os
import sys
import ROOT
from fnmatch import fnmatch


if __name__=='__main__':

    # read command line args
    inputfile = sys.argv[1]
    bfilter = None
    if len(sys.argv)>=3: bfilter = sys.argv[2]

    # read input file
    print(f'Reading input file {inputfile}...')
    f = ROOT.TFile.Open(inputfile)
    tree = f.Get("Events")

    # print branches
    for branch in tree.GetListOfBranches():
        name = branch.GetName()
        if bfilter is not None and not fnmatch(name, bfilter): continue
        print(f' - {name}: {branch.GetClassName()}')
