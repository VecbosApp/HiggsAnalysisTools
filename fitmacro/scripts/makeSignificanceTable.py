#!/usr/bin/env python

import sys
import os
import ROOT
from ROOT import *
import HWWAnalysis.Misc.odict as odict
import hwwlimits
import optparse
import re


 
tagname = 'comb_shape'
basepath = os.getcwd()+'/significance/'

samplename = os.listdir(basepath)
samplename.sort()

points = ['mean', 'observed']

def openTFile(path, option=''):
    f =  ROOT.TFile.Open(path,option)
    if not f.__nonzero__() or not f.IsOpen():
        raise NameError('File '+path+' not open')
    return f

def getTree( file, tree ):
    t = file.Get(tree)
    if not t.__nonzero__():
        raise NameError('Tree '+str(tree)+' doesn\'t exist in '+file.GetName())
    return t


def getValue(file, q):

    cut = 'limit!=0 && limit<100'
    if q is 'mean':
        cut += ' && iToy>0'
    if q is 'observed':
        cut += ' && iToy==0'
    tree = getTree(file,'limit')
    command = 'limit>>h'
    tree.Draw(command,cut,'goff')
    hist = gDirectory.Get("h")
    value = hist.GetMean()
    hist.Delete()
    return value

def printTable(file, table):
    ## latex table
    print >> file, r'\documentclass[a4paper]{article}'
    print >> file, r'\begin{document}'
    print >> file, r'\pagestyle{empty}'
    print >> file, r'\begin{tabular}{|c|cc|}'
    print >> file, r'\hline'
    print >> file, r'Higgs mass [$\mathrm{GeV/c^2}$] & Observed significance & Expected Significance \\'
    print >> file, r'\hline'
    print >> file, r'\hline'
    for mass in table:
        print >> file, r' {mass} & {observed:.2f} & {mean:.2f} \\'.format(mass=mass,**(table[mass]))
    print >> file, r'\hline'
    print >> file, r'\end{tabular}'
    print >> file, r'\end{document}'  

def main():
    usage = 'usage: %prog [dir] [cmd]'
    parser = optparse.OptionParser(usage)
    (opt, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('One and only one datacard tag at the time')

    tag = args[0]

    if tag not in hwwlimits.dcnames['all']:
        parser.error('Wrong tag: '+', '.join(sorted(hwwlimits.dcnames['all'])))

    tagname = tag+'_shape'

    print tagname

    table = odict.OrderedDict()
    for sample in samplename:
        if not '.root' in sample:
            continue
        if not tagname in sample:
            continue
        if '123456' in sample:
            continue
        line = {}
        print sample

        for point in points:
            if 'observed' in point:
                path = basepath+'/'+sample
            if 'mean' in point:
                reBase = re.match('(.+)\.root',sample)
                if reBase:
                    path = basepath+'/'+reBase.group(1)+'.123456.root'
            f = openTFile(path)
            reMass = re.compile('.+\.mH(\d+)\.(.*)root')

            value = getValue(f,point)
            line[point] = value

            m = reMass.match(sample)
            if not m: 
                raise ValueError('Mass label not found in '+sample)
            
            mass = int(m.group(1))

            table[mass] = line

    latex = open(basepath+'/'+tagname+'.tex','w')
    printTable(latex,table)
    latex.close()


    ## summary file
    summary = open(basepath+'/'+tagname+'.summary', 'w')
    for mass in table:
        summary.write('{mass} {observed:.3f} 99 {mean:.3f} \n'.format( mass=mass, **(table[mass]) ) )

    summary.close()

if __name__ == '__main__':
    main()







## hww-1.55fb.mH120.comb.txt
## Observed Limit: r < 3.5430
## Expected  2.5%: r < 1.2811
## Expected 16.0%: r < 1.7041
## Expected 50.0%: r < 2.3613
## Expected 84.0%: r < 3.2799
## Expected 97.5%: r < 4.3577


## hww-1.55fb.mH130.comb.txt
## Observed Limit: r < 1.7859
## Expected  2.5%: r < 0.6184
## Expected 16.0%: r < 0.8226
## Expected 50.0%: r < 1.1398
## Expected 84.0%: r < 1.5832
## Expected 97.5%: r < 2.1035
