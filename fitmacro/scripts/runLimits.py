#!/usr/bin/env python

import sys
import os.path
import optparse
import fnmatch
import subprocess
import hwwlimits

indir  = os.getcwd()+'/datacards/'
outdir = os.getcwd()+'/limits/'
def main():

    usage = '''usage: %prog [opts] dctag'''
    parser = optparse.OptionParser(usage)
    parser.add_option('-s','--stepping',dest='stepping',help='Switch stepping on ', action='store_true', default=False)
    parser.add_option('-1',dest='minuit1',help='Minuit ', action='store_true', default=False)
    parser.add_option('-n',dest='dryrun',help='Dry run ', action='store_true', default=False)
    parser.add_option('-o',dest='observed',help='Observed only', action='store_true', default=False)
    parser.add_option('-S','--significance',dest='significance',help='Compute the significance instead of the limit ', action='store_true', default=False)
    parser.add_option('--prefix','-p',dest='prefix',help='prefix',default=None)
    parser.add_option('-l', '--lumi'     , dest='lumi'        , help='Luminosity'                            , default=None   , type='float'   )

    (opt, args) = parser.parse_args()
    print 'Running with lumi',opt.lumi
    
    constraints = {
        '*':'--rMin=0.01'
    }

    if len(args) != 1:
        parser.error('One and only one datacard tag at the time')

    tag = args[0]

    if tag not in hwwlimits.dcnames['all']:
        parser.error('Wrong tag: '+', '.join(sorted(hwwlimits.dcnames['all'])))

    tmpl = 'hww-{lumi:.1f}fb.mH{mass}.{tag}_shape.txt'
    masses = [115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 180]

    if opt.prefix:
        os.chdir(opt.prefix)
        print 'Running in ',os.getcwd()

    if not os.path.exists(indir):
        print 'Error: datacard directory',indir,'not found'
        sys.exit(-1)

    allcards = [(mass,os.path.join(indir,tmpl.format(mass=mass, tag=tag, lumi=opt.lumi))) for mass in masses] 
    for mass,card in allcards:
        if not os.path.exists(card):
            print 'Error: missing datacard: '+card
            sys.exit(-1)

    os.system('mkdir -p '+outdir)

    tagname = 'HWW_'+tag+'_shape'
    for mass,card in allcards:
        exe  = 'combine '
        flags = ' -n %s -m %s %s'%(tagname,mass,card)
        if opt.stepping:
            flags = ' --minosAlgo stepping'+flags
        if opt.minuit1:
            flags = ' --minimizerAlgo Minuit'+flags
        if opt.observed:
            flags = ' --run observed'+flags
        if opt.significance:
            flags = ' -M ProfileLikelihood --significance --expectSignal=1 -t 100 '+flags
        else:
            flags = ' -M Asymptotic '+flags
        if not opt.significance:
            for c,flag in constraints.iteritems():
                if fnmatch.fnmatch(str(mass),c):
                    flags = ' '+flag+flags

        command = exe+flags

        print '-'*50
        print command
        if opt.dryrun: continue
        code = os.system(command)
        if opt.significance:
            move = 'mv higgsCombine%s.ProfileLikelihood.mH%d.*.root %s' % (tagname,mass,outdir)
        else:
            move = 'mv higgsCombine%s.Asymptotic.mH%d.root %s' % (tagname,mass,outdir)
        print move
        os.system(move)
        
        if code: sys.exit(code)

if __name__ == '__main__':
    main()

