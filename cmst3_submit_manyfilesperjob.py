#! /usr/bin/env python
import os
import sys
import time
import commands
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 7:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName queue prefix"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
inputlist = "cmst3_35X/MC/"+process+"/"+dataset+".list"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
queue = sys.argv[5]
ijobmax = int(sys.argv[3])
application = sys.argv[4]
prefix = sys.argv[6]
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/m/meridian/VecBos21X/OutputSelectionInvertIsolation/"
diskoutputdir = "/cmsrm/pc18/crovelli/data/Higgs.3.5.X/"+prefix+"/"
outputmain = castordir+"/"+process+"/"+output
diskoutputmain = diskoutputdir+"/"+prefix+"/"+process+"/"+output
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir -p "+prefix+"/"+process+"/"+output)
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/log/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/input/")
os.system("mkdir -p "+prefix+"/"+process+"/"+output+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm18 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

#copy the configuration in the actual run directory
os.system("cp -r config "+prefix)

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+prefix+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    # prepare the script to run
    outputname = prefix+"/"+process+"/"+output+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp '+pwd+'/pdfs_MC.root $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+prefix+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+"\n")

#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm18:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+prefix+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+prefix+"/"+process+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    ijob = ijob+1
    time.sleep(5)
    continue
