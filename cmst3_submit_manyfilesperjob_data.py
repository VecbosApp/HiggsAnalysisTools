#! /usr/bin/env python
import os, sys, time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py process dataset nfileperjob applicationName queue dirname
#######################################
if len(sys.argv) != 7:
    print "usage cmst3_submit_manyfilesperjob.py process dataset nfileperjob applicationName queue dirname"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
inputlist = "cmst3_38X/"+process+"/"+dataset+".list"
#settingfile = "config/RSZZsettings.txt"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[5]
nfileperjob = int(sys.argv[3])
#application = "VecbosApp"
application = sys.argv[4]
dirname = sys.argv[6]
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/e/emanuele/Higgs2010/"+dirname+"/"
diskoutputdir = "/cmsrm/pc21_2/emanuele/data/Higgs3.9.X/"+dirname+"/"
outputmain = castordir+"/"+process+"/"+output
diskoutputmain = diskoutputdir+"/"+process+"/"+output
# prepare job to write on the cmst3 cluster disks
################################################
os.system("rm -rf "+process+"/"+output+"/"+dirname)
os.system("mkdir -p "+process+"/"+output+"/"+dirname)
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/log/")
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/input/")
os.system("mkdir -p "+process+"/"+output+"/"+dirname+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
#    os.system("rfrm -r "+outputmain)
    os.system("rfmkdir -p "+castordir)
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+castordir)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 rm -rf "+diskoutputmain)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputdir)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
#print inputlist
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0
while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+process+"/"+output+"/"+dirname+"/"+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(nfileperjob,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = process+"/"+output+"/"+dirname+"/"+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp '+pwd+'/pdfs_MC.root $WORKDIR\n')
    outputfile.write('cp -r '+pwd+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+"\n")
#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
    outputfile.write('ls *.root | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/"+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+process+"/"+output+"/"+dirname+"/"+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    time.sleep(5)
    ijob = ijob+1
    continue
