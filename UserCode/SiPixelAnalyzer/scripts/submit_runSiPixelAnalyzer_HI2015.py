#!/usr/bin/env python
import os,sys
import optparse
import commands
import time
import re

#command line configuration
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,dest='queue'  ,help='batch queue'          ,default='cmscaf1nd')
parser.add_option('-j', '--jobs'       ,dest='jobs'   ,help='number of jobs'       ,default=1,    type=int)
parser.add_option('-i', '--inputF'     ,dest='inputF',help='input file list'       ,default='run262699.list', type='string')
parser.add_option('-n', '--nevts'      ,dest='nevts'  ,help='number of evetns/job' ,default=2003,  type=int)
parser.add_option(      '--proxy'      ,dest='proxy'  ,help='proxy to be used'     ,default=None, type='string')
parser.add_option('-o', '--output'     ,dest='output' ,help='output directory'     ,default='/store/group/phys_heavyions/mverweij/pixel/tree/v1')

(opt, args) = parser.parse_args()

#prepare working directory
cmsswBase=os.environ['CMSSW_BASE']
workBase=os.getcwd()
jobsBase='%s/FARM%s'%(workBase,time.time())
os.system('mkdir -p %s'%jobsBase)

#init a new proxy if none has been passed
if opt.proxy is None:
    print 'Initiating a new proxy'
    os.system('voms-proxy-init --voms cms --valid 72:00')
    os.system('cp /tmp/x509up_u`id -u \`whoami\`` %s/proxyforprod' % workBase)
    print 'Production proxy is now available in %s/proxyforprod (valid for 72h)' % workBase

#loop over the required number of jobs
inputfile='%s'% (opt.inputF)
newfile='new_%s' % (opt.inputF)
f = open(inputfile,'r')
newf = open(newfile,'w')
jobCounter=0
m = re.compile(r"000/\d\d\d/\d\d\d") #used to find run number
for line in f:
    if not line.find('root://eoscms//eos/cms') :
        if not line.startswith('#') :
            print line
            newf.write(line.replace('root://eoscms//eos/cms', '# root://eoscms//eos/cms'))

            found = m.findall(line)
            outdir = '%s/%s' % (opt.output,found[0])
                       
            #create bash script to be submitted
            scriptFile = open('%s/runJob_%d.sh'%(jobsBase,jobCounter), 'w')
            scriptFile.write('#!/bin/bash\n')
            scriptFile.write('export X509_USER_PROXY=%s/proxyforprod\n' % workBase)
            scriptFile.write('cd %s/src\n'%cmsswBase)
            scriptFile.write('eval `scram r -sh`\n')
            scriptFile.write('cd -\n')
            scriptFile.write('cp %s/runSiPixelAnalyzerAndHLT.py %s \n' % (workBase,jobsBase))
            scriptFile.write('cp %s/runSiPixelAnalyzerAndHLT.py .\n' % workBase)
            scriptFile.write('cmsRun runSiPixelAnalyzerAndHLT.py outputFile=histos_%d.root maxEvents=-1 inputFiles=%s\n' % (jobCounter,line) )
            scriptFile.write('cmsMkdir %s\n' % outdir)
            scriptFile.write('ls\n')
            scriptFile.write('cmsStage histos_%dHLT.root %s/histos_%dHLT.root\n' % (jobCounter,outdir,jobCounter) )
            scriptFile.write('cmsStage histos_%dSiPixelAnalyzer.root %s/histos_%dSiPixelAnalyzer.root\n' % (jobCounter,outdir,jobCounter) )
            scriptFile.write('cmsStage tree_%dSiPixelAnalyzer.root %s/tree_%dSiPixelAnalyzer.root\n' % (jobCounter,outdir,jobCounter) )
            scriptFile.write('rm *root\n')
            #scriptFile.write('rm histos_%d.root\n' % (jobCounter))
            scriptFile.close()

            #preare to run it
            os.system('chmod u+rwx %s/runJob_%d.sh' % (jobsBase,jobCounter))

            #submit it to the batch or run it locally
            if opt.queue=='':
                print 'Job #%d will run locally' % jobCounter
                os.system('%s/runJob_%d.sh' % (jobsBase,jobCounter) )
            else:
                print 'Job #%d will run remotely' % jobCounter
                #os.system("echo bsub -q %s -R \"swp>1000 && pool>30000\" -J FOREST_%d \'%s/runJob_%d.sh\'" % (opt.queue,jobCounter,jobsBase,jobCounter) )
                os.system("bsub -q %s -R \"swp>1000 && pool>30000\" -J PIXEL_%d \'%s/runJob_%d.sh\'" % (opt.queue,jobCounter,jobsBase,jobCounter) )

            jobCounter+=1

