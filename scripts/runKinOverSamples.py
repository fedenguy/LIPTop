#!/usr/bin/env python

import os,sys
import json
import ROOT
import getopt
import os.path
import commands

#print usage
def usage(msg='') :
    print msg
    print ' '
    print 'runKinOverSamples.py [options]'
    print '  -s : submit or not to batch'
    print '  -j : json file containing the samples'
    print '  -e : number of events per job'
    print '  -i : initial event (-1 by default)'
    print '  -n : number of events to run per sample'
    print '  -d : directory containing the event summaries'
    print '  -t : select tag to submit'
    print '  -p : extra parameters to pass to KIN'
    print ' '
    exit(-1)

#parse the options
try:
    # retrive command line options
    shortopts  = "s:j:e:n:p:d:i:t:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
                                                              
#open the file which describes the sample
subToBatch=""
samplesDB = ''
evPerJob=-1
initEvent=-1
evToRun=-1
extraParams=''
summaryDir=''
tagSel=''
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-s'): subToBatch=a
    elif o in('-j'): samplesDB = a
    elif o in('-t'): tagSel=a
    elif o in('-e'): evPerJob=int(a)
    elif o in('-n'): evToRun=int(a)
    elif o in('-i'): initEvent=int(a)
    elif o in('-p'): extraParams = a
    elif o in('-d'): summaryDir = a

if(len(summaryDir)==0):
    usage('Summary dir is missing')
    exit(-1)
    
jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()
scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapKinAnalysisRun.sh')

#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :
        
        #run over items in process
        data = desc['data']
        for d in data :
            dtag = d['dtag']
            if(len(tagSel)>0):
                if(dtag.find(tagSel)<0) : continue
                
            fileName=summaryDir + '/' + dtag + '.root'
            if(fileName.find('/castor/cern.ch')==0) : fileName = 'rfio://' + fileName
            elif(fileName.find('/store/cmst3')==0)  : fileName = commands.getstatusoutput('cmsPfn ' + fileName)[1]

            print "***** " + fileName + " ******"

            #check number of jobs
            njobs=1
            if(evPerJob>0) :
                fin = ROOT.TFile.Open(fileName)
                data=fin.Get("evAnalyzer/data")
                nentries=data.GetEntriesFast()
                if(evToRun>0) : nentries=evToRun
                fin.Close()
                njobs=nentries/evPerJob+1
                print "Total entries=" + str(nentries)

            fileName=fileName.replace("&","\&")
            #submit jobs    
            for ijob in range(njobs) :
                evStart= ijob*evPerJob
                evEnd= (ijob+1)*evPerJob
                if(evStart<initEvent) : continue
                params = '-src=' + fileName + ' -f=' + str(evStart) + ' -e=' + str(evEnd) + ' ' + extraParams
                print params
                if(len(subToBatch)==0) :
                    os.system(scriptFile + ' '  + params)
                else :
                    os.system('submit2batch.sh -q' + subToBatch + ' ' + scriptFile + ' ' + params)
