#!/usr/bin/env python

import os,sys
import json
import ROOT
                
if(len(sys.argv)<3):
    print 'runKinOverSamples.py subToBatch samples.json eventsPerJob'
    exit(-1)

#open the file which describes the sample
subToBatch=int(sys.argv[1])
samplesDB = sys.argv[2]
jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()
evPerJob=int(sys.argv[3])
scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapKinAnalysisRun.sh')

params=''
if(len(sys.argv)>4) :
    for i in xrange(5,len(sys.argv)) :
        params += sys.argv[i] + ' '

#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :
        
        #run over items in process
        data = desc['data']
        for d in data :
            dtag = d['dtag']
            try:
            	dir = d['summarydir']
            except:
		continue
	    fileName=dir + '/' + dtag + '.root'
            print "*****"

            #check number of jobs
            njobs=1
            if(evPerJob>0) :
                fin = ROOT.TFile.Open(fileName)
                data=fin.Get("evAnalyzer/data")
                nentries=data.GetEntriesFast()
                fin.Close()
                njobs=nentries/evPerJob+1
                print "Nentries=" + str(nentries)

            #submit jobs    
            for ijob in range(njobs) :
                evStart= ijob*evPerJob
                evEnd= (ijob+1)*evPerJob
                params = '-src=' + fileName + ' -f=' + str(evStart) + ' -e=' + str(evEnd)
                if(subToBatch==0) :
                    os.system(scriptFile + ' '  + params)
                else :
                    os.system('submit2batch.sh ' + scriptFile + ' ' + params)
