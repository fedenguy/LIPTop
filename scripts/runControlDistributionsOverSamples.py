#!/usr/bin/env python
import os,sys
import json
import ROOT
                
if(len(sys.argv)<3):
    print 'runControlDistributionsOverSamples.py samples.json outdir [lumi=1]'
    exit(-1)

#open the file which describes the sample
samplesDB = sys.argv[1]
jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#configure
outdir=sys.argv[2]
lumi=1
if(len(sys.argv)>3) :
    lumi=int(sys.argv[3])

#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :
        
        #run over items in process
        isdata=desc['isdata']
        data = desc['data']
        for d in data :
            dtag = d['dtag']
            try :
                dir = d['summarydir']
            except :
                continue
            eventsFile=dir + '/' + dtag + '.root'
            os.system('showControlDistributions ' + eventsFile + ' ' + outdir + ' ' + str(int(not isdata)) )

#run plotter over results
os.system('mkdir -p ' + outdir + '/plots')
os.system('runPlotterOverSamples.py ' + samplesDB + ' ' + str(lumi) + ' ' + outdir + ' ' + outdir + '/ctrl ctrlAnalyzer')
