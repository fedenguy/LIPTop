#!/usr/bin/env python

import os,sys
import json
import getopt
from math import sqrt,pow

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal
    
#print usage
def usage() :
    print ' '
    print 'runMassMeasurement.py [options]'
    print '  -j : json file containing the samples'
    print '  -i : event summary file'
    print '  -p : fit parameters file'
    exit(-1)

#parse the options
try:
    # retrive command line options
    shortopts  = "j:i:p:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

#configure
samplesDB=''
ifile=''
fitParsFile='MassParFile_exclusive.txt'
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-j'): samplesDB = a
    elif o in('-i'): ifile=a
    elif o in('-p'): fitParsFile=a

    
# load macros
import ROOT
ROOT.gSystem.Load('${CMSSW_BASE}/lib/${SCRAM_ARCH}/libLIPTop.so')
from ROOT import MassMeasurement, EnsembleMeasurement_t, eventHandlerFactory

inF=ROOT.TFile.Open(ifile)

jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#run over sample
evHandler       = eventHandlerFactory()
massFitter      = MassMeasurement(fitParsFile,'data')
ensemble=EnsembleMeasurement_t()
ensemble.nEvents=0
ensemble.status=False
ensemble.mass=0
ensemble.err=0
evtCtr=0
for proc in procList :
    for desc in proc[1] :
        isdata = getByLabel(desc,'isdata',False)
        if(not isdata): continue
        data = desc['data']
        for d in data :

            tag = getByLabel(d,'dtag','')
            t=inF.Get(tag+'/data')
            try :
                t.GetEntriesFast()
            except:
                continue
            
            attResult=evHandler.attachToTree(t)
            nevtsSel = evHandler.getEntries()
            if(attResult is False) : continue
                    
            for ievt in xrange(0,nevtsSel) :
                evHandler.getEntry(ievt)
                mtop=evHandler.evSummary_.evmeasurements[0]
                if(mtop==0) : continue
                nbtags=evHandler.evSummary_.evmeasurements[5]
                if(nbtags==0) : continue
                if(nbtags>2) : nbtags=2
                isZcand=evHandler.evSummary_.evmeasurements[2]
                if(isZcand) : continue
                catToFill=nbtags-1
                if(evHandler.evSummary_.cat==3) : catToFill+=2
                ensemble.evMasses.__setitem__(evtCtr,mtop)
                ensemble.evCategories.__setitem__(evtCtr,catToFill)
                evtCtr+=1                                                     
ensemble.nEvents=evtCtr
ensembleMeasurement = massFitter.DoMassFit(ensemble,True)
sys.exit(-1)


