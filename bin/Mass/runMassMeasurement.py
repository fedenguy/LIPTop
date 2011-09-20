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
    print '  -l : luminosity (pb)'
    print '  -j : json file containing the samples'
    print '  -i : event summary file'
    print '  -n : number of ensembles to test'
    print '  -p : fit parameters file'
    exit(-1)

#parse the options
try:
    # retrive command line options
    shortopts  = "l:j:i:n:p:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

#configure
lumi=200
samplesDB=''
ifile=''
nensemble=1
fitParsFile='MassParFile_bmult.txt'
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-l'): lumi=float(a)
    elif o in('-j'): samplesDB = a
    elif o in('-i'): ifile=a
    elif o in('-n'): nensemble=int(a)
    elif o in('-p'): fitParsFile=a

    
# load macros
import ROOT
ROOT.gSystem.Load('${CMSSW_BASE}/lib/${SCRAM_ARCH}/libLIPTop.so')
from ROOT import MassMeasurement, EventSummaryHandler, getNewCanvas, showPlotsAndMCtoDataComparison, setStyle, formatForCmsPublic, formatPlot

inF=ROOT.TFile.Open(ifile)

jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#run over sample
evHandler       = EventSummaryHandler()
massFitter      = MassMeasurement(fitParsFile)

ensembleInfo = '<big><b>Summary of ensemble tests (data/MC)</b></big>'
for ipe in xrange(0,nensemble+1) :
        
    #progress bar
    print '.',
    sys.stdout.flush()
        
    ensembleHandler = EventSummaryHandler()
    ensembleInfo  += '<table>'
    ensembleInfo += '<tr><th><b>Info for '
    if(ipe==0) : ensembleInfo += 'data'
    else       : ensembleInfo += ' ensemble #' + str(ipe)
    ensembleInfo += '</b></th></tr>'
        
    for proc in procList :

        #run over processes
        id=0
        for desc in proc[1] :
            isdata = getByLabel(desc,'isdata',False)
                
            if(ipe>0 and isdata) : continue
            if(ipe==0 and not isdata) : continue
            
            #run over items in process
            data = desc['data']
            for d in data :
                tag = getByLabel(d,'dtag','')
 
                #get tree of events from file
                t=inF.Get(tag+'/data')
                    
                try :
                    t.GetEntriesFast()
                except:
                    continue
            
                attResult=evHandler.attachToTree(t)
                nevtsSel = evHandler.getEntries()
                if(attResult is False) : continue
                    
                #clone (will use the same address as the original tree)
                id=id+1
                if(id==1):
                    ROOT.gROOT.cd()
                    ensembleHandler.initTree(t.CloneTree(0), False)
                    ensembleHandler.getTree().SetDirectory(0)
                    

                #generate number of events for ensemble    
                nevts=nevtsSel
                nevtsExpected=nevts
                if(ipe>0):
                    evHandler.getEntry(0)
                    nevtsExpected=lumi*(evHandler.evSummary_.weight)*nevtsSel
                    nevts = int(ROOT.gRandom.Poisson( nevtsExpected ))

                genEvts=[]
                for ievt in xrange(0,nevts) :
                    rndEvt=ievt
                    if(ipe>0):
                        rndEvt = int(ROOT.gRandom.Uniform(0,nevtsSel))
                        if( rndEvt in genEvts ) : continue
                    evHandler.getEntry(rndEvt)
                    ensembleHandler.fillTree()
                    genEvts.append(rndEvt)

                ensembleInfo += '<tr><td><b>' + tag + '<b></td></tr>'
                if(ipe>0) :
                    ensembleInfo += '<tr><td>&lt;N<sub>expected</sub>&gt;=' + str(nevtsExpected) + '</td></tr>'
                    ensembleInfo += '<tr><td>&lt;N<sub>generated</sub>&gt;=' + str(nevts) + '</td></tr>'
                else :
                    ensembleInfo += '<tr><td>N<sub>events</sub>=' + str(nevtsExpected) + '</td></tr>'
                
        #take control of the filled tree now
        print ensembleHandler.getTree().GetEntriesFast()
        ensembleHandler.attachToTree( ensembleHandler.getTree() )
        ensembleMeasurement = massFitter.DoMassFit(ensembleHandler,True)

        raw_input('any key to continue')
        #save ensemble info
        #ensembleInfo += '<tr><td><i>' + cat + ' events</i></td></tr>'
        #ensembleInfo += '<tr><td>k-factor=' + str(kNorm) + '</td></tr>'
        #ensembleInfo += '<tr><td>f<sub>correct</sub>=' + str(fcorrect) + "+/-" + str(fcorrectErr)
        #if(ipe>0) : ensembleInfo += ' (MC truth=' + str(fcorrectTrue)+ ')'
        #ensembleInfo += '</td><tr>'
        #ensembleInfo += '<tr><td></td></tr>'

        ensembleHandler.getTree().Delete("all")
        


#save ensemble info to file
fout = open('ensemble_info.html', 'w')
fout.write(ensembleInfo)
fout.close()
                


