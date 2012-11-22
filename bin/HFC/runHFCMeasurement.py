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
    print 'runHFCMeasurement.py [options]'
    print '  -l : luminosity (pb)'
    print '  -j : json file containing the samples'
    print '  -i : event summary file'
    print '  -t : event summary type (0-trees, 1-histos)'
    print '  -f : fit type'
    print '  -a : algo=TCHEL'
    print '  -c : override discriminator cut for algo'
    print '  -m : override mistag rate for algo'
    exit(-1)


#parse the options
try:
    # retrive command line options
    shortopts  = "b:l:j:i:n:p:f:t:a:c:m:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)


#configure
btagWP='TCHEL'
lumi=100
samplesDB=''
ifile=''
fileType=0
jetbin=0
fitType=-1
btagWPcut=-1
mistagRate=-1
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-l'): lumi=float(a)
    elif o in('-j'): samplesDB = a
    elif o in('-i'): ifile=a
    elif o in('-t'): fileType=int(a)
    elif o in('-f'): fitType=int(a)
    elif o in('-b'): jetbin=int(a)
    elif o in('-a'): btagWP=a
    elif o in('-c'): btagWPcut=float(a)
    elif o in('-m'): mistagRate=float(a)
    
# load macros
import ROOT
ROOT.gSystem.Load('${CMSSW_BASE}/lib/${SCRAM_ARCH}/libLIPTop.so')
from ROOT import HFCMeasurement, eventHandlerFactory


inF=ROOT.TFile.Open(ifile)
jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#create the fitter and configure fit from file
hfcFitter = HFCMeasurement(fitType)

#fitParamsFile = open('hfcFitter_cfg.json','r')
fitParamsFile = open('hfcFitter_2012_cfg.json','r')
fitParams=json.load(fitParamsFile,encoding='utf-8')

vareb=0
vareq=-1
varfcorrect=0
varttbar=0
varst=0

btagAlgos=fitParams['btagalgos']
if(btagWPcut<0) : btagWPcut=btagAlgos[btagWP]['cut']
hfcFitter.configureBtagAlgo   (btagWP,btagWPcut)
hfcFitter.setBtagEfficiency   (btagAlgos[btagWP]['effb'][0], btagAlgos[btagWP]['sfb'][0]+vareb*btagAlgos[btagWP]['sfb'][1], btagAlgos[btagWP]['sfb'][1] )
mistagAlgo=btagAlgos[btagWP]['effq'][0]
if(mistagRate>0) : mistagAlgo=mistagRate
hfcFitter.setMistagEfficiency (mistagAlgo, btagAlgos[btagWP]['sfq'][0]+vareq*btagAlgos[btagWP]['sfq'][1] , btagAlgos[btagWP]['sfq'][1] )


for ch in ['emu','mumu','ee'] :
    catParams=fitParams[ch]
    for jbin in [2,3] :
        key=str(jbin)+'jets'

        localEffb=btagAlgos[btagWP]['effb'][0]
        try :
            locallEffb=catParams[key][btagWP+'effb']
        except:
            print 'Assuming global eff_b'

        localEffq=btagAlgos[btagWP]['effq'][0]
        try :
            locallEffq=catParams[key][btagWP+'effq']
        except:
            print 'Assuming global eff_q'


        hfcFitter.setParametersForCategory(catParams[key]['fcorrect'][0]  +varfcorrect*catParams[key]['fcorrect'][1],     catParams[key]['fcorrect'][1],
                                           catParams[key]['fttbar'][0]    +varttbar*catParams[key]['fttbar'][1],       catParams[key]['fttbar'][1],
                                           catParams[key]['fsingletop'][0]+varst*catParams[key]['fsingletop'][1],   catParams[key]['fsingletop'][1],
                                           localEffb/btagAlgos[btagWP]['effb'][0], localEffq/btagAlgos[btagWP]['effq'][0],
                                           jbin, ch)

fitParamsFile.close()

print "****************************************************"
print "Using " + btagWP + " with cut " + str(btagWPcut)

#run over sample
if(fileType==0) :
    
    evHandler       = eventHandlerFactory()
    ensembleHandler = eventHandlerFactory()
    for proc in procList :
        
        #run over processes
        id=0
        for desc in proc[1] :
            isdata = getByLabel(desc,'isdata',False)
            if(not isdata) : continue
        
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
                for ievt in xrange(0,nevtsSel) :
                    evHandler.getEntry(ievt)
                    ensembleHandler.fillTree()
                    
    #take control of the filled tree now
    ensembleHandler.attachToTree( ensembleHandler.getTree() )
    #run the fitter
    hfcFitter.fitHFCtoEnsemble( ensembleHandler,1,True) #True )

#get the histo and fit it
else:
    
    btagHisto = inF.Get('data/'+btagWP+'btagsextended')
    hfcFitter.fitHFCtoMeasurement(btagHisto,1,True)


print "Dilepton exclusive"
for icat in xrange(0,6):
    print str(icat) + ' ' + str(hfcFitter.model.rFit[icat]) + ' +' + str(hfcFitter.model.rFitAsymmErrHi[icat]) + ' ' + str(hfcFitter.model.rFitAsymmErrLo[icat])
print "Dilepton combined"
for icat in xrange(0,3):
    print str(icat) + ' ' + str(hfcFitter.model.rCombFit[icat*2]) + ' +' + str(hfcFitter.model.rCombFitAsymmErrHi[icat*2]) + ' ' + str(hfcFitter.model.rCombFitAsymmErrLo[icat*2])
print "All combined"
print str(hfcFitter.model.rFitResult) + ' +' + str(hfcFitter.model.rFitResultAsymmErrHi) + ' ' + str(hfcFitter.model.rFitResultAsymmErrLo)
print "Feldman-Cousins: " +str(hfcFitter.model.rFitLowerLimit) + ' - ' + str(hfcFitter.model.rFitUpperLimit)

raw_input('*********************************')
#ensembleHandler.getTree().Delete("all")
#inF.Close()



