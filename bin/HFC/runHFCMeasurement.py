
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
    print '  -n : number of ensembles to test'
    print '  -b : jet bin (2=2 jets exclusive, 3=3 jets inclusive,4=4 jets inclusive, otherwise all inclusive)'
    print '  -f : fit type'
    print '  -p : plot only'
    exit(-1)


#parse the options
try:
    # retrive command line options
    shortopts  = "b:l:j:i:n:p:f:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

#configure
lumi=100
samplesDB=''
ifile=''
nensemble=1
plotOnly=False
jetbin=0
fitType=-1
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-l'): lumi=float(a)
    elif o in('-j'): samplesDB = a
    elif o in('-i'): ifile=a
    elif o in('-n'): nensemble=int(a)
    elif o in('-f'): fitType=int(a)
    elif o in('-p'): plotOnly = (a=='True')
    elif o in('-b'): jetbin=int(a)

# load macros
import ROOT
ROOT.gSystem.Load('${CMSSW_BASE}/lib/${SCRAM_ARCH}/libLIPTop.so')
from ROOT import HFCMeasurement, MisassignmentMeasurement, EventSummaryHandler, getNewCanvas, showPlotsAndMCtoDataComparison, setStyle, formatForCmsPublic, formatPlot

if(not plotOnly) :
    inF=ROOT.TFile.Open(ifile)

    jsonFile = open(samplesDB,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()

    #run over sample
    evHandler       = EventSummaryHandler()
    misMeasurement  = MisassignmentMeasurement()

    misMeasurement.setBiasCorrections("all",-0.0188118558)
    misMeasurement.setBiasCorrections("emu",0.0088410018)
    misMeasurement.setBiasCorrections("ll",-0.0334519059)
    if(jetbin==2) :
        misMeasurement.setBiasCorrections("all",0.0208030449)
        misMeasurement.setBiasCorrections("emu",-0.0345264636)
        misMeasurement.setBiasCorrections("ll",-0.018515963)
    elif(jetbin==3) :
        misMeasurement.setBiasCorrections("all",-0.032760782)
        misMeasurement.setBiasCorrections("emu",-0.0129281559)
        misMeasurement.setBiasCorrections("ll",-0.0758820633)
    elif(jetbin==4):
        misMeasurement.setBiasCorrections("all",-0.0606490967)
        misMeasurement.setBiasCorrections("emu",-0.0220872477)
        misMeasurement.setBiasCorrections("ll",0.0904430383)
        #    else :
    misMeasurement.setBiasCorrections("all",0)
    misMeasurement.setBiasCorrections("emu",0)
    misMeasurement.setBiasCorrections("ll",0)
        
    hfcFitter=None
    if(fitType>=0):
        hfcFitter = HFCMeasurement(3,fitType)
    
    ensembleInfo = '<big><b>Summary of ensemble tests (data/MC) for events with'
    if(jetbin==2):   ensembleInfo += '= 2 jets'
    elif(jetbin==3): ensembleInfo += '= 3 jets'
    elif(jetbin==4): ensembleInfo += '= 4 jets'
    else :           ensembleInfo += '&gt;= 2 jets'
    ensembleInfo += '</b></big>'
    print '[HFCmeasurement]',
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
                    if(ipe>0) : ensembleInfo += '<tr><td>&lt;N<sub>expected</sub>&gt;=' + str(nevtsExpected) + '</td></tr>'
                    ensembleInfo += '<tr><td>&lt;N<sub>generated</sub>&gt;=' + str(nevts) + '</td></tr>'
                    #if(ipe>0) : ensembleInfo += '<tr><td><small>' + str(genEvts) + '</small></td></tr>'
                
            #take control of the filled tree now
            ensembleHandler.attachToTree( ensembleHandler.getTree() )
            misMeasurement.measureMisassignments( ensembleHandler, 180, 40, (ipe==0), jetbin )

            #print info
            categories=['all', 'emu','ll']
            for cat in categories : 
                fcorrectTrue  = misMeasurement.getTrueCorrectPairsFraction(cat)
                fcorrectEst   = misMeasurement.getCorrectPairsFraction(cat)
                fcorrect      = fcorrectEst[0]
                fcorrectErr   = fcorrectEst[1]
                kNorm         = misMeasurement.getNorm(cat)

                if(hfcFitter is not None and cat =='all'):

                    #configure fit from file
                    fitParamsFile = open('hfcFitter_cfg.json','r')
                    fitParams=json.load(fitParamsFile,encoding='utf-8')

                    btagWP='TCHEL'

                    btagAlgos=fitParams['btagalgos']
                    hfcFitter.configureBtagAlgo   (btagWP,btagAlgos[btagWP]['cut'])
                    hfcFitter.setBtagEfficiency   (btagAlgos[btagWP]['effb'][0], btagAlgos[btagWP]['sfb'][0], btagAlgos[btagWP]['sfb'][1] )
                    hfcFitter.setMistagEfficiency (btagAlgos[btagWP]['effq'][0], btagAlgos[btagWP]['sfq'][0], btagAlgos[btagWP]['sfq'][1] )

                    catParams=fitParams[cat]
                    hfcFitter.setSelectionFractions (catParams['inclusive']['fcorrect'][0],catParams['inclusive']['fcorrect'][1],
                                                     catParams['inclusive']['fttbar'][0],catParams['inclusive']['fttbar'][1],
                                                     catParams['inclusive']['fsingletop'][0],catParams['inclusive']['fsingletop'][1],0)
                    hfcFitter.setSelectionFractions (catParams['2jets']['fcorrect'][0],catParams['2jets']['fcorrect'][1],
                                                     catParams['2jets']['fttbar'][0],catParams['2jets']['fttbar'][1],
                                                     catParams['2jets']['fsingletop'][0],catParams['2jets']['fsingletop'][1],2)
                    hfcFitter.setSelectionFractions (catParams['3jets']['fcorrect'][0],catParams['3jets']['fcorrect'][1],
                                                     catParams['3jets']['fttbar'][0],catParams['3jets']['fttbar'][1],
                                                     catParams['3jets']['fsingletop'][0],catParams['3jets']['fsingletop'][1],3)
                    fitParamsFile.close()

                    #run the fitter
                    hfcFitter.fitHFCtoEnsemble( ensembleHandler, cat )
                    

                    raw_input('any key to continue')
                #save ensemble info
                ensembleInfo += '<tr><td><i>' + cat + ' events</i></td></tr>'
                ensembleInfo += '<tr><td>k-factor=' + str(kNorm) + '</td></tr>'
                ensembleInfo += '<tr><td>f<sub>correct</sub>=' + str(fcorrect) + "+/-" + str(fcorrectErr)
                if(ipe>0) : ensembleInfo += ' (MC truth=' + str(fcorrectTrue)+ ')'
                ensembleInfo += '</td><tr>'
                ensembleInfo += '<tr><td></td></tr>'
            
            ensembleHandler.getTree().Delete("all")
    misMeasurement.saveMonitoringHistograms()
    inF.Close()

    #save ensemble info to file
    fout = open('ensemble_info.html', 'w')
    fout.write(ensembleInfo)
    fout.close()
                

#display results
plotF = ROOT.TFile.Open("MisassignmentMeasurement.root")
cats=['all','emu','ll']
setStyle()
cnv = getNewCanvas('misc','misc',False)
cnv.SetWindowSize(1500,500)
cnv.Divide(3,1)
cnv2 =getNewCanvas('misresc','misresc',False)
cnv2.SetWindowSize(1500,500)
cnv2.Divide(3,1)
cnv3 =getNewCanvas('fractfitc','fracfitc',False)
cnv3.SetWindowSize(1500,500)
cnv3.Divide(3,1)
fracFitter=[]
icnv=0
for c in cats:
    prefix=''
    if(c!='all') : prefix=c+'_'
    icnv = icnv+1

    #inclusive distributions
    cnv.cd(icnv)
    spimpose = ROOT.TList()
    mcmodelH=plotF.Get('localAnalysis/'+c+'/'+prefix+'avgwrongmodelmlj')
    mcmodelH.SetTitle('MC model')
    formatPlot(mcmodelH,8,9,2,1,0,True,False,8,8,8)
    datamodelH =plotF.Get('localAnalysis/'+c+'/'+prefix+'datawrongmodelmlj')
    datamodelH.SetTitle('data model')
    formatPlot(datamodelH,1,9,2,1,0,True,False,1,1,1)
    spimpose.Add(mcmodelH)
    spimpose.Add(datamodelH)
    
    stack    = ROOT.TList()
    mcWrongH=plotF.Get('localAnalysis/'+c+'/'+prefix+'avgwrongmlj')
    formatPlot(mcWrongH,809,1,1,0,1001,True,False,1,809,809)
    mcWrongH.SetTitle('Wrong assignments')
    mcCorrectH=plotF.Get('localAnalysis/'+c+'/'+prefix+'avgcorrectmlj')
    formatPlot(mcCorrectH,614,1,1,0,1001,True,False,1,614,614)
    mcCorrectH.SetTitle('Correct assignments')
    stack   .Add(mcWrongH)
    stack   .Add(mcCorrectH)

    data    = ROOT.TList()
    dataH=plotF.Get('localAnalysis/'+c+'/'+prefix+'datainclusivemlj')
    dataH.SetDirectory(0)
    dataH.SetTitle('data')
    data    .Add(dataH)
    
    pad=cnv.cd(icnv)
    leg=showPlotsAndMCtoDataComparison(pad,stack,spimpose,data)
    subpad1=pad.cd(1)    
    subpad1.SetLogx()
    subpad2=pad.cd(2)
    subpad2.SetLogx()
    if(icnv==1) : formatForCmsPublic(pad.cd(1),leg,'CMS preliminary, #sqrt{s}=7 TeV, #int L=%3.0f pb^{-1}' % lumi ,2)
    else        : leg.Delete()

    #subtracted distributions
    cnv2.cd(icnv)

    spimpose2    = ROOT.TList()
    mcSubtractedH=plotF.Get('localAnalysis/'+c+'/'+prefix+'avginclusivemlj').Clone('mcres')
    mcSubtractedH.Add(mcmodelH,-1)
    mcSubtractedH.SetTitle('MC residual')
    formatPlot(mcSubtractedH,8,9,2,1,0,True,False,8,8,8)
    spimpose2.Add(mcSubtractedH)

    data2 = ROOT.TList()
    dataSubtractedH = dataH.Clone('datares')
    dataSubtractedH.Add(datamodelH,-1)
    dataSubtractedH.SetTitle('data residual')
    dataSubtractedH.SetDirectory(0)
    data2.Add(dataSubtractedH)
    
    stack2 = ROOT.TList()
    stack2.Add(mcCorrectH)
    
    pad2=cnv2.cd(icnv)
    leg=showPlotsAndMCtoDataComparison(pad2,stack2,spimpose2,data2)
    subpad12=pad2.cd(1)
    subpad12.SetLogx()
    subpad22=pad2.cd(2)
    subpad22.SetLogx()
    if(icnv==1) : formatForCmsPublic(pad2.cd(1),leg,'CMS preliminary, #sqrt{s}=7 TeV, #int L=%3.0f pb^{-1}' % lumi ,2)
    else        : leg.Delete()

    #fraction fit
    pad3=cnv3.cd(icnv)

    mcArray = ROOT.TObjArray(3);
    mcArray.Add(mcCorrectH)
    mcArray.Add(mcWrongH)
    mfracFitter=ROOT.TFractionFitter(dataH, mcArray)
    fracFitter.append( mfracFitter );
    mfracFitter.Constrain(1,0.0,1.0);      
    #fracFitter.SetRangeX(1,15);  #bins to use
    fracFitterStatus = mfracFitter.Fit()
    sigval=ROOT.Double(0.0)
    sigvalerr=ROOT.Double(0.0)
    bkgval=ROOT.Double(0.0)
    bkgvalerr=ROOT.Double(0.0)
    if (fracFitterStatus == 0) :
        resultH = mfracFitter.GetPlot()
        resultH.SetDirectory(0)
        mfracFitter.GetResult(0, sigval, sigvalerr);
        mfracFitter.GetResult(1, bkgval, bkgvalerr)
        resultH.SetLineColor(mcCorrectH.GetLineColor())
        resultH.SetMarkerColor(mcCorrectH.GetMarkerColor())
        resultH.SetFillColor(mcCorrectH.GetFillColor())
        resultH.SetLineStyle(mcCorrectH.GetLineStyle())
        resultH.SetMarkerStyle(mcCorrectH.GetMarkerStyle())
        resultH.SetFillStyle(mcCorrectH.GetFillStyle())
        resultH.Draw("hist")
        mcWrongH.DrawNormalized("histsame",resultH.Integral()*bkgval);
        dataH.Draw("e1psame")

    if(icnv==1) : formatForCmsPublic(pad3,leg,'CMS preliminary, #sqrt{s}=7 TeV, #int L=%3.0f pb^{-1}' % lumi ,2)
    #else        : leg.Delete()

    fcorrsfH=plotF.Get('localAnalysis/'+c+'/'+prefix+'fcorrsf')
    fcorrH=plotF.Get('localAnalysis/'+c+'/'+prefix+'fcorr')
    truefCorrH=plotF.Get('localAnalysis/'+c+'/'+prefix+'truefcorr')
    avgFcorr=fcorrH.GetMean()
    fCorrErr=fcorrH.GetRMS()
    avgFcorrTrue=truefCorrH.GetMean()
    fCorrTrueErr=truefCorrH.GetRMS()
    biasH=plotF.Get('localAnalysis/'+c+'/'+prefix+'bias')
    biasFit=biasH.GetFunction('gaus')
    fCorrBias=biasFit.GetParameter(1)
    fCorrBiasErr=biasFit.GetParError(1)

    #debug
    print c,
    #print ' f_{correct}=' + str(avgFcorr) + ' +/- ' + str(fCorrErr),
    print ' <f_{correct}^{MC}>=' + str(avgFcorrTrue) + ' +/- ' + str(fCorrTrueErr),
    print ' bias=' + str(fCorrBias) + ' +/- ' + str(fCorrBiasErr),
    print ' f_{correct}^{fraction fitter}= ' + str(sigval) + ' +/- ' + str(sigvalerr)
    # print ' SF= ' + str(fcorrsfH.GetMean()) + ' +/- ' + str(fcorrsfH.GetRMS())
    
cnv.cd()
cnv.Draw()
cnv.Modified()
cnv.Update()
cnv.SaveAs('mljwithmodel.C')

cnv2.cd()
cnv2.Draw()
cnv2.Modified()
cnv2.Update()
cnv2.SaveAs('mljsubtracted.C')

    
raw_input(' *** Any key to end')    
plotF.Close()

