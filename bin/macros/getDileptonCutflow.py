#!/usr/bin/env python

import ROOT
import os,sys
import getopt
from math import sqrt,pow,fabs
myscriptpath=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}')
sys.path.append( myscriptpath )
from rounding import toLatexRounded

def usage() :
    print ' '
    print 'getEventYields.py [options]'
    print '  -i : input file'
    print '  -c : counting histo (=evtflow by default)'
    print '  -b : bin to compute systematic variations for (=6 by default)'
    print '  -s : scale factors and uncertainties for DY emu/mumu/ee e.g. -s  1.2:0.26/1.34:0.38/1.38:0.39 to be applied to the summary yields'
    print '  -o : overall trigger efficiencies and uncertainties for all MC emu/umu/ee e.g. -o 0.946:0.008/0.920:0.003/0.966:0.008 to be applied to all yields'
    print '  '

"""
prints histogram to table row
""" 
def getInTableFormat(tit,h,isData=False,addErrors=True):
    tableRow='$'+tit+'$'
    tableRow=tableRow.replace('#','\\')
    tableRow=tableRow.replace(' ','~')
    for ibin in xrange(1,h.GetXaxis().GetNbins()+1) :
        val = h.GetBinContent(ibin)
        valerr= h.GetBinError(ibin)
        try:
            tableRow += ' & '
            fmtValue = str(int(val))
            if(not isData) :
                #fmtValue = '%3.0f' % val
                #fmtValue += '$\\pm$'
                #fmtValue += '%3.0f' % valerr
                fmtValue = '%3.2f' % val
                if(addErrors) :
                    fmtValue = toLatexRounded(val,valerr)
                    fmtValue = fmtValue.replace("[","(")
                    fmtValue = fmtValue.replace("]",")")
            tableRow += fmtValue
        except:
            tableRow += ' & '

    return tableRow

    
#parse the options
try:
    # retrive command line options
    shortopts  = "i:c:b:s:o:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)


evCats=['','emu','mumu','ee']
evTitles=['$Inclusive$','$e\\mu$','$\\mu\\mu$','$ee$']

inputFile='plotter.root'
countHisto='evtflow'
ibinSyst=6
dySfactors={
    "emu":[1.0,0.0],
    "mumu":[1.0,0.0],
    "ee":[1.0,0.0]
    }
trigEffs={
    "emu":[1.0,0.0],
    "mumu":[1.0,0.0],
    "ee":[1.0,0.0]
    }
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(1)
    elif o in('-i'): inputFile = a
    elif o in('-c'): countHisto = a
    elif o in('-b'): ibinSyst = int(a)
    elif o in('-s'):
        sfactors=a.split('/')
        for i in xrange(0,len(sfactors)) :
            sf=sfactors[i].split(':')
            val=float(sf[0])
            err=float(sf[1])
            dySfactors[evCats[i+1]]=[val,err]
    elif o in ('-o'):
        treffs=a.split('/')
        for i in xrange(0,len(treffs)) :
            te=treffs[i].split(':')
            val=float(te[0])
            err=float(te[1])
            trigEffs[evCats[i+1]]=[val,err]
            
# standard analysis
procs=[
    "Di-bosons",
    "Single top",
    "other t#bar{t}",
    "W+jets",
    "Z-#gamma^{*}+jets#rightarrow ll",
#    "t#bar{t} dileptons",
    "data"
    ]
systs     = ['jer', 'jesup','jesdown','puup', 'pudown']
singleSyst= [True,  False,   False,    False,  False]

# charged higgs
#procs=[
#    "WH80", "HH80",
#    "WH90", "HH90",
#    "WH100","HH100",
#    "WH120","HH120",
#    "WH140","HH140",
#    "WH150","HH150",
#    "WH155","HH155",
#    "WH160","HH160"
#    ]
#systs     = ['jer', 'jesup','jesdown','puup', 'pudown']
#singleSyst= [True,  False,   False,    False,  False]

# top systematics
#procs=["t#bar{t} matching down",
#       "t#bar{t} matching up",
#       "t#bar{t} scale down",
#       "t#bar{t} scale up",
#       "Single top (DS)"
#       ]
#systs=[]
#singleSyst=[]

f = ROOT.TFile.Open(inputFile)
href=f.Get(procs[0]+'/'+countHisto)
print href
hsummary    = ROOT.TH2F("summaryyields","summaryyields",len(evCats),0,len(evCats),len(procs),0,len(procs))

#create tables for each event category
iec=0
for ec in evCats:

    #init table
    colfmt='l'
    colnames=''
    for ibin in xrange(1,href.GetXaxis().GetNbins()+1) :
        colfmt += 'c'
        binlabel='$'+href.GetXaxis().GetBinLabel(ibin)+'$'
        binlabel=binlabel.replace('#','\\')
        binlabel=binlabel.replace(' ','~')
        colnames += ' & ' + binlabel
        
    prefix=ec
    if(len(prefix)>0): prefix = prefix+'_'

    #prepare the total expected from simulation
    hcatsum = href.Clone(prefix+'totalmc')
    hcatsum.SetTitle('Total expected')
    hcatsum.Reset('ICE')

    #prepare the variations histogram
    systVarH = ROOT.TH2F(prefix+'systvar',prefix+'systvar',len(procs),0,len(procs),len(systs),0,len(systs)) 

    #print to table
    tabtex =  '\\begin{table}[htp]\n'
    tabtex += '\\begin{center}\n'
    tabtex += '\\caption{$'+prefix+'-$yields}\n'
    tabtex += '\\label{tab:'+prefix+'yields}\n'
    tabtex += '\\tiny\n'
    tabtex += '\\begin{tabular}{'+colfmt+'} \\hline\n'
    tabtex += 'Process ' + colnames + '\\\\ \\hline\\hline\n'
    procCtr=0
    for p in procs:
        h=f.Get(p + "/" +prefix+countHisto)

        sf=1.0
        sferr=0.0
        if(p.find("Z-#gamma^{*}+jets#rightarrow ll")>=0) :
            if(len(ec)>0):
                sf=dySfactors[ec][0]
                sferr=dySfactors[ec][1]
                print 'Will apply DY-SF to ' + ec + ' ' +  str(sf) + ' +/- ' + str(sferr) 
        te=1.0
        teerr=0.0
        if(p!='data') :
            if(len(ec)>0):
                te=trigEffs[ec][0]
                teerr=trigEffs[ec][1]
                print p + ': Will apply trigger efficiency to ' + ec + ' ' + str(te) + ' +/- ' + str(teerr)
        hsummary.SetBinContent(iec+1,procCtr+1, h.GetBinContent(ibinSyst)*sf*te)
        tempErr=sqrt(pow(h.GetBinError(ibinSyst)*sf,2)+pow(h.GetBinContent(ibinSyst)*sferr,2))
        if(p!='data') :
            hsummary.SetBinError(iec+1,procCtr+1, sqrt(pow(tempErr*te,2)+pow(h.GetBinContent(ibinSyst)*teerr,2)) )
        else :
            hsummary.SetBinError(iec+1,procCtr+1,tempErr)
        hsummary.GetXaxis().SetBinLabel(iec+1,ec)
        hsummary.GetYaxis().SetBinLabel(procCtr+1,p)
        
        if(p!='data') :
            hcatsum.Add(h)
            systCtr=0
            for isyst in systs:
                hvar=f.Get(p + "/" +prefix+countHisto+isyst)
                if(hvar==None):continue
                yieldDiff=0
                if(h.GetBinContent(ibinSyst)>0):
                    yieldDiff=(hvar.GetBinContent(ibinSyst)-h.GetBinContent(ibinSyst))
                    if(singleSyst[systCtr]) : yieldDiff = fabs(yieldDiff)*0.5
                    yieldDiff=100*yieldDiff/h.GetBinContent(ibinSyst)
                systVarH.Fill(procCtr,systCtr,yieldDiff)
                systCtr+=1
        else :
            tabtex += '\\hline\n'
            tabtex +=getInTableFormat(hcatsum.GetTitle(),hcatsum) + '\\\\\\hline\n'

        if(p=='data') : tabtex += getInTableFormat(p,h,True) + '\\\\\\hline\\hline\n'            
        else          : tabtex += getInTableFormat(p,h) + '\\\\\n'

        procCtr+=1
        
    tabtex += '\\end{tabular}\n'
    tabtex += '\\end{center}\n'
    tabtex += '\\end{table}\n'

    #now print systematic variations
    systTabtex=''
    if(len(procs)>1):
        colnames=''
        colfmt='l'
        for isyst in systs: 
            colfmt +='c'
            colnames += ' & ' + isyst
        systTabtex =  '\\begin{table}[htp]\n'
        systTabtex += '\\begin{center}\n'
        systTabtex += '\\caption{$'+prefix+ '-$systematics on yields}\n'
        systTabtex += '\\label{tab:'+prefix+'systuncyields}\n'
        systTabtex += '\\begin{tabular}{'+colfmt+'} \\hline\n'
        systTabtex += 'Process' + colnames + '\\\\ \\hline\\hline\n'
        for iy in xrange(1,systVarH.GetYaxis().GetNbins()+2):
            projyH=systVarH.ProjectionY('_py',iy,iy)
            #        if(projyH.Integral()>0):
            systTabtex += getInTableFormat(procs[iy-1],projyH,False,False)
            systTabtex += '\\\\ \\hline\n'
        systTabtex += '\\hline\\hline\n'

        systTabtex += '\\end{tabular}\n'
        systTabtex += '\\end{center}\n'
        systTabtex += '\\end{table}\n'

    #save to file
    fileObj = open(prefix+'yields.tex',"w")
    fileObj.write(tabtex)
    fileObj.write('\n\n')
    fileObj.write(systTabtex)
    fileObj.close()
    iec +=1

#dump the summary table

colfmt='l'
colnames=''
for ec in evTitles:
    colfmt += 'c' 
    colnames +=' & ' + ec

tabtex =  '\\begin{table}[htp]\n'
tabtex += '\\begin{center}\n'
tabtex += '\\caption{}\n'
tabtex += '\\label{tab:summaryyields}\n'
tabtex += '\\begin{tabular}{'+colfmt+'} \\hline\n'
tabtex += 'Channel '+colnames+'\\\\ \\hline\\hline\n'

sumFinal = hsummary.ProjectionY('_py',2,len(evCats))
totFinalErr=ROOT.Double(0.0)
totFinal=sumFinal.IntegralAndError(1,len(procs)-1,totFinalErr)
hsummaryfinal    = ROOT.TH1F("summaryfinalyields","Total expected",len(evCats),0,len(evCats))
for ip in xrange(1,hsummary.GetYaxis().GetNbins()+1):
    projxH=hsummary.ProjectionX('_px',ip,ip)
    tit=hsummary.GetYaxis().GetBinLabel(ip)
    if(tit=='data') : 
        tabtex += '\\hline\n'
        tabtex += getInTableFormat(hsummaryfinal.GetTitle(),hsummaryfinal)
        tabtex += '\\\\\n\\hline\n'
    else :
        hsummaryfinal.Add(projxH)
    tabtex += getInTableFormat(tit,projxH,(tit=="data"))
    tabtex += '\\\\\n'
tabtex += '\\hline\\hline\n'

tabtex += '\\end{tabular}\n'
tabtex += '\\end{center}\n'
tabtex += '\\end{table}\n'


#save to file
fileObj = open('summaryyields.tex',"w")
fileObj.write(tabtex)
fileObj.close()
          


f.Close()


    
