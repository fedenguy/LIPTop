#!/usr/bin/env python

import ROOT
import os,sys
import getopt
from math import sqrt,pow,fabs
myscriptpath=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}')
sys.path.append( myscriptpath )
from rounding import toLatexRounded

"""
help 
"""
def usage() :
    print ' '
    print 'getEventYields.py [options]'
    print '  -i : input file'
    print '  -c : counting histo (=evtflow by default)'
    print '  -b : bin to summarize'
    print '  -m : mode=simple,std,chhiggs,syst'
    print '  -s : scale factors for DY emu/mumu/ee e.g. -s  1.0/1.9/2.1 to be applied'
    print '  -o : overall trigger efficiencies for emu/umu/ee e.g -o 0.946/0.920/0.966 to be applied'
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
                fmtValue = '%3.0f' % val
                #if(addErrors) :
                #    fmtValue += '$\\pm$'
                #    fmtValue += '%3.0f' % valerr
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
    shortopts  = "i:c:b:m:s:o:h?"
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

inputFile='plotter.root'
evCats   = ['emu','mumu','ee']
mode='std'
countHisto='evtflow'
ibinSummary=5
dySfactors = { "emu":1, "mumu":1, "ee":1 }
trigEffs   = { "emu":1, "mumu":1, "ee":1 }
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(1)
    elif o in('-i'): inputFile = a
    elif o in('-m'): mode=a
    elif o in('-c'): countHisto = a
    elif o in('-b'): ibinSummary = int(a)
    elif o in('-s'):
        sfactors=a.split('/')
        for i in xrange(0,len(sfactors)-1) : dySfactors[evCats[i+1]]=float(i)
    elif o in ('-o'):
        treffs=a.split('/')
        for i in xrange(0,len(treffs)-1) : trigEffs[evCats[i+1]]=float(i)

#analysis specific configuration
procs=[]
systs=[]
singleSyst=[]
if(mode=='std' or mode=='simple') :
    procs=[ "Di-bosons", "Single top", "W+jets", "Z-#gamma^{*}+jets#rightarrow ll", "other t#bar{t}", "t#bar{t} dileptons", "data" ]
    if(mode=='std') :
        systs     = ['jer', 'jesup','jesdown','puup', 'pudown']
        singleSyst= [True,  False,   False,    False,  False]
elif(mode=='chhiggs') :
    procs = [ "WH80", "HH80", "WH90", "HH90", "WH100","HH100", "WH120","HH120", "WH140","HH140", "WH150","HH150", "WH155","HH155", "WH160","HH160" ]
    systs = ['jer', 'jesup','jesdown','puup', 'pudown']
    singleSyst = [True,  False,   False,    False,  False]
elif(mode=='syst'):
    procs=[ "t#bar{t} matching down", "t#bar{t} matching up", "t#bar{t} scale down", "t#bar{t} scale up" ]

if(len(procs)==0) :
    usage()
    sys.exit(1)

#now get plots from file
f    = ROOT.TFile.Open(inputFile)
href = f.Get(procs[0]+'/'+countHisto)

evTitles = ['$e\\mu$','$\\mu\\mu$','$ee$']
hsummary = ROOT.TH2F("summaryyields","summaryyields",len(evCats),0,len(evCats),len(procs),0,len(procs))

#create tables for each event category
for iec in xrange(0,len(evCats)):
    ec=evCats[iec]
    prefix=ec+'_'

    #prepare the total expected from simulation and the variations histogram
    hcatsum = href.Clone(prefix+'totalmc')
    hcatsum.SetTitle('Total expected')
    hcatsum.Reset('ICE')
    systVarH = ROOT.TH2F(prefix+'systvar',prefix+'systvar',len(procs),0,len(procs),len(systs),0,len(systs)) 

    #init tex table
    tabtex = '\\documentclass[]{beamer}\n'
    tabtex += '\\usepackage[utf8]{inputenc}\n'
    tabtex += '\\begin{document}\n'
    tabtex += '\\begin{table}[htp]\n'
    tabtex += '\\caption{'+evTitles[iec]+' yields for an integrated luminosity of $2165^{-1}$pb}\n'
    tabtex += '\\label{tab:'+ec+'yields}\n'
    tabtex += '\\hspace*{-0.9cm}\n'
    tabtex += '\\scalebox{0.45}{\n'
    colfmt='l'
    colnames=''
    for ibin in xrange(1,href.GetXaxis().GetNbins()+1) :
        colfmt += 'c'
        binlabel='$'+href.GetXaxis().GetBinLabel(ibin)+'$'
        binlabel=binlabel.replace('#','\\')
        binlabel=binlabel.replace(' ','~')
        colnames += ' & ' + binlabel
    tabtex += '\\begin{tabular}{'+colfmt+'} \\hline\n'
    tabtex += 'Process ' + colnames + '\\\\ \\hline\\hline\n'

    for procCtr in xrange(0,len(procs)) :
        p = procs[procCtr]
        h=f.Get(p + "/" +prefix+countHisto)

        #correct the yields for the summary bin
        sf=1.0
        if(p.find("Z-#gamma^{*}+jets#rightarrow ll")>=0 and len(ec)>0) : sf=dySfactors[ec]
        te=1.0
        if(p!='data' and len(ec)>0): te=trigEffs[ec]
        h.SetBinContent(ibinSummary, h.GetBinContent(ibinSummary)*sf*te)
        h.SetBinError(ibinSummary,   h.GetBinError(ibinSummary)*sf*te)
        hsummary.SetBinContent(iec+1, procCtr+1,  h.GetBinContent(ibinSummary))
        hsummary.SetBinError(iec+1,   procCtr+1,  h.GetBinError(ibinSummary))
        hsummary.GetXaxis().SetBinLabel(iec+1,ec)
        hsummary.GetYaxis().SetBinLabel(procCtr+1,p)

        #increment MC expectations until data is found
        if(p!='data') :
            hcatsum.Add(h)
            systCtr=0
            for isyst in systs:
                hvar=f.Get(p + "/" +prefix+countHisto+isyst)
                if(hvar==None):continue
                yieldDiff=0
                if(h.GetBinContent(ibinSummary)>0):
                    yieldDiff=(hvar.GetBinContent(ibinSummary)*sf*te-h.GetBinContent(ibinSummary))
                    if(singleSyst[systCtr]) : yieldDiff = fabs(yieldDiff)*0.5
                    yieldDiff=100*yieldDiff/h.GetBinContent(ibinSummary)
                systVarH.Fill(procCtr,systCtr,yieldDiff)
                systCtr+=1
        else :
            tabtex += '\\hline\n'
            tabtex +=getInTableFormat(hcatsum.GetTitle(),hcatsum) + '\\\\\\hline\n'

        if(p=='data') : tabtex += getInTableFormat(p,h,True) + '\\\\\\hline\\hline\n'            
        else          : tabtex += getInTableFormat(p,h) + '\\\\\n'

    tabtex += '\\end{tabular}\n'
    tabtex += '}n'
    tabtex += '\\end{table}\n'
    tabtex += '\\end{document}\n'

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
    if(len(systs)>0) :
        fileObj.write('\n\n')
        fileObj.write(systTabtex)
    fileObj.close()

#dump the summary table
tabtex ='\\documentclass[]{beamer}\n'
tabtex +='\\usepackage[utf8]{inputenc}\n'
tabtex +='\\begin{document}\n'
tabtex +='\\begin{table}[htp]\n'
bintitle=href.GetXaxis().GetBinLabel(ibinSummary)
tabtex +='\\caption{$'+bintitle+'$ summary events for $2165^{-1}$pb}\n'
tabtex +='\\label{tab:summaryyields}\n'
tabtex +='\\hspace*{-0.9cm}\n'
tabtex +='\\scalebox{0.70}{\n'
colfmt='lc'
for ec in evTitles:
    colfmt += 'c' 
    colnames +=' & ' + ec
colnames=' & Total'
tabtex += '\\begin{tabular}{'+colfmt+'} \\hline\n'
tabtex += 'Channel '+colnames+'\\\\ \\hline\\hline\n'
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
tabtex +='\\end{tabular}\n'
tabtex +='}\n'
tabtex +='\\end{table}\n'
tabtex +='\\end{document}\n'

#save to file
fileObj = open('summaryyields.tex',"w")
fileObj.write(tabtex)
fileObj.close()
 


f.Close()


    
