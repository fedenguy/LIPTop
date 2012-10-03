#/usr/bin/python

fInUrl='crossSectionSummary.dat'
f=open(fInUrl, 'r').readlines()

plotTitle  = [f[x].split(':') for x, line in enumerate(f) if "plotTitle" in line ][0][1].strip()
valLabel   = [f[x].split(':') for x, line in enumerate(f) if "valLabel" in line ][0][1].strip()
channels   = [f[x].split() for x, line in enumerate(f) if "Channel" in line ][0]
categories = [f[x].split() for x, line in enumerate(f) if "Category" in line ]
values     = [f[x+1].split() for x, line in enumerate(f) if "Category" in line ]
statUnc    = [f[x+2].split() for x, line in enumerate(f) if "Category" in line ]
totalUnc   = [f[x+3].split() for x, line in enumerate(f) if "Category" in line ]

import ROOT
totalUncGr=[ROOT.TGraphErrors() for i in range(1,len(channels))]
totalStatGr=[ROOT.TGraphErrors() for i in range(1,len(channels))]
ptLegends=[]
lLines=[]

for i in xrange(0,len(categories)):
    cat=categories[i][1].replace('eq','=')

    np=totalUncGr[0].GetN()
    jx=np*(len(categories)+1)
    pt = ROOT.TPaveText(0.75,jx+0.8,0.8,jx+1.2,"br")
    pt.SetBorderSize(0)
    pt.SetFillColor(0)
    pt.SetFillStyle(0)
    pt.SetTextFont(42);
    pt.SetTextSize(0.04);
    pt.AddText(cat)
    ptLegends.append(pt)

    l = ROOT.TLine(0.75,jx+0.8,1.15,jx+0.8)
    l.SetLineColor(17)
    l.SetLineStyle(9)
    lLines.append(l)

    for j in xrange(1,len(channels)) :
        np=totalUncGr[j-1].GetN()
        jx=np*(len(categories)+1)+j
        totalUncGr[j-1].SetPoint(np,float(values[i][j]),jx)
        totalUncGr[j-1].SetPointError(np,float(totalUnc[i][j]),0.08)
        totalStatGr[j-1].SetPoint(np,float(values[i][j]),jx)
        totalStatGr[j-1].SetPointError(np,float(statUnc[i][j]),0.)
        

ROOT.gSystem.Load('${CMSSW_BASE}/lib/${SCRAM_ARCH}/libCMGToolsHtoZZ2l2nu.so')
ROOT.setStyle()
c=ROOT.TCanvas("c1", "c1",600,600)
c.SetTickx(1);
c.SetTicky(1);
   
multigraphTotal=ROOT.TMultiGraph();
multigraphTotal.SetName("syst");
multigraphTotal.SetTitle("syst");
for item in totalUncGr :
    multigraphTotal.Add(item,'e2')
    item.SetFillColor(17)
    item.SetLineColor(2)
    item.SetLineWidth(2)
    item.SetMarkerColor(2)
    item.SetMarkerStyle(2)

multigraphStat=ROOT.TMultiGraph();
multigraphStat.SetName("stat");
multigraphStat.SetTitle("stat");
iitem=0
for item in totalStatGr :
    multigraphStat.Add(item,'p')
    item.SetFillStyle(0)
    if(iitem==0) : item.SetMarkerStyle(20); item.SetLineWidth(2)
    else :         item.SetMarkerStyle(24+iitem)
    iitem=iitem+1
    
multigraphTotal.Draw('ap')
multigraphStat.Draw('p')

multigraphTotal.GetXaxis().SetTitle(valLabel)
multigraphTotal.GetYaxis().SetNdivisions(0)

for l in lLines :
    l.SetX1(multigraphTotal.GetXaxis().GetXmin())
    l.SetX2(multigraphTotal.GetXaxis().GetXmax())
    l.Draw()
for p in ptLegends :
    delta=(multigraphTotal.GetXaxis().GetXmax()-multigraphTotal.GetXaxis().GetXmin())/10
    p.SetX1(multigraphTotal.GetXaxis().GetXmin()+delta/4)
    p.SetX2(multigraphTotal.GetXaxis().GetXmin()+5*delta/4)
    p.Draw()

leg=ROOT.TLegend(0.18,0.93,0.35,0.75)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetTextFont(42)
for iitem in xrange(0,len(totalStatGr)):
    label=channels[iitem+1].replace('mu','#mu')
    leg.AddEntry(totalStatGr[iitem],label,"lp")
leg.Draw()

pt=ROOT.TPaveText(0.1,0.97,0.8,0.98,"brNDC");
pt.SetBorderSize(0);
pt.SetFillColor(19);
pt.SetFillStyle(0);
#pt.SetTextAlign(11);
pt.SetTextFont(42);
pt.SetTextSize(0.04);
pt.AddText(plotTitle)
pt.Draw();

c.Modified()
c.Update()

raw_input('')
