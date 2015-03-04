import glob, os, sys
import commands
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *
from glob import *
import time
import argparse

ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasUtils.C")
SetAtlasStyle()
gStyle.SetTextFont(43)
gStyle.SetTextSize(18)
gStyle.SetLabelFont(43)
gStyle.SetLabelSize(18)

files=[]

parser = argparse.ArgumentParser(description='Create NP plot.')
parser.add_argument('fileList', metavar='S', type=string, nargs='+',
                   help='list of input txt files')
parser.add_argument('--outFile', metavar='O', type=string, 
                   help='output file')
args = parser.parse_args()

out = ""

files = args.fileList
out   = args.outFile

print "Out file = ",out
print "In files = ",files

for file in files:
  print file

# "xcheckResults/TextFileFitResult/GlobalFit_fitres_unconditionnal_mu0.txt"]
# files=["dummy.txt"]


# --- OPTIONS --- #

xmin = -2.9
xmax = 2.9
max = 0

npToExclude = ["SigXsecOverSM","gamma_"]
brazilian = True
grayLines = False

# --- --- --- --- #


g = []

for file in files:
  g.append(TGraphErrors())

idx=0

Names=[]
Cval=[]
Err=[]

b = []

i = 0
for file in files:
  print "Opening file",file
  f     = open(str(file))
  print "Opened"
  lines = f.readlines()
  for line in lines:
    line2 = line.split("\n")[0]

    if "NUISANCE" in line2:
      continue
    if line2=="":
      continue
    if "Luminosity" in line2:
      continue
    if "CORRELATION" in line2:
      break

    ii = 0
    name = line2.split(" ")[ii]
    ii += 1
    val  = line2.split(" ")[ii]
    ii += 1
    while val == "":
      val = line2.split(" ")[ii]
      ii += 1
    err  = line2.split(" ")[ii]
    ii += 1
    while err == "":
      err = line2.split(" ")[ii]
      ii += 1

    #print name,val,err

    if any(np in line for np in npToExclude):
      continue
    
    g[i].SetPoint(idx,float(val),float(idx+0.5) )
    g[i].SetPointError(idx,float(err),float(0) )

    Names.append(name)

    if grayLines:
      b.append(TBox(xmin,float(idx),xmax,float(idx+1)))
      if (idx % 2) == 0:
        b[idx].SetFillColor(kGray)
      else:
        b[idx].SetFillColor(kWhite)
    
    idx+=1

  if idx > max:
    max = idx


# start drawing...

#defHeight = 1200
#defNpar = 48
lineHeight = 20
offsetUp = 10
offsetDown = 40
offset = offsetUp+offsetDown
newHeight = offset+int(max*lineHeight)
c = TCanvas("c","c",600,newHeight)
gPad.SetLeftMargin(0.05)
gPad.SetRightMargin(0.33)
# absolute top margin should be 0.02*defHeight -> rel = 0.02
# so relative new will be 0.02*defHeight/newHeight ( so that (...)*newHeight = 0.02*defHeight)
#gPad.SetTopMargin(0.02*defHeight/newHeight)
#gPad.SetBottomMargin(0.07*defHeight/newHeight)
gPad.SetTopMargin(1.*offsetUp/newHeight)
gPad.SetBottomMargin(1.*offsetDown/newHeight)
    
h_dummy = TH1F("h_dummy","h_dummy",10,xmin,xmax)
h_dummy.SetMaximum(max)
h_dummy.SetLineWidth(0)
h_dummy.SetFillStyle(0)
h_dummy.SetLineColor(kWhite)
h_dummy.SetFillColor(kWhite)
h_dummy.SetMinimum(0.)
h_dummy.GetYaxis().SetLabelSize(0)
h_dummy.Draw()

for idx in range(0,max):
  if grayLines:
    b[idx].Draw("SAME")
  systs=TLatex()
  #systs.SetTextSize(0.02/(max/48.))
  systs.DrawLatex(3.,idx+0.3,Names[idx])

#h_dummy.GetXaxis().SetLabelSize(h_dummy.GetXaxis().GetLabelSize()*0.5/(max/48.))
c.SetTicks(1,0)
h_dummy.GetYaxis().SetNdivisions(0)

if grayLines:
  l0 = TLine(0,0,0,max)
  l0.SetLineStyle(7)
  l0.SetLineColor(kBlue)
  l0.Draw("SAME")
  l1 = TLine(-1,0,-1,max)
  l1.SetLineStyle(2)
  l1.SetLineColor(kRed)
  l1.Draw("SAME")
  l2 = TLine(1,0,1,max)
  l2.SetLineStyle(2)
  l2.SetLineColor(kRed)
  l2.Draw("SAME")

if brazilian:
  l0 = TLine(0,0,0,max)
  l0.SetLineStyle(7)
  l0.SetLineColor(kBlack)
  b1 = TBox(-1,0,1,max)
  b2 = TBox(-2,0,2,max)
  b1.SetFillColor(kGreen)
  b2.SetFillColor(kYellow)
  b2.Draw("SAME")
  b1.Draw("SAME")
  l0.Draw("SAME")

for i in range(0,len(files)):
  g[i].Draw("psame")

gPad.RedrawAxis();

c.Update()
##time.sleep(20)

if out == "":
  c.SaveAs("NuisPar_compare.png")
else:
  c.SaveAs(str(out))
