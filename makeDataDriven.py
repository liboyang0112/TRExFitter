import ROOT

def makeHistograms(yCR1, yCR2, ySR1, ySR2, statUncertainty=0.02):
    statU = statUncertainty
    # control region
    histCR = ROOT.TH1F("xCR", "xCR", 2, 0, 2)
    histCR.SetBinContent(1, yCR1)
    histCR.SetBinError(1, statU*yCR1)
    histCR.SetBinContent(2, yCR2)
    histCR.SetBinError(2, statU*yCR2)    
    # signal region
    histSR = ROOT.TH1F("xSR", "xSR", 2, 0, 2)
    histSR.SetBinContent(1, ySR1)
    histSR.SetBinError(1, statU*ySR1)
    histSR.SetBinContent(2, ySR2)
    histSR.SetBinError(2, statU*ySR2)
    return histCR, histSR

def makeSample(name, yCR1, yCR2, ySR1, ySR2, exDir='.'):
    fileTmp = ROOT.TFile(exDir + '/' + name + '.root', "RECREATE")
    histCR, histSR = makeHistograms(yCR1, yCR2, ySR1, ySR2)
    fileTmp.Write()
    fileTmp.Close()
    
def makeDataDrivenWithSig(exampleDir='exampleDataDriven'):
    # two bin example
    # data in signal region SR and control region CR
    # background bkg1 is MC based
    # background bkg2 is data-driven
    # one root file per sample
    ROOT.gSystem.mkdir(exampleDir)
    # data
    makeSample('hist_data', 150, 70, 90, 110, exDir=exampleDir) 
    # signal given by MC
    makeSample('hist_signal', 5, 40, 10, 80, exDir=exampleDir)     
    # bkg1 MC
    makeSample('hist_bkg1', 10, 10, 20, 20, exDir=exampleDir)
    # bkg2 data-driven
    # assume that extrapolation
    makeSample('hist_bkg2', 120, 10, 60, 5, exDir=exampleDir)

def makeDataDrivenNoSig(exampleDir='exampleDataDriven'):
    # two bin example
    # data in signal region SR and control region CR
    # background bkg1 is MC based
    # background bkg2 is data-driven
    # one root file per sample
    ROOT.gSystem.mkdir(exampleDir)
    # data
    makeSample('hist_data', 150, 70, 90, 110, exDir=exampleDir) 
    # signal given by MC
    makeSample('hist_signal', 5, 40, 10, 80, exDir=exampleDir)     
    # bkg1 MC
    makeSample('hist_bkg1', 10, 10, 20, 70, exDir=exampleDir)
    # bkg2 data-driven
    # assume that extrapolation
    makeSample('hist_bkg2', 138, 55, 77, 30, exDir=exampleDir)
    


makeDataDrivenNoSig()
