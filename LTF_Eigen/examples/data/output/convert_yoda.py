#! /usr/bin/env python


############################################################################################################
#                                                                                                          #
#  python script to convert yoda file to root                                                              #
#  Takes the (compressed) yoda file as input                                                               #
#  Script modified from                                                                                    #
#  https://gitlab.cern.ch/atlas/athena/-/blob/main/Generators/Rivet_i/examples/convert2root                #
#  To run the conversion do                                                                                #
#  $ python convert_yoda.py WbWb_singlelepton.yoda.gz                                                      #
#                                                                                                          #
############################################################################################################

from array import array
import ROOT as rt
import yoda, sys

fName = str(sys.argv[1])
yodaAOs = yoda.read(fName) # creates dictionary holding all the hists

rtFile = rt.TFile(fName[:fName.find('.yoda.gz')] + '.root', 'recreate')

for name in yodaAOs:
  print("Now processing "+name)
  yodaAO = yodaAOs[name];  rtAO = None
  if 'Histo1D' in str(yodaAO):
    rtAO = rt.TH1D(name, '', yodaAO.numBins(), array('d', yodaAO.xEdges()))
    rtAO.Sumw2(); rtErrs = rtAO.GetSumw2()
    for i in range(rtAO.GetNbinsX()):
      rtAO.SetBinContent(i + 1, yodaAO.bin(i+1).sumW())
      rtErrs.AddAt(yodaAO.bin(i+1).sumW2(), i+1)
  elif 'Histo2D' in str(yodaAO):
    if type(yodaAO.xEdges()[0]) is str:
      continue
    rtAO = rt.TH2D(name, '', yodaAO.numBinsX(), array('d', yodaAO.xEdges()), yodaAO.numBinsY(), array('d', yodaAO.yEdges()))
    for j in range(len(yodaAO.yEdges()) + 1):
      for i in range(len(yodaAO.xEdges()) + 1 ):
        bin_weight = yodaAO.bin(i,j).sumW()
        bin_weight2 = yodaAO.bin(i,j).sumW2()
        rtAO.SetBinContent(i,j, bin_weight)
        rtAO.SetBinError(i,j, m.sqrt(bin_weight2))
### For objects of the estimate class you need to change the name (drop "/RAW")
#  elif 'Estimate1D' in str(yodaAO):
#    if type(yodaAO.xEdges()[0]) is str:
#      continue
#    rtAO = rt.TH1D(name, '', yodaAO.numBins(), array('d', yodaAO.xEdges()))
#    values = yodaAO.vals()
#    err = yodaAO.mkScatter().errs(1)
#    for i in range(rtAO.GetNbinsX()):
#      rtAO.SetBinContent(i + 1, values[i])
#      rtAO.SetBinError(i + 1, err[i][0])
#
#  elif 'Estimate2D' in str(yodaAO):
#    if type(yodaAO.xEdges()[0]) is str:
#      continue
#    rtAO = rt.TH2D(name, '', yodaAO.numBinsX(), array('d', yodaAO.xEdges()), yodaAO.numBinsY(), array('d', yodaAO.yEdges()))
#    content = yodaAO.vals()
#    #err = yodaAO.mkScatter().errs(2)                                                                                                                                                                       
#    for i in range(len(yodaAO.xEdges()) - 1):
#      for j in range(len(yodaAO.yEdges()) - 1):
#        rtAO.SetBinContent(i+1,j+1, content[i*(len(yodaAO.yEdges())-1)+j])

  ### This is not tested yet
  elif 'Scatter2D' in str(yodaAO):
    rtAO = rt.TGraphAsymmErrors(yodaAO.numPoints())
    for i in range(yodaAO.numPoints()):
      x = yodaAO.point(i).x(); y = yodaAO.point(i).y()
      xLo, xHi = yodaAO.point(i).xErrs()
      yLo, yHi = yodaAO.point(i).yErrs()
      rtAO.SetPoint(i, x, y)
      rtAO.SetPointError(i, xLo, xHi, yLo, yHi)
  else:
    continue
  rtAO.Write(str(name.replace("/WbWb_singlelepton/","")))
  print("Write out histogram " + name.replace("/WbWb_singlelepton/",""))
rtFile.Close()
