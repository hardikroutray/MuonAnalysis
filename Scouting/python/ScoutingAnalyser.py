# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inFile", help="text file containing list of ntuples")
parser.add_argument("outFile", help="name of ROOT file to write")
args = parser.parse_args()


from ROOT import TH1F, TH2F, TFile, TLorentzVector, TChain

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events


data_chain = TChain('scoutingntuplizer')
data_list = open(args.inFile, 'r')
data_files = data_list.readlines()
for i in data_files:
    data_chain.Add(i.strip())


#disable branches we don't want
data_chain.SetBranchStatus("*", 0)
data_chain.SetBranchStatus("muon_pt", 1)
data_chain.SetBranchStatus("muon_eta", 1)
data_chain.SetBranchStatus("muon_phi", 1)

#book histos

outfile = TFile(args.outFile, "recreate")

mu_pt = TH1F("pt", "pt", 100, 0, 20)
mu_eta = TH1F("eta", "eta", 100, -2, 2)
mu_phi = TH1F("phi", "phi", 90, -3.14, 3.14)
mumuM = TH1F("mumuM", "Mass", 1200, 0, 120)
mumuMlow = TH1F("mumuMlow", "Mass", 250, 0, 25)
mumuMverylow = TH1F("mumuMverylow", "Mass", 1000, 0, 5.0)
mvpt = TH2F("mvpt","mvpt",50,0,100,60,0,15)
mumuMDeltaCutZero = TH1F("mumuMDeltaCutZero", "Mass", 500, 0, 5.0)
mumuMDeltaCut4 = TH1F("mumuMDeltaCut4", "Mass", 1000, 0, 5.0)
mumuMDeltaCut6 = TH1F("mumuMDeltaCut6", "Mass", 1000, 0, 5.0)
mumuMDeltaCut8 = TH1F("mumuMDeltaCut8", "Mass", 1000, 0, 5.0)
mumuMDeltaCut10 = TH1F("mumuMDeltaCut10", "Mass", 1000, 0, 5.0)
mumuMDeltaCut12 = TH1F("mumuMDeltaCut12", "Mass", 1000, 0, 5.0)
mumuMDeltaCut14 = TH1F("mumuMDeltaCut14", "Mass", 1000, 0, 5.0)
mumuMDeltaCut16 = TH1F("mumuMDeltaCut16", "Mass", 1000, 0, 5.0)
mumuMDeltaCut18 = TH1F("mumuMDeltaCut18", "Mass", 1000, 0, 5.0)
mumuMDeltaCut20 = TH1F("mumuMDeltaCut20", "Mass", 1000, 0, 5.0)

mpipi = TH1F("mpipi", "Masspipi", 4000, 0,20)
mKK = TH1F("mKK", "MassKK", 4000, 0,20)

nevents = data_chain.GetEntries()
for i in range(nevents):       
    data_chain.GetEntry(i)

    p1 = TLorentzVector(0, 0, 0, 0)
    p2 = TLorentzVector(0, 0, 0, 0)


    p1.SetPtEtaPhiM(data_chain.muon_pt[0], data_chain.muon_eta[0], data_chain.muon_phi[0], 0.1056583745)
    p2.SetPtEtaPhiM(data_chain.muon_pt[1], data_chain.muon_eta[1], data_chain.muon_phi[1], 0.1056583745)


    p1pi = TLorentzVector(0, 0, 0, 0)
    p2pi = TLorentzVector(0, 0, 0, 0)

    p1pi.SetPtEtaPhiM(data_chain.muon_pt[0], data_chain.muon_eta[0], data_chain.muon_phi[0], 0.13957018)
    p2pi.SetPtEtaPhiM(data_chain.muon_pt[1], data_chain.muon_eta[1], data_chain.muon_phi[1], 0.13957018)

    
    p1K = TLorentzVector(0, 0, 0, 0)
    p2K = TLorentzVector(0, 0, 0, 0)
    
    p1K.SetPtEtaPhiM(data_chain.muon_pt[0], data_chain.muon_eta[0], data_chain.muon_phi[0], 0.493677)
    p2K.SetPtEtaPhiM(data_chain.muon_pt[1], data_chain.muon_eta[1], data_chain.muon_phi[1], 0.493677)


    mumuM.Fill((p1 + p2).M())
    mumuMlow.Fill((p1 + p2).M())
    mumuMverylow.Fill((p1 + p2).M())
    mvpt.Fill(abs(p1.Pt())+abs(p2.Pt()), (p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 0.0 :
        mumuMDeltaCutZero.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 4 :
        mumuMDeltaCut4.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 6 :
        mumuMDeltaCut6.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 8 :
        mumuMDeltaCut8.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 10 :
        mumuMDeltaCut10.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 12 :
        mumuMDeltaCut12.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 14 :
        mumuMDeltaCut14.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 16 :
        mumuMDeltaCut16.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 18 :
        mumuMDeltaCut18.Fill((p1 + p2).M())

    if  (p1 + p2).M() < abs(p1.Pt())+abs(p2.Pt()) - 20 :
        mumuMDeltaCut20.Fill((p1 + p2).M())


    mpipi.Fill((p1pi+p2pi).M())
    mKK.Fill((p1K+p2K).M())



    if( i % 200000 == 0):  
      print "At 2:", i, (p1 + p2).M() 
      print " -------------------------------"


outfile.cd()
mu_pt.Write()
mu_eta.Write()
mu_phi.Write()
mumuM.Write()
mumuMlow.Write()
mumuMverylow.Write()
mvpt.Write()
mumuMDeltaCutZero.Write()
mumuMDeltaCut4.Write()
mumuMDeltaCut6.Write()
mumuMDeltaCut8.Write()
mumuMDeltaCut10.Write()
mumuMDeltaCut12.Write()
mumuMDeltaCut14.Write()
mumuMDeltaCut16.Write()
mumuMDeltaCut18.Write()
mumuMDeltaCut20.Write()
mpipi.Write()
mKK.Write()
outfile.Close()
