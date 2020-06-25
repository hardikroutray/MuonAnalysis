import sys
from array import array                                                                                                             

import ROOT
import numpy
 
from ROOT import TTree, TH1F, TFile, TLorentzVector, TCanvas, TVector3, TH2F
                                                             
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

file_mu = ROOT.TFile('file:./Wmunu_output/ntuple_WMuNu_MC_1.root')
tree_mu = file_mu.Get('ttree')

outfile = ROOT.TFile('result_WMuNu_MC_1.root', 'recreate')
outfile.cd()

leading_muon_pt = ROOT.TH1F("leading_muon_pt", "leading_muon_pt", 50, 0, 50)
leading_muon_eta = ROOT.TH1F("leading_muon_eta", "leading_muon_eta", 200, -10, 10)
dimuon_mass = ROOT.TH1F("dimuon_mass", "dimuon_mass", 120, 0, 120)
munu_transversemass = ROOT.TH1F("munu_transversemass", "munu_transversemass", 120, 0, 120)

#Print Branches of the tree
tree_mu.Print()

#######Start the main code#########
numberofevents = 0

mu1p = TLorentzVector(0, 0, 0, 0)
mu2p = TLorentzVector(0, 0, 0, 0)

mup = TLorentzVector(0, 0, 0, 0)
nup = TLorentzVector(0, 0, 0, 0)
#Print and Fill histograms with the values from the flat tree
for ev, event in enumerate(tree_mu):
    
    numberofevents += 1
    
    print "In CMS event", event.event_num, "and CMS run", event.run, "and CMS lumi", event.lumi, "number of muons =", event.num_muons

    if event.num_muons != 0:

        leading_muon_pt.Fill(event.muon_pt[0])
        leading_muon_eta.Fill(event.muon_eta[0])

#####################dimuon mass calculation######################## 

    if event.num_muons > 2:
        
        if event.muon_q[0]*event.muon_q[1] < 0:

            mu1p.SetPtEtaPhiM(event.muon_pt[0], event.muon_eta[0], event.muon_phi[0], 0.1056583745)
            mu2p.SetPtEtaPhiM(event.muon_pt[1], event.muon_eta[1], event.muon_phi[1], 0.1056583745)

            print "dimuon mass", (mu1p+mu2p).M() 

            dimuon_mass.Fill((mu1p+mu2p).M())


#####################munu mass calculation######################## 

    if event.num_muons == 1:

        if event.muon_pt[0] > 25 and event.met_pt > 0:

            mup.SetPtEtaPhiM(event.muon_pt[0], event.muon_eta[0], event.muon_phi[0], 0.1056583745)
            nup.SetPtEtaPhiM(event.met_pt, event.met_eta, event.met_phi, 0)

            print "munu transverse mass", (mup+nup).Mt()

            munu_transversemass.Fill((mup+nup).Mt())


################################################################



print "Number of events processed", numberofevents

#outfile.cd()

leading_muon_pt.Write()
leading_muon_eta.Write()
dimuon_mass.Write()
munu_transversemass.Write()

outfile.Write()
outfile.Close()
                        
                        
