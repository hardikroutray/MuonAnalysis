import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from array import array                                                                                                             

from DataFormats.FWLite import Handle, Events
import ROOT
import numpy

#from ROOT import TH1F, TFile, TLorentzVector, TCanvas                                                                               
from ROOT import TTree, TH1F, TFile, TLorentzVector, TCanvas, TVector3, TH2F
                                                             
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

muons, muonLabel = Handle("vector<pat::Muon>"), "slimmedMuons"
mets, metLabel = Handle("vector<pat::MET>"), "slimmedMETs"
jets, jetLabel = Handle("vector<pat::Jet>"), "slimmedJets"

files = ['DYJetsToLL.root']
#files = [#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/80000/009229C7-D9A9-E811-8084-001E67E6F7CE.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/80000/5A653255-2AA9-E811-AD8F-A4BF0112BE48.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/80000/BE60FBA0-74A9-E811-8CDA-001E67792422.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/80000/A83F9F22-9FA9-E811-B9F1-A4BF01158940.root'
#'root://cmsxrootd.fnal.gov//store/data/Run2017B/SingleMuon/MINIAOD/22Jun2017-v1/10000/DE62F0AB-A359-E711-84C1-3417EBE2EEC6.root',
#'root://cmsxrootd.fnal.gov//store/data/Run2017B/SingleMuon/MINIAOD/23Jun2017-v1/120000/AE7C4E13-A959-E711-9183-0025904C6624.root'
#]
events = Events(files)

# print "start running on file"

# lines=[line.rstrip('\n') for line in open('SingleMuon.txt', 'r')]

# print "started running"

# if (int(sys.argv[1])>=len(lines)):
#       sys.exit('File number out of range')

# print "Running on :", lines[int(sys.argv[1])]

# events = Events(lines[int(sys.argv[1])])


# outfile = TFile('file:./Wmunu_output/ntuple_SingleMuon_data_wjets_{0}.root'.format(sys.argv[1]), 'recreate')


#outfile = TFile('ntuple_data_2.root', 'recreate')
outfile.cd()

ttree = TTree('ttree', 'ttree with WMuNu MC')

maxn = 100

run = array("i", [0])
lumi = array("i", [0])
event_num = array("f", [0.])
met_pt = array("f", [0.0])
met_eta = array("f", [0.0])
met_phi = array("f", [0.0])
num_muons = array("i", [0])
muon_pt = array("f", maxn*[0.0])
muon_eta = array("f", maxn*[0.0])
muon_phi = array("f", maxn*[0.0])
muon_q = array("f", maxn*[0.0])
num_jets = array("i", [0])
jet_pt = array("f", maxn*[0.0])
jet_eta = array("f", maxn*[0.0])
jet_phi = array("f", maxn*[0.0])
jet_mass = array("f", maxn*[0.0])
jet_csv = array("f", maxn*[0.0])
#munu_transversemass = array("f", [0.0])
#muon1_pt = array("f", [0.0])
#muon1_eta = array("f", [0.0])
#muon1_phi = array("f", [0.0]) 
#muon1_q = array("f", [0.0])
#muon2_pt = array("f", [0.0])
#muon2_eta = array("f", [0.0])
#muon2_phi = array("f", [0.0])
#muon2_q = array("f", [0.0])
#dimuon_mass = array("f", [0.0])


ttree.Branch("run", run, "run/I")
ttree.Branch("lumi", lumi, "lumi/I")
ttree.Branch("event_num", event_num, "event_num/I")
ttree.Branch("met_pt", met_pt, "met_pt/F")
ttree.Branch("met_eta", met_eta, "met_eta/F")
ttree.Branch("met_phi", met_pt, "met_phi/F")
ttree.Branch("num_muons", num_muons, "num_muons/I")
ttree.Branch("muon_pt", muon_pt, "muon_pt[num_muons]/F")
ttree.Branch("muon_eta", muon_eta, "muon_eta[num_muons]/F")
ttree.Branch("muon_phi", muon_phi, "muon_phi[num_muons]/F")
ttree.Branch("muon_q", muon_q, "muon_q[num_muons]/F")
ttree.Branch("num_jets", num_jets, "num_jets/I")
ttree.Branch("jet_pt", jet_pt, "jet_pt[num_jets]/F")
ttree.Branch("jet_eta", jet_eta, "jet_eta[num_jets]/F")
ttree.Branch("jet_phi", jet_phi, "jet_phi[num_jets]/F")
ttree.Branch("jet_mass", jet_mass, "jet_mass[num_jets]/F")
ttree.Branch("jet_csv", jet_csv, "jet_csv[num_jets]/F")
#ttree.Branch("munu_transversemass", munu_transversemass, "munu_transversemass/F")

#ttree.Branch("muon1_pt", muon1_pt, "muon1_pt/F")
#ttree.Branch("muon1_eta", muon1_eta, "muon1_eta/F")
#ttree.Branch("muon1_phi", muon1_phi, "muon1_phi/F")
#ttree.Branch("muon1_q", muon1_q, "muon1_q/F")
#ttree.Branch("muon2_pt", muon2_pt, "muon2_pt/F")
#ttree.Branch("muon2_eta", muon2_eta, "muon2_eta/F")
#ttree.Branch("muon2_phi", muon2_phi, "muon2_phi/F")
#ttree.Branch("muon2_q", muon2_q, "muon2_q/F")
#ttree.Branch("dimuon_mass", dimuon_mass, "dimuon_mass/F")

#leading_muon_pt = ROOT.TH1F("leading_muon_pt", "leading_muon_pt", 100, -10, 10)
#leading_muon_eta = ROOT.TH1F("leading_muon_eta", "leading_muon_eta", 100, -10, 10)
#dimuon_mass_hist = ROOT.TH1F("dimuon_mass_hist", "dimuon_mass_hist", 1200, 0, 120)

#munu_mass_hist = ROOT.TH1F("munu_mass_hist", "munu_mass_hist", 120, 0, 120)


i = 0
numberofmuons = 0
numberofmets = 0
numberofjets = 0
#txtfile = open("runlumieventmask_data_4.txt", "w")
#txtfile = open("eventinfo_datafile4.txt", "w")

for event in events:
      i+= 1
      
#      if i > 100:
#            break

#    print "Event", i                                                                                                               
  
      mu1p = TLorentzVector(0, 0, 0, 0)
      mu2p = TLorentzVector(0, 0, 0, 0)

#      mup = TLorentzVector(0, 0, 0, 0)
#      nup = TLorentzVector(0, 0, 0, 0)

#      print "run number", event.eventAuxiliary().run()
#      print "event number", event.eventAuxiliary().event()
#      print "lumisection", event.eventAuxiliary().luminosityBlock()

      run[0] = event.eventAuxiliary().run()
      event_num[0] = event.eventAuxiliary().event()
      lumi[0] = event.eventAuxiliary().luminosityBlock()

      event.getByLabel(muonLabel, muons)
      event.getByLabel(metLabel, mets)
      event.getByLabel(jetLabel, jets)

      numberofmuons+= muons.product().size()
      numberofmets+= mets.product().size()
      numberofjets+= jets.product().size()

      num_muons[0] = muons.product().size() 
      num_jets[0] = jets.product().size()


#      print "met number", mets.product().size()

#      if muons.product().size() == 1 and mets.product().size() == 1:

#            for n1, neutrino1 in enumerate(mets.product()):

#                  if neutrino1.pt() > 40:

#                        nup.SetPtEtaPhiM(neutrino1.pt(), neutrino1.eta(), neutrino1.phi(), 0)

#                  print "met pt", neutrino1.pt(), "met_eta", neutrino1.eta(), "met_phi", neutrino1.phi()

#            for m1, muon1 in enumerate(muons.product()):

#                  if muon1.isGlobalMuon():


#                        muon_pt[0] = muon1.pt()                                                                                    
#                        muon_eta[0] = muon1.eta()                                                                                  
#                        muon_phi[0] = muon1.phi()                                                                                  
#                        muon_q[0] = muon1.charge() 

#                        mup.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), 0.1056583745)
                  
#                        munu_transversemass[0] = (mup+nup).Mt()

#                        if muon1.pt() > 40 and numpy.abs(muon1.eta()) < 2.1 and ((muon1.trackIso() + muon1.caloIso())/muon1.pt()) < 0.15: 

#                              munu_mass_hist.Fill((mup+nup).Mt())


#                        print "Event number", event.eventAuxiliary().event()
#                        print "Lumi number", event.eventAuxiliary().luminosityBlock()
#                        print "run number", event.eventAuxiliary().run()



#                        txtfile.write(str(int(event.eventAuxiliary().run())) + ":" + str(int(event.eventAuxiliary().luminosityBlock())) + ":" + str(int(event.eventAuxiliary().event())) + ",")

#      if muons.product().size() > 4:                                                                                               
#             continue                

      if mets.product().size() != 1:

            print "the no of mets in the event", mets.product().size()
#      print "the no of muons in the event", muons.product().size()

      for met1, met in enumerate(mets.product()): 
            
            met_pt[0] = met.pt()
            met_eta[0] = met.eta()
            met_phi[0] = met.phi()



      for m, muon in enumerate(muons.product()):

#            print m, "muon_pt", muon.pt()
            
            muon_pt[m] = muon.pt()                                                                                     
            muon_eta[m] = muon.eta()                                                                                   
            muon_phi[m] = muon.phi()                                                                                   
            muon_q[m] = muon.charge()



      for j, jet in enumerate(jets.product()):

#            print j, "jet_pt", jet.pt(), "jet_eta", jet.eta(), "jet_phi", jet.phi(), "jet_mass", jet.mass(),
#            print j, "jetcsv", jet.bDiscriminator("pfJetBProbabilityBJetTags")

            jet_pt[j] = jet.pt()
            jet_eta[j] = jet.eta()
            jet_phi[j] = jet.phi()
            jet_mass[j] = jet.mass()
            jet_csv[j] = jet.bDiscriminator("pfJetBProbabilityBJetTags") 

      '''
      
      if muons.product().size() != 2:
             continue

      for m1, muon1 in enumerate(muons.product()):

            for m2, muon2 in enumerate(muons.product()):

                  if m2 <= m1:
                        continue

                  if muon1.isGlobalMuon() and muon2.isGlobalMuon():

#                  print "In Event", i , "muon1 pt is", muon1.pt()

                        mu1p.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), 0.1056583745)
                        mu2p.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), 0.1056583745)
                        motherp = mu1p + mu2p


                        leading_muon_pt.Fill(muon1.pt())
                        leading_muon_eta.Fill(muon1.eta())
                        dimuon_mass_hist.Fill(motherp.M())

                        muon1_pt[0] = muon1.pt()
                        muon1_eta[0] = muon1.eta()
                        muon1_phi[0] = muon1.phi()
                        muon1_q[0] = muon1.charge()
                        muon2_pt[0] = muon2.pt()
                        muon2_eta[0] = muon2.eta()
                        muon2_phi[0] = muon2.phi()
                        muon2_q[0] = muon2.charge()
                        dimuon_mass[0] = motherp.M()

#                        if numpy.abs(muon1.eta()) < 0.1 and numpy.abs(muon2.eta()) < 0.1:

#                              print "Event number", event.eventAuxiliary().event()
#                              print "Lumi number", event.eventAuxiliary().luminosityBlock()
#                              print "run number", event.eventAuxiliary().run()

#                              txtfile.write(str(int(event.eventAuxiliary().run())) + ":" + str(int(event.eventAuxiliary().luminosityBlock())) + ":" + str(int(event.eventAuxiliary().event())) + ",")

#                              txtfile.write(str(int(event.eventAuxiliary().run())) + " " + str(int(event.eventAuxiliary().luminosityBlock())) + " " + str(int(event.eventAuxiliary().event())) + " " + str(float(muon1.pt())) + " " + str(float(muon1.eta())) + " " + str(float(muon1.phi())) + " " + str(float(muon1.charge())) +  " " + str(float(muon2.pt())) + " " +  str(float(muon2.eta())) + " " + str(float(muon2.phi())) + " " + str(float(muon2.charge())) + " " + str(float(motherp.M())) + "\n")

      '''
                        
      ttree.Fill()

#      print "next event"


print "The number of muons processed in the file", numberofmuons


#leading_muon_pt.Write()
#leading_muon_eta.Write()
#dimuon_mass_hist.Write()
#munu_mass_hist.Write()
outfile.Write()
outfile.Close()

                        
                        
