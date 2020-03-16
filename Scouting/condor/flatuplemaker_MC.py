# import ROOT in batch mode                                                                                                          
import os
import sys

#import sys                                                                                                                          
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ROOT.gROOT.LoadMacro("propagation_utils.cc")

import numpy as np
from array import array

from ROOT import TH1F, TH1D, TH2D, TFile, TLorentzVector, TVector3, TChain, TProfile, TTree

# load FWLite C++ libraries                                                                                                          
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries                                                                                                       
from DataFormats.FWLite import Handle, Events


#file_mu = ROOT.TFile("./scouting_ntuple_test.root")                                                                                
file_mu = ROOT.TFile("./Phi_mass20_ct2_trkparam.root") 
tree_mu = file_mu.Get('scoutingntuplizer')

#file_mu = ROOT.TFile('/cms/routray/crab_output/muontuples_10percent/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017C_v1/200131_190828/0000/scouting_ntuple_{0}.root'.format(sys.argv[1]))
#tree_mu = file_mu.Get('scoutingntuplizer')

outfile = ROOT.TFile('flat_dimuon_tree_Phi_mass20_ct2_baseline_noiso_trkparam.root', 'recreate')
outfile.cd()

#outfile = TFile("flatdimuontrees_baselineselections/C/flatdimuontree_{0}.root".format(sys.argv[1]), "recreate")
#outfile.cd()

tree = ROOT.TTree('events', 'flat tree with dimuon info')

run = array("f", [0.0])
lumi = array("f", [0.0])
event_num = array("d", [0.0])
dimuon_pt = array("f", [0.0])
dimuon_mass = array("f", [0.0])
dimuon_mass_uncorr = array("f", [0.0])
lxy = array("f", [0.0])
dphidimudv = array("f", [0.0])
detadimudv = array("f", [0.0])
ctau = array("f", [0.0])
muon1_dxy = array("f", [0.0])
muon2_dxy = array("f", [0.0])
muon1_edxy = array("f", [0.0])
muon2_edxy = array("f", [0.0])
muon1_dxycorr = array("f", [0.0])
muon2_dxycorr = array("f", [0.0])
muon1_trkiso = array("f", [0.0])
muon2_trkiso = array("f", [0.0])
muon1_chi2overndof = array("f", [0.0])
muon2_chi2overndof = array("f", [0.0])
muon1_pt = array("f", [0.0])
muon2_pt = array("f", [0.0])
muon1_eta = array("f", [0.0])
muon2_eta = array("f", [0.0])
muon1_phi = array("f", [0.0])
muon2_phi = array("f", [0.0])
muon1_phicorr = array("f", [0.0])
muon2_phicorr = array("f", [0.0])
DV_chi2overndof = array("f", [0.0])
PVx = array("f", [0.0])
PVy = array("f", [0.0])
PVz = array("f", [0.0])
BSx = array("f", [0.0])
BSy = array("f", [0.0])
BSz = array("f", [0.0])
DVx = array("f", [0.0])
DVy = array("f", [0.0])
DVz = array("f", [0.0])
DVxerr = array("f", [0.0])
DVyerr = array("f", [0.0])
DVzerr = array("f", [0.0])
detamumu = array("f", [0.0])
dphimumu = array("f", [0.0])
dRmumu = array("f", [0.0])
dRmuon1jet = array("f", [0.0]) 
dRmuon2jet = array("f", [0.0])


#tree.Branch('muon1_pt', muon1_pt, 'muon1_pt/F')
tree.Branch('run', run, 'run/F')
tree.Branch('lumi', lumi, 'lumi/F')
tree.Branch('event_num', event_num, 'event_num/D')
tree.Branch('dimuon_pt', dimuon_pt, 'dimuon_pt/F')
tree.Branch('dimuon_mass', dimuon_mass, 'dimuon_mass/F')
tree.Branch('dimuon_mass_uncorr', dimuon_mass_uncorr, 'dimuon_mass_uncorr/F')
#tree.Branch('Lxy', Lxy, 'Lxy/F')
tree.Branch('lxy', lxy, 'lxy/F')
#tree.Branch('rho', rho, 'rho/F')
#tree.Branch('theta', theta, 'theta/F')
tree.Branch('dphidimudv', dphidimudv, 'dphidimudv/F')
tree.Branch('detadimudv', detadimudv, 'detadimudv/F')
tree.Branch('ctau', ctau, 'ctau/F')
tree.Branch('muon1_dxy', muon1_dxy, 'muon1_dxy/F')
tree.Branch('muon2_dxy', muon2_dxy, 'muon2_dxy/F')
tree.Branch('muon1_edxy', muon1_edxy, 'muon1_edxy/F')
tree.Branch('muon2_edxy', muon2_edxy, 'muon2_edxy/F')
tree.Branch('muon1_dxycorr', muon1_dxycorr, 'muon1_dxycorr/F')
tree.Branch('muon2_dxycorr', muon2_dxycorr, 'muon2_dxycorr/F')
tree.Branch('muon1_trkiso', muon1_trkiso, 'muon1_trkiso/F')
tree.Branch('muon2_trkiso', muon2_trkiso, 'muon2_trkiso/F')
tree.Branch('muon1_chi2overndof', muon1_chi2overndof, 'muon1_chi2overndof/F')
tree.Branch('muon2_chi2overndof', muon2_chi2overndof, 'muon2_chi2overndof/F')
tree.Branch('muon1_pt', muon1_pt, 'muon1_pt/F')
tree.Branch('muon2_pt', muon2_pt, 'muon2_pt/F')
tree.Branch('muon1_eta', muon1_eta, 'muon1_eta/F')
tree.Branch('muon2_eta', muon2_eta, 'muon2_eta/F')
tree.Branch('muon1_phi', muon1_phi, 'muon1_phi/F')
tree.Branch('muon2_phi', muon2_phi, 'muon2_phi/F')
tree.Branch('muon1_phicorr', muon1_phicorr, 'muon1_phicorr/F')
tree.Branch('muon2_phicorr', muon2_phicorr, 'muon2_phicorr/F')
tree.Branch('DV_chi2overndof', DV_chi2overndof, 'DV_chi2overndof/F')
tree.Branch('PVx', PVx, 'PVx/F')
tree.Branch('PVy', PVy, 'PVy/F')
tree.Branch('PVz', PVz, 'PVz/F')
tree.Branch('BSx', BSx, 'BSx/F')
tree.Branch('BSy', BSy, 'BSy/F')
tree.Branch('BSz', BSz, 'BSz/F')
tree.Branch('DVx', DVx, 'DVx/F')
tree.Branch('DVy', DVy, 'DVy/F')
tree.Branch('DVz', DVz, 'DVz/F')
tree.Branch('DVxerr', DVxerr, 'DVxerr/F')
tree.Branch('DVyerr', DVyerr, 'DVyerr/F')
tree.Branch('DVzerr', DVzerr, 'DVzerr/F')
tree.Branch('detamumu', detamumu, 'detamumu/F')
tree.Branch('dphimumu', dphimumu, 'dphimumu/F')
tree.Branch('dRmumu', dRmumu, 'dRmumu/F')
tree.Branch('dRmuon1jet', dRmuon1jet, 'dRmuon1jet/F')
tree.Branch('dRmuon2jet', dRmuon2jet, 'dRmuon2jet/F')


#DeltaRmumu = ROOT.TH1F("DeltaRmumu", "DeltaRmumu", 100, -10, 10)
#muchi2overndof = ROOT.TH1F("muchi2overndof", "muchi2overndof", 100, -10, 10)
#ntrackerlayerswithmeasurement = ROOT.TH1F("ntrackerlayerswithmeasurement", "ntrackerlayerswithmeasurement", 100, -10, 10)
#nvalidmuonhits = ROOT.TH1F("nvalidmuonhits", "nvalidmuonhits", 100, -10, 10)
#nvalidpixelhits = ROOT.TH1F("nvalidpixelhits", "nvalidpixelhits", 100, -10, 10)
#nvalidstriphits = ROOT.TH1F("nvalidstriphits", "nvalidstriphits", 100, -10, 10)

numberofevents = 0

for i, event in enumerate(tree_mu):


    numberofevents += 1
    
    mu1p = TLorentzVector(0, 0, 0, 0)
    mu2p = TLorentzVector(0, 0, 0, 0)
    mu1p_uncorr = TLorentzVector(0, 0, 0, 0)
    mu2p_uncorr = TLorentzVector(0, 0, 0, 0)
    motherp = TLorentzVector(0, 0, 0, 0)
    PVvector2D = TVector3(0.,0.,0.)
    BSvector2D = TVector3(0.,0.,0.)
    DVvector2D  = TVector3(0.,0.,0.)
    PVvector3D = TVector3(0.,0.,0.)
    DVvector3D  = TVector3(0.,0.,0.)
    DVPVvector2D = TVector3(0.,0.,0.)
    DVBSvector2D = TVector3(0.,0.,0.)
    DVPVvector3D = TVector3(0.,0.,0.)
    motherv2D = TVector3(0.,0.,0.)
    motherv3D = TVector3(0.,0.,0.)
    jetp = TLorentzVector(0, 0, 0, 0)

    
    passL1_DoubleMu4_SQ_OS_dR_Max1p2 = 0
    passL1_DoubleMu_15_7 = 0
    passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 0  
    passL1 = 0
    passDST_DoubleMu3_noVtx_CaloScouting_v = 0
    passHLT = 0

    print "-----In Run" ,tree_mu.Run , "---Lumi", tree_mu.Lumi, "----Event", tree_mu.Event, "-----------"

    for j, _ in enumerate(event.l1bitmap):

        if event.l1bitmap[j].first == "L1_DoubleMu4_SQ_OS_dR_Max1p2" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu4_SQ_OS_dR_Max1p2 = 1

        if event.l1bitmap[j].first == "L1_DoubleMu_15_7" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu_15_7 = 1

        if event.l1bitmap[j].first == "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 1

    passL1 = passL1_DoubleMu4_SQ_OS_dR_Max1p2 + passL1_DoubleMu_15_7 + passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4

    # for j, _ in enumerate(event.hltbitmap):

    #     if event.hltbitmap[j].first == "DST_DoubleMu3_noVtx_CaloScouting_v" and event.hltbitmap[j].second == 1:
    #         passDST_DoubleMu3_noVtx_CaloScouting_v = 1

    # passHLT = passDST_DoubleMu3_noVtx_CaloScouting_v

#    print "In event", i , "passL1?", passL1, "passHLT?", passHLT

        
    if passL1 == 0:
        continue

    # if passHLT == 0:
    #     continue

    if len(event.muon_pt) < 2:
        continue

    if event.dispvertex_num < 1:
        continue


    print "Number of muons", len(event.muon_pt), "Number of vertices", event.dispvertex_num


    mindverror = 9999
    bestdv = -1

    
    for dv in range(event.dispvertex_num):

        if event.dispvertex_chi2[dv]/event.dispvertex_ndof[dv] >= 5:
            continue
        if event.dispvertex_ex[dv] >= 0.05:
            continue
        if event.dispvertex_ey[dv] >= 0.05:
            continue
        if event.dispvertex_ez[dv] >= 0.1:
            continue
        if np.sqrt((event.dispvertex_x[dv] - event.privertex_x[0])**2 + (event.dispvertex_y[dv] - event.privertex_y[0])**2) >= 11:
            continue


        if np.maximum(event.dispvertex_ex[dv],event.dispvertex_ey[dv]) < mindverror:
            mindverror = np.maximum(event.dispvertex_ex[dv],event.dispvertex_ey[dv])            
            bestdv = dv

    
    print "best dv index(if -1 no best DV)", bestdv

    if bestdv == -1:
        continue

    numberofmuons_invtx = 0
    highestpt = -1
    lowestpt = 9999
    sublead = -1
    lead = -1


    for index in range(len(event.vertex_index)):
        for subindex in event.vertex_index[index]:
                
            print subindex, "muon",index
                
            if subindex == bestdv:

                numberofmuons_invtx += 1
                if numberofmuons_invtx > 2:
                    continue
                    
                        
                if event.muon_pt[index] > highestpt:
                    lead = index

                if event.muon_pt[index] < lowestpt:
                    sublead = index

                highestpt = event.muon_pt[index]
                lowestpt = event.muon_pt[index]
        

    print "Number of muons in vtx", numberofmuons_invtx
    print "leading muon index",lead
    print "subleading muon index", sublead
        

    if numberofmuons_invtx != 2:
        continue
 
        
    # a = np.array(event.muon_pt)
    # a = np.argsort(a)[::-1]

    print "DV chi2overndof", event.dispvertex_chi2[bestdv]/event.dispvertex_ndof[bestdv], "DV_xerr", event.dispvertex_ex[bestdv], "DVyerr", event.dispvertex_ey[bestdv], "DVzerr", event.dispvertex_ez[bestdv]


    # if event.dispvertex_chi2[bestdv]/event.dispvertex_ndof[bestdv] >= 5:
    #     continue
    # if event.dispvertex_ex[bestdv] >= 0.05:
    #     continue
    # if event.dispvertex_ey[bestdv] >= 0.05:
    #     continue
    # if event.dispvertex_ez[bestdv] >= 0.1:
    #     continue


    print "muon1 chi2overndof", event.muon_chi2[lead]/event.muon_ndof[lead], "muon2 chi2overndof", event.muon_chi2[sublead]/event.muon_ndof[sublead]


    if event.muon_chi2[lead]/event.muon_ndof[lead] >= 3:
        continue
    if event.muon_chi2[sublead]/event.muon_ndof[sublead] >= 3:
        continue


    print "muon1 pt", event.muon_pt[lead], "muon2 pt", event.muon_pt[sublead]

    if event.muon_pt[lead] <= 3:
        continue
    if event.muon_pt[sublead] <= 3:
        continue


    print "muon1 eta", event.muon_eta[lead], "muon2 eta", event.muon_eta[sublead] 

    if np.abs(event.muon_eta[lead]) >= 2.4:
        continue
    if np.abs(event.muon_eta[sublead]) >= 2.4:
        continue


    print "muon 1 trackerlayerwm", event.ntrackerlayerswithmeasurement[lead], "muon2 trackerlayerwm", event.ntrackerlayerswithmeasurement[sublead]
        

    if event.ntrackerlayerswithmeasurement[lead] <= 5:                           
        continue  
    if event.ntrackerlayerswithmeasurement[sublead] <= 5:
        continue


    print "muon1 track iso", event.muon_trackIso[lead], "muon2 track iso", event.muon_trackIso[sublead]


    # if event.muon_trackIso[lead] >= 0.1:
    #     continue
    # if event.muon_trackIso[sublead] >= 0.1:
    #     continue


    print "* sign of charges", event.muon_q[lead]*event.muon_q[sublead]


    if event.muon_q[lead]*event.muon_q[sublead] > 0:
        continue

    mu1p_uncorr.SetPtEtaPhiM(event.muon_pt[lead], event.muon_eta[lead], event.muon_phi[lead], 0.1056583745)
    mu2p_uncorr.SetPtEtaPhiM(event.muon_pt[sublead], event.muon_eta[sublead], event.muon_phi[sublead], 0.1056583745)

    dxy1corr = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(event.muon_phi[lead]) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(event.muon_phi[lead])
    dxy2corr = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(event.muon_phi[sublead]) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(event.muon_phi[sublead])

    muon1_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[lead], event.muon_trketa[lead], event.muon_trkphi[lead], event.muon_trkdsz[lead], event.muon_dz[lead], event.muon_trklambda[lead], dxy1corr)
    muon1_phicorrected = ROOT.recalculate_phi_at_DV(muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], mu1p_uncorr.Px(), mu1p_uncorr.Py(), mu1p_uncorr.Pz(), event.muon_q[lead], event.dispvertex_x[bestdv], event.dispvertex_y[bestdv] )
    print "muon1 PCA", muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], "muon1_phicorr", muon1_phicorrected

    muon2_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[sublead], event.muon_trketa[sublead], event.muon_trkphi[sublead], event.muon_trkdsz[sublead], event.muon_dz[sublead], event.muon_trklambda[sublead], dxy2corr)
    muon2_phicorrected = ROOT.recalculate_phi_at_DV(muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], mu2p_uncorr.Px(), mu2p_uncorr.Py(), mu2p_uncorr.Pz(), event.muon_q[sublead], event.dispvertex_x[bestdv], event.dispvertex_y[bestdv] )
    print "muon2 PCA", muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], "muon2_phicorr", muon2_phicorrected


    mu1p.SetPtEtaPhiM(event.muon_pt[lead], event.muon_eta[lead], muon1_phicorrected, 0.1056583745)
    mu2p.SetPtEtaPhiM(event.muon_pt[sublead], event.muon_eta[sublead], muon2_phicorrected, 0.1056583745)
    motherp = mu1p + mu2p
    
    dRmuon1jetmin = 9999
    dRmuon2jetmin = 9999
    
    for jet in range(event.jet_num):

        jetp.SetPtEtaPhiM(event.jet_pt[jet], event.jet_eta[jet], event.jet_phi[jet], event.jet_m[jet])

        if mu1p.DeltaR(jetp) < dRmuon1jetmin: 
            dRmuon1jetmin = mu1p.DeltaR(jetp)

        if mu2p.DeltaR(jetp) < dRmuon1jetmin:
            dRmuon2jetmin = mu2p.DeltaR(jetp)


    print "muon1 min dR jet", dRmuon1jetmin, "muon2 min dR jet", dRmuon2jetmin


    # if dRmuon1jetmin <= 0.3:
    #     continue
    # if dRmuon2jetmin <= 0.3:
    #     continue
            
    motherv2D.SetXYZ(motherp.Px() , motherp.Py(), 0)
    motherv3D.SetXYZ(motherp.Px() , motherp.Py(), motherp.Pz())
    PVvector2D.SetXYZ(event.privertex_x[0], event.privertex_y[0], 0)
    BSvector2D.SetXYZ(event.BS_x, event.BS_y, 0)
    PVvector3D.SetXYZ(event.privertex_x[0], event.privertex_y[0], event.privertex_z[0])
    DVvector3D.SetXYZ(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv])
    DVvector2D.SetXYZ(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], 0)
    DVPVvector2D = DVvector2D - PVvector2D
    DVPVvector3D = DVvector3D - PVvector3D
    DVBSvector2D = DVvector2D - BSvector2D

    decaylength2Dunsigned = np.abs((DVPVvector2D.Dot(motherv2D))/(motherv2D.Mag())) 
    decaylength2Dsigned = (DVPVvector2D.Dot(motherv2D))/(motherv2D.Mag())
    decaylength2Dmod = np.sqrt(DVPVvector2D.Dot(DVPVvector2D))
    ctausigned = decaylength2Dsigned*(motherp.M()/motherp.Pt())
    angle = np.arccos((DVPVvector2D.Dot(motherv2D))/(motherv2D.Mag()*DVPVvector2D.Mag()))
    angle1 = motherv3D.DeltaPhi(DVPVvector3D)
    angle2 = motherv3D.Eta() - DVPVvector3D.Eta()


    print "lxy", decaylength2Dmod, "dphi DVPV dimupt", angle1, "cos(angle1)", np.cos(angle1), "dphi muon1 muon2", mu1p.DeltaPhi(mu2p)


    if decaylength2Dmod >= 11:
        continue
    if np.abs(angle1) >= 0.02:
       continue
    if np.cos(angle1) <= 0:
        continue
    if np.abs(mu1p.DeltaPhi(mu2p)) >= 2.8:
        continue


    print "I have passed and am getting filled in the tree woohoo"


    run[0] = tree_mu.Run
    lumi[0] = tree_mu.Lumi
    event_num[0] = tree_mu.Event
    dimuon_pt[0] = ((mu1p + mu2p).Pt())
    dimuon_mass[0] = ((mu1p + mu2p).M())
    dimuon_mass_uncorr[0] = ((mu1p_uncorr + mu2p_uncorr).M())
    lxy[0] = decaylength2Dmod
    dphidimudv[0] = angle1
    detadimudv[0] = angle2
    ctau[0] = ctausigned
    muon1_dxy[0] = event.muon_dxy[lead]
    muon2_dxy[0] = event.muon_dxy[sublead]
    muon1_edxy[0] = event.muon_edxy[lead]
    muon2_edxy[0] = event.muon_edxy[sublead]
    muon1_dxycorr[0] = dxy1corr
    muon2_dxycorr[0] = dxy2corr
#    muon1_dxycorr[0] = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(mu1p.Phi()) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(mu1p.Phi()) 
#    muon2_dxycorr[0] = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(mu2p.Phi()) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(mu2p.Phi())
    muon1_trkiso[0] = event.muon_trackIso[lead]
    muon2_trkiso[0] = event.muon_trackIso[sublead]
    muon1_chi2overndof[0] = event.muon_chi2[lead]/event.muon_ndof[lead]
    muon2_chi2overndof[0] = event.muon_chi2[sublead]/event.muon_ndof[sublead]
    muon1_pt[0] = mu1p.Pt()
    muon2_pt[0] = mu2p.Pt()
    muon1_eta[0] = mu1p.Eta()
    muon2_eta[0] = mu2p.Eta()
    muon1_phi[0] = mu1p.Phi()
    muon2_phi[0] = mu2p.Phi()
    muon1_phicorr[0] = muon1_phicorrected
    muon2_phicorr[0] = muon2_phicorrected
    PVx[0] = event.privertex_x[0]
    PVy[0] = event.privertex_y[0]
    PVz[0] = event.privertex_z[0]
    BSx[0] = event.BS_x
    BSy[0] = event.BS_y
    BSz[0] = event.BS_z
    DVx[0] = event.dispvertex_x[bestdv]
    DVy[0] = event.dispvertex_y[bestdv]
    DVz[0] = event.dispvertex_z[bestdv]
    DVxerr[0] = event.dispvertex_ex[bestdv]
    DVyerr[0] = event.dispvertex_ey[bestdv]
    DVzerr[0] = event.dispvertex_ez[bestdv]
    DV_chi2overndof[0] = event.dispvertex_chi2[bestdv]/event.dispvertex_ndof[bestdv]
    detamumu[0] = mu1p.Eta() - mu2p.Eta()
    dphimumu[0] = mu1p.DeltaPhi(mu2p)
    dRmumu[0] = mu1p.DeltaR(mu2p)
    dRmuon1jet[0] = dRmuon1jetmin
    dRmuon2jet[0] = dRmuon2jetmin

    tree.Fill()

print "number of events in the file", numberofevents

#DeltaRmumu.Write()
#muchi2overndof.Write()
#ntrackerlayerswithmeasurement.Write()
#nvalidmuonhits.Write()
#nvalidpixelhits.Write()
#nvalidstriphits.Write()

outfile.Write()
outfile.Close()


