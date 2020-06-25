# import ROOT in batch mode                                                                                                          
import os
import sys

#import sys                                                                                                                          
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ROOT.gROOT.LoadMacro("calculate_pixel.cc")
ROOT.gROOT.LoadMacro("propagation_utils.cc")

# print ROOT.point_in_which_module(0., 0., 0.)

import numpy as np
from array import array

from ROOT import TH1F, TH1D, TH2D, TFile, TLorentzVector, TVector3, TChain, TProfile, TTree

# load FWLite C++ libraries                                                                                                          
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries                                                                                                       
from DataFormats.FWLite import Handle, Events


# file_mu = ROOT.TFile("/cms/scoutingmuon/hardik/condor_output/ggPhi_ntuples/ggPhi_m2_ct100_ntuple.root")                           
# tree_mu = file_mu.Get('scoutingntuplizer')

file_mu = ROOT.TFile('/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017D_v1/200307_225553/0000/scouting_ntuple_{0}.root'.format(sys.argv[1]))
tree_mu = file_mu.Get('scoutingntuplizer')

# outfile = ROOT.TFile('histos.root', 'recreate')
# outfile.cd()


print "Running on file", file_mu

outfile = TFile("output/2017/D/flatdimuontree_{0}.root".format(sys.argv[1]), "recreate")
outfile.cd()

# ttree = ROOT.TTree('events', 'flat tree with muon, jet, met info')

# run = array("f", [0.0])
# lumi = array("f", [0.0])
# event_num = array("d", [0.0])
# met_pt = array("f", [0.0])
# met_eta = array("f", [0.0])
# met_phi = array("f", [0.0])
# num_muons = array("i", [0])
# muon_pt = array("f", maxn*[0.0])
# muon_eta = array("f", maxn*[0.0])
# muon_phi = array("f", maxn*[0.0])
# muon_phicorr = array("f", maxn*[0.0])
# muon_q = array("f", maxn*[0.0])
# muon_excesspixelhits = array("f", maxn*[0.0])
# muon_chi2overndof = array("f", maxn*[0.0])
# muon_dxy = array("f", maxn*[0.0])
# muon_dxycorr = array("f", maxn*[0.0])
# num_jets = array("i", [0])
# jet_pt = array("f", maxn*[0.0])
# jet_eta = array("f", maxn*[0.0])
# jet_phi = array("f", maxn*[0.0])
# jet_mass = array("f", maxn*[0.0])
# jet_csv = array("f", maxn*[0.0])
# PVx = array("f", [0.0])
# PVy = array("f", [0.0])
# PVz = array("f", [0.0])
# num_DV = array("i", [0])
# DVx = array("f", maxn*[0.0])
# DVy = array("f", maxn*[0.0])
# DVz = array("f", maxn*[0.0])
# DVxerr = array("f", maxn*[0.0])
# DVyerr = array("f", maxn*[0.0])
# DVzerr = array("f", maxn*[0.0])
# DV_chi2overndof = array("f", maxn*[0.0])

# ttree.Branch("run", run, "run/I")
# ttree.Branch("lumi", lumi, "lumi/I")
# ttree.Branch("event_num", event_num, "event_num/I")
# ttree.Branch("met_pt", met_pt, "met_pt/F")
# ttree.Branch("met_eta", met_eta, "met_eta/F")
# ttree.Branch("met_phi", met_pt, "met_phi/F")
# ttree.Branch("num_muons", num_muons, "num_muons/I")
# ttree.Branch("muon_pt", muon_pt, "muon_pt[num_muons]/F")
# ttree.Branch("muon_eta", muon_eta, "muon_eta[num_muons]/F")
# ttree.Branch("muon_phi", muon_phi, "muon_phi[num_muons]/F")
# ttree.Branch("muon_phicorr", muon_phicorr, "muon_phicorr[num_muons]/F")
# ttree.Branch("muon_excesspixelhits", muon_excesspixelhits, "muon_excesspixelhits[num_muons]/F")
# ttree.Branch("muon_q", muon_q, "muon_q[num_muons]/F")
# ttree.Branch("muon_chi2overndof", muon_chi2overndof, "muon_chi2overndof[num_muons]/F")
# ttree.Branch("muon_dxy", muon_dxy, "muon_dxy[num_muons]/F")
# ttree.Branch("muon_dxycorr", muon_dxycorr, "muon_dxycorr[num_muons]/F")
# ttree.Branch("num_jets", num_jets, "num_jets/I")
# ttree.Branch("jet_pt", jet_pt, "jet_pt[num_jets]/F")
# ttree.Branch("jet_eta", jet_eta, "jet_eta[num_jets]/F")
# ttree.Branch("jet_phi", jet_phi, "jet_phi[num_jets]/F")
# ttree.Branch("jet_mass", jet_mass, "jet_mass[num_jets]/F")
# ttree.Branch("jet_csv", jet_csv, "jet_csv[num_jets]/F")
# ttree.Branch("PVx", PVx, "PVx/F")
# ttree.Branch("PVy", PVy, "PVy/F")
# ttree.Branch("PVz", PVz, "PVz/F")
# ttree.Branch("num_DV", num_DV, "num_DV/I")
# ttree.Branch("DVx", DVx, "DVx[num_DV]/F")
# ttree.Branch("DVy", DVy, "DVy[num_DV]/F")
# ttree.Branch("DVz", DVz, "DVz[num_DV]/F")
# ttree.Branch("DVxerr", DVxerr, "DVxerr[num_DV]/F")
# ttree.Branch("DVyerr", DVyerr, "DVyerr[num_DV]/F")
# ttree.Branch("DVzerr", DVzerr, "DVzerr[num_DV]/F")
# ttree.Branch("DV_chi2overndof", DV_chi2overndof, "DV_chi2overndof[num_DV]/F")
# ttree.Branch("dRmuon_closestjet", dRmuon_closestjet, "dRmuon_closestjet[num_muons]/F")



############ Declare Histos #####################


munu_transversemass = ROOT.TH1F("munu_transversemass", "munu_transversemass", 1000, 0, 10)

dimuon_mass_corr = ROOT.TH1F("dimuon_mass_corr", "dimuon_mass_corr", 1000, 0, 10) 
dimuon_mass_uncorr = ROOT.TH1F("dimuon_mass_uncorr", "dimuon_mass_uncorr", 1000, 0, 10)


#DeltaRmumu = ROOT.TH1F("DeltaRmumu", "DeltaRmumu", 100, -10, 10)
#muchi2overndof = ROOT.TH1F("muchi2overndof", "muchi2overndof", 100, -10, 10)
#ntrackerlayerswithmeasurement = ROOT.TH1F("ntrackerlayerswithmeasurement", "ntrackerlayerswithmeasurement", 100, -10, 10)
#nvalidmuonhits = ROOT.TH1F("nvalidmuonhits", "nvalidmuonhits", 100, -10, 10)
#nvalidpixelhits = ROOT.TH1F("nvalidpixelhits", "nvalidpixelhits", 100, -10, 10)
#nvalidstriphits = ROOT.TH1F("nvalidstriphits", "nvalidstriphits", 100, -10, 10)



numberofevents = 0

for i, event in enumerate(tree_mu):

    numberofevents += 1

    # if numberofevents > 10000:
    #     break

    ################Initialize four vectors , three vectors, variables etc ####################

    mup = TLorentzVector(0, 0, 0, 0)
    nup = TLorentzVector(0, 0, 0, 0)

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
    passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 0
    passL1_DoubleMu_15_7 = 0
    passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 0  
    passL1 = 0
    passDST_DoubleMu3_noVtx_CaloScouting_v = 0
    passDST_DoubleMu1_noVtx_CaloScouting_v = 0
    passHLT = 0



    # print "-----In Run" ,tree_mu.Run , "---Lumi", tree_mu.Lumi, "----Event", tree_mu.Event, "-----------"

    ############### Select L1 and HLT seeds #####################


    for j, _ in enumerate(event.l1bitmap):

        # if event.l1bitmap[j].first == "L1_DoubleMu4_SQ_OS_dR_Max1p2" and event.l1bitmap[j].second == 1:
        #     passL1_DoubleMu4_SQ_OS_dR_Max1p2 = 1

        if event.l1bitmap[j].first == "L1_DoubleMu4p5_SQ_OS_dR_Max1p2" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 1

        if event.l1bitmap[j].first == "L1_DoubleMu_15_7" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu_15_7 = 1

        if event.l1bitmap[j].first == "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 1

    # passL1 = passL1_DoubleMu4_SQ_OS_dR_Max1p2 + passL1_DoubleMu_15_7 + passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4
    passL1 = passL1_DoubleMu4_SQ_OS_dR_Max1p2 + passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 + passL1_DoubleMu_15_7 + passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4

    for j, _ in enumerate(event.hltbitmap):

        if event.hltbitmap[j].first == "DST_DoubleMu3_noVtx_CaloScouting_v" and event.hltbitmap[j].second == 1:
            passDST_DoubleMu3_noVtx_CaloScouting_v = 1

    for j, _ in enumerate(event.hltbitmap):

        if event.hltbitmap[j].first == "DST_DoubleMu1_noVtx_CaloScouting_v" and event.hltbitmap[j].second == 1:
            passDST_DoubleMu1_noVtx_CaloScouting_v = 1

    passHLT = passDST_DoubleMu3_noVtx_CaloScouting_v + passDST_DoubleMu1_noVtx_CaloScouting_v

#    print "In event", i , "passL1?", passL1, "passHLT?", passHLT
        
    # if passL1 == 0:
    #     continue

    if passHLT == 0:
        continue

    if len(event.muon_pt) < 1:
        continue

    # if event.dispvertex_num == 0:
    #     continue

    # print "Number of muons", len(event.muon_pt), "Number of vertices", event.dispvertex_num


######################Start munu Analysis####################################################################


    if len(event.muon_pt) == 1:

        if event.muon_pt[0] > 25 and event.met_pt > 0:

            mup.SetPtEtaPhiM(event.muon_pt[0], event.muon_eta[0], event.muon_phi[0], 0.1056583745)
            nup.SetPtEtaPhiM(event.met_pt, event.met_eta, event.met_phi, 0)

            print "munu transverse mass", (mup+nup).Mt()

            munu_transversemass.Fill((mup+nup).Mt())


#####################End munu Analysis########################################################################


#########################Start Displaced Dimuon Analysis#######################################################

    if len(event.muon_pt) < 2:
        continue

    if event.dispvertex_num < 1:                                                                                                  
        continue 

    ###### Decay Vertex Quality Cut ########

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
        # if event.dispvertex_tracksize[dv] != 2:
        #     print "track size in DV--", dv, "----is----", event.dispvertex_tracksize[dv]
        #     continue
        if np.maximum(event.dispvertex_ex[dv],event.dispvertex_ey[dv]) < mindverror:
            mindverror = np.maximum(event.dispvertex_ex[dv],event.dispvertex_ey[dv])            
            bestdv = dv
    
    # print "best dv index(if -1 no best DV)", bestdv

    if bestdv == -1:
        continue

    ########Excess Pixel Hits Calculation#################

    # NModules = 1856
    # for i in range(1856)
    # print ROOT.point_in_which_module(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv])
    # print ROOT.dist_to_imodule_plane(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv], 1878)     
         
    # if ROOT.point_in_which_module(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv]) != -1:
    #     continue


    # distPixelmin = 9999

    for module in range(1856):
        imodule = ROOT.point_in_which_module(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv])
    # print imodule
    distPixelmin = ROOT.dist_to_imodule_plane(event.dispvertex_x[bestdv], event.dispvertex_y[bestdv], event.dispvertex_z[bestdv], imodule)
    # print distPixelmin

    ########Associating best DV to chosen muon pair###################


    numberofmuons_invtx = 0
    highestpt = -1
    lowestpt = 9999
    sublead = -1
    lead = -1


    for index in range(len(event.vertex_index)):
        for subindex in event.vertex_index[index]:
                
            # print subindex, "muon",index
                
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
        

    # print "Number of muons in vtx", numberofmuons_invtx
    # print "leading muon index", lead
    # print "subleading muon index", sublead
        

    if numberofmuons_invtx != 2:
        continue
         


    #######Muon Quality Cuts#########################



    # print "DV chi2overndof", event.dispvertex_chi2[bestdv]/event.dispvertex_ndof[bestdv], "DV_xerr", event.dispvertex_ex[bestdv], "DVyerr", event.dispvertex_ey[bestdv], "DVzerr", event.dispvertex_ez[bestdv]
 

    # print "muon1 chi2overndof", event.muon_chi2[lead]/event.muon_ndof[lead], "muon2 chi2overndof", event.muon_chi2[sublead]/event.muon_ndof[sublead]


    if event.muon_chi2[lead]/event.muon_ndof[lead] >= 3:
        continue
    if event.muon_chi2[sublead]/event.muon_ndof[sublead] >= 3:
        continue


    # print "muon1 pt", event.muon_pt[lead], "muon2 pt", event.muon_pt[sublead]

    if event.muon_pt[lead] <= 3:
        continue
    if event.muon_pt[sublead] <= 3:
        continue


    # print "muon1 eta", event.muon_eta[lead], "muon2 eta", event.muon_eta[sublead] 

    if np.abs(event.muon_eta[lead]) >= 2.4:
        continue
    if np.abs(event.muon_eta[sublead]) >= 2.4:
        continue


    # print "muon 1 trackerlayerwm", event.ntrackerlayerswithmeasurement[lead], "muon2 trackerlayerwm", event.ntrackerlayerswithmeasurement[sublead]
        

    if event.ntrackerlayerswithmeasurement[lead] <= 5:                           
        continue  
    if event.ntrackerlayerswithmeasurement[sublead] <= 5:
        continue


    ############Muon track iso cut##############################3

    # print "muon1 track iso", event.muon_trackIso[lead], "muon2 track iso", event.muon_trackIso[sublead]


    # if event.muon_trackIso[lead] >= 0.1:
    #     continue
    # if event.muon_trackIso[sublead] >= 0.1:
    #     continue


    #######Main part of Dimuon Analysis with extra cuts and phi correction at DV##############


    # print "* sign of charges", event.muon_q[lead]*event.muon_q[sublead]


    if event.muon_q[lead]*event.muon_q[sublead] > 0:
        continue

    mu1p_uncorr.SetPtEtaPhiM(event.muon_pt[lead], event.muon_eta[lead], event.muon_phi[lead], 0.1056583745)
    mu2p_uncorr.SetPtEtaPhiM(event.muon_pt[sublead], event.muon_eta[sublead], event.muon_phi[sublead], 0.1056583745)

    dxy1corr = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(event.muon_phi[lead]) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(event.muon_phi[lead])
    dxy2corr = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(event.muon_phi[sublead]) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(event.muon_phi[sublead])

    muon1_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[lead], event.muon_trketa[lead], event.muon_trkphi[lead], event.muon_trkdsz[lead], event.muon_dz[lead], event.muon_trklambda[lead], dxy1corr)
    muon1_phicorrected = ROOT.recalculate_phi_at_DV(muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], mu1p_uncorr.Px(), mu1p_uncorr.Py(), mu1p_uncorr.Pz(), event.muon_q[lead], event.dispvertex_x[bestdv], event.dispvertex_y[bestdv] )
    # print "muon1 PCA", muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], "muon1_phicorr", muon1_phicorrected

    muon2_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[sublead], event.muon_trketa[sublead], event.muon_trkphi[sublead], event.muon_trkdsz[sublead], event.muon_dz[sublead], event.muon_trklambda[sublead], dxy2corr)
    muon2_phicorrected = ROOT.recalculate_phi_at_DV(muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], mu2p_uncorr.Px(), mu2p_uncorr.Py(), mu2p_uncorr.Pz(), event.muon_q[sublead], event.dispvertex_x[bestdv], event.dispvertex_y[bestdv] )
    # print "muon2 PCA", muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], "muon2_phicorr", muon2_phicorrected

    mu1p.SetPtEtaPhiM(event.muon_pt[lead], event.muon_eta[lead], muon1_phicorrected, 0.1056583745)
    mu2p.SetPtEtaPhiM(event.muon_pt[sublead], event.muon_eta[sublead], muon2_phicorrected, 0.1056583745)
    motherp = mu1p + mu2p


    ##################Muon dR(mu,closest jet) Isolation cut##############################

    dRmuon1jetmin = 9999
    dRmuon2jetmin = 9999

    highestjetpt = -1
    leadjet = -1
    btagvalue = -99
    mvavalue = -99
    
    for jet in range(event.jet_num):

        jetp.SetPtEtaPhiM(event.jet_pt[jet], event.jet_eta[jet], event.jet_phi[jet], event.jet_m[jet])

        if mu1p.DeltaR(jetp) < dRmuon1jetmin: 
            dRmuon1jetmin = mu1p.DeltaR(jetp)

        if mu2p.DeltaR(jetp) < dRmuon1jetmin:
            dRmuon2jetmin = mu2p.DeltaR(jetp)

        if event.jet_pt[jet] > highestjetpt:
            highestjetpt = event.jet_pt[jet]
            leadjet = jet

    if leadjet != -1:
        btagvalue = event.jet_btagDiscriminator[leadjet]
        mvavalue = event.jet_mvaDiscriminator[leadjet]
    # print "muon1 min dR jet", dRmuon1jetmin, "muon2 min dR jet", dRmuon2jetmin

    # if dRmuon1jetmin <= 0.3:
    #     continue
    # if dRmuon2jetmin <= 0.3:
    #     continue
            
    ############Calculation of decay length and other dimuon variables#########33

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

    # print "lxy", decaylength2Dmod, "dphi DVPV dimupt", angle1, "cos(angle1)", np.cos(angle1), "dphi muon1 muon2", mu1p.DeltaPhi(mu2p)

    #############More cuts################################33

    if decaylength2Dmod >= 11:
        continue
    # if np.abs(angle1) >= 0.02:
    #     continue
    if np.cos(angle1) <= 0:
        continue
    if np.abs(mu1p.DeltaPhi(mu2p)) >= 2.8:
        continue

    passexcesshits = 0

    if decaylength2Dmod < 3.5:
        passexcesshits = 1
    elif (event.nvalidpixelhits[lead] - event.nexpectedhitsmultiple[lead]) <= 0 and (event.nvalidpixelhits[sublead] - event.nexpectedhitsmultiple[sublead]) <= 0:
        passexcesshits = 1
        
    if passexcesshits != 1:
        continue


    ########Fill Dimuon Histos##############

    dimuon_mass_corr.Fill((mu1p + mu2p).M())
    dimuon_mass_uncorr.Fill((mu1p_uncorr + mu2p_uncorr).M())
    
    #######End Dimuon Histos################
    
    # print "I have passed and am getting filled in the flat tree woohoo"


    # run[0] = tree_mu.Run
    # lumi[0] = tree_mu.Lumi
    # event_num[0] = tree_mu.Event
    # dimuon_pt[0] = ((mu1p + mu2p).Pt())
    # dimuon_mass[0] = ((mu1p + mu2p).M())
    # dimuon_mass_uncorr[0] = ((mu1p_uncorr + mu2p_uncorr).M())
    # lxy[0] = decaylength2Dmod
    # dphidimudv[0] = angle1
    # detadimudv[0] = angle2
    # ctau[0] = ctausigned
    # muon1_dxy[0] = event.muon_dxy[lead]
    # muon2_dxy[0] = event.muon_dxy[sublead]
    # muon1_edxy[0] = event.muon_edxy[lead]
    # muon2_edxy[0] = event.muon_edxy[sublead]
    # muon1_dxycorr[0] = dxy1corr
    # muon2_dxycorr[0] = dxy2corr
    # # muon1_dxycorr[0] = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(mu1p.Phi()) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(mu1p.Phi()) 
    # # muon2_dxycorr[0] = -1*(event.dispvertex_x[bestdv] - event.privertex_x[0])*np.sin(mu2p.Phi()) + (event.dispvertex_y[bestdv] - event.privertex_y[0])*np.cos(mu2p.Phi())
    # muon1_trkiso[0] = event.muon_trackIso[lead]
    # muon2_trkiso[0] = event.muon_trackIso[sublead]
    # muon1_chi2overndof[0] = event.muon_chi2[lead]/event.muon_ndof[lead]
    # muon2_chi2overndof[0] = event.muon_chi2[sublead]/event.muon_ndof[sublead]
    # muon1_pt[0] = mu1p.Pt()
    # muon2_pt[0] = mu2p.Pt()
    # muon1_eta[0] = mu1p.Eta()
    # muon2_eta[0] = mu2p.Eta()
    # muon1_phi[0] = mu1p.Phi()
    # muon2_phi[0] = mu2p.Phi()
    # muon1_phicorr[0] = muon1_phicorrected
    # muon2_phicorr[0] = muon2_phicorrected
    # muon1_excesspixelhits[0] = event.nvalidpixelhits[lead] - event.nexpectedhitsmultiple[lead]
    # muon2_excesspixelhits[0] = event.nvalidpixelhits[sublead] - event.nexpectedhitsmultiple[sublead]
    # PVx[0] = event.privertex_x[0]
    # PVy[0] = event.privertex_y[0]
    # PVz[0] = event.privertex_z[0]
    # BSx[0] = event.BS_x
    # BSy[0] = event.BS_y
    # BSz[0] = event.BS_z
    # DVx[0] = event.dispvertex_x[bestdv]
    # DVy[0] = event.dispvertex_y[bestdv]
    # DVz[0] = event.dispvertex_z[bestdv]
    # DVxerr[0] = event.dispvertex_ex[bestdv]
    # DVyerr[0] = event.dispvertex_ey[bestdv]
    # DVzerr[0] = event.dispvertex_ez[bestdv]
    # DV_chi2overndof[0] = event.dispvertex_chi2[bestdv]/event.dispvertex_ndof[bestdv]
    # distPixel[0] = distPixelmin
    # detamumu[0] = mu1p.Eta() - mu2p.Eta()
    # dphimumu[0] = mu1p.DeltaPhi(mu2p)
    # dRmumu[0] = mu1p.DeltaR(mu2p)
    # dRmuon1jet[0] = dRmuon1jetmin
    # dRmuon2jet[0] = dRmuon2jetmin
    # njet[0] = event.jet_num
    # btag[0] = btagvalue
    # mva[0] = mvavalue

    # ttree.Fill()


#########################################End Dimuon Analysis###################################################


print "number of events in the file", numberofevents

munu_transversemass.Write()
dimuon_mass_corr.Write()
dimuon_mass_uncorr.Write()

#DeltaRmumu.Write()
#muchi2overndof.Write()
#ntrackerlayerswithmeasurement.Write()
#nvalidmuonhits.Write()
#nvalidpixelhits.Write()
#nvalidstriphits.Write()

outfile.Write()
outfile.Close()


