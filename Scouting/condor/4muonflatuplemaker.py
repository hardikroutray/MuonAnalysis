# import ROOT in batch mode                                                  
                                                       
import os
import sys

#import sys                                                                  
                                                       
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--MC17", default=False, action="store_true")
parser.add_argument("--MC18", default=False, action="store_true")
parser.add_argument("--Data17", default=False, action="store_true")
parser.add_argument("--Data18", default=False, action="store_true")

parser.add_argument("--storetreepreobjsel", default=False, action="store_true")
parser.add_argument("--storeDV", default=False, action="store_true")
# parser.add_argument("--dobaselinenoisocuts", default=False, action="store_true")
# parser.add_argument("--dobaselineisocuts", default=False, action="store_true")
# parser.add_argument("--doextracuts", default=False, action="store_true")
# parser.add_argument("--dodxycuts", default=False, action="store_true")

parser.add_argument("inputpath", help="path where input files are located",type=str)
parser.add_argument("file_number", help="file to run on",type=int)

args = parser.parse_args()

if args.MC17:
    print "running on", "2017", "MC file", args.file_number
elif args.MC18:
    print "running on", "2018", "MC file", args.file_number
elif args.Data17:
    print "running on", args.inputpath, "Data file", args.file_number
elif args.Data18:
    print "running on", args.inputpath, "Data file", args.file_number
else:
    print "Not 2017 or 2018 MC/Data"
    exit()

ROOT.gROOT.LoadMacro("calculate_pixel.cc")
ROOT.gROOT.LoadMacro("vertexing_utils.cc")

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

'''                             
file_mu = ROOT.TFile("/cms/scoutingmuon/hardik/condor_output/BTo4mu_ntuples/B0To4mu_mphi1_full.root")
tree_mu = file_mu.Get('scoutingntuplizer')
outfile = ROOT.TFile('flat_dimuon_tree_B0To4mu_mphi1_flattree_vf1_test.root', 'recreate')
outfile.cd()
'''

#'''
if args.inputpath == "2017-C":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017C_v1/200307_225537/0000/'
elif args.inputpath == "2017-D":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017D_v1/200307_225553/0000/'
elif args.inputpath == "2017-E":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017E_v1/200307_225609/0000/'
elif args.inputpath == "2017-F":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2017F_v1/200308_000417/0000/'
elif args.inputpath == "2018-A":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2018A_v1/200308_060230/0000/'
elif args.inputpath == "2018-B":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2018B_v1/200308_060242/0000/' 
elif args.inputpath == "2018-C":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2018C_v1/200308_060301/0000/'
elif args.inputpath == "2018-D":
    path = '/cms/routray/crab_output/muontuples_10percent_v2/ScoutingCaloMuon/ScoutingCaloMuon_Ntuples_Run2018D_v1/200308_060314/0000/'

if os.path.isfile(path + '/scouting_ntuple_{}.root'.format(args.file_number)) != 1:
    print "file does not exist"
    exit()

file_mu = ROOT.TFile(path + 'scouting_ntuple_{}.root'.format(args.file_number))
tree_mu = file_mu.Get('scoutingntuplizer')

if args.storeDV:
    outfile = TFile("flat4muontrees_withDV/{}/{}/flatdimuontree_{}.root".format(args.inputpath.rsplit("-")[0], args.inputpath.split("-")[1], args.file_number), "recreate")
else:
    outfile = TFile("flat4muontrees/{}/{}/flatdimuontree_{}.root".format(args.inputpath.rsplit("-")[0], args.inputpath.split("-")[1], args.file_number), "recreate")

outfile.cd()
#'''

print file_mu
print outfile

if args.storetreepreobjsel:

    tree0 = ROOT.TTree('events0', 'flat tree with dimuon info-preobjsel')

    run0 = array("f", [0.0])
    lumi0 = array("f", [0.0])
    event_num0 = array("d", [0.0])
    L10 = array("f", [0.0])
    HLT0 = array("f", [0.0])
    # passobjsel = array("f", [0.0])

    tree0.Branch('run0', run0, 'run0/F')
    tree0.Branch('lumi0', lumi0, 'lumi0/F')
    tree0.Branch('event_num0', event_num0, 'event_num0/D')
    tree0.Branch('L10', L10, 'L10/F')
    tree0.Branch('HLT0', HLT0, 'HLT0/F')
    # tree0.Branch('passobjsel', passobjsel, 'passobjsel/F')

    for i0, event0 in enumerate(tree_mu):

        pass0L1_DoubleMu4_SQ_OS_dR_Max1p2 = 0
        pass0L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 0
        pass0L1_DoubleMu_15_7 = 0
        pass0L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 0
        pass0L1 = 0
        pass0DST_DoubleMu3_noVtx_CaloScouting_v = 0
        pass0DST_DoubleMu1_noVtx_CaloScouting_v = 0
        pass0HLT = 0

        for j0, _0 in enumerate(event0.l1bitmap):

            if args.MC17 or args.Data17:

                if event0.l1bitmap[j0].first == "L1_DoubleMu4_SQ_OS_dR_Max1p2" and event0.l1bitmap[j0].second == 1:
                    pass0L1_DoubleMu4_SQ_OS_dR_Max1p2 = 1
                else:
                    pass

            elif args.MC18 or args.Data18:

                if event0.l1bitmap[j0].first == "L1_DoubleMu4p5_SQ_OS_dR_Max1p2" and event0.l1bitmap[j0].second == 1:
                    pass0L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 1
                else:
                    pass

            else:
                pass

            if event0.l1bitmap[j0].first == "L1_DoubleMu_15_7" and event0.l1bitmap[j0].second == 1:
                pass0L1_DoubleMu_15_7 = 1
            else:
                pass

            if event0.l1bitmap[j0].first == "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4" and event0.l1bitmap[j0].second == 1:
                pass0L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 1
            else:
                pass

        pass0L1 = pass0L1_DoubleMu4_SQ_OS_dR_Max1p2 + pass0L1_DoubleMu4p5_SQ_OS_dR_Max1p2 + pass0L1_DoubleMu_15_7 + pass0L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4

        for j0, _0 in enumerate(event0.hltbitmap):

            if event0.hltbitmap[j0].first == "DST_DoubleMu3_noVtx_CaloScouting_v" and event0.hltbitmap[j0].second == 1:
                pass0DST_DoubleMu3_noVtx_CaloScouting_v = 1

        for j0, _ in enumerate(event0.hltbitmap):

            if event0.hltbitmap[j0].first == "DST_DoubleMu1_noVtx_CaloScouting_v" and event0.hltbitmap[j0].second == 1:
                pass0DST_DoubleMu1_noVtx_CaloScouting_v = 1

        pass0HLT = pass0DST_DoubleMu3_noVtx_CaloScouting_v + pass0DST_DoubleMu1_noVtx_CaloScouting_v

        run0[0] = tree_mu.Run
        lumi0[0] = tree_mu.Lumi
        event_num0[0] = tree_mu.Event
        L10[0] = pass0L1
        HLT0[0] = pass0HLT

        tree0.Fill()

else:
    pass

tree = ROOT.TTree('events', 'flat tree with dimuon info-L1HLTobjsel')

run = array("f", [0.0])
lumi = array("f", [0.0])
event_num = array("d", [0.0])
L1 = array("f", [0.0])
HLT = array("f", [0.0])

fourmuonmass = array("f", [0.0])

M01 = array("f", [0.0])
M02 = array("f", [0.0])
M03 = array("f", [0.0])
M12 = array("f", [0.0])
M13 = array("f", [0.0])
M23 = array("f", [0.0])

pt01 = array("f", [0.0])
pt02 = array("f", [0.0])
pt03 = array("f", [0.0])
pt12 = array("f", [0.0])
pt13 = array("f", [0.0])
pt23 = array("f", [0.0])

pair1_mass = array("f", [0.0])
pair2_mass = array("f", [0.0])
pair1_pt = array("f", [0.0])
pair2_pt = array("f", [0.0])
pass_masspair = array("i", [0])

jpsipair_mass = array("f", [0.0])
phipair_mass = array("f", [0.0])
jpsipair_pt = array("f", [0.0])
phipair_pt = array("f", [0.0])
pass_jpsiphi = array("i", [0])

pass_ZZ = array("i", [0])

mu0q = array("i", [0])
mu1q = array("i", [0])
mu2q = array("i", [0])
mu3q = array("i", [0])

mu0pt = array("f", [0.0])
mu1pt = array("f", [0.0])
mu2pt = array("f", [0.0])
mu3pt = array("f", [0.0])

mu0eta = array("f", [0.0])
mu1eta = array("f", [0.0])
mu2eta = array("f", [0.0])
mu3eta = array("f", [0.0])

mu0phi = array("f", [0.0])
mu1phi = array("f", [0.0])
mu2phi = array("f", [0.0])
mu3phi = array("f", [0.0])

mu0trkiso = array("f", [0.0])
mu1trkiso = array("f", [0.0])
mu2trkiso = array("f", [0.0])
mu3trkiso = array("f", [0.0])

mu0dxy = array("f", [0.0])
mu1dxy = array("f", [0.0])
mu2dxy = array("f", [0.0])
mu3dxy = array("f", [0.0])
mu0edxy = array("f", [0.0])
mu1edxy = array("f", [0.0])
mu2edxy = array("f", [0.0])
mu3edxy = array("f", [0.0])

PVx = array("f", [0.0])
PVy = array("f", [0.0])
PVz = array("f", [0.0])

scoutingDVx = array("f", [0.0])
scoutingDVy = array("f", [0.0])
scoutingDVz = array("f", [0.0])
scoutingDVlen = array("i", [0])

reDVx = array("f", [0.0])
reDVy = array("f", [0.0])
reDVz = array("f", [0.0])
reDCA = array("f", [0.0])
remaxDCA = array("f", [0.0])

if args.MC17 or args.MC18:
    genDVx = array("f", [0.0])
    genDVy = array("f", [0.0])
    genDVz = array("f", [0.0])

#tree.Branch('muon1_pt', muon1_pt, 'muon1_pt/F')
tree.Branch('run', run, 'run/F')
tree.Branch('lumi', lumi, 'lumi/F')
tree.Branch('event_num', event_num, 'event_num/D')
tree.Branch('L1', L1, 'L1/F')
tree.Branch('HLT', HLT, 'HLT/F')

tree.Branch('fourmuonmass', fourmuonmass, 'fourmuonmass/F')

tree.Branch('M01', M01, 'M01/F')
tree.Branch('M02', M02, 'M02/F')
tree.Branch('M03', M03, 'M03/F')
tree.Branch('M12', M12, 'M12/F')
tree.Branch('M13', M13, 'M13/F')
tree.Branch('M23', M23, 'M23/F')

tree.Branch('pt01', pt01, 'pt01/F')
tree.Branch('pt02', pt02, 'pt02/F')
tree.Branch('pt03', pt03, 'pt03/F')
tree.Branch('pt12', pt12, 'pt12/F')
tree.Branch('pt13', pt13, 'pt13/F')
tree.Branch('pt23', pt23, 'pt23/F')

tree.Branch('pair1_mass', pair1_mass, 'pair1_mass/F')
tree.Branch('pair2_mass', pair2_mass, 'pair2_mass/F')
tree.Branch('pair1_pt', pair1_pt, 'pair1_pt/F')
tree.Branch('pair2_pt', pair2_pt, 'pair2_pt/F')
tree.Branch('pass_masspair', pass_masspair, 'pass_masspair/I')

tree.Branch('jpsipair_mass', jpsipair_mass, 'jpsipair_mass/F')
tree.Branch('phipair_mass', phipair_mass, 'phipair_mass/F')
tree.Branch('jpsipair_pt', jpsipair_pt, 'jpsipair_pt/F')
tree.Branch('phipair_pt', phipair_pt, 'phipair_pt/F')
tree.Branch('pass_jpsiphi', pass_jpsiphi, 'pass_jpsiphi/I')

tree.Branch('pass_ZZ', pass_ZZ, 'pass_ZZ/I')

tree.Branch('mu0q', mu0q, 'mu0q/I')
tree.Branch('mu1q', mu1q, 'mu1q/I')
tree.Branch('mu2q', mu2q, 'mu2q/I')
tree.Branch('mu3q', mu3q, 'mu3q/I')

tree.Branch('mu0pt', mu0pt, 'mu0pt/F')
tree.Branch('mu1pt', mu1pt, 'mu1pt/F')
tree.Branch('mu2pt', mu2pt, 'mu2pt/F')
tree.Branch('mu3pt', mu3pt, 'mu3pt/F')

tree.Branch('mu0eta', mu0eta, 'mu0eta/F')
tree.Branch('mu1eta', mu1eta, 'mu1eta/F')
tree.Branch('mu2eta', mu2eta, 'mu2eta/F')
tree.Branch('mu3eta', mu3eta, 'mu3eta/F')

tree.Branch('mu0phi', mu0phi, 'mu0phi/F')
tree.Branch('mu1phi', mu1phi, 'mu1phi/F')
tree.Branch('mu2phi', mu2phi, 'mu2phi/F')
tree.Branch('mu3phi', mu3phi, 'mu3phi/F')

tree.Branch('mu0trkiso', mu0trkiso, 'mu0trkiso/F')
tree.Branch('mu1trkiso', mu1trkiso, 'mu1trkiso/F')
tree.Branch('mu2trkiso', mu2trkiso, 'mu2trkiso/F')
tree.Branch('mu3trkiso', mu3trkiso, 'mu3trkiso/F')

tree.Branch('mu0dxy', mu0dxy, 'mu0dxy/F')
tree.Branch('mu1dxy', mu1dxy, 'mu1dxy/F')
tree.Branch('mu2dxy', mu2dxy, 'mu2dxy/F')
tree.Branch('mu3dxy', mu3dxy, 'mu3dxy/F')
tree.Branch('mu0edxy', mu0edxy, 'mu0edxy/F')
tree.Branch('mu1edxy', mu1edxy, 'mu1edxy/F')
tree.Branch('mu2edxy', mu2edxy, 'mu2edxy/F')
tree.Branch('mu3edxy', mu3edxy, 'mu3edxy/F')

tree.Branch('PVx', PVx, 'PVx/F')
tree.Branch('PVy', PVy, 'PVy/F')
tree.Branch('PVz', PVz, 'PVz/F')

tree.Branch('scoutingDVx', scoutingDVx, 'scoutingDVx/F')
tree.Branch('scoutingDVy', scoutingDVy, 'scoutingDVy/F')
tree.Branch('scoutingDVz', scoutingDVz, 'scoutingDVz/F')
tree.Branch('scoutingDVlen', scoutingDVlen, 'scoutingDVlen/I')

tree.Branch('reDVx', reDVx, 'reDVx/F')
tree.Branch('reDVy', reDVy, 'reDVy/F')
tree.Branch('reDVz', reDVz, 'reDVz/F')
tree.Branch('reDCA', reDCA, 'reDCA/F')
tree.Branch('remaxDCA', remaxDCA, 'remaxDCA/F')

if args.MC17 or args.MC18:
    tree.Branch('genDVx', genDVx, 'genDVx/F')
    tree.Branch('genDVy', genDVy, 'genDVy/F')
    tree.Branch('genDVz', genDVz, 'genDVz/F')

M4mu = ROOT.TH1F("M4mu", "M4mu", 2000, 0, 200)
muon0pt = ROOT.TH1F("muon0pt", "muon0pt", 1000, 0, 100)  
muon0eta = ROOT.TH1F("muon0eta", "muon0eta", 200, -10, 10)  
muon0phi = ROOT.TH1F("muon0phi", "muon0phi", 200, -10, 10)  
muon0trackiso = ROOT.TH1F("muon0trackiso", "muon0trackiso", 110, -1, 10)  

numberofevents = 0

for i, event in enumerate(tree_mu):

    numberofevents += 1
    
    mu0p = TLorentzVector(0., 0., 0., 0.)
    mu1p = TLorentzVector(0., 0., 0., 0.)
    mu2p = TLorentzVector(0., 0., 0., 0.)
    mu3p = TLorentzVector(0., 0., 0., 0.)

    mu0p_ = TLorentzVector(0., 0., 0., 0.)
    mu1p_ = TLorentzVector(0., 0., 0., 0.)
    mu2p_ = TLorentzVector(0., 0., 0., 0.)
    mu3p_ = TLorentzVector(0., 0., 0., 0.)

    passL1_DoubleMu4_SQ_OS_dR_Max1p2 = 0
    passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 0
    passL1_DoubleMu_15_7 = 0
    passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 0  
    passL1 = 0
    passDST_DoubleMu3_noVtx_CaloScouting_v = 0
    passDST_DoubleMu1_noVtx_CaloScouting_v = 0
    passHLT = 0

    # print "-----In Run" ,tree_mu.Run , "---Lumi", tree_mu.Lumi, "----Event", tree_mu.Event, "-----------"

    for j, _ in enumerate(event.l1bitmap):

        if args.MC17 or args.Data17:

            if event.l1bitmap[j].first == "L1_DoubleMu4_SQ_OS_dR_Max1p2" and event.l1bitmap[j].second == 1:
                passL1_DoubleMu4_SQ_OS_dR_Max1p2 = 1
            else:
                pass

        elif args.MC18 or args.Data18:

            if event.l1bitmap[j].first == "L1_DoubleMu4p5_SQ_OS_dR_Max1p2" and event.l1bitmap[j].second == 1:
                passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 = 1
            else:
                pass

        else:
            pass

        if event.l1bitmap[j].first == "L1_DoubleMu_15_7" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu_15_7 = 1
        else:
            pass

        if event.l1bitmap[j].first == "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4" and event.l1bitmap[j].second == 1:
            passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = 1
        else:
            pass

    # passL1 = passL1_DoubleMu4_SQ_OS_dR_Max1p2 + passL1_DoubleMu_15_7 + passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4
    passL1 = passL1_DoubleMu4_SQ_OS_dR_Max1p2 + passL1_DoubleMu4p5_SQ_OS_dR_Max1p2 + passL1_DoubleMu_15_7 + passL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4

    for j, _ in enumerate(event.hltbitmap):

        if event.hltbitmap[j].first == "DST_DoubleMu3_noVtx_CaloScouting_v" and event.hltbitmap[j].second == 1:
            passDST_DoubleMu3_noVtx_CaloScouting_v = 1
        else:
            pass
            
    for j, _ in enumerate(event.hltbitmap):

        if event.hltbitmap[j].first == "DST_DoubleMu1_noVtx_CaloScouting_v" and event.hltbitmap[j].second == 1:
            passDST_DoubleMu1_noVtx_CaloScouting_v = 1
        else:
            pass

    passHLT = passDST_DoubleMu3_noVtx_CaloScouting_v + passDST_DoubleMu1_noVtx_CaloScouting_v

#    print "In event", i , "passL1?", passL1, "passHLT?", passHLT

    # if passL1 == 0:
    #     continue

    if passHLT == 0:
        continue

    # print "Number of muons", len(event.muon_pt), "Number of vertices", event.dispvertex_num                                                     

    if len(event.muon_pt) != 4:
        continue

    if event.dispvertex_num !=0 and args.storeDV:

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
            else:
                pass

        numberofmuons_invtx = 0
        for index in range(len(event.vertex_index)):
            for subindex in event.vertex_index[index]:
                # print subindex, "muon",index                           
                if subindex == bestdv:
                    numberofmuons_invtx += 1
                    # if numberofmuons_invtx > 4:
                    #     continue
                else:
                    pass

        if numberofmuons_invtx != 4:                         
            continue

        print "The best dv index is ", bestdv, "and number of muons in the dv is", numberofmuons_invtx    

    else:
        pass


    a = np.array(event.muon_pt)
    a = np.argsort(a)[::-1]

    # print a

    if args.MC17 or args.MC18:
        nphi = 0
        for gen in range(len(event.gid)):
            if event.gid[gen] == 6000211:
                nphi+=1
        if nphi != 2:
            continue

    if event.muon_q[a[0]]+event.muon_q[a[1]]+event.muon_q[a[2]]+event.muon_q[a[3]] != 0:
        continue

    mu0p.SetPtEtaPhiM(event.muon_pt[a[0]], event.muon_eta[a[0]], event.muon_phi[a[0]], 0.1056583745)
    mu1p.SetPtEtaPhiM(event.muon_pt[a[1]], event.muon_eta[a[1]], event.muon_phi[a[1]], 0.1056583745)
    mu2p.SetPtEtaPhiM(event.muon_pt[a[2]], event.muon_eta[a[2]], event.muon_phi[a[2]], 0.1056583745)
    mu3p.SetPtEtaPhiM(event.muon_pt[a[3]], event.muon_eta[a[3]], event.muon_phi[a[3]], 0.1056583745)

    muon0_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[a[0]], event.muon_trketa[a[0]], event.muon_trkphi[a[0]], event.muon_trkdsz[a[0]], event.muon_dz[a[0]], event.muon_trklambda[a[0]], event.muon_dxy[a[0]])

    muon1_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[a[1]], event.muon_trketa[a[1]], event.muon_trkphi[a[1]], event.muon_trkdsz[a[1]], event.muon_dz[a[1]], event.muon_trklambda[a[1]], event.muon_dxy[a[1]])

    muon2_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[a[2]], event.muon_trketa[a[2]], event.muon_trkphi[a[2]], event.muon_trkdsz[a[2]], event.muon_dz[a[2]], event.muon_trklambda[a[2]], event.muon_dxy[a[2]])

    muon3_PCA =  ROOT.calculate_track_reference_point(event.muon_trkpt[a[3]], event.muon_trketa[a[3]], event.muon_trkphi[a[3]], event.muon_trkdsz[a[3]], event.muon_dz[a[3]], event.muon_trklambda[a[3]], event.muon_dxy[a[3]])

    reDV01 =  ROOT.recalculate_DV(muon0_PCA[0], muon0_PCA[1], muon0_PCA[2], muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], mu0p.Px(), mu0p.Py(), mu0p.Pz(), mu1p.Px(), mu1p.Py(), mu1p.Pz(), event.muon_q[a[0]], event.muon_q[a[1]])

    reDV02 =  ROOT.recalculate_DV(muon0_PCA[0], muon0_PCA[1], muon0_PCA[2], muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], mu0p.Px(), mu0p.Py(), mu0p.Pz(), mu2p.Px(), mu2p.Py(), mu2p.Pz(), event.muon_q[a[0]], event.muon_q[a[2]])

    reDV03 =  ROOT.recalculate_DV(muon0_PCA[0], muon0_PCA[1], muon0_PCA[2], muon3_PCA[0], muon3_PCA[1], muon3_PCA[2], mu0p.Px(), mu0p.Py(), mu0p.Pz(), mu3p.Px(), mu3p.Py(), mu3p.Pz(), event.muon_q[a[0]], event.muon_q[a[3]])

    reDV12 =  ROOT.recalculate_DV(muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], mu1p.Px(), mu1p.Py(), mu1p.Pz(), mu2p.Px(), mu2p.Py(), mu2p.Pz(), event.muon_q[a[1]], event.muon_q[a[2]])

    reDV13 =  ROOT.recalculate_DV(muon1_PCA[0], muon1_PCA[1], muon1_PCA[2], muon3_PCA[0], muon3_PCA[1], muon3_PCA[2], mu1p.Px(), mu1p.Py(), mu1p.Pz(), mu3p.Px(), mu3p.Py(), mu3p.Pz(), event.muon_q[a[1]], event.muon_q[a[3]])

    reDV23 =  ROOT.recalculate_DV(muon2_PCA[0], muon2_PCA[1], muon2_PCA[2], muon3_PCA[0], muon3_PCA[1], muon3_PCA[2], mu2p.Px(), mu2p.Py(), mu2p.Pz(), mu3p.Px(), mu3p.Py(), mu3p.Pz(), event.muon_q[a[2]], event.muon_q[a[3]])

    reDV = [(reDV01[0] + reDV02[0] + reDV03[0] + reDV12[0] + reDV13[0] + reDV23[0])/6, (reDV01[1] + reDV02[1] + reDV03[1] + reDV12[1] + reDV13[1] + reDV23[1])/6, (reDV01[2] + reDV02[2] + reDV03[2] + reDV12[2] + reDV13[2] + reDV23[2])/6, (reDV01[3] + reDV02[3] + reDV03[3] + reDV12[3] + reDV13[3] + reDV23[3])/6]

    remaxDCA_ = max(reDV01[3], reDV02[3], reDV03[3], reDV12[3], reDV13[3], reDV23[3])

    scoutingDVx_ = 0.
    scoutingDVy_ = 0.
    scoutingDVz_ = 0.
    for dv in range(event.dispvertex_num):
        
        scoutingDVx_ += event.dispvertex_x[dv]
        scoutingDVy_ += event.dispvertex_y[dv]
        scoutingDVz_ += event.dispvertex_z[dv]

        # print "scouting {}: DVx ".format(dv), event.dispvertex_x[dv], " DVy ", event.dispvertex_y[dv], " DVz ", event.dispvertex_z[dv], " chi2/ndf ", event.dispvertex_chi2[dv]/event.dispvertex_ndof[dv]  

    try:
        scoutingDV = [scoutingDVx_/event.dispvertex_num, scoutingDVy_/event.dispvertex_num, scoutingDVz_/event.dispvertex_num]
    except:
        print "weird: No Scouting DV"
        scoutingDV = [9999.,9999.,9999.]
    # print "scouting : DVx ", scoutingDV[0], " DVy ", scoutingDV[1], " DVz ", scoutingDV[2]

    # print "recalculated 01: DVx ", reDV01[0], " DVy ", reDV01[1], " DVz ", reDV01[2], " reDCA ", reDV01[3]
    # print "recalculated 02: DVx ", reDV02[0], " DVy ", reDV02[1], " DVz ", reDV02[2], " reDCA ", reDV02[3]
    # print "recalculated 03: DVx ", reDV03[0], " DVy ", reDV03[1], " DVz ", reDV03[2], " reDCA ", reDV03[3]
    # print "recalculated 12: DVx ", reDV12[0], " DVy ", reDV12[1], " DVz ", reDV12[2], " reDCA ", reDV12[3]
    # print "recalculated 13: DVx ", reDV13[0], " DVy ", reDV13[1], " DVz ", reDV13[2], " reDCA ", reDV13[3]
    # print "recalculated 23: DVx ", reDV23[0], " DVy ", reDV23[1], " DVz ", reDV23[2], " reDCA ", reDV23[3]
    # print "recalculated : DVx ", reDV[0], " DVy ", reDV[1], " DVz ", reDV[2], " reDCA ", reDV[3]

    if args.MC17 or args.MC18:
        for gen in range(len(event.gid)):
            # if event.gid[gen] == 13 and event.mothergid[gen] == 6000211:                       
            #     print "generated: DVx ", event.gendispvertex[gen][0], " DVy ", event.gendispvertex[gen][1], " DVz ", event.gendispvertex[gen][2]                                                
            if event.gid[gen] == 6000211:
                genDV = [event.gendispvertex[gen][0], event.gendispvertex[gen][1], event.gendispvertex[gen][2]]
    #             print "generated: DVx ", event.gendispvertex[gen][0], " DVy ", event.gendispvertex[gen][1], " DVz ", event.gendispvertex[gen][2]

    # print "      "


    M4mu.Fill((mu0p+mu1p+mu2p+mu3p).M())
    muon0pt.Fill(mu0p.Pt())
    muon0eta.Fill(mu0p.Eta())
    muon0phi.Fill(mu0p.Phi())
    muon0trackiso.Fill(event.muon_trackIso[a[0]])

    run[0] = tree_mu.Run
    lumi[0] = tree_mu.Lumi
    event_num[0] = tree_mu.Event
    L1[0] = passL1
    HLT[0] = passHLT

    fourmuonmass[0] = (mu0p+mu1p+mu2p+mu3p).M()

    mu0q[0] = event.muon_q[a[0]]
    mu1q[0] = event.muon_q[a[1]]
    mu2q[0] = event.muon_q[a[2]]
    mu3q[0] = event.muon_q[a[3]]

    mu0pt[0] = mu0p.Pt()
    mu1pt[0] = mu1p.Pt()
    mu2pt[0] = mu2p.Pt()
    mu3pt[0] = mu3p.Pt()

    mu0eta[0] = mu0p.Eta()
    mu1eta[0] = mu1p.Eta()
    mu2eta[0] = mu2p.Eta()
    mu3eta[0] = mu3p.Eta()
    
    mu0phi[0] = mu0p.Phi()
    mu1phi[0] = mu1p.Phi()
    mu2phi[0] = mu2p.Phi()
    mu3phi[0] = mu3p.Phi()

    mu0trkiso[0] = event.muon_trackIso[a[0]]
    mu1trkiso[0] = event.muon_trackIso[a[1]]
    mu2trkiso[0] = event.muon_trackIso[a[2]]
    mu3trkiso[0] = event.muon_trackIso[a[3]]

    mu0dxy[0] = event.muon_dxy[a[0]]
    mu1dxy[0] = event.muon_dxy[a[1]]
    mu2dxy[0] = event.muon_dxy[a[2]]
    mu3dxy[0] = event.muon_dxy[a[3]]
    mu0edxy[0] = event.muon_edxy[a[0]]
    mu1edxy[0] = event.muon_edxy[a[1]]
    mu2edxy[0] = event.muon_edxy[a[2]]
    mu3edxy[0] = event.muon_edxy[a[3]]

    PVx[0] = event.privertex_x[0]
    PVy[0] = event.privertex_y[0]
    PVz[0] = event.privertex_z[0]

    scoutingDVx[0] = scoutingDV[0]
    scoutingDVy[0] = scoutingDV[1]
    scoutingDVz[0] = scoutingDV[2]
    scoutingDVlen[0] = event.dispvertex_num

    reDVx[0] = reDV[0]
    reDVy[0] = reDV[1]
    reDVz[0] = reDV[2]
    reDCA[0] = reDV[3]
    remaxDCA[0] = remaxDCA_ 
    
    if args.MC17 or args.MC18:
        genDVx[0] = genDV[0]
        genDVy[0] = genDV[1]
        genDVz[0] = genDV[2]

    b = np.array(event.muon_q)
    b = np.argsort(b)[::-1]
    
    # print(event.muon_q[b[0]], event.muon_q[b[1]], event.muon_q[b[2]], event.muon_q[b[3]])
    mu0p_.SetPtEtaPhiM(event.muon_pt[b[0]], event.muon_eta[b[0]], event.muon_phi[b[0]], 0.1056583745)
    mu1p_.SetPtEtaPhiM(event.muon_pt[b[1]], event.muon_eta[b[1]], event.muon_phi[b[1]], 0.1056583745)
    mu2p_.SetPtEtaPhiM(event.muon_pt[b[2]], event.muon_eta[b[2]], event.muon_phi[b[2]], 0.1056583745)
    mu3p_.SetPtEtaPhiM(event.muon_pt[b[3]], event.muon_eta[b[3]], event.muon_phi[b[3]], 0.1056583745)


    # print "M02", (mu0p_ + mu2p_).M(), "M13", (mu1p_ +mu3p_).M(), "M03", (mu0p_ +mu3p_).M(), "M12", (mu1p_ +mu2p_).M()
        
    pair1_mass_ = -1.
    pair2_mass_ = -1.
    pair1_pt_ = -1.
    pair2_pt_ = -1.
    pass_masspair_ = 0

    if 0.95 < (mu0p_ + mu2p_).M() /  (mu1p_ +mu3p_).M() < 1.05 and not (0.95 < (mu0p_ + mu3p_).M() /  (mu1p_ +mu2p_).M() < 1.05):
        pass_masspair_ += 1
        if (mu0p_ + mu2p_).Pt() > (mu1p_ +mu3p_).Pt():
            pair1_mass_ = (mu0p_ + mu2p_).M()
            pair2_mass_ = (mu1p_ +mu3p_).M()
            pair1_pt_ = (mu0p_ + mu2p_).Pt()
            pair2_pt_ = (mu1p_ +mu3p_).Pt()
        else:
            pair1_mass_ = (mu1p_ +mu3p_).M()
            pair2_mass_ = (mu0p_ + mu2p_).M()
            pair1_pt_ = (mu1p_ + mu3p_).Pt()
            pair2_pt_ = (mu0p_ + mu2p_).Pt()
    
    elif 0.95 < (mu0p_ + mu3p_).M() /  (mu1p_ +mu2p_).M() < 1.05 and not (0.95 < (mu0p_ + mu2p_).M() /  (mu1p_ +mu3p_).M() < 1.05):
        pass_masspair_ += 1
        if (mu0p_ + mu3p_).Pt() > (mu1p_ +mu2p_).Pt():
            pair1_mass_ = (mu0p_ + mu3p_).M()
            pair2_mass_ = (mu1p_ +mu2p_).M()
            pair1_pt_ = (mu0p_ + mu3p_).Pt()
            pair2_pt_ = (mu1p_ +mu2p_).Pt()
        else:
            pair1_mass_ = (mu1p_ +mu2p_).M()
            pair2_mass_ = (mu0p_ + mu3p_).M()
            pair1_pt_ = (mu1p_ + mu2p_).Pt()
            pair2_pt_ = (mu0p_ +mu3p_).Pt()
    elif 0.95 < (mu0p_ + mu3p_).M() /  (mu1p_ +mu2p_).M() < 1.05 and (0.95 < (mu0p_ + mu2p_).M() /  (mu1p_ +mu3p_).M() < 1.05):
        pass_masspair_ += 1
        print "Both pair of muon pairs in this event ", i ," has mass within 5% of each other"
        # if (mu0p_ + mu3p_).Pt() > (mu1p_ +mu2p_).Pt():
        #     pair1_mass_ = (mu0p_ + mu3p_).M()
        #     pair2_mass_ = (mu1p_ +mu2p_).M()
        #     pair1_pt_ = (mu0p_ + mu3p_).Pt()
        #     pair2_pt_ = (mu1p_ +mu2p_).Pt()
        # else:
        #     pair1_mass_ = (mu1p_ +mu2p_).M()
        #     pair2_mass_ = (mu0p_ + mu3p_).M()
        #     pair1_pt_ = (mu1p_ + mu2p_).Pt()
        #     pair2_pt_ = (mu0p_ +mu3p_).Pt()

        if (mu0p_ + mu2p_).Pt() > (mu1p_ +mu3p_).Pt():
            pair1_mass_ = (mu0p_ + mu2p_).M()
            pair2_mass_ = (mu1p_ +mu3p_).M()
            pair1_pt_ = (mu0p_ + mu2p_).Pt()
            pair2_pt_ = (mu1p_ +mu3p_).Pt()
        else:
            pair1_mass_ = (mu1p_ +mu3p_).M()
            pair2_mass_ = (mu0p_ + mu2p_).M()
            pair1_pt_ = (mu1p_ + mu3p_).Pt()
            pair2_pt_ = (mu0p_ + mu2p_).Pt()

    else:
        # print "No pair of muon pairs in this event ", i , "has mass within 5% of each other"
        pass
    
    # print "pair1 mass", pair1_mass_, "pair2 mass", pair2_mass_, "pair1 pt", pair1_pt_, "pair2 pt", pair2_pt_, "pass masspaircut", pass_masspair_
    # print " "

    pair1_mass[0] = pair1_mass_ 
    pair2_mass[0] = pair2_mass_
    pair1_pt[0] = pair1_pt_
    pair2_pt[0] = pair2_pt_
    pass_masspair[0] = pass_masspair_


    jpsipair_mass_ = -1
    phipair_mass_ = -1
    jpsipair_pt_ = -1
    phipair_pt_ = -1
    pass_jpsiphi_ = 0

    if (2.9 < (mu0p_ + mu2p_).M() < 3.1) and (0.9 < (mu1p_ +mu3p_).M() < 1.1):
        pass_jpsiphi_ += 1
        jpsipair_mass_ = (mu0p_ + mu2p_).M()
        phipair_mass_ = (mu1p_ +mu3p_).M()
        jpsipair_pt_ = (mu0p_ + mu2p_).Pt()
        phipair_pt_ = (mu1p_ +mu3p_).Pt()
    elif (0.9 < (mu0p_ + mu2p_).M() < 1.1) and (2.9 < (mu1p_ +mu3p_).M() < 3.1):
        pass_jpsiphi_ += 1
        jpsipair_mass_ = (mu1p_ +mu3p_).M()
        phipair_mass_ = (mu0p_ + mu2p_).M()
        jpsipair_pt_ = (mu1p_ +mu3p_).Pt()
        phipair_pt_ = (mu0p_ + mu2p_).Pt()
    elif (2.9 < (mu0p_ + mu3p_).M() < 3.1) and (0.9 < (mu1p_ +mu2p_).M() < 1.1):
        pass_jpsiphi_ += 1
        jpsipair_mass_ = (mu0p_ +mu3p_).M()
        phipair_mass_ = (mu1p_ + mu2p_).M()
        jpsipair_pt_ = (mu0p_ +mu3p_).Pt()
        phipair_pt_ = (mu1p_ + mu2p_).Pt()
    elif (0.9 < (mu0p_ + mu3p_).M() < 1.1) and (2.9 < (mu1p_ +mu2p_).M() < 3.1):
        pass_jpsiphi_ += 1
        jpsipair_mass_ = (mu1p_ +mu2p_).M()
        phipair_mass_ = (mu0p_ + mu3p_).M()
        jpsipair_pt_ = (mu1p_ +mu2p_).Pt()
        phipair_pt_ = (mu0p_ + mu3p_).Pt()
    else:
        # print "No pair in this event having mass in both jpsi and phi mass range found"
        pass
    
    jpsipair_mass[0] = jpsipair_mass_
    phipair_mass[0] = phipair_mass_
    jpsipair_pt[0] = jpsipair_pt_
    phipair_pt[0] = phipair_pt_
    pass_jpsiphi[0] = pass_jpsiphi_


    pass_ZZ_ = 0

    if (88 < (mu0p_ + mu2p_).M() < 92) and (88 < (mu1p_ +mu3p_).M() < 92):
        pass_ZZ_ += 1
    elif (88 < (mu0p_ + mu3p_).M() < 92) and (88 < (mu1p_ +mu2p_).M() < 92):
        pass_ZZ_ += 1
    else:
        # print "No pairs in this event having both masses in Z mass range found"
        pass
    
    pass_ZZ[0] = pass_ZZ_


    if event.muon_q[a[0]]*event.muon_q[a[1]] == -1:
        M01[0] = (mu0p + mu1p).M()
        M23[0] = (mu2p +mu3p).M()
        pt01[0] = (mu0p + mu1p).Pt()
        pt23[0] = (mu2p +mu3p).Pt()
    else:
        M01[0] = -1
        M23[0] = -1
        pt01[0] = -1
        pt23[0] = -1

    if event.muon_q[a[0]]*event.muon_q[a[2]] == -1:
        M02[0] = (mu0p +mu2p).M()
        M13[0] = (mu1p +mu3p).M()
        pt02[0] = (mu0p +mu2p).Pt()
        pt13[0] = (mu1p +mu3p).Pt()
    else:
        M02[0] = -1
        M13[0] = -1
        pt02[0] = -1
        pt13[0] = -1

    if event.muon_q[a[0]]*event.muon_q[a[3]] == -1:
        M03[0] = (mu0p +mu3p).M()
        M12[0] = (mu1p +mu2p).M()
        pt03[0] = (mu0p +mu3p).Pt()
        pt12[0] = (mu1p +mu2p).Pt()
    else:
        M03[0] = -1
        M12[0] = -1
        pt03[0] = -1
        pt12[0] = -1

    tree.Fill()


    # print "-----In Run" ,tree_mu.Run , "---Lumi", tree_mu.Lumi, "----Event", tree_mu.Event, "-----------"
    # print "Number of muons", len(event.muon_pt), "Number of vertices", event.dispvertex_num                                                     

    # print "The charges of muons are", event.muon_q[a[0]], event.muon_q[a[1]], event.muon_q[a[2]], event.muon_q[a[3]]
    # print "The pt of muons are", mu0p.Pt(), mu1p.Pt(), mu2p.Pt(), mu3p.Pt()
    # print "The eta of muons are", mu0p.Eta(), mu1p.Eta(), mu2p.Eta(), mu3p.Eta()
    # print "The phi of muons are", mu0p.Phi(), mu1p.Phi(), mu2p.Phi(), mu3p.Phi()
    # print "The trackiso of the muons are", event.muon_trackIso[a[0]], event.muon_trackIso[a[1]], event.muon_trackIso[a[2]], event.muon_trackIso[a[3]] 

print "number of events in the file", numberofevents

# M4mu.Write()
# muon0pt.Write()
# muon0eta.Write()
# muon0phi.Write()
# muon0trackiso.Write()

outfile.Write()
outfile.Close()


