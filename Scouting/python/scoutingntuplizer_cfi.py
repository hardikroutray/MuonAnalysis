import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')

options.register('is2017data',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Determines whether sample is 2017 data or not.")
options.register('is2018data',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Determines whether sample is 2018 data or not.")
options.register('is2017MC',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Determines whether sample is 2017 MC or not.")
options.register('is2018MC',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Determines whether sample is 2018 MC or not.")
options.register('condInputFolder',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "folder containing list of files for condor")
options.register('condInputNum',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "the file to run on")
options.register('condInputtxt',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "txt file containing list of root files for condor")
options.parseArguments()


scoutingntuplizer = cms.EDAnalyzer(
   'ScoutingNtuplizer',
   hltProcess=cms.string("HLT"),
   bits          = cms.InputTag("TriggerResults", "", "HLT"),
   bsCollection = cms.InputTag("offlineBeamSpot"),
   #triggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
   AlgInputTag             = cms.InputTag("gtStage2Digis"),
   l1tAlgBlkInputTag       = cms.InputTag("gtStage2Digis"),
   l1tExtBlkInputTag       = cms.InputTag("gtStage2Digis"),
   jet_collection      = cms.InputTag('hltScoutingCaloPacker'),
   rho                 = cms.InputTag('hltScoutingCaloPacker:rho'),
   muon_collection     = cms.InputTag('hltScoutingMuonPackerCalo'),
   primary_vertex_collection = cms.InputTag('hltScoutingPrimaryVertexPackerCaloMuon:primaryVtx'),
   displaced_vertex_collection   = cms.InputTag('hltScoutingMuonPackerCalo:displacedVtx'),
   MET_pt              = cms.InputTag('hltScoutingCaloPacker:caloMetPt'),
   MET_phi             = cms.InputTag('hltScoutingCaloPacker:caloMetPhi'),
   output_file_name    = cms.string('scouting_ntuple.root'),
   mu_min              = cms.int32(2),
   hltseeds            = cms.vstring("DST_DoubleMu3_noVtx_CaloScouting_v",
                                     "DST_DoubleMu3_noVtx_CaloScouting_Monitoring_v",
                                     "DST_DoubleMu1_noVtx_CaloScouting_v",
                                     "DST_HT250_CaloScouting_v",
                                     "HLT_PFJet40_v",
                                     "DST_HT250_CaloScouting_v",
                                     "DST_HT410_PFScouting_v",
                                     "DST_ZeroBias_CaloScouting_PFScouting_v",
                                     "DST_ZeroBias_v",
                                     "DST_CaloJet40_CaloBTagScouting_v",
                                     "DST_DoubleMu3_noVtx_Mass10_PFScouting_v",
                                     "DST_CaloJet40_CaloScouting_PFScouting_v"
                                      ),
   l1seeds             = cms.vstring("L1_DoubleMu_12_5",
                                     "L1_DoubleMu_12_8",
                                     "L1_DoubleMu_13_6",
                                     "L1_DoubleMu_15_5",
                                     "L1_DoubleMu_15_7",
                                     "L1_DoubleMu18er2p1",
                                     "L1_DoubleMu22er2p1",
                                     "L1_TripleMu_4_4_4",
                                     "L1_TripleMu_5_0_0",
                                     "L1_TripleMu_5_3_3",
                                     "L1_TripleMu_5_5_3",
                                     "L1_QuadMu0",
                                     "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
                                     "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18",
                                     "L1_DoubleMu4_SQ_OS_dR_Max1p2",
                                     "L1_DoubleMu5_SQ_OS_Mass7to18",
                                     "L1_DoubleMu_20_2_SQ_Mass_Max20",
                                     "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
                                     "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
                                     "L1_DoubleMu6_SQ_OS",
                                     "L1_DoubleMu0er1p5_SQ_dR_Max1p4",
                                     "L1_DoubleMu0er2_SQ_dR_Max1p4",
                                     "L1_DoubleMu0_SQ"),
    is2017data                 = cms.bool(options.is2017data),
    is2018data                 = cms.bool(options.is2018data),
    is2017MC                 = cms.bool(options.is2017MC),
    is2018MC                 = cms.bool(options.is2018MC)

)
