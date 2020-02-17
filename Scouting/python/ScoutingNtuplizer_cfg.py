import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList



process = cms.Process("ScoutingNtuplizer")



import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
import EventFilter.L1TRawToDigi.gtStage2Digis_cfi
import EventFilter.L1TRawToDigi.caloStage2Digis_cfi
import EventFilter.L1TRawToDigi.gmtStage2Digis_cfi

# standard unpacking sequence                                                                                                       
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

# L1 Trigger sequences                                                                                                              

# l1tMonitor and l1tMonitorEndPathSeq                                                                                               
process.load("DQM.L1TMonitor.L1TMonitor_cff")

# L1 trigger synchronization module - it uses also HltHighLevel filter                                                              
process.load("DQM.L1TMonitor.L1TSync_cff")


# l1tMonitorClient and l1tMonitorClientEndPathSeq                                                                                   
process.load("DQM.L1TMonitorClient.L1TMonitorClient_cff")

process.gtDigis.DaqGtInputTag = cms.InputTag("hltFEDSelectorL1")

process.gtStage2Digis = EventFilter.L1TRawToDigi.gtStage2Digis_cfi.gtStage2Digis.clone(InputLabel = cms.InputTag("hltFEDSelectorL1"))
process.caloStage2Digis = EventFilter.L1TRawToDigi.caloStage2Digis_cfi.caloStage2Digis.clone(InputLabel = cms.InputTag("hltFEDSelectorL1"))
process.gmtStage2Digis = EventFilter.L1TRawToDigi.gmtStage2Digis_cfi.gmtStage2Digis.clone(InputLabel = cms.InputTag("hltFEDSelectorL1"))

process.l1DigiSeq = cms.Sequence(process.gtDigis*process.gtStage2Digis*process.caloStage2Digis*process.gmtStage2Digis)



process.load("FWCore.MessageService.MessageLogger_cfi")


options = VarParsing.VarParsing('analysis')
options.outputFile = 'scouting_ntuple.root'
#options.inputFiles = 'root://xrootd-cms.infn.it///store/data/Run2017D/ScoutingCaloMuon/RAW/v1/000/302/031/00000/301F9B4E-648D-E711-8480-02163E012748.root' 
#'root://cmsxrootd.fnal.gov//store/data/Run2017E/ScoutingCaloMuon/RAW/v1/000/303/832/00000/DAF942AC-CAA1-E711-BD1B-02163E01A69F.root'
options.inputFiles = '/store/data/Run2017D/ScoutingCaloMuon/RAW/v1/000/302/033/00000/9C0FCC26-8B8D-E711-8A1F-02163E01273D.root'
#'/store/data/Run2017C/ScoutingCaloMuon/RAW/v1/000/302/026/00000/F815A335-528D-E711-AEF7-02163E014129.root'
#'/store/data/Run2017D/ScoutingCaloMuon/RAW/v1/000/302/663/00000/DCB0847B-1498-E711-B26A-02163E01A3FB.root' 
#'/store/data/Run2017E/ScoutingCaloMuon/RAW/v1/000/304/204/00000/447CB8DC-0AA7-E711-A550-02163E01A362.root' 
#'root://cmsxrootd.fnal.gov//store/data/Run2017F/ScoutingCaloMuon/RAW/v1/000/305/377/00000/180B769B-0CB7-E711-9825-02163E01A4CE.root'

options.maxEvents = -1
#options.maxEvents = 100
options.register('reportEvery',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process before reporting progress.")
options.parseArguments()



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(
    options.reportEvery)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)


process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.autoCond import autoCond
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_HLT_v7', '')

#process.source.lumisToProcess = LumiList.LumiList(filename = './JSON.txt').getVLuminosityBlockRange()
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.load('MuonAnalysis.Scouting.scoutingntuplizer_cfi')

#process.scoutingntuplizer = cms.EDAnalyzer(
#   'ScoutingNtuplizer',
#   jet_collection      = cms.InputTag('hltScoutingCaloPacker'),
#   rho                 = cms.InputTag('hltScoutingCaloPacker:rho'),
#   muon_collection     = cms.InputTag('hltScoutingMuonPackerCalo'),
#   vertex_collection   = cms.InputTag('hltScoutingMuonPackerCalo:displacedVtx'),
#   MET_pt              = cms.InputTag('hltScoutingCaloPacker:caloMetPt'),
#   MET_phi             = cms.InputTag('hltScoutingCaloPacker:caloMetPhi'),
#   output_file_name    = cms.string('scouting_ntuple.root'),
#   mu_min              = cms.int32(2)
#)



process.scoutingntuplizer.output_file_name = cms.string(options.outputFile)


process.p = cms.Path(process.l1DigiSeq*process.offlineBeamSpot*process.scoutingntuplizer)