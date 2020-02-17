import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TheoryMiniTuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")


options = VarParsing.VarParsing('analysis')
options.outputFile = 'theory_minituple.root'
options.inputFiles = 'root://cmsxrootd.fnal.gov//store/data/Run2017E/ScoutingCaloMuon/RAW/v1/000/303/832/00000/DAF942AC-CAA1-E711-BD1B-02163E01A69F.root'
options.maxEvents = -1
options.register('reportEvery',
                 1000000,
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

process.load('MuonAnalysis.Scouting.theoryminituplizer_cfi')


process.theoryminituplizer.output_file_name = cms.string(options.outputFile)


process.p = cms.Path(process.theoryminituplizer)
