import FWCore.ParameterSet.Config as cms

theoryminituplizer = cms.EDAnalyzer(
   'TheoryMiniTuplizer',
   jet_collection      = cms.InputTag('hltScoutingCaloPacker'),
   rho                 = cms.InputTag('hltScoutingCaloPacker:rho'),
   muon_collection     = cms.InputTag('hltScoutingMuonPackerCalo'),
   MET_pt              = cms.InputTag('hltScoutingCaloPacker:caloMetPt'),
   MET_phi             = cms.InputTag('hltScoutingCaloPacker:caloMetPhi'),
   output_file_name    = cms.string('theory_minituple.root')
)
