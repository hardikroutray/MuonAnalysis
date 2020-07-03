# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/BPH-RunIIFall17GS-00016-fragment.py --fileout file:BPH-RunIIFall17GS-00016.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 93X_mc2017_realistic_v3 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN,SIM --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename BPH-RunIIFall17GS-00016_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring
import FWCore.ParameterSet.Config as cms
import os
import sys

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2017Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring('file:./testlhe.root'),
    fileNames = cms.untracked.vstring('file:./LHEROOT_output/LHEROOT_6000211_100events_{0}.root'.format(sys.argv[2])),
    inputCommands = cms.untracked.vstring('keep *',
        'drop LHEXMLStringProduct_*_*_*'),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
)

process.options = cms.untracked.PSet(

)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/BPH-RunIIFall17GS-00016-fragment.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    #fileName = cms.untracked.string('file:./test.root'),
    fileName = cms.untracked.string('file:./GENSIM_output/gg_llH_mumu_GENSIM_{0}.root'.format(sys.argv[2])),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_mc2017_realistic_v3', '')



'''
process.decayfilter = cms.EDFilter("PythiaDauVFilter",
    DaughterIDs = cms.untracked.vint32(13, -13),
    MaxEta = cms.untracked.vdouble(9999.0, 9999.0),
    MinEta = cms.untracked.vdouble(-9999.0, -9999.0),
    MinPt = cms.untracked.vdouble(-99.0, -99.0),
    NumberDaughters = cms.untracked.int32(2),
    ParticleID = cms.untracked.int32(511),
    verbose = cms.untracked.int32(1)
)


process.bfilter = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(9999.0),
    MinEta = cms.untracked.double(-9999.0),
    ParticleID = cms.untracked.int32(511)
)

'''

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2',
            'Main:timesAllowErrors = 10000',
            'Check:epTolErr = 0.01',
            'Beams:setProductionScalesFromLHEF = off',
            'SLHA:keepSM = on',
            'SLHA:minMassSM = 1000.',
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tau0Max = 1000',
            'ParticleDecays:allowPhotonRadiation = on'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14',
            'Tune:ee 7',
            'MultipartonInteractions:pT0Ref=2.4024',
            'MultipartonInteractions:ecmPow=0.25208',
            'MultipartonInteractions:expPow=1.6'),
        parameterSets = cms.vstring('pythia8CommonSettings',
            'pythia8CUEP8M1Settings','processParameters'),
        processParameters = cms.vstring(
                                        'JetMatching:setMad = off',
					'JetMatching:scheme = 1',
					'JetMatching:merge = on',
					'JetMatching:jetAlgorithm = 2',
					'JetMatching:etaJetMax = 999.',
					'JetMatching:coneRadius = 1.',
					'JetMatching:slowJetPower = 1',
					'JetMatching:qCut = 3.', #this is the actual merging scale
					'JetMatching:doFxFx = off',
					'JetMatching:qCutME = 1.',#this must match the ptj cut in the lhe generation step
					'JetMatching:nQmatch = 4', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
					'JetMatching:nJetMax = 1', #number of partons in born matrix element for highest
					'JetMatching:doShowerKt = off',  #off for MLM matching, turn on for shower-kT matching


			                'LesHouches:setLifetime = 2',
			                '6000211:new',
                                        '6000211:name = higgs_bsm',
                                        '6000211:spinType = 0',
                                        '6000211:chargeType = 0',
                                        '6000211:colType = 0',
                                        '6000211:m0 = 1.7',
                                        '6000211:tau0 = 5',
                                        #'46:7:bRatio  = 1.0',                                                                      
                                        #'46:7:meMode  = 100'                                                                       
                                        '6000211:addChannel = 1 1.0 0 13 -13'
                                        #'25:addChannel = 1 0.5 0 111 111 111'                                                      
                                        )
    )
)




# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
# filter all path with the production filter sequence
#for path in process.paths:
#	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 
for path in process.paths:
        getattr(process,path)._seq = process.generator * getattr(process,path)._seq


# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
