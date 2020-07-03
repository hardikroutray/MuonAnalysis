
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: lhe --mc --no_exec -n 10 --filein file:unweighted_events.lhe --fileout file:LHEROOT.root --python_filename LHEROOT_cfg.py --eventcontent LHE --datatier LHE --step NONE --conditions MCRUN2_71_V1::All --customise Configuration/DataProcessing/Utils.addMonitoring
import FWCore.ParameterSet.Config as cms
import os
import sys

INPUTs = int("{0}00".format(sys.argv[2]))

print INPUTs 

process = cms.Process('LHE')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('file:./unweighted_events_6000211_100000events_mass1p7.lhe'),
#    fileNames = cms.untracked.vstring('file:./unweighted_events_100000_mass125.lhe'),
    skipEvents=cms.untracked.uint32(INPUTs)
                            
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('lhe nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.LHEoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.LHEEventContent.outputCommands,
#    fileName = cms.untracked.string('file:./test_lheroot.root'),
    fileName = cms.untracked.string('file:./LHEROOT_output/LHEROOT_6000211_100events_{0}.root'.format(sys.argv[2])),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('LHE')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_mc2017_realistic_v3', '')

# Path and EndPath definitions
process.LHEoutput_step = cms.EndPath(process.LHEoutput)

# Schedule definition
process.schedule = cms.Schedule(process.LHEoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
