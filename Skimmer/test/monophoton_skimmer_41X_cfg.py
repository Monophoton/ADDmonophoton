process = cms.Process("ADDSkim")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#======

## global tag for data

process.GlobalTag.globaltag = cms.string('GR_P_V14::All')


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/853/2E97032B-5453-E011-A481-0030487A3C9A.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/827/5A9A2B84-8A53-E011-BC97-0030487CD7E0.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/819/EA548277-8753-E011-8EA5-001D09F2B30B.root',
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

#closes files after code is done running on that file
#process.options = cms.untracked.PSet(
#	fileMode = cms.untracked.string('NOMERGE')
#)
#process.source.skipEvents = cms.untracked.uint32(400)
# let it run
process.monoSkim = cms.EDFilter("MonoPhotonSkimmer",
  phoTag = cms.InputTag("photons::RECO"),
  selectEE = cms.bool(True),
  selectTrack = cms.bool(False),
  ecalisoOffsetEB = cms.double(5000.),
  ecalisoSlopeEB = cms.double(0.15),
  hcalisoOffsetEB = cms.double(7000.0),
  hcalisoSlopeEB = cms.double(0.),
  hadoveremEB = cms.double(0.05),
  minPhoEtEB = cms.double(30.),
  scHighEtThreshEB = cms.double(100.),
  ecalisoOffsetEE = cms.double(5000.),
  ecalisoSlopeEE = cms.double(0.15),
  hcalisoOffsetEE = cms.double(10000.),
  hcalisoSlopeEE = cms.double(0.),
  hadoveremEE = cms.double(0.05),
  minPhoEtEE = cms.double(30.),
  scHighEtThreshEE = cms.double(100.),
 )

process.p = cms.Path(process.monoSkim)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.options.wantSummary = True

# define output event selection to be that which satisfies 'p'
process.patEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)
# setup event content: drop everything before 
process.patEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *')
)


# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
    process.patEventSelection,
    process.patEventContent,
    #verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string('monopho.root')
)

process.outpath = cms.EndPath(process.out)

#process.schedule=cms.Schedule(process.p)
try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")



