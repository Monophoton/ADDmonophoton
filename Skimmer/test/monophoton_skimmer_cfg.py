process = cms.Process("ADDSkim")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#spike cleaning
process.load("EGamma.EGammaSkims.promptRecoTrackCorrections_cff")
process.load("EGamma.EGammaSkims.cleanReRecoSequence_cff")
#======

process.ecalCleanClustering = cms.Sequence(process.cleanedEcalClusters*process.cleanedEgammaSkimReco)
## global tag for data

process.GlobalTag.globaltag = cms.string('GR_R_36X_V12B::All')


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_1.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_2.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_3x.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_4.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_5.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_6.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_7.root',
'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_8.root'
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

#closes files after code is done running on that file
process.options = cms.untracked.PSet(
	fileMode = cms.untracked.string('NOMERGE')
)
#process.source.skipEvents = cms.untracked.uint32(400)
# let it run
process.monoSkimC = cms.EDFilter("MonoPhotonSkimmer",
  phoTag = cms.InputTag("photons"),
  selectEE = cms.bool(True),
  selectTrack = cms.bool(True),
  ecalisoOffsetEB = cms.double(5.),
  ecalisoSlopeEB = cms.double(0.15),
  hcalisoOffsetEB = cms.double(7.0),
  hcalisoSlopeEB = cms.double(0.),
  hadoveremEB = cms.double(0.5),
  minPhoEtEB = cms.double(25.),
  scHighEtThreshEB = cms.double(100.),
  ecalisoOffsetEE = cms.double(5.),
  ecalisoSlopeEE = cms.double(0.15),
  hcalisoOffsetEE = cms.double(10.),
  hcalisoSlopeEE = cms.double(0.),
  hadoveremEE = cms.double(0.5),
  minPhoEtEE = cms.double(25.),
  scHighEtThreshEE = cms.double(100.),
 )
print
print "============== Warning =============="
print "technical trigger filter:    ENABLED"
print "physics declare bit filter:  ENABLED"
print "primary vertex filter:       ENABLED"

process.p = cms.Path(
   process.ecalCleanClustering*monoSkimC

    )


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
#process.options.wantSummary = True


process.schedule=cms.Schedule(process.p)
try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")



