# This is an example PAT configuration showing the usage of PAT on minbias data

process = cms.Process("PAT")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#spike correction
process.load("EGamma.EGammaSkims.promptRecoTrackCorrections_cff")
process.load("EGamma.EGammaSkims.cleanReRecoSequence_cff")
#======
### prevent bias in pixel vertex
#process.pixelVertices.useBeamConstraint = False

### pixelTracks
#---- new parameters ----
#process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = 4.06
#process.pixelTracks.RegionFactoryPSet.RegionPSet.originHalfLength = cms.double(40.6)

### 0th step of iterative tracking
#---- new parameters ----
#process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(
#4.06)
#process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 40.6

### 2nd step of iterative tracking
#---- new parameters ----
process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)
process.secTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 44.7

## Primary Vertex
process.offlinePrimaryVerticesWithBS.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVerticesWithBS.TkClusParameters.TkGapClusParameters.zSeparation = 1
process.offlinePrimaryVertices.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVertices.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVertices.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVertices.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVertices.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVertices.TkClusParameters.TkGapClusParameters.zSeparation= 1

#=========
process.ecalCleanClustering = cms.Sequence(process.cleanedEcalClusters*process.cleanedEgammaSkimReco)
## global tag for data
#process.GlobalTag.globaltag = cms.String('Start3x_v27::All')
process.GlobalTag.globaltag = cms.string('GR_R_36X_V10::All')
#process.GlobalTag.globaltag = cms.string('MC_3XY_V27:All')
# turn off MC matching for the process
#removeMCMatching(process, ['All'])
from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])
addPfMET(process, 'PF')
addTcMET(process,"TC")


# get the 900 GeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332")

# run ak5 gen jets
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#run33xOnReRecoMC( process, "ak5GenJets")

# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
#switchOnTrigger( process )

# Select jets
process.selectedPatJets.cut = cms.string('pt > 10 & abs(eta) < 3.0')
process.selectedPatJetsAK5PF.cut = cms.string('pt > 8 & abs(eta) < 3.0')

#process.hltHighLevel = cms.EDFilter("HLTHighLevel",
#                            TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#                            HLTPaths = cms.vstring("HLT_Photon20_Cleaned_L1R"),
#                            eventSetupPathsKey = cms.string(''),
#                            andOr = cms.bool(False),
#                            throw = cms.bool(True)
#                            )

# newly added
#import FWCore.ParameterSet.Config as cms
#
#CosMuCountFilter = cms.EDFilter("CandViewCountFilter",
#    src = cms.InputTag("muonsFromCosmics"),
#    minNumber = cms.uint32(1)
#)
#
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

print
print "============== Warning =============="
print "technical trigger filter:    ENABLED"
print "physics declare bit filter:  ENABLED"
print "primary vertex filter:       ENABLED"
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.demo = cms.EDAnalyzer('Analyzer',
                              electronTag      = cms.untracked.InputTag("selectedPatElectrons"),
                              tauTag           = cms.untracked.InputTag("selectedPatTaus"),
                              muonTag          = cms.untracked.InputTag("selectedPatMuons"),
                              jetTag           = cms.untracked.InputTag("selectedPatJets"),
                              photonTag        = cms.untracked.InputTag("selectedPatPhotons"),
                              rechitBTag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEE"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag           = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag           = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              outFile          = cms.untracked.string("/uscms/home/miceli/scratch/data/2010JulyData/monophotonNtuples/MuPR4_fromSkim_test2final.root"),
                              runphotons       = cms.untracked.bool(True),
                              runHErechit      = cms.untracked.bool(True),
                              runrechit        = cms.untracked.bool(True),
                              runmet           = cms.untracked.bool(True),
                              rungenmet        = cms.untracked.bool(False),
                              runPFmet         = cms.untracked.bool(True),
                              runTCmet         = cms.untracked.bool(True),
                              runjets          = cms.untracked.bool(True),
                              runelectrons     = cms.untracked.bool(True),
                              runtaus          = cms.untracked.bool(True),
                              runmuons         = cms.untracked.bool(True),
                              runcosmicmuons   = cms.untracked.bool(True),
                              rungenParticleCandidates = cms.untracked.bool(False),
                              runHLT           = cms.untracked.bool(True),
                              runL1            = cms.untracked.bool(True),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True)
                              )
#process.out.fileName = "DROPPED"

process.p = cms.Path(
   # process.hltHighLevel*
   process.ecalCleanClustering*
 #  process.hltHighLevel*
   process.patDefaultSequence*
   process.demo
    )

# rename output file
#process.out.fileName = cms.untracked.string('/tmp/sandhya/pat.root')
#process.out.fileName = cms.untracked.string('5K_patb1.root')
# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.options.wantSummary = True

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
#from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
#process.out.outputCommands = patEventContentNoCleaning
#process.out.outputCommands += patExtraAodEventContent
#process.out.outputCommands += patTriggerEventContent
#process.out.outputCommands += ['keep recoPFCandidates_particleFlow_*_*']
#process.out.outputCommands += ['keep *']
process.schedule=cms.Schedule(process.p)
try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")
