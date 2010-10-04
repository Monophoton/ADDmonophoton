import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDtuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoEcal.Configuration.RecoEcal_cff")
from Configuration.StandardSequences.Reconstruction_cff import *
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#spike cleaning
#process.load("EGamma.EGammaSkims.promptRecoTrackCorrections_cff")
#process.load("EGamma.EGammaSkims.cleanReRecoSequence_cff")
#======

#process.ecalCleanClustering = cms.Sequence(process.cleanedEcalClusters*process.cleanedEgammaSkimReco)
## global tag for data

process.GlobalTag.globaltag = cms.string('GR_R_38X_V13A::All')

from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])
addPfMET(process, 'PF')
addTcMET(process,"TC")


# get the 900 GeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Summer09_7TeV_ReReco332")

# run ak5 gen jets
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *


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

# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
'dcache:/pnfs/cms/WAX/11/store/user/askew/HighPT/Jul16th/HighETSkim_9_1.root',
'dcache:/pnfs/cms/WAX/11/store/user/askew/HighPT/Jul16th/HighETSkim_99_1.root'
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_1.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_2.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_3x.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_4.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_5.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_6.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_7.root',
#'file:/uscms_data/d2/shruti/EventEGv4/EventsR01/EG_EVENTS_8.root'
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
                              outFile          = cms.untracked.string("output_00.root"),
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

##This part is only needed if you're going to redo
##any part of the clustering, e.g. if you're
##running 3_8 over 3_6 data
########################
    ## TRACKING:
    ## Skip events with HV off
process.newSeedFromTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.newSeedFromPairs.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
process.secTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.fifthSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 20000
process.fourthPLSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
process.thTripletsA.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000
process.thTripletsB.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000

    ###### FIXES TRIPLETS FOR LARGE BS DISPLACEMENT ######

    ### prevent bias in pixel vertex
process.pixelVertices.useBeamConstraint = False
    
    ### pixelTracks
    #---- new parameters ----
process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = 4.06
process.pixelTracks.RegionFactoryPSet.RegionPSet.originHalfLength = cms.double(40.6)

    ### 0th step of iterative tracking
    #---- new parameters ----
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(4.06)  
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 40.6

    ### 2nd step of iterative tracking
    #---- new parameters ----
process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)  
process.secTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 44.7

    ## ECAL 
process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]


    ## particle flow HF cleaning
process.particleFlowRecHitHCAL.LongShortFibre_Cut = 30.
process.particleFlowRecHitHCAL.ApplyPulseDPG = True

    ## HF cleaning for data only
process.hcalRecAlgos.SeverityLevels[3].RecHitFlags.remove("HFDigiTime")
process.hcalRecAlgos.SeverityLevels[4].RecHitFlags.append("HFDigiTime")

    ##beam-halo-id for data only
process.CSCHaloData.ExpectedBX = cms.int32(3)

    ## hcal hit flagging
process.hfreco.PETstat.flagsToSkip  = 2
process.hfreco.S8S1stat.flagsToSkip = 18
process.hfreco.S9S1stat.flagsToSkip = 26

process.cleanedEcalClusters = cms.Sequence(ecalClusters) 

process.cleanedLocalreco = cms.Sequence(siPixelRecHits*ckftracks*particleFlowCluster)

process.cleanedElectronseeding = cms.Sequence(trackerDrivenElectronSeeds*ecalDrivenElectronSeeds*electronMergedSeeds)

process.cleanedElectrontracking = cms.Sequence(electronCkfTrackCandidates*electronGsfTracks)

process.cleanedElectronreco = cms.Sequence(particleFlowReco*pfElectronTranslator*gsfElectronSequence)

process.cleanedEgammaSkimReco = cms.Sequence(process.cleanedLocalreco*process.cleanedElectronseeding*process.cleanedElectrontracking*process.cleanedElectronreco*process.conversionSequence*process.trackerOnlyConversionSequence*process.photonSequence*process.interestingEgammaIsoDetIds*process.photonIDSequence*process.eIdSequence)
process.ecalCleanClustering = cms.Sequence(process.cleanedEcalClusters*process.cleanedEgammaSkimReco)
################End extra stuff.
#Comment out process.ecalCleanClustering if you don't
#need to redo clustering.
process.p = cms.Path(
   # process.hltHighLevel*
   process.ecalCleanClustering*
 #  process.hltHighLevel*
   process.patDefaultSequence*
   process.demo
    )


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.options.wantSummary = True


process.schedule=cms.Schedule(process.p)
try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")
