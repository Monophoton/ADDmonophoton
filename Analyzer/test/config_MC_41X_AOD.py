import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDtuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
#---Needed to Reconsctruct on the fly from uncleaned SCs without timing cut for slpikes
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('RecoEgamma.EgammaPhotonProducers.conversionTracks_cff')

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoEcal.Configuration.RecoEcal_cff")
from Configuration.StandardSequences.Reconstruction_cff import *
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## global tag for data
process.GlobalTag.globaltag = cms.string('MC_311_V2::All')


from PhysicsTools.PatAlgos.tools.metTools import *                       
                                                                         
addPfMET(process,'PF')
addTcMET(process,"TC")       

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger

# Select calo jets
process.patJetCorrFactors.levels = cms.vstring(['L2Relative','L3Absolute','L5Flavor', 'L7Parton'])
process.selectedPatJets.cut = cms.string('pt > 10 & abs(eta) < 3.0')

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute','L5Flavor','L7Parton'])),
                 doType1MET    = True,
                 doL1Cleaning  = True,
                 doL1Counters  = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID       = True,
                 jetIdLabel    = "ak5"
                )
process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')

# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
#MC
        '/store/mc/Spring11/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/AODSIM/PU_S1_START311_V1G1-v1/0015/803351A7-375D-E011-BE83-002618943885.root',
        '/store/mc/Spring11/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/AODSIM/PU_S1_START311_V1G1-v1/0015/6A546083-375D-E011-9800-002618943875.root',
        '/store/mc/Spring11/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/AODSIM/PU_S1_START311_V1G1-v1/0015/20A3BE89-375D-E011-AD77-001A92810A92.root'
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

#closes files after code is done running on that file
process.options = cms.untracked.PSet(
	fileMode = cms.untracked.string('NOMERGE')
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.demo = cms.EDAnalyzer('Analyzer',
                              electronTag      = cms.untracked.InputTag("selectedPatElectrons"),
                              tauTag           = cms.untracked.InputTag("selectedPatTaus"),
                              muonTag          = cms.untracked.InputTag("selectedPatMuons"),
                              cosMuonTag       = cms.untracked.InputTag("muonsFromCosmics"),
                              jetTag           = cms.untracked.InputTag("selectedPatJets"),
                              pfjetTag         = cms.untracked.InputTag("selectedPatJetsAK5PF"),
                              genjetTag        = cms.untracked.InputTag("ak5GenJets"),
                              photonTag        = cms.untracked.InputTag("selectedPatPhotons"),
                              cscTag           = cms.untracked.InputTag("cscSegments"),
                              rpcTag           = cms.untracked.InputTag("rpcRecHits"),
                              rechitBTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                              hcalrechitTag    = cms.untracked.InputTag("reducedHcalRecHits:hbhereco"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag         = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag         = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","",""),
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              BeamHaloSummary  = cms.untracked.InputTag("BeamHaloSummary"),
                              pileup           = cms.untracked.InputTag("PileUpInfo"),
                              outFile          = cms.untracked.string("Histo_MC_AOD.root"),
                              runphotons       = cms.untracked.bool(True),
                              runHErechit      = cms.untracked.bool(True),
                              runrechit        = cms.untracked.bool(True),
                              runmet           = cms.untracked.bool(True),
                              rungenmet        = cms.untracked.bool(True),
                              runPFmet         = cms.untracked.bool(True),
                              runTCmet         = cms.untracked.bool(True),
                              runjets          = cms.untracked.bool(True),
                              runpfjets        = cms.untracked.bool(True),
                              rungenjets        = cms.untracked.bool(True),
                              runelectrons     = cms.untracked.bool(True),
                              runtaus          = cms.untracked.bool(True),
                              runmuons         = cms.untracked.bool(True),
                              runcosmicmuons   = cms.untracked.bool(True),
                              rungenParticleCandidates = cms.untracked.bool(True),
                              runHLT           = cms.untracked.bool(True),
                              runL1            = cms.untracked.bool(True),
                              runscraping      = cms.untracked.bool(True),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True),
                              ##OFF for AOD---
                              runCSCseg        = cms.untracked.bool(False),
                              runRPChit        = cms.untracked.bool(False),
                              #---------------
                              runBeamHaloSummary= cms.untracked.bool(True),
                              runPileUp         = cms.untracked.bool(True),
                              isAOD             = cms.untracked.bool(True),
                              debug            = cms.untracked.bool(False)
                              )




#process.monoSkim = cms.EDFilter("MonoPhotonSkimmer",
#  phoTag = cms.InputTag("photons"),
#  selectEE = cms.bool(True),
#  selectTrack = cms.bool(False),
#  ecalisoOffsetEB = cms.double(5000.),
#  ecalisoSlopeEB = cms.double(0.15),
#  hcalisoOffsetEB = cms.double(7000.0),
#  hcalisoSlopeEB = cms.double(0.),
#  hadoveremEB = cms.double(0.05),
#  minPhoEtEB = cms.double(30.),
#  scHighEtThreshEB = cms.double(100.),
#  ecalisoOffsetEE = cms.double(5000.),
#  ecalisoSlopeEE = cms.double(0.15),
#  hcalisoOffsetEE = cms.double(10000.),
#  hcalisoSlopeEE = cms.double(0.),
#  hadoveremEE = cms.double(0.05),
#  minPhoEtEE = cms.double(30.),
#  scHighEtThreshEE = cms.double(100.),
# )



#Remove the time severity flag from cleanedHybridSuperCluster
#process.cleanedHybridSuperClusters.RecHitSeverityToBeExcluded= cms.vint32(4,5)
#For Uncleaned ones
#process.uncleanedHybridSuperClusters.RecHitSeverityToBeExcluded= cms.vint32(999)


#Reconstrucition of photon from SC without the timing cuts
#process.photonReReco = cms.Sequence(
#        process.ecalClusters*
#	ckfTracksFromConversions*
#	process.conversionSequence*
#	process.photonSequence*
#	process.photonIDSequence)

#All paths are here
process.p = cms.Path(
#   process.photonReReco*
   process.patDefaultSequence*
   process.demo
   )



# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.options.wantSummary = True
#process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.schedule=cms.Schedule(process.p)
try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")

#print process.dumpPython()
