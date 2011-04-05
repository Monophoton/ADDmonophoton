import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDtuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
#---Needed to Reconsctruct on the fly from uncleaned SCs without timing cut for slpikes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoEgamma.EgammaPhotonProducers.conversionTracks_cff')

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
process.GlobalTag.globaltag = cms.string('GR_R_39X_V6::All')


from PhysicsTools.PatAlgos.tools.metTools import *                       
removeMCMatching(process, ['All'],
          outputInProcess = False )                                       
                                                                         
addPfMET(process,'PF')
addTcMET(process,"TC")       

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger

# Select jets
process.selectedPatJets.cut = cms.string('pt > 10 & abs(eta) < 3.0')

# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
#'dcache:/pnfs/cms/WAX/11/store/user/lpcgg/miceli/Photon/PhotonSkim_Run2010B-Nov4ReReco_v1_146240-147829/3c009990a403ffefbef1d67824f3333a/phoskim_88_2_Faz.root'
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E8139119-2110-E011-855E-003048C69318.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E80A4E9E-1F10-E011-B243-003048D439A8.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E6F8E442-2010-E011-B25D-0030487D5EB1.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E6D83A9A-1F10-E011-9F55-0030487D5EB1.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E4698D2B-2010-E011-89CF-003048C693C8.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E299636E-1F10-E011-A71F-0030487D5D7B.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E25FD4C8-1F10-E011-98DA-0030487E4EC5.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E05C8E7D-1F10-E011-9B13-0030487F1665.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/E0088F3F-1F10-E011-85E6-003048C693BA.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/DE243F7E-2010-E011-8A58-003048C6903C.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/DC8D550F-2110-E011-80B0-003048C69292.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/DC20AA28-2010-E011-B7BD-0030487D5EB1.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/DAE6744F-1F10-E011-A3F6-0030487D5DA5.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/DA222625-2010-E011-A7DC-003048C69318.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/D8DB462E-2110-E011-A12E-003048C69318.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/D8BBDFC5-1F10-E011-80A1-003048C693C8.root',
        '/store/data/Run2010A/EG/RECO/Dec22ReReco_v1/0027/D8484310-2110-E011-8AD5-003048C693D0.root'

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
                              photonTag        = cms.untracked.InputTag("selectedPatPhotons"),
                              cscTag           = cms.untracked.InputTag("cscSegments"),
                              rpcTag           = cms.untracked.InputTag("rpcRecHits"),
                              rechitBTag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEE"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag           = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag           = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              outFile          = cms.untracked.string("Histo_Data_A.root"),
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
                              runscraping      = cms.untracked.bool(True),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True),
                              runCSCseg        = cms.untracked.bool(True),
                              runRPChit        = cms.untracked.bool(True),
                              debug            = cms.untracked.bool(False)
                              )




process.monoSkim = cms.EDFilter("MonoPhotonSkimmer",
  phoTag = cms.InputTag("photons"),   # photons::RECO will look for original collection
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

#Remove the time severity flag from cleanedHybridSuperCluster
process.cleanedHybridSuperClusters.RecHitSeverityToBeExcluded= cms.vint32(4,5)

#For Uncleaned ones
process.uncleanedHybridSuperClusters.RecHitSeverityToBeExcluded= cms.vint32(999)


#Reconstrucition of photon from SC without the timing cuts
process.photonReReco = cms.Sequence(
        process.ecalClusters*
	ckfTracksFromConversions*
	process.conversionSequence*
	process.photonSequence*
	process.photonIDSequence)

#All paths are here
process.p = cms.Path(
   process.photonReReco*
   process.monoSkim*
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

