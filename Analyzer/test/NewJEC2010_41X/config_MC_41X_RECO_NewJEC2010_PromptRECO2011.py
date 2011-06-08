import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDtuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
#---Needed to Reconsctruct on the fly from uncleaned SCs without timing cut for slpikes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
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
process.GlobalTag.globaltag = cms.string('GR_R_41_V0::All')


from PhysicsTools.PatAlgos.tools.metTools import *                       
                                                                         
addPfMET(process,'PF')
addTcMET(process,"TC")       

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger

# Select calo jets
process.patJetCorrFactors.levels = cms.vstring(['L1Offset','L2Relative','L3Absolute'])
process.selectedPatJets.cut = cms.string('pt > 10 & abs(eta) < 3.0')


# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute'])),
                 doType1MET    = True,
                 doL1Cleaning  = True,
                 doL1Counters  = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID       = True,
                 jetIdLabel    = "ak5"
                )
process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')


# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# Modify defaults setting to avoid an over-efficiency in the presence of OFT PU
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
 rocess.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)


#---For FastJet JEC for 41X--
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)
 
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Ghost_EtaMax = cms.double(5.0)
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)                                                                                                                                     
##----NEW PAT  PHOTN COLLLECTION
# make a new patPhotons
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
process.uncleanpatPhotons = patPhotons.clone(
   photonSource = cms.InputTag("photons::ADDSkim")
) 
# make a new selectedPatCandidates
from PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff import *
process.selecteduncleanPatPhotons = selectedPatPhotons.clone(
    src = cms.InputTag("uncleanpatPhotons"),                                                                                                                 
    cut = cms.string('') 
)                 

# make cleanPatCandidates
from PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff import *
process.uncleanPatPhotons = cleanPatPhotons.clone(            
    src = cms.InputTag("selecteduncleanPatPhotons")                                                                                                                 
)      

##--default pick new one so explictely spefifying RECO here
process.patPhotons.photonSource = cms.InputTag("photons::RECO")
process.patPhotons.photonIDSources = cms.PSet(
        PhotonCutBasedIDTight = cms.InputTag("PhotonIDProd","PhotonCutBasedIDTight","RECO"),
        PhotonCutBasedIDLoose = cms.InputTag("PhotonIDProd","PhotonCutBasedIDLoose","RECO")
    )



# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
        '/store/mc/Spring11/QCD_MuPt5EtaFilter_7TeV-pythia6/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0002/E4C98778-D75B-E011-887C-00304867C034.root',
        '/store/mc/Spring11/QCD_MuPt5EtaFilter_7TeV-pythia6/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0002/E2714E05-CA5B-E011-9811-0018F3D09634.root',
        '/store/mc/Spring11/QCD_MuPt5EtaFilter_7TeV-pythia6/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0002/E2188350-025C-E011-B0DC-001A92810A98.root'
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

                              uncleanphotonTag = cms.untracked.InputTag("selecteduncleanPatPhotons"),
                              cscTag           = cms.untracked.InputTag("cscSegments"),
                              rpcTag           = cms.untracked.InputTag("rpcRecHits"),
                              rechitBTag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("ecalRecHit:EcalRecHitsEE"),
                              hcalrechitTag    = cms.untracked.InputTag("hbhereco"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag           = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag           = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),#check what you need here HLT or REDIGI3..
                              triggerEventTag  = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),#check what you need here HLT or REDIGI3..
                              hltlabel          = cms.untracked.string("HLT"),#check what you need here HLT or REDIGI3..
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              BeamHaloSummary  = cms.untracked.InputTag("BeamHaloSummary"),
                              pileup           = cms.untracked.InputTag("PileUpInfo"),
                              outFile          = cms.untracked.string("Histo_MC_A_RECO.root"),
                              runphotons       = cms.untracked.bool(True),
                              rununcleanphotons= cms.untracked.bool(True),
                              runHErechit      = cms.untracked.bool(True),
                              runrechit        = cms.untracked.bool(True),
                              runmet           = cms.untracked.bool(True),
                              rungenmet        = cms.untracked.bool(True),
                              runPFmet         = cms.untracked.bool(True),
                              runTCmet         = cms.untracked.bool(True),
                              runjets          = cms.untracked.bool(True),
                              runpfjets        = cms.untracked.bool(True),
                              rungenjets       = cms.untracked.bool(True),
                              runelectrons     = cms.untracked.bool(True),
                              runtaus          = cms.untracked.bool(True),
                              runDetailTauInfo = cms.untracked.bool(True),
                              runmuons         = cms.untracked.bool(True),
                              runcosmicmuons   = cms.untracked.bool(True),
                              rungenParticleCandidates = cms.untracked.bool(True),
                              runHLT           = cms.untracked.bool(True),
                              runL1            = cms.untracked.bool(True),
                              runscraping      = cms.untracked.bool(False),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True),
                              runCSCseg        = cms.untracked.bool(True),
                              runBeamHaloSummary= cms.untracked.bool(True),
                              runPileUp         = cms.untracked.bool(True),
                              runRPChit        = cms.untracked.bool(True),
                              isAOD            = cms.untracked.bool(False),
                              debug            = cms.untracked.bool(False)
                              )



process.NewPatPhotons = cms.Sequence(
               process.uncleanpatPhotons*
               process.selecteduncleanPatPhotons
               #process.uncleanPatPhotons
               )

#All paths are here
process.p = cms.Path(
   process.HBHENoiseFilter*
   process.kt6PFJets*
   process.ak5PFJets*
   process.NewPatPhotons*
   process.patDefaultSequence*
   process.demo
   )


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
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

#print process.dumpPython()
