import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDSkim")
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

from Configuration.StandardSequences.Reconstruction_cff import *
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi import *
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## global tag for data
process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')

# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
         '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/3A1103FE-817C-E011-B3C4-001A92971B3A.root'  
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)


#--------For Uncleaned Photon-----
process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi") 
process.uncleanSCRecovered.cleanScCollection = cms.InputTag("correctedHybridSuperClusters")
process.photonCore.scHybridBarrelProducer = cms.InputTag("uncleanSCRecovered:uncleanHybridSuperClusters")
photons.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
photons.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")
 
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *
newisolationSumsCalculator = isolationSumsCalculator.clone()
newisolationSumsCalculator.barrelEcalRecHitProducer = cms.string('reducedEcalRecHitsEB')
newisolationSumsCalculator.endcapEcalRecHitProducer = cms.string('reducedEcalRecHitsEE')
newisolationSumsCalculator.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
newisolationSumsCalculator.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
photons.isolationSumsCalculatorSet = newisolationSumsCalculator



# let it run
process.monoSkim = cms.EDFilter("MonoPhotonSkimmer",
  phoTag = cms.InputTag("photons::RECO"),
  selectEE = cms.bool(True),
  selectTrack = cms.bool(False),
  ecalisoOffsetEB = cms.double(7000.),
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


#Sequence for ucleaned photon
process.uncleanPhotons = cms.Sequence(
               process.uncleanSCRecovered*
               process.photonSequence*
               process.photonIDSequence
               )            
                            
# let it run
process.monoUncleanSkim = cms.EDFilter("MonoPhotonSkimmer",
  phoTag = cms.InputTag("photons::ADDSkim"),
  selectEE = cms.bool(True),
  selectTrack = cms.bool(False),
  ecalisoOffsetEB = cms.double(7000.),
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



#All paths are here
process.p = cms.Path(process.monoSkim)

process.pp = cms.Path(process.uncleanPhotons*
                     process.monoUncleanSkim)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
#process.options.wantSummary = True


# define output event selection to be that which satisfies 'p'
process.patEventSelection = cms.PSet(
        SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p','pp')
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
    fileName = cms.untracked.string('phoskim.root')
)


process.outpath = cms.EndPath(process.out)


try:
   import readline
except ImportError:
   print "Module readline not available."
else:
   import rlcompleter
   readline.parse_and_bind("tab: complete")
