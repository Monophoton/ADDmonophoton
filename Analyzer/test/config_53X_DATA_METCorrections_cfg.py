import FWCore.ParameterSet.Config as cms

process = cms.Process("ADDtuple")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## global tag for data
process.GlobalTag.globaltag = cms.string('GR_P_V40_AN1::All')

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'],
                 outputModules = [])


# Add MET collection for PAT
from PhysicsTools.PatAlgos.tools.metTools import *                       
addPfMET(process,'PF')
addTcMET(process,"TC")       


#---------Type MET Type-1,0 and phi Corrections---------------------------
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
process.patType1CorrectedPFMet.srcType1Corrections = cms.VInputTag(
       cms.InputTag('patPFJetMETtype1p2Corr','type1'),
       cms.InputTag('patPFMETtype0Corr'),
       cms.InputTag('pfMEtSysShiftCorr') 
)      

#make sure it is 2012 one
process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
#tell which  pf-jet otherwise it will complaint that jet is not tuype of PF!!
process.selectedPatJetsForMETtype1p2Corr.src = cms.InputTag('patJetsAK5PF')
process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L2L3Residual") # NOTE: use "L3Absolute" for MC 


# Select calo jets
process.patJetCorrFactors.levels = cms.vstring(['L1Offset','L2Relative','L3Absolute','L2L3Residual'])
process.selectedPatJets.cut = cms.string('pt > 10 & abs(eta) < 3.0')

# Add PF jets this also apply MET-1 corrections
from PhysicsTools.PatAlgos.tools.jetTools import *

addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])),
                 doType1MET    = True,
                 doL1Cleaning  = True,
                 doL1Counters  = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID       = True,
                 jetIdLabel    = "ak5"
                )

process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
##Need this for valumap to know which jet 
process.puJetId.jets=cms.InputTag("selectedPatJetsAK5PF")
process.puJetMva.jets=cms.InputTag("selectedPatJetsAK5PF")
process.pileupJetIdProducer.jets=cms.InputTag("selectedPatJetsAK5PF")

#---------Fast Rho calculation-------------------------
#Rho for eta= 2.5 
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets25.Ghost_EtaMax = cms.double(2.5)
process.kt6PFJets25.doAreaFastjet = True 
process.kt6PFJets25.voronoiRfact = 0.9
process.fastjetSequence25 = cms.Sequence( process.kt6PFJets25 )


#Fast jet correction
process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = True 


process.load('EGamma.EGammaAnalysisTools.photonIsoProducer_cfi')
process.phoPFIso.verbose = False
process.phoPFIso.photonTag= 'photons::RECO'





#-------------ALL MET and other Cleaning Filters are here------------
# HB + HE noise filtering 
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#hcalLaserFilter
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
#Bad EE SC filter, not needed but goot to have them 
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
#Tracking Failure filter
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.goodVertices = cms.EDFilter("VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
 )
##ECAL dead cell filter 
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
#CSC Halo Filters
#process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

process.AllMETFilters= cms.Sequence( process.HBHENoiseFilter
                                    *process.hcalLaserEventFilter
                                     *process.eeBadScFilter
                                     *(process.goodVertices*process.trackingFailureFilter)
                                     *process.EcalDeadCellTriggerPrimitiveFilter
                                     #*process.CSCTightHaloFilter 
                                     )

process.patPhotons.photonSource = cms.InputTag("photons::RECO")
process.patPhotons.photonIDSources = cms.PSet(
    PhotonCutBasedIDTight =
    cms.InputTag("PhotonIDProd","PhotonCutBasedIDTight","RECO"),
    PhotonCutBasedIDLoose =
    cms.InputTag("PhotonIDProd","PhotonCutBasedIDLoose","RECO")
    )


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles.extend( [
    'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/weinberg/monoPhoton/skims/singlePhoton_run2012C_promptReco_v1/monoPhotonSkim_10_2_fmK.root'


    
    ] );

process.source = cms.Source("PoolSource",
    fileNames = readFiles
)

#closes files after code is done running on that file
process.options = cms.untracked.PSet(
	fileMode = cms.untracked.string('NOMERGE')
)



##---input to analyzer
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
                              caloTowerTag   = cms.untracked.InputTag("towerMaker"),
                              cscTag           = cms.untracked.InputTag("cscSegments"),
                              rpcTag           = cms.untracked.InputTag("rpcRecHits"),
                              rechitBTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                              hcalrechitTag    = cms.untracked.InputTag("reducedHcalRecHits:hbhereco"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag         = cms.untracked.InputTag("patType1CorrectedPFMet"),
                              TCmetTag         = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),#check what you need here HLT or REDIGI3..
                              triggerEventTag  = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),#check what you need here HLT or REDIGI3..
                              hltlabel          = cms.untracked.string("HLT"),   #check what you need here HLT or REDIGI3..
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              BeamHaloSummary  = cms.untracked.InputTag("BeamHaloSummary"),
                              pileup           = cms.untracked.InputTag("addPileupInfo"),
                              rhoLabel         = cms.untracked.InputTag("kt6PFJets", "rho"),
                              sigmaLabel       = cms.untracked.InputTag("kt6PFJets", "sigma"),
                              rhoLabel25      = cms.untracked.InputTag("kt6PFJets25", "rho"),
                              sigmaLabel25     = cms.untracked.InputTag("kt6PFJets25", "sigma"),
                              runphotons       = cms.untracked.bool(True),
                              rununcleanphotons= cms.untracked.bool(False),
                              runHErechit      = cms.untracked.bool(True),
                              runrechit        = cms.untracked.bool(True),
                              runmet           = cms.untracked.bool(True),
                              rungenmet        = cms.untracked.bool(False),
                              runPFmet         = cms.untracked.bool(True),
                              runTCmet         = cms.untracked.bool(True),
                              runjets          = cms.untracked.bool(True),
                              runpfjets        = cms.untracked.bool(True),
                              rungenjets        = cms.untracked.bool(False),
                              runelectrons     = cms.untracked.bool(True),
                              runtaus          = cms.untracked.bool(True),
                              runDetailTauInfo = cms.untracked.bool(True),
                              runmuons         = cms.untracked.bool(True),
                              runcosmicmuons   = cms.untracked.bool(True),
                              rungenParticleCandidates = cms.untracked.bool(False),
                              runHLT           = cms.untracked.bool(True),
                              runL1            = cms.untracked.bool(True),
                              runscraping      = cms.untracked.bool(True),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True),
                              ##OFF for AOD---
                              runCSCseg        = cms.untracked.bool(False),
                              runRPChit        = cms.untracked.bool(False),
                              #---------------
                              runcaloTower     = cms.untracked.bool(True),
                              runBeamHaloSummary= cms.untracked.bool(True),
                              runPileUp         = cms.untracked.bool(False),
                              isAOD             = cms.untracked.bool(True),
                              debug            = cms.untracked.bool(False),
                              outFile          = cms.untracked.string('Histo_AOD.root'),
                              Photons = cms.InputTag('photons::RECO'),
                              IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                                                           cms.InputTag('phoPFIso:phIsoForGsfEle'),
                                                           cms.InputTag('phoPFIso:nhIsoForGsfEle'),
                                                           ),
                              
                              
                              )




#All paths are here
process.p = cms.Path(
    process.AllMETFilters *
    process.fastjetSequence25 *
    process.ak5PFJets *
    process.phoPFIso *
    process.pfMEtSysShiftCorrSequence *
    process.patDefaultSequence *
    process.producePatPFMETCorrections *    # Produce MET corrections
    process.puJetIdSqeuence *               # pileup based jet id
    process.demo
  )


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

# process all the events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )


