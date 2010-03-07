# This is an example PAT configuration showing the usage of PAT on minbias data

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("Analysis.Analyzer.EXOSinglePhoOct09_mod_cff")
## global tag for data
process.GlobalTag.globaltag = cms.string('GR09_P_V8_34X::All')

# turn off MC matching for the process
removeMCMatching(process, 'All')

# get the 900 GeV jet corrections
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "900GeV")

# remove the tag infos
process.allLayer1Jets.addTagInfos = False
# require jet pt > 20 (L2+L3 corrected)
#process.selectedLayer1Jets.cut = cms.string('pt > 10')
# look for at least one jet
#process.countLayer1Jets.minNumber = 0

from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process,"TC")
addPfMET(process,"PF")


# Add the files for runs 123575 and 123596
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0006/5A5773A3-B9ED-DE11-A99A-0026189437F8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/SD_InterestingEvents-Dec19thSkim_341_v1/0006/CA156CCC-B9ED-DE11-B022-00261894390B.root'
        ] );
process.source.fileNames = readFiles

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41)')

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.demo = cms.EDAnalyzer('Analyzer',
                                                            electronTag      = cms.untracked.InputTag("selectedLayer1Electrons"),
                                                            tauTag           = cms.untracked.InputTag("selectedLayer1Taus"),
                                                            muonTag          = cms.untracked.InputTag("selectedLayer1Muons"),
                                                            jetTag           = cms.untracked.InputTag("selectedLayer1Jets"),
                                                            photonTag        = cms.untracked.InputTag("selectedLayer1Photons"),
                                                            rechitBTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                                                            rechitETag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                                                            metTag           = cms.untracked.InputTag("layer1METs"),
                                                            PFmetTag           = cms.untracked.InputTag("layer1METsPF"),
                                                            TCmetTag           = cms.untracked.InputTag("layer1METsTC"),
                                                            HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                                            Tracks           = cms.untracked.InputTag("generalTracks"),
                                                            Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","","RECO"),
                                                            outFile          = cms.untracked.string("Histo.root"),
                                                            runphotons       = cms.untracked.bool(True),
                                                            runrechit        = cms.untracked.bool(True),
                                                            runmet           = cms.untracked.bool(True),
                                                            rungenmet        = cms.untracked.bool(False),
                                                            runPFmet         = cms.untracked.bool(True),
                                                            runTCmet         = cms.untracked.bool(True),
                                                            runjets          = cms.untracked.bool(True),
                                                            runelectrons     = cms.untracked.bool(True),
                                                            runtaus          = cms.untracked.bool(True),
                                                            runmuons         = cms.untracked.bool(True),
                                                            rungenParticleCandidates = cms.untracked.bool(False),
                                                            runHLT           = cms.untracked.bool(True),
                                                            runL1            = cms.untracked.bool(True),
                                                            runtracks        = cms.untracked.bool(True),
                                                            runvertex        = cms.untracked.bool(True)
                                                            )

# let it run
process.p = cms.Path(
    process.hltLevel1GTSeed+
    process.exoticaRecoSinglePhoFilter+
    process.patDefaultSequence+
    process.demo
    )

# rename output file
#process.out.fileName = cms.untracked.string('/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/RealData_SinglePhoSkim/firstdata_pat.root')
process.out.fileName = cms.untracked.string('MinBias_pat.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = 1000
process.options.wantSummary = True


from Configuration.EventContent.EventContent_cff import RECOEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += (['keep *'])
process.out.outputCommands += RECOEventContent.outputCommands
process.out.outputCommands += patEventContent
process.out.outputCommands += patEventContentNoLayer1Cleaning
