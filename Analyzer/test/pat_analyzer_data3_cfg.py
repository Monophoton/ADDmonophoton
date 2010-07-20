# This is an example PAT configuration showing the usage of PAT on minbias data

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *
## global tag for data
#process.GlobalTag.globaltag = cms.String('Start3x_v27::All')
process.GlobalTag.globaltag = cms.string('GR_R_36X_V10::All')
# turn off MC matching for the process
#removeMCMatching(process, ['All'])
# add pf met
#INPUTFILELIST="Analysis.Analyzer.Files_DataSkimmed_cff"
#INPUTFILELIST="Analysis.Analyzer.Files_PromptReco_cff"
#INPUTFILELIST="Analysis.Analyzer.Files_DataSkimmed22Apr_1_cff"
#process.load(INPUTFILELIST)
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


process.hltHighLevel = cms.EDFilter("HLTHighLevel",
                                    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                    HLTPaths = cms.vstring("HLT_Photon20_L1R"),
                                    eventSetupPathsKey = cms.string(''),
                                    andOr = cms.bool(False),
                                    throw = cms.bool(True)
                                    )


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
#    'rfio:/castor/cern.ch/user/b/berzano/Skim/7TeV/MonoPhoton/MC/MonoPhotonSkim_1.root',
    #'/store/data/Run2010A/Cosmics/RECO/v2/000/137/243/F019D09A-3173-DF11-97E1-001D09F24489.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/137/436/0CF7FFDA-8F74-DF11-8D44-0030487A195C.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/137/433/7AAAE230-7C74-DF11-8575-0030487C8CB6.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/137/415/EAE5A512-8674-DF11-AEFB-0030487C912E.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/137/411/82D19F49-7774-DF11-8821-0030487CD7EE.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/137/411/3613EABF-7F74-DF11-924D-0030487CD7B4.root'
    # '/store/data/Run2010A/EGMonitor/RECO/v2/000/136/088/381B68C5-7267-DF11-B2B5-0030487A18F2.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/136/088/269975FB-6D67-DF11-9797-0030487CD718.root',
    #'/store/data/Run2010A/EGMonitor/RECO/v2/000/136/087/BEBA3AF3-7467-DF11-9CAB-00304879EDEA.root'
    # '/store/data/Run2010A/EG/RECO/v2/000/136/088/D887BA58-7167-DF11-8DF0-0030487C6062.root',
    #       '/store/data/Run2010A/EG/RECO/v2/000/136/088/B6CC1FC5-7267-DF11-B838-0030487CD6DA.root',
    #      '/store/data/Run2010A/EG/RECO/v2/000/136/088/887017DA-6D67-DF11-BDCA-0030487CD716.root',
    #     '/store/data/Run2010A/EG/RECO/v2/000/136/088/0EEF4E97-8367-DF11-8D7E-0030487CD162.root',
    #    '/store/data/Run2010A/EG/RECO/v2/000/136/087/D62BE917-6967-DF11-A82E-0030487CD718.root',
    '/store/data/Run2010A/EG/RECO/v2/000/136/087/C086B14B-7D67-DF11-B58C-0030487CD6DA.root'
    
    # '/store/data/Run2010A/Mu/RECO/v2/000/136/100/B640EC18-ED67-DF11-A369-001D09F23A84.root',
    #     '/store/data/Run2010A/Mu/RECO/v2/000/136/100/B0E0A03F-0468-DF11-8331-000423D98930.root',
    #  '/store/data/Run2010A/Mu/RECO/v2/000/136/100/AEAD3EDD-D567-DF11-A4B0-0030487C90D4.root',
    #   '/store/data/Run2010A/Mu/RECO/v2/000/136/100/A0EEBCD1-DA67-DF11-8C30-000423D991D4.root',
    #  '/store/data/Run2010A/Mu/RECO/v2/000/136/100/9A8E77CA-D867-DF11-B1E9-0019B9F70607.root'
    
    ] );
process.source.fileNames = readFiles

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
                              rechitBTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                              metTag           = cms.untracked.InputTag("patMETs"),
                              PFmetTag           = cms.untracked.InputTag("patMETsPF"),
                              TCmetTag           = cms.untracked.InputTag("patMETsTC"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","",""),
                              outFile          = cms.untracked.string("outputEGmon.root"),
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
#process.out.fileName = "DROPPED"

process.p = cms.Path(
    # process.hltHighLevel*
    process.patDefaultSequence*
    process.demo
    )

# rename output file
#process.out.fileName = cms.untracked.string('rfio:/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/RealData/ExpressPhysics/132440_Files/pat.root')
#process.out.fileName = cms.untracked.string('/tmp/sandhya/pat.root')
#process.out.fileName = cms.untracked.string('5K_patb1.root')
# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# process all the events
process.maxEvents.input =6687
process.options.wantSummary = True

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
#from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
#process.out.outputCommands = patEventContentNoCleaning
#process.out.outputCommands += patExtraAodEventContent
#process.out.outputCommands += patTriggerEventContent
#process.out.outputCommands += ['keep recoPFCandidates_particleFlow_*_*']
#process.out.outputCommands += ['keep *']

