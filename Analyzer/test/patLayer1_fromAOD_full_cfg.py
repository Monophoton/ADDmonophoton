# This is an example PAT configuration showing the usage of PAT on full sim samples

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# note that you can use a bunch of core tools of PAT 
# to taylor your PAT configuration; for a few examples
# uncomment the following lines

from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, 'Muons')
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Photons', 'Electrons', 'Muons', 'Taus'])
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.demo = cms.EDAnalyzer('Analyzer',
                              electronTag      = cms.untracked.InputTag("cleanLayer1Electrons"),
                              tauTag           = cms.untracked.InputTag("cleanLayer1Taus"),
                              muonTag          = cms.untracked.InputTag("cleanLayer1Muons"),
                              jetTag           = cms.untracked.InputTag("cleanLayer1Jets"),
                              photonTag        = cms.untracked.InputTag("cleanLayer1Photons"),
                              rechitBTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                              rechitETag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                              metTag           = cms.untracked.InputTag("layer1METs"),
                              HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                              Tracks           = cms.untracked.InputTag("generalTracks"),
                              Vertices         = cms.untracked.InputTag("offlinePrimaryVertices","","RECO"),
                              outFile          = cms.untracked.string("Histo.root"),
                              runphotons       = cms.untracked.bool(True),
                              runrechit        = cms.untracked.bool(True),
                              runmet           = cms.untracked.bool(True),
                              runjets          = cms.untracked.bool(True),
                              runelectrons     = cms.untracked.bool(True),
                              runtaus          = cms.untracked.bool(True),
                              runmuons         = cms.untracked.bool(True),
                              rungenParticleCandidates = cms.untracked.bool(True),
                              runHLT           = cms.untracked.bool(True),
                              runtracks        = cms.untracked.bool(True),
                              runvertex        = cms.untracked.bool(True)
                              )

# let it run
process.outpath = cms.EndPath(
    process.patDefaultSequence
#    +process.content
    +process.demo
)

# In addition you usually want to change the following parameters:
#
process.GlobalTag.globaltag = 'MC_31X_V5::All'   ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
process.maxEvents.input = 1         ##  (e.g. -1 to run on all events)
process.out.outputCommands = ['keep *']  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#process.out.fileName = '/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/3_1_2/Pat.root'            ##  (e.g. 'myTuple.root')
process.source.fileNames = ['/store/relval/CMSSW_3_1_2/RelValSingleGammaPt35/GEN-SIM-RECO/MC_31X_V3-v1/0007/123C78F5-9078-DE11-8BAD-001D09F23A61.root','/store/relval/CMSSW_3_1_2/RelValSingleGammaPt35/GEN-SIM-RECO/MC_31X_V3-v1/0006/E6C7ED95-4878-DE11-B082-000423D98BE8.root']
process.options.wantSummary = False       ##  (to suppress the long output at the end of the job)    
