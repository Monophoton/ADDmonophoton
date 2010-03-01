# This is an example PAT configuration showing the usage of PAT on full sim samples

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF') 

# note that you can use a bunch of core tools of PAT 
# to taylor your PAT configuration; for a few examples
# uncomment the following lines

from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, 'Muons')
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Photons', 'Electrons', 'Muons', 'Taus'])

process.source = cms.Source("PoolSource",
           fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/MC_Generation/ZInvGamma/PYTHIA6_ZInvgamma_7TeV_STARTUP_RECO.root',
                                             'rfio:/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/MC_Generation/ZInvGamma/PYTHIA6_ZInvgamma_7TeV_STARTUP_RECO_1.root'),
           duplicateCheckMode=cms.untracked.string('noDuplicateCheck')           
           )
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
                              rungenmet           = cms.untracked.bool(True),
                              runPFmet           = cms.untracked.bool(True),
                              runTCmet           = cms.untracked.bool(True),
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
process.maxEvents.input = 5         ##  (e.g. -1 to run on all events)
process.out.outputCommands = ['keep *']  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#process.out.fileName = '/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/3_1_2/Pat.root'            ##  (e.g. 'myTuple.root')
#process.source.fileNames = ['/store/relval/CMSSW_3_1_2/RelValSingleGammaPt35/GEN-SIM-RECO/MC_31X_V3-v1/0006/E6C7ED95-4878-DE11-B082-000423D98BE8.root']
#process.source.fileNames = ['rfio:/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/MC_Generation/ZInvGamma/PYTHIA6_ZInvgamma_7TeV_1E31_RECO.root','rfio:/castor/cern.ch/user/s/sandhya/Monophoton_3_X_X/MC_Generation/ZInvGamma/PYTHIA6_ZInvgamma_7TeV_1E31_RECO_1.root']
process.source.fileNames = ['rfio:/castor/cern.ch/user/b/berzano/ADD/Monophoton_3_X_X/MC/314/monophoton_MD1d2-sherpa-MC_31X_V5-GEN-SIM-RECO_1.root','rfio:/castor/cern.ch/user/b/berzano/ADD/Monophoton_3_X_X/MC/314/monophoton_MD1d2-sherpa-MC_31X_V5-GEN-SIM-RECO_2.root']
#process.source.fileNames = ['rfio:/castor/cern.ch/user/b/berzano/ADD/Monophoton_3_X_X/MC/314/monophoton_MD1d2_pt30-sherpa-MC_31X_V5-GEN-SIM-RECO_1.root','rfio:/castor/cern.ch/user/b/berzano/ADD/Monophoton_3_X_X/MC/314/monophoton_MD1d2_pt30-sherpa-MC_31X_V5-GEN-SIM-RECO_2.root']
process.options.wantSummary = False       ##  (to suppress the long output at the end of the job)    
#process.cleanLayer1Jets.checkOverlaps.photons.requireNoOverlaps = True
