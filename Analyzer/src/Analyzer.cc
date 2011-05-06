
// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc Analysis/Analyzer/src/Analyzer.cc
  
 
 Implementation:
 <Notes on implementation>
*/
//
// Original Author:  Sandhya Jain
//         Created:  Fri Apr 17 11:00:06 CEST 2009
// $Id: Analyzer.cc,v 1.48 2011/05/05 12:38:38 schauhan Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"//(phi1,phi2) does phi1-phi2
#include "DataFormats/Math/interface/deltaR.h"//(eta1,phi1,eta2,phi2)
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

#include <map>
#include <string>

#include "ADDmonophoton/Analyzer/interface/Analyzer.h"

#include <iostream>
#include <fstream>
using namespace std;
using namespace ROOT::Math::VectorUtil ;

//utility function prototypes
float correct_phi(float phi);
float Theta(float eta);
float Pl(float P,float Pt);


//myEcalSeverityLevelAlgo.cc function prototypes ---- for calculating e2e9                                                                            
float recHitE( const DetId id, const EcalRecHitCollection &recHits );
float recHitE( const DetId id, const EcalRecHitCollection & recHits, int dEta, int dPhi );
float recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits );


//
// class decleration
//
class genPho{
public:
  genPho(){};  ~genPho(){};
  float pt,px,py,pz,eta,phi,E,motherPt,motherEta,motherPhi , GrandmotherPt,GrandmotherEta,GrandmotherPhi;
  int motherID,GrandmotherID, status, motherStatus, GrandmotherStatus;
};   

class PtSortCriterium{
public:
  bool operator() (genPho p1,genPho p2){
    return p1.pt >= p2.pt;
  }
};


class PtSortCriterium3{
public:
  bool operator() (reco::Track p1,reco::Track p2 ){
    return p1.pt() > p2.pt();
  }
};

class EnergySortCriterium{
public:
  bool operator() (CrystalInfo p1,CrystalInfo p2 ){
    return p1.energy > p2.energy;
  }
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig):
//histocontainer_(),
  eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  cosMuoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("cosMuonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  pfjetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfjetTag")),
  genjetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("genjetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  PFmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFmetTag")),
  TCmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("TCmetTag")),
  phoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag")),
  cscLabel_(iConfig.getUntrackedParameter<edm::InputTag>("cscTag")),
  rpcLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rpcTag")),
  rechitBLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitBTag")),
  rechitELabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitETag")),
  hcalrechitLabel_(iConfig.getUntrackedParameter<edm::InputTag>("hcalrechitTag")),
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults")),
  Tracks_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks")),
  BeamHaloSummaryLabel_(iConfig.getUntrackedParameter<edm::InputTag>("BeamHaloSummary")),
  Vertices_(iConfig.getUntrackedParameter<edm::InputTag>("Vertices")),
  pileupLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pileup")), 
  outFile_(iConfig.getUntrackedParameter<string>("outFile")),
  rungenParticleCandidates_(iConfig.getUntrackedParameter<bool>("rungenParticleCandidates")),
  runphotons_(iConfig.getUntrackedParameter<bool>("runphotons")),
  runmet_(iConfig.getUntrackedParameter<bool>("runmet")),
  rungenmet_(iConfig.getUntrackedParameter<bool>("rungenmet")),
  runPFmet_(iConfig.getUntrackedParameter<bool>("runPFmet")),
  runTCmet_(iConfig.getUntrackedParameter<bool>("runTCmet")),
  runelectrons_(iConfig.getUntrackedParameter<bool>("runelectrons")),
  runmuons_(iConfig.getUntrackedParameter<bool>("runmuons")),
  runcosmicmuons_(iConfig.getUntrackedParameter<bool>("runcosmicmuons")),
  runjets_(iConfig.getUntrackedParameter<bool>("runjets")),
  runpfjets_(iConfig.getUntrackedParameter<bool>("runpfjets")),
  rungenjets_(iConfig.getUntrackedParameter<bool>("rungenjets")),
  runtaus_(iConfig.getUntrackedParameter<bool>("runtaus")),
  runHLT_(iConfig.getUntrackedParameter<bool>("runHLT")),
  runL1_(iConfig.getUntrackedParameter<bool>("runL1")),
  runscraping_(iConfig.getUntrackedParameter<bool>("runscraping")),
  runtracks_(iConfig.getUntrackedParameter<bool>("runtracks")),
  runrechit_(iConfig.getUntrackedParameter<bool>("runrechit")),
  runHErechit_(iConfig.getUntrackedParameter<bool>("runHErechit")),
  runvertex_(iConfig.getUntrackedParameter<bool>("runvertex")),
  runCSCseg_(iConfig.getUntrackedParameter<bool>("runCSCseg")),
  runBeamHaloSummary_(iConfig.getUntrackedParameter<bool>("runBeamHaloSummary")),
  runRPChit_(iConfig.getUntrackedParameter<bool>("runRPChit")),
  isAOD_(iConfig.getUntrackedParameter<bool>("isAOD")),
  debug_(iConfig.getUntrackedParameter<bool>("debug")),
  init_(false)
{
  //now do what ever initialization is needed
  all_triggers.clear();
  all_triggerprescales.clear();
  all_ifTriggerpassed.clear();
  nevents =0;
  ntriggers = 0;
  n_signal_events = 0; 
  n_Z_events   = 0; 
  n_W_events    = 0 ;
  n_Zelec_events =  0; 
  n_Zmu_events = 0; 
  n_Ztau_events = 0; 
  n_Znunu_events = 0;
  n_Welec_events =  0; 
  n_Wmu_events = 0; 
  n_Wtau_events = 0;
  n_diphoton_events=0;  
  n_SingleHardPhoton_events = 0;
  ngenphotons  = 0; //for every event, how many stable photons 
  nhardphotons = 0; //for every event, how many hard photons 
  Vertex_n         = 0;
  Track_n          = 0;
  Photon_n         = 0;
  HERecHit_subset_n = 0;
  Jet_n            = 0;
  pfJet_n          = 0;
  EBRecHit_size    = 0;
  EERecHit_size    = 0;  
  Electron_n       = 0; 
  Muon_n           = 0;
  CosmicMuon_n	   = 0;
  Tau_n            = 0;
  genJet_n         = 0;
}


Analyzer::~Analyzer(){
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
  f->Close();
  delete f;
}


//
// member functions
//
void 
Analyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){
  if(debug_){
    cout<<"\n----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"BEGIN NEW RUN: "<<iRun.run()<<endl;
  }

   MaxN=200;

  bool changed(true);
  if(hltConfig_.init(iRun,iSetup,"HLT",changed)){  // HLTlabel_ is grabbed from the cfg file
    // if init returns TRUE, initialisation has succeeded!
    if(debug_) std::cout << "Initalizing HLTConfigProvider"  << std::endl;
    if(changed){
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histogr ams or do anything else dependent on the revised HLT config
      if(debug_) cout<< "HLT config has changed wrt the previous Run"  << std::endl;
      photon_triggers_in_run.clear();
      met_triggers_in_run.clear();
      cosmic_triggers_in_run.clear();
      halo_triggers_in_run.clear();
      if(debug_){
        cout<<" Trigger Table : "<<hltConfig_.tableName()<<endl;
        cout<<" Available photon, met, cosmic, and halo triggers : "<<endl;
      }
      unsigned int ntriggers = hltConfig_.size();
      
      // Loop over all available triggers
      for(unsigned int t=0;t<ntriggers;++t){
        std::string hltname(hltConfig_.triggerName(t));
        string string_search ("HLT_Photon");
        string string_search1 ("HLT_MET");
        string string_search_cosmic("cosmic");
        string string_search_Cosmic("Cosmic");
        string string_search_halo("halo");
        string string_search_Halo("Halo");
        
        //search the trigger name for string_search. 
        size_t found = hltname.find(string_search);
        size_t found1 = hltname.find(string_search1);
        size_t found_cosmic = hltname.find(string_search_cosmic);
        size_t found_Cosmic = hltname.find(string_search_Cosmic);
        size_t found_halo = hltname.find(string_search_halo);
        size_t found_Halo = hltname.find(string_search_Halo);
        
        //cout<<"Does "<<hltname<<" contain "<<string_search<<" ? found returns "<<int(found)<<endl;
  
        if(found!=string::npos ){
          if(debug_) cout<<"  "<<t<<" Photon HLT path in this run\t"<<hltname<<endl;
          photon_triggers_in_run.push_back(hltname);
        }
	
        if(found1!=string::npos ){
          if(debug_)cout<<"  "<<t<<" MET HLT path in this run\t"<<hltname<<endl;
            met_triggers_in_run.push_back(hltname);    
        }
        
        if(found_cosmic!=string::npos ){
          if(debug_) cout<<"  "<<t<<" cosmic trigger found in this run\t"<<hltname<<endl;
          cosmic_triggers_in_run.push_back(hltname);
        }
        
        if(found_Cosmic!=string::npos ){
          if(debug_) cout<<"  "<<t<<" Cosmic trigger found in this run\t"<<hltname<<endl;
          cosmic_triggers_in_run.push_back(hltname);
        }
        
        if(found_halo!=string::npos ){
          if(debug_) cout<<"  "<<t<<" halo trigger found in this run\t"<<hltname<<endl;
          halo_triggers_in_run.push_back(hltname);
        }
        
        if(found_Halo!=string::npos ){
          if(debug_) cout<<"  "<<t<<" Halo trigger found in this run\t"<<hltname<<endl;
          halo_triggers_in_run.push_back(hltname);
        }
        
      }//loop over ntriggers

      //SAVE PHOTON TRIGGER INFO
      for(int x = 0; x< (int)photon_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==photon_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found)
          all_triggers.push_back(photon_triggers_in_run[x]);
      }//loop photon triggers
      
      
      //SAVE MET TRIGGER INFO
      for(int x= 0; x< (int)met_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==met_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found) all_triggers.push_back(met_triggers_in_run[x]);
      }//loop met triggers
      
      
      //SAVE COSMIC TRIGGER INFO
      for(int x= 0; x< (int)cosmic_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==cosmic_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found) all_triggers.push_back(cosmic_triggers_in_run[x]);
      }//loop cosmic triggers
      
      
      //SAVE HALO TRIGGER INFO
      for(int x= 0; x< (int)halo_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==halo_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found) all_triggers.push_back(halo_triggers_in_run[x]);
      }//loop halo triggers
      
      
    }//loop over all available triggers
  }else{
      // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
      // with the file and/or code and needs to be investigated!
    std::cout << " HLT config extraction failure with name " << hlTriggerResults_ << std::endl;
      // In this case, all access methods will return empty values!
  }
    
  ntriggers = (int)all_triggers.size();
  if(debug_){
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"the triggers in HLT list till now are:"<<endl;
    for(int i = 0; i< ntriggers;i++) cout<<"\t"<<all_triggers[i]<<endl;
  }  
}




// ------------ method called to for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  if(debug_) cout<<"DEBUG: new event"<<endl;
  using namespace edm;
  using namespace reco;
	
  RunNumber   = iEvent.id().run();
  EventNumber = iEvent.id().event();
  LumiNumber  = iEvent.id().luminosityBlock();
  BXNumber = iEvent.bunchCrossing();

  // get ConditionsInLumiBlock
  const edm::LuminosityBlock& iLumi = iEvent.getLuminosityBlock();
  edm::Handle<edm::ConditionsInLumiBlock> condInLumiBlock;
  iLumi.getByLabel("conditionsInEdm", condInLumiBlock);

  //initialize to -1 in case isValid fails
  totalIntensityBeam1 = -1;
  totalIntensityBeam2 = -1;
  if(condInLumiBlock.isValid()){//check pointer condInLumiBlock is not null
    totalIntensityBeam1 = condInLumiBlock->totalIntensityBeam1;
    totalIntensityBeam2 = condInLumiBlock->totalIntensityBeam2;
  }
	
  // get LumiSummary
  edm::Handle<LumiSummary> lumiSummary;
  iLumi.getByLabel("lumiProducer", lumiSummary);
	
	//initialize to -1 in case isValid fails
  avgInsDelLumi = -1.;
  avgInsDelLumiErr = -1.;
  avgInsRecLumi = -1.;
  avgInsRecLumiErr = -1.;
  if(lumiSummary.isValid()){//check pointer lumiSummary is not null
    if(lumiSummary->isValid()){//data are valid only if run exists from all sources lumi,trg ,hlt
      avgInsDelLumi = lumiSummary->avgInsDelLumi();
      avgInsDelLumiErr = lumiSummary->avgInsDelLumiErr();
      avgInsRecLumi = lumiSummary->avgInsRecLumi();
      avgInsRecLumiErr = lumiSummary->avgInsRecLumiErr();
    }
  }
  if(debug_) cout<<"DEBUG: start saving stuff"<<endl;
  nevents++;
 

 if(runPileUp_){
   //pile up info from MC only 
   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   iEvent.getByLabel(pileupLabel_, PupInfo);
                    
   std::vector<PileupSummaryInfo>::const_iterator PVI;   
    ootnpuVertices =0;
    npuVertices    =0;            
                                                                                 
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing()== 0){npuVertices += PVI->getPU_NumInteractions();
       }
       else{
             ootnpuVertices += PVI->getPU_NumInteractions();
           }

  if(debug_)std::cout << " Pileup Information: bunchXing, nvtx,: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions()<< std::endl;
                    
   }//loop over pileup infor               
  }//Run over Pileup

  
  //getting handle to generator level information
  if( rungenParticleCandidates_ ){
	  
    Handle<GenEventInfoProduct > genEvent;
    iEvent.getByLabel( "generator", genEvent );
    gen_pthat = (genEvent->hasBinningValues() ? (genEvent->binningValues())[0] : 0.0); 
	  
    ngenphotons  = 0;
    nhardphotons = 0;
    is_signal_event = false; 
    is_Z_event   = false; 
    is_W_event    = false; 
    is_Zelec_event  = false; 
    is_Zmu_event = false; 
    is_Ztau_event = false; 
    is_Znunu_event =false  ;
    is_Welec_event  = false; 
    is_Wmu_event = false; 
    is_Wtau_event = false;
    is_SingleHardPhoton_event=false;  
    is_diphoton_event=false;
    is_isr_photon_event = false;
    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles); 
    std::vector<genPho>  mygenphoton_container;
    mygenphoton_container.clear();
    genPho genphoton;
    int ii =0;
    for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++) {
      //getting information from hard scattered Graviton 
      if(debug_) cout<<"start gen particle collection event:"<<nevents<<" particle:"<<ii<<endl;
      if (genparticle->pdgId()==39 && genparticle->status()==3) { 
	is_signal_event = true;
	n_signal_events++;
	gen_graviton_pt  = genparticle->pt();
	gen_graviton_px  = genparticle->px();
	gen_graviton_py  = genparticle->py();
	gen_graviton_pz  = genparticle->pz();
	gen_graviton_phi = correct_phi(genparticle->phi());
	gen_graviton_eta = genparticle->eta();
	gen_graviton_E   = genparticle->energy();
      }
      //getting information from Z
      if (genparticle->pdgId()==23 && genparticle->status()==3){ 
	is_Z_event = true;
	n_Z_events++;
	//if(debug_) cout<<" getting information from Z now"<<endl;
	gen_Zboson_pt  = genparticle->pt();
	gen_Zboson_px  = genparticle->px();
	gen_Zboson_py  = genparticle->py();
	gen_Zboson_pz  = genparticle->pz();
	gen_Zboson_phi = correct_phi(genparticle->phi());
	gen_Zboson_eta = genparticle->eta();
	gen_Zboson_E   = genparticle->energy();
	int daughters  = genparticle->numberOfDaughters();
	int iDaughter=0;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  //cout<<"genparticle->daughter(i)"<<genparticle->daughter(i)<<std::endl;
	  if(abs(daughter->pdgId())==12||abs(daughter->pdgId())==14||abs(daughter->pdgId())==16){is_Znunu_event=true; n_Znunu_events++;}
	  if(daughter->pdgId()==11) { is_Zelec_event=true; n_Zelec_events++;}
	  if(daughter->pdgId()==13) { is_Zmu_event=true  ; n_Zmu_events++  ;}
	  if(daughter->pdgId()==15) { is_Ztau_event=true ; n_Ztau_events++  ;}
	  if(daughter->pdgId()!=23) {
	    gen_Zdaughter_pt[iDaughter]    = daughter->pt();
	    gen_Zdaughter_px[iDaughter]    = daughter->px();
	    gen_Zdaughter_py[iDaughter]    = daughter->py();
	    gen_Zdaughter_pz[iDaughter]    = daughter->pz();
	    gen_Zdaughter_phi[iDaughter]   = correct_phi(daughter->phi());
	    gen_Zdaughter_eta[iDaughter]   = daughter->eta();
	    gen_Zdaughter_E[iDaughter]     = daughter->energy();
	    gen_Zdaughter_charge[iDaughter]= daughter->charge();
	    gen_Zdaughter_ID[iDaughter]    = daughter->pdgId();
	    iDaughter++;
	  }//end of if(daughter->pdgId()!=23)
	}//end of for loop for daughters
      }//end for loop of Z information    
      //getting information from W+/W-

      if (abs(genparticle->pdgId())==24 && genparticle->status()==3) { 
	is_W_event = true;
	if(debug_) cout<<" W motherID:" << genparticle->mother()->pdgId()<<endl; 
	n_W_events++;
	gen_Wboson_pt      = genparticle->pt();
	gen_Wboson_px      = genparticle->px();
	gen_Wboson_py      = genparticle->py();
	gen_Wboson_pz      = genparticle->pz();
	gen_Wboson_phi     = correct_phi(genparticle->phi());
	gen_Wboson_eta     = genparticle->eta();
	gen_Wboson_E       = genparticle->energy();
	gen_Wboson_charge  = genparticle->charge();
	gen_Wboson_ID      = genparticle->pdgId();
	int daughters      = genparticle->numberOfDaughters();
	int iDaughter =0;
	if(debug_) cout<<"W pt:"<<genparticle->pt()<<endl;
	if(debug_) cout<<"W daughters:"<<endl;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  if(abs(daughter->pdgId())==11) {is_Welec_event=true; n_Welec_events++;}
	  if(abs(daughter->pdgId())==13) {is_Wmu_event=true  ; n_Wmu_events++  ;}
	  if(abs(daughter->pdgId())==15) {is_Wtau_event=true ; n_Wtau_events++ ;}  
	  if(debug_) cout<<"ID, Status,Pt:"<<abs(daughter->pdgId())<<"   "<<daughter->status()<<"   "<<daughter->pt()<<endl;
	  //getting leptons decaying from W
	  if(abs(daughter->pdgId())!=24) {
	    gen_Wdaughter_pt[iDaughter]     = daughter->pt();
	    gen_Wdaughter_px[iDaughter]     = daughter->px();
	    gen_Wdaughter_py[iDaughter]     = daughter->py();
	    gen_Wdaughter_pz[iDaughter]     = daughter->pz();
	    gen_Wdaughter_phi[iDaughter]    = correct_phi(daughter->phi());
	    gen_Wdaughter_eta[iDaughter]    = daughter->eta();
	    gen_Wdaughter_E[iDaughter]      = daughter->energy();
	    gen_Wdaughter_charge[iDaughter] = daughter->charge();
	    gen_Wdaughter_ID[iDaughter]     = daughter->pdgId();
	    iDaughter++;
	  }//if(abs(daughter->pdgId())!=24)
	}//end of for loop for daughters
      }//end for loop of W information

		//getting info from decay of muons 
      if( abs(genparticle->pdgId())==13 ){
	if(debug_) cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	int daughters   = genparticle->numberOfDaughters();
	int iDaughter=0;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  if(debug_) cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
	  gen_Muon_ID[iDaughter]     = genparticle->pdgId();
	  gen_Muon_Status[iDaughter] = genparticle->status();
	  gen_Muon_Pt[iDaughter]     = genparticle->pt();
	  gen_MuonDaughter_pt[iDaughter]           = daughter->pt();
	  gen_MuonDaughter_px[iDaughter]           = daughter->px();
	  gen_MuonDaughter_py[iDaughter]           = daughter->py();
	  gen_MuonDaughter_pz[iDaughter]           = daughter->pz();
	  gen_MuonDaughter_phi[iDaughter]          = correct_phi(daughter->phi());
	  gen_MuonDaughter_eta[iDaughter]          = daughter->eta();
	  gen_MuonDaughter_E[iDaughter]            = daughter->energy();
	  gen_MuonDaughter_ID[iDaughter]           = daughter->pdgId();
	  gen_MuonDaughter_status[iDaughter]       = daughter->status();
	  gen_MuonDaughter_charge[iDaughter]       = daughter->charge();
	  iDaughter++;
	}//for(int i = 0;i<daughters;i++)
      }//if((abs(genparticle->pdgId())==13)

		//getting info from decay of taus 
      if( abs(genparticle->pdgId())==15 ){
	if(debug_) cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	int daughters   = genparticle->numberOfDaughters();
	int iDaughter=0;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  //if(debug_) cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
	  gen_tau_ID[iDaughter]           = genparticle->pdgId();
	  gen_tau_Status[iDaughter]       = genparticle->status();
	  gen_tau_Pt[iDaughter]           = genparticle->pt();
	  gen_tauDaughter_pt[iDaughter]           = daughter->pt();
	  gen_tauDaughter_px[iDaughter]           = daughter->px();
	  gen_tauDaughter_py[iDaughter]           = daughter->py();
	  gen_tauDaughter_pz[iDaughter]           = daughter->pz();
	  gen_tauDaughter_phi[iDaughter]          = correct_phi(daughter->phi());
	  gen_tauDaughter_eta[iDaughter]          = daughter->eta();
	  gen_tauDaughter_E[iDaughter]            = daughter->energy();
	  gen_tauDaughter_ID[iDaughter]           = daughter->pdgId();
	  gen_tauDaughter_status[iDaughter]       = daughter->status();
	  gen_tauDaughter_charge[iDaughter]       = daughter->charge();
	  iDaughter++;
	}//for(int i = 0;i<daughters;i++)
      }//if((abs(genparticle->pdgId())==15)
      if(debug_) cout<<"DEBUG 6"<<endl;
      if(debug_) cout<<"getting gen photon information now"<<endl;
      //getting information from all photons
      //if (genparticle->pdgId()==22 && genparticle->status()==1 && genparticle->pt()>5.)
      if (genparticle->pdgId()==22){ 
	//doing it this way, as I want to sort it in pt after filling everything in the container
        if(genparticle->mother()!=NULL){
	  const reco::Candidate *mom = genparticle->mother();
	  genphoton.motherID = mom->pdgId();
	  genphoton.motherStatus = mom->status();
	  if(debug_) cout<<"DEBUG mom pdgid is "<<mom->pdgId()<<" status is "<<mom->status()<<" pt "<<mom->pt()<< " e "<<mom->energy() <<endl;
	  genphoton.motherPt = mom->pt();
	  genphoton.motherEta = mom->eta();
	  genphoton.motherPhi = correct_phi(mom->phi());
	  if(genparticle->status()==1 && mom->pdgId()==2212 && mom->mother()==NULL)
	    is_isr_photon_event = true;
	  if(mom->mother()!=NULL){	  
	    if(debug_) cout<<"DEBUG gramma is "<<mom->mother()->pdgId()<<endl;
	    genphoton.GrandmotherID = mom->mother()->pdgId();
	    genphoton.GrandmotherStatus = mom->mother()->status();
	    genphoton.GrandmotherPt = mom->mother()->pt();
	    genphoton.GrandmotherEta = mom->mother()->eta();
	    genphoton.GrandmotherPhi = correct_phi(mom->mother()->phi());
            if(genparticle->status()==1 && genparticle->mother()->mother()->pdgId()==2212)
              is_isr_photon_event = true;
          } else{//if grandma is null but mom ok, initialize grandma stuff
	    genphoton.GrandmotherID = 0;
	    genphoton.GrandmotherStatus = -99;
	    genphoton.GrandmotherPt = -99.;
	    genphoton.GrandmotherEta = -99.;
	    genphoton.GrandmotherPhi = -99.;
	  }//end grandma else
        } else {//if mother is null, initialize mom and grandma stuff
	  genphoton.motherID = 0;
	  genphoton.motherStatus = -99;
	  genphoton.motherPt = -99.;
	  genphoton.motherEta = -99.;
	  genphoton.motherPhi = -99.;
	  genphoton.GrandmotherID = 0;
	  genphoton.GrandmotherStatus = -99;
	  genphoton.GrandmotherPt = -99.;
	  genphoton.GrandmotherEta = -99.;
	  genphoton.GrandmotherPhi = -99.;
	}
	genphoton.pt  = genparticle->pt();
	genphoton.px  = genparticle->px();
	genphoton.py  = genparticle->py();
	genphoton.pz  = genparticle->pz();
	genphoton.phi = correct_phi(genparticle->phi());
	genphoton.eta = genparticle->eta();
	genphoton.E   = genparticle->energy();
	genphoton.status = genparticle->status();
	mygenphoton_container.push_back(genphoton); ngenphotons++;
	if(debug_) cout<<"DEBUG done with pho"<<endl;
      }//end of if (genparticle->pdgId()==22 && genparticle->status()==1)
      if(debug_) cout<<"DEBUG before checking mothers event:"<< nevents<<" particle:"<<ii<<endl;  
      if(genparticle->mother()!=NULL && genparticle->mother()->mother()!=NULL)
        if(genparticle->pdgId()==22 && genparticle->status()==1 && genparticle->mother()->mother()->pdgId()==2212)
          is_isr_photon_event = true;
      if (genparticle->pdgId()==22 && genparticle->status()==3) {
	gen_Hpho_pt[nhardphotons]  = genparticle->pt();
	gen_Hpho_px[nhardphotons]  = genparticle->px();
	gen_Hpho_py[nhardphotons]  = genparticle->py();
	gen_Hpho_pz[nhardphotons]  = genparticle->pz();
	gen_Hpho_phi[nhardphotons] = correct_phi(genparticle->phi());
	gen_Hpho_eta[nhardphotons] = genparticle->eta();
	gen_Hpho_E[nhardphotons]   = genparticle->energy();
	nhardphotons++;
      }//end of if (genparticle->pdgId()==22 && genparticle->status()==3)
      if(debug_) cout<<"DEBUG end loop over gen particles"<<endl;
      ii++;
    }//end of for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++)
    if(debug_) cout<<"DEBUG: finish gen"<<endl;	
    if (nhardphotons==1) {
      n_SingleHardPhoton_events++;
      is_SingleHardPhoton_event=true;
    }
    if (nhardphotons==2){
      n_diphoton_events++; 
      is_diphoton_event=true;
    }
		
    if(mygenphoton_container.size()!=0){
      std::sort(mygenphoton_container.begin(),mygenphoton_container.end(),PtSortCriterium());
      for(unsigned int x=0;x < mygenphoton_container.size(); x++){
	//std::if(debug_) cout<<"genphoton motherID:"<<mygenphoton_container[x].motherID<<endl;
	gen_pho_motherPt[x]   = mygenphoton_container[x].motherPt;
	gen_pho_motherEta[x]  = mygenphoton_container[x].motherEta;
	gen_pho_motherPhi[x]  = mygenphoton_container[x].motherPhi;
	gen_pho_motherID[x]   = mygenphoton_container[x].motherID;
	gen_pho_motherStatus[x]   = mygenphoton_container[x].motherStatus;
	gen_pho_GrandmotherPt[x]   = mygenphoton_container[x].GrandmotherPt;
	gen_pho_GrandmotherEta[x]  = mygenphoton_container[x].GrandmotherEta;
	gen_pho_GrandmotherPhi[x]  = mygenphoton_container[x].GrandmotherPhi;
	gen_pho_GrandmotherID[x]   = mygenphoton_container[x].GrandmotherID;
	gen_pho_GrandmotherStatus[x]   = mygenphoton_container[x].GrandmotherStatus;
	gen_pho_pt[x]         = mygenphoton_container[x].pt;
	gen_pho_px[x]         = mygenphoton_container[x].px;
	gen_pho_py[x]         = mygenphoton_container[x].py;
	gen_pho_pz[x]         = mygenphoton_container[x].pz;
	gen_pho_phi[x]        = mygenphoton_container[x].phi;
	gen_pho_eta[x]        = mygenphoton_container[x].eta;
	gen_pho_E[x]          = mygenphoton_container[x].E;
	gen_pho_status[x]          = mygenphoton_container[x].status;
	/*if( gen_pho_pt[x] > 100. && (abs(gen_pho_motherID[x])<=6||abs(gen_pho_motherID[x])==11||abs(gen_pho_motherID[x])==9|| abs(gen_pho_motherID[x])==21 )) 
	  if(debug_) cout<<"[x]:"<<x<< "pho_pt:" << gen_pho_pt[x]<<" pho_status:"<< gen_pho_status[x]<<" motherID:" << gen_pho_motherID[x]<<" motherStatus: "<< gen_pho_motherStatus[x]<< " GrandmotherID: "<< gen_pho_GrandmotherID[x]<< " GrandmotherStatus:"<< gen_pho_GrandmotherStatus[x]<< endl;
        */
	// if(debug_) cout<<"got the photon info right"<<endl;
      }//end of for loop
    }//end of if((mygenphoton_container.size()!=0)
    //if(debug_) cout<<"mygenphoton_container loop ended"<<std::endl; 
  }//end of if(rungenParticleCandidates_)
  if(debug_) cout<<"DEBUG: finish saving gen"<<endl;
  
  
  ///// L1
  if(runL1_){
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
    const L1GtTriggerMenu* menu = menuRcd.product();
    
    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
    const DecisionWord dWord = gtRecord->decisionWord();  // this will get the decision word *before* masking disabled bits
    
    //Initialize and fill the L1 flags
    for(map<string,int>::iterator iter = L1_chosen.begin(); iter != L1_chosen.end();iter++){
      iter->second = 0;
      string L1name =  iter->first;
      if ( menu->gtAlgorithmResult( L1name, dWord) ) L1_chosen[L1name] = 1;
      else L1_chosen[ L1name ]=0;
    } 
  }  
 
  if(runHLT_){
    Handle<TriggerResults> HLTR;
    iEvent.getByLabel(hlTriggerResults_,HLTR);
    if (HLTR.isValid()) {
      all_triggerprescales.clear();
      all_ifTriggerpassed.clear();
      const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*HLTR);
      vector<int> idx;
      for(int i = 0; i< ntriggers;i++){
        all_triggerprescales.push_back(0);
        all_ifTriggerpassed.push_back(0);
        idx.push_back(triggerNames_.triggerIndex(all_triggers[i]));
      }
      Int_t hsize = Int_t(HLTR->size());
      for(int i=0;i<ntriggers;i++){
        if(idx[i] < hsize){
          all_ifTriggerpassed[i]=HLTR->accept(idx[i]);
          all_triggerprescales[i]=hltConfig_.prescaleValue( iEvent, iSetup, all_triggers[i]);
	}
      }
      if(debug_){
        for(int i=0;i<ntriggers;i++){
          cout<<"prescale for "<<all_triggers[i]<<" is: "<< all_triggerprescales[i]<<endl;
          cout<<"if triggger passed for "<<all_triggers[i]<<" : "<<all_ifTriggerpassed[i]<<endl;
        }//loop over ntriggers
      }//debug
    }//if HLTR is Valid
  }//runHLT_
   
  if(runvertex_){
    Handle<reco::VertexCollection> recVtxs;
     iEvent.getByLabel(Vertices_, recVtxs);
     vector<Vertex> my_vertices;
	   
     my_vertices.clear();
     for(reco::VertexCollection::const_iterator v=recVtxs->begin();v!=recVtxs->end(); ++v){
       my_vertices.push_back(*v);
     }
     Vertex_n = 0;  
     for (unsigned int y = 0; y < min(my_vertices.size(),MaxN);y++){
       vx[y]     = my_vertices[y].x();
       vy[y]     = my_vertices[y].y();
       vz[y]     = my_vertices[y].z();
       chi2[y]   = my_vertices[y].chi2();
       vtracksize[y] = my_vertices[y].tracksSize();
       vndof[y] = my_vertices[y].ndof();
       v_isFake[y] = my_vertices[y].isFake();
       v_d0[y] = my_vertices[y].position().rho();
       Vertex_n++;
     }
   }
   
  if(runscraping_){ // taken from DPGAnalysis/Skims/src/FilterOutScraping.cc (CMSSW_3_8_3)
    Scraping_isScrapingEvent = true;
    Scraping_fractionOfGoodTracks = 0;
    Scraping_numOfTracks=0;
 
    // get GeneralTracks collection
    edm::Handle<reco::TrackCollection> tkRef;
    iEvent.getByLabel(Tracks_,tkRef);    
    const reco::TrackCollection* tkColl = tkRef.product();

    int numhighpurity=0;
    reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");
    Scraping_numOfTracks = tkColl->size();

    if(tkColl->size()>10){
      reco::TrackCollection::const_iterator itk = tkColl->begin();
      reco::TrackCollection::const_iterator itk_e = tkColl->end();
      for(;itk!=itk_e;++itk){
        if(itk->quality(_trackQuality)) numhighpurity++;
      }
      Scraping_fractionOfGoodTracks = (float)numhighpurity/(float)tkColl->size();
		
      if(Scraping_fractionOfGoodTracks>0.25) Scraping_isScrapingEvent = false;
    }else{
      //if less than 10 Tracks, mark as not Scraping    
      Scraping_isScrapingEvent = false;
    }
  }	
	
   if(runtracks_){
     Handle<reco::TrackCollection> tracks;
     iEvent.getByLabel(Tracks_,tracks);
     std::vector<reco::Track>  myTrack_container;
     myTrack_container.clear();
     for(reco::TrackCollection::const_iterator Track_iter = tracks->begin();
	 Track_iter != tracks->end(); ++Track_iter) {
       if(Track_iter->pt()>5.){
	 myTrack_container.push_back(*Track_iter);
       }
     }
     
     
     if(myTrack_container.size()>1)
       std::sort(myTrack_container.begin(),myTrack_container.end(),PtSortCriterium3());
     Track_n = 0;
     for(unsigned int x=0;x < min(myTrack_container.size(),MaxN); x++){
       trk_pt[x]  = myTrack_container[x].pt();
       trk_px[x]  = myTrack_container[x].px();
       trk_py[x]  = myTrack_container[x].py();
       trk_pz[x]  = myTrack_container[x].pz();
       trk_vx[x] = myTrack_container[x].vx();
       trk_vy[x] = myTrack_container[x].vy();
       trk_vz[x] = myTrack_container[x].vz();
       trk_phi[x] = correct_phi(myTrack_container[x].phi());
       trk_eta[x] = myTrack_container[x].eta();
       Track_n++;
     }//end of for loop
   }//end of if(runtracks_)

   if(runmuons_){
     edm::Handle<edm::View<pat::Muon> > muonHandle;
     iEvent.getByLabel(muoLabel_,muonHandle);
     vector <pat::Muon> mymuon_container;
     
     const edm::View<pat::Muon> & muons = *muonHandle;   // const ... &, we don't make a copy of it!
     for(edm::View<pat::Muon>::const_iterator muon = muons.begin(); muon!=muons.end(); ++muon){
       //if(muon->isGlobalMuon())
       mymuon_container.push_back(*muon);
     }
     Muon_n = 0;
     for(unsigned int x=0;x < min(mymuon_container.size(),MaxN); x++){
       muon_pt[x]  = mymuon_container[x].pt();
       muon_energy[x]  = mymuon_container[x].energy();
       muon_px[x]  = mymuon_container[x].px();
       muon_py[x]  = mymuon_container[x].py();
       muon_pz[x]  = mymuon_container[x].pz();
       muon_phi[x] = correct_phi(mymuon_container[x].phi());
       muon_eta[x] = mymuon_container[x].eta();
       muon_charge[x] = mymuon_container[x].charge();
       muon_vx[x]  = mymuon_container[x].vx();
       muon_vy[x]  = mymuon_container[x].vy();
       muon_vz[x]  = mymuon_container[x].vz();
       
       //tia's stuff
       muon_isGlobalMuon[x] =     mymuon_container[x].isGlobalMuon();
       muon_isTrackerMuon[x] =    mymuon_container[x].isTrackerMuon();
       muon_isStandAloneMuon[x] = mymuon_container[x].isStandAloneMuon();

       muon_OuterTrack_InnerPoint_x[x] = 0.;
       muon_OuterTrack_InnerPoint_y[x] = 0.;
       muon_OuterTrack_InnerPoint_z[x] = 0.;
       muon_OuterTrack_InnerPoint_px[x] = 0.;
       muon_OuterTrack_InnerPoint_py[x] = 0.;
       muon_OuterTrack_InnerPoint_pz[x] = 0.;
       muon_OuterTrack_OuterPoint_x[x] = 0.;
       muon_OuterTrack_OuterPoint_y[x] = 0.;
       muon_OuterTrack_OuterPoint_z[x] = 0.;
       muon_OuterTrack_OuterPoint_px[x] = 0.;
       muon_OuterTrack_OuterPoint_py[x] = 0.;
       muon_OuterTrack_OuterPoint_pz[x] = 0.;
       muon_InnerTrack_InnerPoint_x[x] = 0.;
       muon_InnerTrack_InnerPoint_y[x] = 0.;
       muon_InnerTrack_InnerPoint_z[x] = 0.;
       muon_InnerTrack_InnerPoint_px[x] = 0.;
       muon_InnerTrack_InnerPoint_py[x] = 0.;
       muon_InnerTrack_InnerPoint_pz[x] = 0.;
       muon_InnerTrack_OuterPoint_x[x] = 0.;
       muon_InnerTrack_OuterPoint_y[x] = 0.;
       muon_InnerTrack_OuterPoint_z[x] = 0.;
       muon_InnerTrack_OuterPoint_px[x] = 0.;
       muon_InnerTrack_OuterPoint_py[x] = 0.;
       muon_InnerTrack_OuterPoint_pz[x] = 0.;

   if(isAOD_){
   //FIXME: Seems needed  top and bottom referene point, but not sure how to do that !!!
      reco::TrackRef moTrkref;
     if((muon_isGlobalMuon[x]) || (muon_isTrackerMuon[x])){
        moTrkref = mymuon_container[x].innerTrack();
        muon_InnerTrack_isNonnull[x] = mymuon_container[x].innerTrack().isNonnull();
        }
        else{   
            moTrkref = mymuon_container[x].outerTrack();
            muon_OuterTrack_isNonnull[x] =   mymuon_container[x].outerTrack().isNonnull();
         }
    //For StandAlone
    muon_OuterPoint_x[x]=0.0;
    muon_OuterPoint_y[x]=0.0;
    muon_OuterPoint_z[x]=0.0;
    //For Global,Tracker
    muon_InnerPoint_x[x]=0.0;
    muon_InnerPoint_y[x]=0.0;
    muon_InnerPoint_z[x]=0.0;

  if((muon_OuterTrack_isNonnull[x])){//stand-alone
    muon_OuterPoint_x[x]= moTrkref->referencePoint().x();
    muon_OuterPoint_y[x]= moTrkref->referencePoint().y();
    muon_OuterPoint_z[x]= moTrkref->referencePoint().z();
    }
   if(muon_InnerTrack_isNonnull[x]){//global,tracker
    muon_InnerPoint_x[x]= moTrkref->referencePoint().x();
    muon_InnerPoint_y[x]= moTrkref->referencePoint().y();
    muon_InnerPoint_z[x]= moTrkref->referencePoint().z();
    }
   }//if(AOD_)
     

if(!isAOD_){
       muon_InnerTrack_isNonnull[x] =   mymuon_container[x].innerTrack().isNonnull();
       muon_OuterTrack_isNonnull[x] =   mymuon_container[x].outerTrack().isNonnull();
       
			
       if(mymuon_container[x].innerTrack().isNonnull()){
	 muon_InnerTrack_InnerPoint_x[x] = mymuon_container[x].innerTrack()->innerPosition().x();
	 muon_InnerTrack_InnerPoint_y[x] = mymuon_container[x].innerTrack()->innerPosition().y();
	 muon_InnerTrack_InnerPoint_z[x] = mymuon_container[x].innerTrack()->innerPosition().z();
	 muon_InnerTrack_InnerPoint_px[x] = mymuon_container[x].innerTrack()->innerMomentum().x();
	 muon_InnerTrack_InnerPoint_py[x] = mymuon_container[x].innerTrack()->innerMomentum().y();
	 muon_InnerTrack_InnerPoint_pz[x] = mymuon_container[x].innerTrack()->innerMomentum().z();
	 muon_InnerTrack_OuterPoint_x[x] = mymuon_container[x].innerTrack()->outerPosition().x();
	 muon_InnerTrack_OuterPoint_y[x] = mymuon_container[x].innerTrack()->outerPosition().y();
	 muon_InnerTrack_OuterPoint_z[x] = mymuon_container[x].innerTrack()->outerPosition().z();
	 muon_InnerTrack_OuterPoint_px[x] = mymuon_container[x].innerTrack()->outerMomentum().x();
	 muon_InnerTrack_OuterPoint_py[x] = mymuon_container[x].innerTrack()->outerMomentum().y();
	 muon_InnerTrack_OuterPoint_pz[x] = mymuon_container[x].innerTrack()->outerMomentum().z();
       }
       if(mymuon_container[x].outerTrack().isNonnull()){
	 muon_OuterTrack_InnerPoint_x[x] = mymuon_container[x].outerTrack()->innerPosition().x();
	 muon_OuterTrack_InnerPoint_y[x] = mymuon_container[x].outerTrack()->innerPosition().y();
	 muon_OuterTrack_InnerPoint_z[x] = mymuon_container[x].outerTrack()->innerPosition().z();
	 muon_OuterTrack_InnerPoint_px[x] = mymuon_container[x].outerTrack()->innerMomentum().x();
	 muon_OuterTrack_InnerPoint_py[x] = mymuon_container[x].outerTrack()->innerMomentum().y();
	 muon_OuterTrack_InnerPoint_pz[x] = mymuon_container[x].outerTrack()->innerMomentum().z();
	 muon_OuterTrack_OuterPoint_x[x] = mymuon_container[x].outerTrack()->outerPosition().x();
	 muon_OuterTrack_OuterPoint_y[x] = mymuon_container[x].outerTrack()->outerPosition().y();
	 muon_OuterTrack_OuterPoint_z[x] = mymuon_container[x].outerTrack()->outerPosition().z();
	 muon_OuterTrack_OuterPoint_px[x] = mymuon_container[x].outerTrack()->outerMomentum().x();
	 muon_OuterTrack_OuterPoint_py[x] = mymuon_container[x].outerTrack()->outerMomentum().y();
	 muon_OuterTrack_OuterPoint_pz[x] = mymuon_container[x].outerTrack()->outerMomentum().z();
       }
 }//if(!isAOD_)
       Muon_n++;
     }//end of for loop
   }// if runmuons_
   
   if(runcosmicmuons_){
     //cosmic muon
     edm::Handle<reco::MuonCollection>cosmicMuonHandle;
     iEvent.getByLabel(cosMuoLabel_,cosmicMuonHandle);
     const reco::MuonCollection & cosmicmuons = *cosmicMuonHandle;
     vector <reco::Muon> mycosmicmuon_container;
     
     for(reco::MuonCollection::const_iterator cosmuon = cosmicmuons.begin(); cosmuon!=cosmicmuons.end(); ++cosmuon){
       mycosmicmuon_container.push_back(*cosmuon);
     }

     CosmicMuon_n = 0;
     for(unsigned int x=0;x < min(mycosmicmuon_container.size(),MaxN);x++){
 
       cosmicmuon_pt[x]  = mycosmicmuon_container[x].pt();
       cosmicmuon_energy[x]  = mycosmicmuon_container[x].energy();
       cosmicmuon_px[x]  = mycosmicmuon_container[x].px();
       cosmicmuon_py[x]  = mycosmicmuon_container[x].py();
       cosmicmuon_pz[x]  = mycosmicmuon_container[x].pz();
       cosmicmuon_phi[x] = correct_phi(mycosmicmuon_container[x].phi());
       cosmicmuon_eta[x] = mycosmicmuon_container[x].eta();
       cosmicmuon_charge[x] = mycosmicmuon_container[x].charge();
       
       //tia's stuff
       cosmicmuon_isGlobalMuon[x] =     mycosmicmuon_container[x].isGlobalMuon();
       cosmicmuon_isTrackerMuon[x] =    mycosmicmuon_container[x].isTrackerMuon();
       cosmicmuon_isStandAloneMuon[x] = mycosmicmuon_container[x].isStandAloneMuon();
       
       cosmicmuon_OuterTrack_InnerPoint_x[x] = 0.;
       cosmicmuon_OuterTrack_InnerPoint_y[x] = 0.;
       cosmicmuon_OuterTrack_InnerPoint_z[x] = 0.;
       cosmicmuon_OuterTrack_InnerPoint_px[x] = 0.;
       cosmicmuon_OuterTrack_InnerPoint_py[x] = 0.;
       cosmicmuon_OuterTrack_InnerPoint_pz[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_x[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_y[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_z[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_px[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_py[x] = 0.;
       cosmicmuon_OuterTrack_OuterPoint_pz[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_x[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_y[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_z[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_px[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_py[x] = 0.;
       cosmicmuon_InnerTrack_InnerPoint_pz[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_x[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_y[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_z[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_px[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_py[x] = 0.;
       cosmicmuon_InnerTrack_OuterPoint_pz[x] = 0.;


   if(isAOD_){
    reco::TrackRef cmoTrkref = mycosmicmuon_container[x].outerTrack();
    cosmicmuon_OuterTrack_isNonnull[x] = mycosmicmuon_container[x].outerTrack().isNonnull();    
 
    cosmicmuon_OuterPoint_x[x]=0.0;
    cosmicmuon_OuterPoint_y[x]=0.0;
    cosmicmuon_OuterPoint_z[x]=0.0;
    // standalone muon variables
    if(cosmicmuon_OuterTrack_isNonnull[x]){
    cosmicmuon_OuterPoint_x[x]= cmoTrkref->referencePoint().x();
    cosmicmuon_OuterPoint_y[x]= cmoTrkref->referencePoint().y();
    cosmicmuon_OuterPoint_z[x]= cmoTrkref->referencePoint().z();}
   }//if(AOD_)
      
  if(!isAOD_){ 
       cosmicmuon_InnerTrack_isNonnull[x] =   mycosmicmuon_container[x].innerTrack().isNonnull();
       cosmicmuon_OuterTrack_isNonnull[x] =   mycosmicmuon_container[x].outerTrack().isNonnull();
       
       if(mycosmicmuon_container[x].innerTrack().isNonnull()){
	 cosmicmuon_InnerTrack_InnerPoint_x[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().x();
	 cosmicmuon_InnerTrack_InnerPoint_y[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().y();
	 cosmicmuon_InnerTrack_InnerPoint_z[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().z();
	 cosmicmuon_InnerTrack_InnerPoint_px[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().x();
	 cosmicmuon_InnerTrack_InnerPoint_py[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().y();
	 cosmicmuon_InnerTrack_InnerPoint_pz[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().z();
	 cosmicmuon_InnerTrack_OuterPoint_x[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().x();
	 cosmicmuon_InnerTrack_OuterPoint_y[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().y();
	 cosmicmuon_InnerTrack_OuterPoint_z[x] = mycosmicmuon_container[x].innerTrack()->outerPosition().z();
	 cosmicmuon_InnerTrack_OuterPoint_px[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().x();
	 cosmicmuon_InnerTrack_OuterPoint_py[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().y();
	 cosmicmuon_InnerTrack_OuterPoint_pz[x] = mycosmicmuon_container[x].innerTrack()->outerMomentum().z();
       }
       if(mycosmicmuon_container[x].outerTrack().isNonnull()){
	 cosmicmuon_OuterTrack_InnerPoint_x[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().x();
	 cosmicmuon_OuterTrack_InnerPoint_y[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().y();
	 cosmicmuon_OuterTrack_InnerPoint_z[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().z();
	 cosmicmuon_OuterTrack_InnerPoint_px[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().x();
	 cosmicmuon_OuterTrack_InnerPoint_py[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().y();
	 cosmicmuon_OuterTrack_InnerPoint_pz[x] = mycosmicmuon_container[x].outerTrack()->innerMomentum().z();
	 cosmicmuon_OuterTrack_OuterPoint_x[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().x();
	 cosmicmuon_OuterTrack_OuterPoint_y[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().y();
	 cosmicmuon_OuterTrack_OuterPoint_z[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().z();
	 cosmicmuon_OuterTrack_OuterPoint_px[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().x();
	 cosmicmuon_OuterTrack_OuterPoint_py[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().y();
	 cosmicmuon_OuterTrack_OuterPoint_pz[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().z();
       }
  }//if(!isAOD_)
       CosmicMuon_n++;
     }//end of for loop
  }//if runcosmicmuons_
   

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloGeometry* caloGeom = geoHandle.product();                                                                                                                         
                            
 
   // declare outside of runphotons_ so that we can use this container inside other "bools_"
   std::vector<pat::Photon> myphoton_container;
   myphoton_container.clear();
   if(runphotons_){
     edm::Handle<edm::View<pat::Photon> > phoHandle;
     iEvent.getByLabel(phoLabel_,phoHandle);
     edm::View<pat::Photon>::const_iterator photon;
     
     set<DetId> HERecHitSet;
     HERecHit_subset_n = 0;
     
     for(photon = phoHandle->begin();photon!=phoHandle->end();++photon){
       myphoton_container.push_back(*photon) ;
     }
     Photon_n = 0;
     if(myphoton_container.size()!=0){
       for(unsigned int x=0; x < min(myphoton_container.size(), MaxN);x++){
	 pho_E[x]                     =  myphoton_container[x].energy();
	 pho_pt[x]                    =  myphoton_container[x].pt();
	 pho_px[x]                    =  myphoton_container[x].px();
	 pho_py[x]                    =  myphoton_container[x].py();
	 pho_pz[x]                    =  myphoton_container[x].pz();
	 pho_vx[x]                    =  myphoton_container[x].vx();
	 pho_vy[x]                    =  myphoton_container[x].vy();
	 pho_vz[x]                    =  myphoton_container[x].vz();
	 pho_et[x]                    =  myphoton_container[x].et();
	 pho_eta[x]                   =  myphoton_container[x].eta();
	 pho_phi[x]                   =  correct_phi(myphoton_container[x].phi());
	 pho_theta[x]                 =  myphoton_container[x].theta();
	 pho_r9[x]                    =  myphoton_container[x].r9();
	 pho_e1x5[x]                  =  myphoton_container[x].e1x5();
	 pho_e2x5[x]                  =  myphoton_container[x].e2x5();
	 pho_e3x3[x]                  =  myphoton_container[x].e3x3();
	 pho_e5x5[x]                  =  myphoton_container[x].e5x5();
	 pho_maxEnergyXtal[x]         =  myphoton_container[x].maxEnergyXtal();
	 pho_SigmaEtaEta[x]           =  myphoton_container[x].sigmaEtaEta();        
	 pho_SigmaIetaIeta[x]         =  myphoton_container[x].sigmaIetaIeta();        
	 pho_r1x5[x]                  =  myphoton_container[x].r1x5();
	 pho_r2x5[x]                  =  myphoton_container[x].r2x5();
	 pho_size[x]                  =  myphoton_container[x].superCluster()->clustersSize();
	 pho_sc_energy[x]             =  myphoton_container[x].superCluster()->energy();
	 pho_sc_eta[x]                =  myphoton_container[x].superCluster()->eta();
	 pho_sc_phi[x]                =  correct_phi(myphoton_container[x].superCluster()->phi());
	 pho_sc_x[x]                  =  myphoton_container[x].superCluster()->x();
         pho_sc_y[x]                  =  myphoton_container[x].superCluster()->y();
         pho_sc_z[x]                  =  myphoton_container[x].superCluster()->z();
	 pho_sc_etaWidth[x]           =  myphoton_container[x].superCluster()->etaWidth();
	 pho_sc_phiWidth[x]           =  myphoton_container[x].superCluster()->phiWidth();
	 pho_HoE[x]                   =  myphoton_container[x].hadronicOverEm();              
	 pho_ecalRecHitSumEtConeDR03[x]      =  myphoton_container[x].ecalRecHitSumEtConeDR03();
	 pho_hcalTowerSumEtConeDR03[x]       =  myphoton_container[x].hcalTowerSumEtConeDR03();
	 pho_trkSumPtHollowConeDR03[x]       =  myphoton_container[x].trkSumPtHollowConeDR03();
	 pho_trkSumPtSolidConeDR03[x]        =  myphoton_container[x].trkSumPtSolidConeDR03();
	 pho_nTrkSolidConeDR03[x]            = myphoton_container[x].nTrkSolidConeDR03();
	 pho_nTrkHollowConeDR03[x]           = myphoton_container[x].nTrkHollowConeDR03();
	 pho_hcalDepth1TowerSumEtConeDR03[x] = myphoton_container[x].hcalDepth1TowerSumEtConeDR03();
	 pho_hcalDepth2TowerSumEtConeDR03[x] = myphoton_container[x].hcalDepth2TowerSumEtConeDR03();
	 pho_ecalRecHitSumEtConeDR04[x]      =  myphoton_container[x].ecalRecHitSumEtConeDR04();
	 pho_hcalTowerSumEtConeDR04[x]       =  myphoton_container[x].hcalTowerSumEtConeDR04();
	 pho_trkSumPtHollowConeDR04[x]       =  myphoton_container[x].trkSumPtHollowConeDR04();
	 pho_trkSumPtSolidConeDR04[x]        =  myphoton_container[x].trkSumPtSolidConeDR04();
	 pho_nTrkSolidConeDR04[x]            = myphoton_container[x].nTrkSolidConeDR04();
	 pho_nTrkHollowConeDR04[x]           = myphoton_container[x].nTrkHollowConeDR04();
	 pho_hcalDepth1TowerSumEtConeDR04[x] = myphoton_container[x].hcalDepth1TowerSumEtConeDR04();
	 pho_hcalDepth2TowerSumEtConeDR04[x] = myphoton_container[x].hcalDepth2TowerSumEtConeDR04();
	 pho_hasPixelSeed[x]                 = myphoton_container[x].hasPixelSeed(); 
	 pho_isEB[x]                         = myphoton_container[x].isEB(); 
	 pho_isEE[x]                         = myphoton_container[x].isEE();
	 pho_isEBGap[x]                      = myphoton_container[x].isEBGap(); 
	 pho_isEEGap[x]                      = myphoton_container[x].isEEGap(); 
	 pho_isEBEEGap[x]                    = myphoton_container[x].isEBEEGap(); 
         pho_hasConvTrk[x]                   = myphoton_container[x].hasConversionTracks();
	 
	 if(myphoton_container[x].genParticleRef().isNonnull()){
	   matchpho_E[x]                =  myphoton_container[x].genPhoton()->energy();
	   matchpho_pt[x]               =  myphoton_container[x].genPhoton()->pt();
	   matchpho_eta[x]              =  myphoton_container[x].genPhoton()->eta();
	   matchpho_phi[x]              =  correct_phi(myphoton_container[x].genPhoton()->phi());
	   matchpho_px[x]               =  myphoton_container[x].genPhoton()->px();
	   matchpho_py[x]               =  myphoton_container[x].genPhoton()->py();
	   matchpho_pz[x]               =  myphoton_container[x].genPhoton()->pz();
	 }
	 else{
	   matchpho_E[x]                = -99.;
	   matchpho_pt[x]               = -99.;
	   matchpho_eta[x]              = -99.;
	   matchpho_phi[x]              = -99.;
	   matchpho_px[x]               = -99.;
	   matchpho_py[x]               = -99.;
	   matchpho_pz[x]               = -99.;
	 }
	 ismatchedpho[x]                =  myphoton_container[x].genParticleRef().isNonnull();

	 pho_nTracks[x]                    = 9999;
	 pho_isConverted[x]                = false;
	 pho_pairInvariantMass[x]          = -99.;
	 pho_pairCotThetaSeparation[x]     = -99.;
	 pho_pairMomentum_x[x]             = -99.;
	 pho_pairMomentum_y[x]             = -99.;
	 pho_pairMomentum_z[x]             = -99.;
	 pho_conv_vx[x]                    = -99.;
	 pho_conv_vy[x]                    = -99.;
	 pho_conv_vz[x]                    = -99.;
	 pho_EoverP[x]                     = -99.;
	 pho_zOfPrimaryVertex[x]           = -99.;
	 pho_distOfMinimumApproach[x]      = -99.;
	 pho_dPhiTracksAtVtx[x]            = -99.;
	 pho_dPhiTracksAtEcal[x]           = -99.;
	 pho_dEtaTracksAtEcal[x]           = -99.;
	 
	 reco::ConversionRefVector conversions   = myphoton_container[x].conversions();
	 for (unsigned int iConv=0; iConv<conversions.size(); iConv++) {
	   reco::ConversionRef aConv=conversions[iConv];
	   if ( aConv->nTracks() <2 ) continue; 
	   if ( aConv->conversionVertex().isValid() ){
	     pho_nTracks[x]                    = aConv->nTracks();
	     pho_isConverted[x]                = aConv->isConverted();
	     pho_pairInvariantMass[x]          = aConv->pairInvariantMass();
	     pho_pairCotThetaSeparation[x]     = aConv->pairCotThetaSeparation();
	     pho_pairMomentum_x[x]             = aConv->pairMomentum().x();
	     pho_pairMomentum_y[x]             = aConv->pairMomentum().y();
	     pho_pairMomentum_z[x]             = aConv->pairMomentum().z();
	     pho_conv_vx[x]                    = aConv->conversionVertex().x();
	     pho_conv_vy[x]                    = aConv->conversionVertex().y();
	     pho_conv_vz[x]                    = aConv->conversionVertex().z();
	     pho_EoverP[x]                     = aConv->EoverP();
	     pho_zOfPrimaryVertex[x]           = aConv->zOfPrimaryVertexFromTracks();
	     pho_distOfMinimumApproach[x]      = aConv->distOfMinimumApproach();
	     pho_dPhiTracksAtVtx[x]            = aConv->dPhiTracksAtVtx();
	     pho_dPhiTracksAtEcal[x]           = aConv->dPhiTracksAtEcal();
	     pho_dEtaTracksAtEcal[x]           = aConv->dEtaTracksAtEcal();
	   }//end of if ( aConv->conversionVertex().isValid() )
	 }//end of for (unsigned int iConv=0; iConv<conversions.size(); iConv++)
	 
	 if(runHErechit_ && pho_isEB[x]){
	   //Store HE rechits
	   edm::Handle<HBHERecHitCollection> hcalRecHitHandle;
	   iEvent.getByLabel(hcalrechitLabel_, hcalRecHitHandle);
	   const HBHERecHitCollection *hbhe =  hcalRecHitHandle.product();
	   
	   for(HBHERecHitCollection::const_iterator hh = hbhe->begin(); hh != hbhe->end() && HERecHit_subset_n<10000; hh++){
	     HcalDetId id(hh->detid());
	     if (id.subdet()==2){
	       edm::ESHandle<CaloGeometry> geoHandle;	       
	       //iSetup.get<CaloGeometryRecord>().get(geoHandle);
	       //const CaloGeometry* caloGeom = geoHandle.product();
	       const CaloCellGeometry *hbhe_cell = caloGeom->getGeometry(hh->id());
	       Global3DPoint hbhe_position = hbhe_cell->getPosition();
	       
	       if(fabs(deltaPhi(pho_sc_phi[x],correct_phi(hbhe_position.phi())) ) < 0.5
		  && hh->energy()>1.){
		 //find the detid in the set
		 set<DetId>::const_iterator HERecHitChecker = HERecHitSet.find(hh->detid());
		 //if detid is not found in the set,(we reached the end), save info!
		 if(HERecHitChecker == HERecHitSet.end()){
		   HERecHitSet.insert(hh->detid());
		   HERecHit_subset_detid[HERecHit_subset_n]  = hh->detid();
		   HERecHit_subset_energy[HERecHit_subset_n] = hh->energy();
		   HERecHit_subset_time[HERecHit_subset_n]   = hh->time();
		   HERecHit_subset_depth[HERecHit_subset_n]  = id.depth();
		   HERecHit_subset_phi[HERecHit_subset_n]    = correct_phi(hbhe_position.phi());
		   HERecHit_subset_eta[HERecHit_subset_n]	  = hbhe_position.eta();
		   HERecHit_subset_x[HERecHit_subset_n]      = hbhe_position.x();
		   HERecHit_subset_y[HERecHit_subset_n]      = hbhe_position.y();
		   HERecHit_subset_z[HERecHit_subset_n]      = hbhe_position.z();
		   HERecHit_subset_n++;
		 }//check to see if hit is already saved
	       }//if delta dphi from photon is small and E>1 try to save
	     }
	   }
	 }// runHErechit_ && pho_isEB
         
	 Photon_n++;
       }//end of for loop over x
     }//if(myphoton_container.size!=0) 
     

     //to get the photon hit information from every crystal of SC
     if(runrechit_){ 
       Handle<EcalRecHitCollection> Brechit;//barrel
       Handle<EcalRecHitCollection> Erechit;//endcap
       iEvent.getByLabel(rechitBLabel_,Brechit);
       iEvent.getByLabel(rechitELabel_,Erechit);

       //this will be needed later for swiss corss
       EcalClusterLazyTools lazyTool(iEvent, iSetup,rechitBLabel_, rechitELabel_ );

       const EcalRecHitCollection* barrelRecHits= Brechit.product();
       const EcalRecHitCollection* endcapRecHits= Erechit.product();


       edm::ESHandle<CaloTopology> pTopology;
       iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
       const CaloTopology *topology = theCaloTopo_.product();

       if(myphoton_container.size()!=0){

	 for(unsigned int x=0; x < myphoton_container.size();x++){   

	   std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = myphoton_container[x].superCluster()->hitsAndFractions();
	   std::vector<CrystalInfo> crystalinfo_container;
	   crystalinfo_container.clear();
	   CrystalInfo crystal;
	   float timing_avg =0.0;
	   int ncrys   = 0;
	   ncrysPhoton[x]= 0;
	   vector< std::pair<DetId, float> >::const_iterator detitr;
	   for(detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr){
	     if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) {
	       EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
	       EcalRecHitCollection::const_iterator thishit;
	       if ( j!= Brechit->end())  thishit = j;
	       if ( j== Brechit->end()){
		 continue;
	       }
	       EBDetId detId  = (EBDetId)((*detitr).first);
	       crystal.rawId  = thishit->id().rawId();
	       crystal.energy = thishit->energy();
	       crystal.time   = thishit->time();
	       crystal.ieta   = detId.ieta();
	       crystal.iphi   = detId.iphi();
	       if(crystal.energy > 0.1){
		 timing_avg  = timing_avg + crystal.time;
		 ncrys++;
	       }  
	     }//end of if ((*detitr).det() == DetId::Ecal && (*detitr).subdetId() == EcalBarrel)
	     else if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap){
	       EcalRecHitCollection::const_iterator j= Erechit->find(((*detitr).first));
	       EcalRecHitCollection::const_iterator thishit;
	       if ( j!= Erechit->end())  thishit = j;
	       if ( j== Erechit->end()){
		 continue;
	       }
	       EEDetId detId  = (EEDetId)((*detitr).first);
	       crystal.energy = thishit->energy();
	       crystal.time   = thishit->time();
	       crystal.rawId  = 999;
	       crystal.ieta   = -99;
	       crystal.iphi   = -99;
	       if(crystal.energy > 0.1){
		 timing_avg  = timing_avg + crystal.time;
		 ncrys++;
	       } 
	     }//end of if EcalEndcap
	     crystalinfo_container.push_back(crystal);  
	   }//End loop over detids
	   std::sort(crystalinfo_container.begin(),crystalinfo_container.end(),EnergySortCriterium());
	   //Without taking into account uncertainty, this time makes no sense.
	   if (ncrys !=0) timing_avg = timing_avg/(float)ncrys;
	   else timing_avg = -99.;
	   ncrysPhoton[x] = crystalinfo_container.size(); 
	   pho_timingavg_xtal[x]      = timing_avg;
	   for (unsigned int y =0; y < 100.;y++){
	     pho_timing_xtal[x][y]         = -99.;
	     pho_energy_xtal[x][y]         = -99.;
	     pho_ieta_xtalEB[x][y]           = -99;
	     pho_iphi_xtalEB[x][y]           = -99;
	   }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++)
	   for (unsigned int y =0; y < crystalinfo_container.size() && y < 100;y++){ 
	     pho_timing_xtal[x][y]         = crystalinfo_container[y].time;
	     pho_energy_xtal[x][y]         = crystalinfo_container[y].energy;
	     pho_ieta_xtalEB[x][y]           = crystalinfo_container[y].ieta;
	     pho_iphi_xtalEB[x][y]           = crystalinfo_container[y].iphi;
	   }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++
           const reco::BasicCluster& seedClus = *(myphoton_container[x].superCluster()->seed());

           //edm::ESHandle<CaloGeometry> geoHandle;	       
	   //iSetup.get<CaloGeometryRecord>().get(geoHandle);
	   //const CaloGeometry* caloGeom = geoHandle.product();

	   if(myphoton_container[x].isEB()){
	     std::vector<float> showershapes_barrel = EcalClusterTools::roundnessBarrelSuperClusters(*(myphoton_container[x].superCluster()),*barrelRecHits,0);
	     pho_roundness[x]    = (float)showershapes_barrel[0];
	     pho_angle[x]        = (float)showershapes_barrel[1];
	     pho_s9[x]           = pho_energy_xtal[x][0]/pho_e3x3[x];

             //-New way to get swiss cross
             const reco::CaloClusterPtr  seed = myphoton_container[x].superCluster()->seed();
             DetId id = lazyTool.getMaximum(*seed).first;
             float swissCross=-99.;

             const EcalRecHitCollection & rechits = ( myphoton_container[x].isEB() ? *Brechit : *Brechit);
             EcalRecHitCollection::const_iterator it = rechits.find( id );

               if( it != rechits.end() ){swissCross = EcalSeverityLevelAlgo::swissCross( id, rechits, 0.08, true);}
               pho_swissCross[x]=swissCross;
               pho_e2e9[x]      = -99.;
               pho_e2e9[x]      = GetE2OverE9(id,rechits);
	     if(debug_ && 1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) 
	       cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl;

            vector<float> stdCov = EcalClusterTools::covariances(seedClus,&(*barrelRecHits),&(*topology),&(*caloGeom));
            vector<float> crysCov = EcalClusterTools::localCovariances(seedClus,&(*barrelRecHits),&(*topology));
            pho_SigmaEtaPhi[x]   = sqrt(stdCov[1]);
            pho_SigmaIetaIphi[x] = sqrt(crysCov[1]);    
            pho_SigmaPhiPhi[x]   = sqrt(stdCov[2]);
            pho_SigmaIphiIphi[x] = sqrt(crysCov[2]);
 
	}//end of if(myphoton_container[x].isEB())
	   else{ 
	     pho_roundness[x]   = -99.;
	     pho_angle[x]       = -99.;
	     pho_s9[x]          = pho_energy_xtal[x][0]/pho_e3x3[x];

             //----New way to get the swiss cross
             const reco::CaloClusterPtr  seed = myphoton_container[x].superCluster()->seed();
             DetId id = lazyTool.getMaximum(*seed).first;
             float swissCross=-99.;

             const EcalRecHitCollection & rechits = ( myphoton_container[x].isEB() ? *Erechit : *Erechit);
             EcalRecHitCollection::const_iterator it = rechits.find( id );

            if( it != rechits.end() ) {swissCross = EcalSeverityLevelAlgo::swissCross( id, rechits, 0.08, true);}

            pho_swissCross[x]= swissCross;
            pho_e2e9[x]      = -99.;
            pho_e2e9[x]      = GetE2OverE9(id,rechits); 

	     if(debug_ && 1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) {
	       cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm." << endl;
	       cout<<"This would be weird since there aren't spikes in the endcap of ECAL"<<endl; 
	     }

             vector<float> stdCov = EcalClusterTools::covariances(seedClus,&(*endcapRecHits),&(*topology),&(*caloGeom));
             vector<float> crysCov = EcalClusterTools::localCovariances(seedClus,&(*endcapRecHits),&(*topology));
             pho_SigmaEtaPhi[x]   = sqrt(stdCov[1]);
             pho_SigmaIetaIphi[x] = sqrt(crysCov[1]);   
             pho_SigmaPhiPhi[x]   = sqrt(stdCov[2]);
             pho_SigmaIphiIphi[x] = sqrt(crysCov[2]);

	   }//end of else (if !EB)i

	 }//end of for loop over x
       }//if(myphoton_container.size!=0)

 
    if(runrechit_){
     //edm::ESHandle<CaloGeometry> geoHandle;
     //iSetup.get<CaloGeometryRecord>().get(geoHandle);
     //const CaloGeometry* caloGeom = geoHandle.product();

     int EBRecHit_n =0;
     int EERecHit_n =0;

     EBRecHit_size= EBRecHit_n;
     EERecHit_size= EERecHit_n;

     for (EcalRecHitCollection::const_iterator it = barrelRecHits->begin();it!=barrelRecHits->end();++it){
       EBDetId dit = it->detid();

       if(EBRecHit_n>=10000)continue; 

       const CaloCellGeometry *this_cell = caloGeom->getGeometry(it->id());
       GlobalPoint position = this_cell->getPosition();
       float ET = (it->energy()*sin(position.theta())); 
       EBRecHit_eta[EBRecHit_n] = position.eta();
       EBRecHit_phi[EBRecHit_n] = position.phi();
       EBRecHit_ieta[EBRecHit_n] = dit.ieta();
       EBRecHit_iphi[EBRecHit_n] = dit.iphi();
       EBRecHit_e[EBRecHit_n] = it->energy();
       EBRecHit_et[EBRecHit_n] = ET;
       EBRecHit_flag[EBRecHit_n] = it->recoFlag();
       EBRecHit_time[EBRecHit_n] = it->time();
       EBRecHit_size = EBRecHit_n;
       EBRecHit_n++;  }

     for (EcalRecHitCollection::const_iterator it = endcapRecHits->begin();it!=endcapRecHits->end();++it){
       EEDetId dit = it->detid();

       if(EERecHit_n>=10000)continue;

       const CaloCellGeometry *this_cell = caloGeom->getGeometry(it->id());
       GlobalPoint position = this_cell->getPosition();

       float ET = (it->energy()*sin(position.theta())); 
       EERecHit_eta[EERecHit_n] = position.eta();
       EERecHit_phi[EERecHit_n] = position.phi();
       EERecHit_ieta[EERecHit_n] = dit.ix();
       EERecHit_iphi[EERecHit_n] = dit.iy();
       EERecHit_e[EERecHit_n] = it->energy();
       EERecHit_et[EERecHit_n] = ET;
       EERecHit_flag[EERecHit_n] = it->recoFlag();
       EERecHit_time[EERecHit_n] = it->time();
       EERecHit_size = EERecHit_n;
       EERecHit_n++;  }
      }//if(runrechit_)    
 
     }//if(runrechit_)
     
   }//if(runphotons_)

   
   if(runCSCseg_){ //Store CSC segments
     // Get the CSC Geometry :
     ESHandle<CSCGeometry> cscGeom;
     iSetup.get<MuonGeometryRecord>().get(cscGeom);
     // get CSC segment collection
     Handle<CSCSegmentCollection> cscSegments;
     iEvent.getByLabel(cscLabel_, cscSegments);
     CSCseg_n = 0;
           
     for(CSCSegmentCollection::const_iterator it=cscSegments->begin(); CSCseg_n<10000 && it != cscSegments->end(); it++){
       CSCDetId id  = (CSCDetId)it->cscDetId();
       LocalVector segDir = it->localDirection();
       LocalPoint localPos = it->localPosition();

       // global transformation
       const CSCChamber* cscchamber = cscGeom->chamber(id);
       if (cscchamber) {
         
         GlobalPoint globalPosition = cscchamber->toGlobal(localPos);
         GlobalVector globalDirection = cscchamber->toGlobal(segDir);
         
         int nhits = 0;
         float timeSum = 0;
         // Get the CSC recHits that contribute to this segment.
         std::vector<CSCRecHit2D> theseRecHits = it->specificRecHits();
         for( vector<CSCRecHit2D>::const_iterator iRH =theseRecHits.begin(); iRH != theseRecHits.end(); iRH++){
           if( !(iRH->isValid()) ) continue;  // only interested in valid hits
           nhits++;
           float rhTime = iRH->tpeak();
           timeSum += rhTime;
         }//end rechit loop

         float segmentTime = -999;
         if (nhits>0)
           segmentTime = timeSum/nhits;
         
         CSCseg_time[CSCseg_n] = segmentTime;
         CSCseg_x[CSCseg_n]    = globalPosition.x();
         CSCseg_y[CSCseg_n]    = globalPosition.y();
         CSCseg_z[CSCseg_n]    = globalPosition.z();
         CSCseg_phi[CSCseg_n]  = globalPosition.phi();
               
         CSCseg_DirectionX[CSCseg_n] = globalDirection.x();
         CSCseg_DirectionY[CSCseg_n] = globalDirection.y();
         CSCseg_DirectionZ[CSCseg_n] = globalDirection.z();
        
         CSCseg_n++;
       }
     }// loop over segs
   }// runCSCseg_
  

 if(runBeamHaloSummary_){
   // Get BeamHaloSummary
   edm::Handle<BeamHaloSummary> TheBeamHaloSummary ;
   iEvent.getByLabel(BeamHaloSummaryLabel_, TheBeamHaloSummary) ;
    
 
  isBeamHaloIDTightPass    = false;
  isBeamHaloIDLoosePass    = false;
     
  isBeamHaloEcalLoosePass   = false;
  isBeamHaloHcalLoosePass   = false;
  isBeamHaloCSCLoosePass    = false;
  isBeamHaloGlobalLoosePass = false;
     
  isBeamHaloEcalTightPass   = false;
  isBeamHaloHcalTightPass   = false;
  isBeamHaloCSCTightPass    = false;
  isBeamHaloGlobalTightPass = false;
    
  
  isSmellsLikeHalo_Tag  = false;
  isLooseHalo_Tag       = false;
  isTightHalo_Tag       = false;
  isExtremeTightHalo_Tag = false;
 
     
  if(TheBeamHaloSummary.isValid()) {
   const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
     
        if(TheSummary.EcalLooseHaloId())isBeamHaloEcalLoosePass   = true;
        if(TheSummary.HcalLooseHaloId())isBeamHaloHcalLoosePass   = true;
     
        if(TheSummary.EcalTightHaloId())isBeamHaloEcalTightPass   = true;
        if(TheSummary.HcalTightHaloId())isBeamHaloHcalTightPass   = true;
     
        if(TheSummary.CSCLooseHaloId())isBeamHaloCSCLoosePass    = true;
        if(TheSummary.CSCTightHaloId())isBeamHaloCSCTightPass    = true;
     
        if(TheSummary.GlobalLooseHaloId())isBeamHaloGlobalLoosePass = true;
        if(TheSummary.GlobalTightHaloId())isBeamHaloGlobalTightPass = true;
    
        if(TheSummary.EventSmellsLikeHalo())isSmellsLikeHalo_Tag = true;
        if(TheSummary.LooseId())isLooseHalo_Tag = true;
        if(TheSummary.TightId())isTightHalo_Tag = true;
        if(TheSummary.ExtremeTightId())isExtremeTightHalo_Tag = true;
   
  
   if( TheSummary.EcalLooseHaloId()  && TheSummary.HcalLooseHaloId() &&
       TheSummary.CSCLooseHaloId()   && TheSummary.GlobalLooseHaloId() )
        isBeamHaloIDLoosePass = true;
     
   if( TheSummary.EcalTightHaloId()  && TheSummary.HcalTightHaloId() &&
       TheSummary.CSCTightHaloId()   && TheSummary.GlobalTightHaloId() )
        isBeamHaloIDTightPass = true;
     
   }//not empty
     
     
 }//if(runBeamHaloSummary_)
     
 
  if(runRPChit_){
     /*
     // Get the RPC Geometry  
     edm::ESHandle<RPCGeometry> rpcGeo;
     iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  
     //get RPC hit collection
     edm::Handle<RPCRecHitCollection> rpcHits;
     iEvent.getByLabel("cscSegments",rpcHits);
     
     for(RPCRecHitCollection::const_iterator RecHitsIt = rpcRecHitRange->begin(); RecHitsIt!=rpcRecHitRange->end(); RecHitsIt++){
     
     
     }//loop over rpc hits
     */
     
    //get geometry of tracking stuff
    //edm::ESHandle<GlobalTrackingGeometry> geometry;
    edm::ESHandle<RPCGeometry> geometry;
    //iSetup.get<GlobalTrackingGeometryRecord>().get(geometry);
    iSetup.get<MuonGeometryRecord>().get(geometry);
    //get RPC hit collection
     edm::Handle<RPCRecHitCollection> rpcHits;
     iEvent.getByLabel(rpcLabel_,rpcHits);
     
     RPChit_n=0;
     for(RPCRecHitCollection::const_iterator RecHitsIt = rpcHits->begin(); RPChit_n<10000 && RecHitsIt!=rpcHits->end(); RecHitsIt++){
       const GeomDet *det = geometry->idToDet(RecHitsIt->geographicalId());
       GlobalPoint pos = det->toGlobal(RecHitsIt->localPosition());
       RPChit_x[RPChit_n] = pos.x();
       RPChit_y[RPChit_n] = pos.y();
       RPChit_z[RPChit_n] = pos.z();
       //BunchX is in units of 25ns. BunchX=0 is the expected value for collision muons.
       RPChit_BunchX[RPChit_n] = RecHitsIt->BunchX();
       RPChit_n++;
     }//loop over rpc hits

   }//runRPChit_
   
   //calomet variables  
   if(runmet_){
     edm::Handle<edm::View<pat::MET> > metHandle;
     iEvent.getByLabel(metLabel_,metHandle);
     if ( metHandle.isValid() ){
      const edm::View<pat::MET> & met = *metHandle;
      //edm::View<pat::MET>::const_iterator met;
       CaloMetSig                            = met[0].mEtSig();
       CaloEtFractionHadronic                = met[0].etFractionHadronic();
       CaloEmEtFraction                      = met[0].emEtFraction();
       CaloHadEtInHB                         = met[0].hadEtInHB();
       CaloHadEtInHO                         = met[0].hadEtInHO();
       CaloHadEtInHE                         = met[0].hadEtInHE();
       CaloHadEtInHF                         = met[0].hadEtInHF();
       CaloEmEtInEB                          = met[0].emEtInEB();
       CaloEmEtInEE                          = met[0].emEtInEE();
       CaloEmEtInHF                          = met[0].emEtInHF();
       CaloMetEz                             = met[0].e_longitudinal();
       CaloMaxEtInEmTowers                   = met[0].maxEtInEmTowers();
       CaloMaxEtInHadTowers                  = met[0].maxEtInHadTowers();
       // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU
       
       CaloMetPt[0]                          = met[0].pt();
       CaloMetPx[0]                          = met[0].px();
       CaloMetPy[0]                          = met[0].py();
       CaloMetPhi[0]                         = correct_phi(met[0].phi());
       CaloMetSumEt[0]                       = met[0].sumEt();
       
       CaloMetPt[1]                          = met[0].uncorrectedPt(pat::MET::uncorrNONE); 
       CaloMetPhi[1]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
       CaloMetPx[1]                          = met[0].corEx(pat::MET::uncorrNONE); 
       CaloMetPy[1]                          = met[0].corEy(pat::MET::uncorrNONE); 
       CaloMetSumEt[1]                       = met[0].corSumEt(pat::MET::uncorrNONE); 
       
       CaloMetPt[2]                          = met[0].uncorrectedPt(pat::MET::uncorrALL); 
       CaloMetPhi[2]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrALL)); 
       CaloMetPx[2]                          = met[0].corEx(pat::MET::uncorrALL); 
       CaloMetPy[2]                          = met[0].corEy(pat::MET::uncorrALL); 
       CaloMetSumEt[2]                       = met[0].corSumEt(pat::MET::uncorrALL); 
       
       CaloMetPt[3]                          = met[0].uncorrectedPt(pat::MET::uncorrJES); 
       CaloMetPhi[3]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrJES)); 
       CaloMetPx[3]                          = met[0].corEx(pat::MET::uncorrJES); 
       CaloMetPy[3]                          = met[0].corEy(pat::MET::uncorrJES); 
       CaloMetSumEt[3]                       = met[0].corSumEt(pat::MET::uncorrJES); 
       
       CaloMetPt[4]                          = met[0].uncorrectedPt(pat::MET::uncorrMUON); 
       CaloMetPhi[4]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
       CaloMetPx[4]                          = met[0].corEx(pat::MET::uncorrMUON); 
       CaloMetPy[4]                          = met[0].corEy(pat::MET::uncorrMUON); 
       CaloMetSumEt[4]                       = met[0].corSumEt(pat::MET::uncorrMUON); 
       
       CaloMetPt[5]                          = met[0].uncorrectedPt(pat::MET::uncorrTAU); 
       CaloMetPhi[5]                         = correct_phi(met[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
       CaloMetPx[5]                          = met[0].corEx(pat::MET::uncorrTAU); 
       CaloMetPy[5]                          = met[0].corEy(pat::MET::uncorrTAU); 
       CaloMetSumEt[5]                       = met[0].corSumEt(pat::MET::uncorrTAU); 
       if(runphotons_==1)
	 if (myphoton_container.size()!=0)
	   Delta_phi                       = fabs(reco::deltaPhi(correct_phi(met[0].phi()),correct_phi(myphoton_container[0].phi())));
       if(rungenmet_){
	 const reco::GenMET *genMet = met[0].genMET();
	 genMetPt     = genMet->et();
	 genMetPhi    = correct_phi(genMet->phi());
	 genMetSumEt  = genMet->sumEt();
	 genMetPx     = genMet->px();
	 genMetPy     = genMet->py();
	 if(runphotons_==1)
	   if (myphoton_container.size()!=0)
	     Delta_phiGEN                        = fabs(reco::deltaPhi(correct_phi(genMet->phi()),correct_phi(myphoton_container[0].phi())));
       }
     }
   }
   if(runPFmet_){
     edm::Handle<edm::View<pat::MET> > metPFHandle;
     iEvent.getByLabel(PFmetLabel_,metPFHandle);
     const edm::View<pat::MET> & metsPF = *metPFHandle;
     if ( metPFHandle.isValid() ){
     // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU
     PFMetPt[0]                          = metsPF[0].pt();
     PFMetPx[0]                          = metsPF[0].px();
     PFMetPy[0]                          = metsPF[0].py();
     PFMetPhi[0]                         = correct_phi(metsPF[0].phi());
     PFMetSumEt[0]                       = metsPF[0].sumEt();
     
     PFMetPt[1]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrNONE); 
     PFMetPhi[1]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
     PFMetPx[1]                          = metsPF[0].corEx(pat::MET::uncorrNONE); 
     PFMetPy[1]                          = metsPF[0].corEy(pat::MET::uncorrNONE); 
     PFMetSumEt[1]                       = metsPF[0].corSumEt(pat::MET::uncorrNONE); 
     
     PFMetPt[2]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrALL); 
     PFMetPhi[2]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrALL)); 
     PFMetPx[2]                          = metsPF[0].corEx(pat::MET::uncorrALL); 
     PFMetPy[2]                          = metsPF[0].corEy(pat::MET::uncorrALL); 
     PFMetSumEt[2]                       = metsPF[0].corSumEt(pat::MET::uncorrALL); 
     
     PFMetPt[3]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrJES); 
     PFMetPhi[3]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrJES)); 
     PFMetPx[3]                          = metsPF[0].corEx(pat::MET::uncorrJES); 
     PFMetPy[3]                          = metsPF[0].corEy(pat::MET::uncorrJES); 
     PFMetSumEt[3]                       = metsPF[0].corSumEt(pat::MET::uncorrJES); 
     
     PFMetPt[4]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrMUON); 
     PFMetPhi[4]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
     PFMetPx[4]                          = metsPF[0].corEx(pat::MET::uncorrMUON); 
     PFMetPy[4]                          = metsPF[0].corEy(pat::MET::uncorrMUON); 
     PFMetSumEt[4]                       = metsPF[0].corSumEt(pat::MET::uncorrMUON); 
     
     PFMetPt[5]                          = metsPF[0].uncorrectedPt(pat::MET::uncorrTAU); 
     PFMetPhi[5]                         = correct_phi(metsPF[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
     PFMetPx[5]                          = metsPF[0].corEx(pat::MET::uncorrTAU); 
     PFMetPy[5]                          = metsPF[0].corEy(pat::MET::uncorrTAU); 
     PFMetSumEt[5]                       = metsPF[0].corSumEt(pat::MET::uncorrTAU); 

     if(runphotons_==1)
       if (myphoton_container.size()!=0)
	 Delta_phiPF  = fabs(reco::deltaPhi(PFMetPhi[0],correct_phi(myphoton_container[0].phi())));
     }
     else{
       LogWarning("METEventSelector") << "No Met results for InputTag " ;
       return;
     }
   }
   if(runTCmet_){
     edm::Handle<edm::View<pat::MET> > metTCHandle;
     iEvent.getByLabel(TCmetLabel_,metTCHandle);
     const edm::View<pat::MET> & metsTC = *metTCHandle;
     if ( metTCHandle.isValid() ){
       // 0=full corr1=uncorrNone  2=uncorrALL 3=uncorrJES  4=uncorrMUON  5=TAU
       TCMetPt[0]                          = metsTC[0].pt();
       TCMetPx[0]                          = metsTC[0].px();
       TCMetPy[0]                          = metsTC[0].py();
       TCMetPhi[0]                         = correct_phi(metsTC[0].phi());
       TCMetSumEt[0]                       = metsTC[0].sumEt();
       
       TCMetPt[1]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrNONE); 
       TCMetPhi[1]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrNONE)); 
       TCMetPx[1]                          = metsTC[0].corEx(pat::MET::uncorrNONE); 
       TCMetPy[1]                          = metsTC[0].corEy(pat::MET::uncorrNONE); 
       TCMetSumEt[1]                       = metsTC[0].corSumEt(pat::MET::uncorrNONE); 
       
       TCMetPt[2]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrALL); 
       TCMetPhi[2]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrALL)); 
       TCMetPx[2]                          = metsTC[0].corEx(pat::MET::uncorrALL); 
       TCMetPy[2]                          = metsTC[0].corEy(pat::MET::uncorrALL); 
       TCMetSumEt[2]                       = metsTC[0].corSumEt(pat::MET::uncorrALL); 
       
       TCMetPt[3]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrJES); 
       TCMetPhi[3]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrJES)); 
       TCMetPx[3]                          = metsTC[0].corEx(pat::MET::uncorrJES); 
       TCMetPy[3]                          = metsTC[0].corEy(pat::MET::uncorrJES); 
       TCMetSumEt[3]                       = metsTC[0].corSumEt(pat::MET::uncorrJES); 
       
       TCMetPt[4]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrMUON); 
       TCMetPhi[4]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrMUON)); 
       TCMetPx[4]                          = metsTC[0].corEx(pat::MET::uncorrMUON); 
       TCMetPy[4]                          = metsTC[0].corEy(pat::MET::uncorrMUON); 
       TCMetSumEt[4]                       = metsTC[0].corSumEt(pat::MET::uncorrMUON); 
       
       TCMetPt[5]                          = metsTC[0].uncorrectedPt(pat::MET::uncorrTAU); 
       TCMetPhi[5]                         = correct_phi(metsTC[0].uncorrectedPhi(pat::MET::uncorrTAU)); 
       TCMetPx[5]                          = metsTC[0].corEx(pat::MET::uncorrTAU); 
       TCMetPy[5]                          = metsTC[0].corEy(pat::MET::uncorrTAU); 
       TCMetSumEt[5]                       = metsTC[0].corSumEt(pat::MET::uncorrTAU); 
       
       if(runphotons_==1)
	 if (myphoton_container.size()!=0)
	   Delta_phiTC  = fabs(reco::deltaPhi(TCMetPhi[0],correct_phi(myphoton_container[0].phi()))); 
     }
     else{
       LogWarning("METEventSelector") << "No Met results for InputTag " ;
       return;
     } 
     
   }//end of if(runTCmet_)


    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;

   if(runjets_){
     //this for jec uncert 
     iSetup.get<JetCorrectionsRecord>().get("AK5Calo",JetCorParColl); 
     JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
     JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
       

     edm::Handle<edm::View<pat::Jet> > jetHandle;
     iEvent.getByLabel(jetLabel_,jetHandle);
     const edm::View<pat::Jet> & jets = *jetHandle;
     
     size_t njetscounter=0;
     std::vector<pat::Jet>  myjet_container;
     myjet_container.clear();
     for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
       if(jet_iter->pt()>30)  njetscounter++;
       myjet_container.push_back(*jet_iter);
     }
     Jet_n = 0;  
    //  std::cout<<"Jet Container Size = "<<myjet_container.size()<<std::endl;
     if(myjet_container.size()!=0){
       for(unsigned int x=0;x < min(myjet_container.size(),MaxN);x++){
     //    if(x==0)std::cout<<"jet pt ="<<myjet_container[x].pt()<<std::endl;
	 jet_pt[x]  = myjet_container[x].pt();
	 jet_px[x]  = myjet_container[x].px();
	 jet_py[x]  = myjet_container[x].py();
         jet_E[x]   = myjet_container[x].energy();
	 jet_pz[x]  = myjet_container[x].pz();
	 jet_vx[x]  = myjet_container[x].vx();
	 jet_vy[x]  = myjet_container[x].vy();
	 jet_vz[x]  = myjet_container[x].vz();
	 jet_phi[x] = correct_phi(myjet_container[x].phi());
	 jet_eta[x] = myjet_container[x].eta();
	 jet_emEnergyFraction[x]= myjet_container[x].emEnergyFraction();
	 jet_energyFractionHadronic[x] = myjet_container[x].energyFractionHadronic();
	 jet_hitsInN90[x]= myjet_container[x].jetID().hitsInN90;
	 jet_n90Hits[x] = myjet_container[x].jetID().n90Hits;
	 jet_fHPD[x] = (float) myjet_container[x].jetID().fHPD;
	 jet_fRBX[x] = (float) myjet_container[x].jetID().fRBX;
	 jet_RHF[x] = (float)(myjet_container[x].jetID().fLong - myjet_container[x].jetID().fShort)/(myjet_container[x].jetID().fLong + myjet_container[x].jetID().fShort);
	 jet_nTowers[x] = myjet_container[x].jetID().nECALTowers + myjet_container[x].jetID().nHCALTowers ;
         //jet energy uncertiany
         jecUnc->setJetEta(jet_eta[x]);
         jecUnc->setJetPt(jet_pt[x]);
         jet_jecUncer[x] = jecUnc->getUncertainty(true);
         if(myjet_container[x].jecFactor("Uncorrected") != 0 )
             {jet_jecCorr[x] = 1./(myjet_container[x].jecFactor("Uncorrected")); 
              }  
               else{jet_jecCorr[x] =0.;}

         //get the uncorrected jet and fill them
         pat::Jet uncjet = myjet_container[x].correctedJet("Uncorrected");
         ucjet_pt[x] = uncjet.pt();
         ucjet_px[x] = uncjet.px();
         ucjet_py[x] = uncjet.py();
         ucjet_pz[x] = uncjet.pz();
         ucjet_E[x]   = uncjet.energy();
         ucjet_eta[x] = uncjet.eta();
         ucjet_phi[x] = uncjet.phi();

 	 Jet_n++;
       }//end of for loop
     }
   }



   if(runpfjets_){
     //for jec uncert
     edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
     iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
     JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
     JetCorrectionUncertainty *pfjecUnc = new JetCorrectionUncertainty(JetCorPar);


     edm::Handle<edm::View<pat::Jet> > pfjetHandle;
     iEvent.getByLabel(pfjetLabel_,pfjetHandle);
     const edm::View<pat::Jet> & pfjets = *pfjetHandle;
  
     size_t npfjetscounter=0;
     std::vector<pat::Jet>  mypfjet_container;
     mypfjet_container.clear();

      for(edm::View<pat::Jet>::const_iterator pfjet_iter = pfjets.begin(); pfjet_iter!=pfjets.end(); ++pfjet_iter){
       if(pfjet_iter->pt()>30)  npfjetscounter++;
       mypfjet_container.push_back(*pfjet_iter);
     }

     pfJet_n = 0;  
     //std::cout<<" PF Jet Container size "<<mypfjet_container.size()<<std::endl;
     if(mypfjet_container.size()!=0){ 
       for(unsigned int x=0;x < min(mypfjet_container.size(), MaxN); x++){
         pfjet_pt[x]  = mypfjet_container[x].pt();
         pfjet_px[x]  = mypfjet_container[x].px();
         pfjet_py[x]  = mypfjet_container[x].py();
         pfjet_E[x]   = mypfjet_container[x].energy();
         pfjet_pz[x]  = mypfjet_container[x].pz();
         pfjet_vx[x]  = mypfjet_container[x].vx();
         pfjet_vy[x]  = mypfjet_container[x].vy();
         pfjet_vz[x]  = mypfjet_container[x].vz();
         pfjet_phi[x] = correct_phi(mypfjet_container[x].phi());
         pfjet_eta[x] = mypfjet_container[x].eta();
         
        if(mypfjet_container[x].jecFactor("Uncorrected")!= 0)
         {pfjet_jecCorr[x]  = (1.0/mypfjet_container[x].jecFactor("Uncorrected")); 
         }
         else{pfjet_jecCorr[x] =0.;}

         /*
         pfjet_emEnergyFraction[x]= mypfjet_container[x].emEnergyFraction();
         pfjet_energyFractionHadronic[x] = mypfjet_container[x].energyFractionHadronic();
         pfjet_hitsInN90[x]= mypfjet_container[x].jetID().hitsInN90;
         pfjet_n90Hits[x] = mypfjet_container[x].jetID().n90Hits;
         pfjet_fHPD[x] = (float) mypfjet_container[x].jetID().fHPD;
         pfjet_fRBX[x] = (float) mypfjet_container[x].jetID().fRBX;
         pfjet_RHF[x] = (float)(mypfjet_container[x].jetID().fLong - mypfjet_container[x].jetID().fShort)/(mypfjet_container[x].jetID().fLong + mypfjet_container[x].jetID().fShort);
         pfjet_nTowers[x] = mypfjet_container[x].jetID().nECALTowers + mypfjet_container[x].jetID().nHCALTowers ;
         */
         //jet energy uncertiany
         pfjecUnc->setJetEta(pfjet_eta[x]);
         pfjecUnc->setJetPt(pfjet_pt[x]);
         pfjet_jecUncer[x] = pfjecUnc->getUncertainty(true);
         //get the uncorrected jet and fill them
         pat::Jet uncpfjet = mypfjet_container[x].correctedJet("Uncorrected");
         ucpfjet_pt[x] = uncpfjet.pt();
         ucpfjet_px[x] = uncpfjet.px();
         ucpfjet_py[x] = uncpfjet.py();
         ucpfjet_pz[x] = uncpfjet.pz();
         ucpfjet_E[x]  = uncpfjet.energy();
         ucpfjet_eta[x]= uncpfjet.eta();
         ucpfjet_phi[x]= uncpfjet.phi();
 
         pfJet_n++;
       }//end of for loop
     }
   }


if(rungenjets_){

         Handle<GenJetCollection>  GenJetHandle;
         iEvent.getByLabel(genjetLabel_,GenJetHandle);
         reco::GenJetCollection::const_iterator genjet_itr;

         std::vector<reco::GenJet> genJet_container;
         genJet_container.clear();

      for(reco::GenJetCollection::const_iterator genjet_itr =GenJetHandle->begin();genjet_itr !=GenJetHandle->end();++genjet_itr){
          genJet_container.push_back(*genjet_itr);
          }

        genJet_n=0; 
        for(unsigned int y=0; y<min(genJet_container.size(),MaxN); y++)
         {
           genjet_pt[y]  =genJet_container[y].pt();
           genjet_E[y]   =genJet_container[y].energy();
           genjet_eta[y] =genJet_container[y].eta();
           genjet_phi[y] =genJet_container[y].phi();
           genjet_px[y]  =genJet_container[y].px();
           genjet_py[y]  =genJet_container[y].py();
           genjet_pz[y]  =genJet_container[y].pz();
           genJet_n++;
         }//for(y=0....)


}//if(rungenJet_)

   if(runelectrons_){
     edm::Handle<edm::View<pat::Electron> > electronHandle;
     iEvent.getByLabel(eleLabel_,electronHandle);
     vector<pat::Electron> myelectron_container;
     
     const edm::View<pat::Electron> & electrons = *electronHandle;   // const ... &, we don't make a copy of it!
     for(edm::View<pat::Electron>::const_iterator electron = electrons.begin(); electron!=electrons.end(); ++electron){
       myelectron_container.push_back(*electron);
     }
     Electron_n = 0;
     for(unsigned int x=0;x < min(myelectron_container.size(),MaxN);x++){
       electron_pt[x]       = myelectron_container[x].pt();
       electron_energy[x]   = myelectron_container[x].energy();
       electron_px[x]       = myelectron_container[x].px();
       electron_py[x]       = myelectron_container[x].py();
       electron_pz[x]       = myelectron_container[x].pz();
       electron_vx[x]       = myelectron_container[x].vx();
       electron_vy[x]       = myelectron_container[x].vy();
       electron_vz[x]       = myelectron_container[x].vz();
       electron_phi[x]      = correct_phi(myelectron_container[x].phi());
       electron_eta[x]      = myelectron_container[x].eta();
       electron_charge[x]       = myelectron_container[x].charge();
       electron_HoE[x]          = myelectron_container[x].hadronicOverEm();
       electron_SigmaIetaIeta[x]= myelectron_container[x].sigmaIetaIeta();               
       electron_trkIso[x]       = myelectron_container[x].dr03TkSumPt() ;
       electron_ecalIso[x]      = myelectron_container[x].dr03EcalRecHitSumEt();
       electron_hcalIso[x]      = myelectron_container[x].dr03HcalTowerSumEt();
       electron_dEtaIn[x]       = myelectron_container[x].deltaEtaSuperClusterTrackAtVtx();
       electron_dPhiIn[x]       = myelectron_container[x].deltaPhiSuperClusterTrackAtVtx();
       electron_sc_energy[x]    = myelectron_container[x].superCluster()->energy();
       electron_sc_eta[x]       = myelectron_container[x].superCluster()->eta();
       electron_sc_phi[x]       = correct_phi(myelectron_container[x].superCluster()->phi());
       Electron_n++;
     }//end of for loop
   }
	
   if(runtaus_){
     edm::Handle<edm::View<pat::Tau> > tauHandle;
     iEvent.getByLabel(tauLabel_,tauHandle);
     vector <pat::Tau> mytau_container;
     
     const edm::View<pat::Tau> & taus = *tauHandle;   // const ... &, we don't make a copy of it!
     for(edm::View<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau){
       mytau_container.push_back(*tau);
     }
     Tau_n = 0;
     for(unsigned int x=0;x < min(mytau_container.size(), MaxN);x++){
       tau_pt[x]  = mytau_container[x].pt();
       tau_energy[x]  = mytau_container[x].energy();
       tau_px[x]  = mytau_container[x].px();
       tau_py[x]  = mytau_container[x].py();
       tau_pz[x]  = mytau_container[x].pz();
       tau_vx[x]  = mytau_container[x].vx();
       tau_vy[x]  = mytau_container[x].vy();
       tau_vz[x]  = mytau_container[x].vz();
       tau_phi[x] = correct_phi(mytau_container[x].phi());
       tau_eta[x] = mytau_container[x].eta();
       tau_charge[x] = mytau_container[x].charge();
       Tau_n++;
     }//end of for loop
   }
   // if (Photon_n>0)
     myEvent->Fill();
   if(debug_) cout<<"DEBUG: analyze loop done"<<endl;
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob(){
  triggernames     = &all_triggers; 
  triggerprescales = &all_triggerprescales; 
  ifTriggerpassed  = &all_ifTriggerpassed; 
  f=new TFile(outFile_.c_str(),"RECREATE");
  myEvent = new TTree("myEvent","a tree with histograms");
  myEvent->Branch("nevents",&nevents,"nevents/I");
  myEvent->Branch("run",&RunNumber,"RunNumber/i");
  myEvent->Branch("event",&EventNumber,"EventNumber/i");
  myEvent->Branch("luminosityBlock",&LumiNumber,"LumiNumber/i");
  myEvent->Branch("beamCrossing",&BXNumber,"BXNumber/i");
  myEvent->Branch("totalIntensityBeam1",&totalIntensityBeam1,"totalIntensityBeam1/i");  
  myEvent->Branch("totalIntensityBeam2",&totalIntensityBeam2,"totalIntensityBeam2/i");
  myEvent->Branch("avgInsDelLumi",&avgInsDelLumi,"avgInsDelLumi/F");
  myEvent->Branch("avgInsDelLumiErr",&avgInsDelLumiErr,"avgInsDelLumiErr/F");
  myEvent->Branch("avgInsRecLumi",&avgInsRecLumi,"avgInsRecLumi/F");
  myEvent->Branch("avgInsRecLumiErr",&avgInsRecLumiErr,"avgInsRecLumiErr/F");
  myEvent->Branch("ntriggers",&ntriggers,"ntriggers/I");
  myEvent->Branch("triggernames","vector<std::string>",&triggernames);

  if(runHLT_){
  myEvent->Branch("triggerprescales","vector<int>",&triggerprescales);
  myEvent->Branch("ifTriggerpassed","vector<bool>",&ifTriggerpassed);
  }
  
  if(runvertex_){
    myEvent->Branch("Vertex_n",&Vertex_n,"Vertex_n/I");
    myEvent->Branch("Vertex_x",vx,"vx[Vertex_n]/F");
    myEvent->Branch("Vertex_y",vy,"vy[Vertex_n]/F");
    myEvent->Branch("Vertex_z",vz,"vz[Vertex_n]/F");
    myEvent->Branch("Vertex_tracksize",vtracksize,"vtracksize[Vertex_n]/I");
    myEvent->Branch("Vertex_ndof",vndof,"vndof[Vertex_n]/I");
    myEvent->Branch("Vertex_chi2",chi2,"chi2[Vertex_n]/F");
    myEvent->Branch("Vertex_d0",v_d0,"v_d0[Vertex_n]/F");
    myEvent->Branch("Vertex_isFake",v_isFake,"v_isFake[Vertex_n]/O");
  }
  
  if(runscraping_){
    myEvent->Branch("Scraping_isScrapingEvent",&Scraping_isScrapingEvent,"Scraping_isScrapingEvent/O");
    myEvent->Branch("Scraping_numOfTracks",&Scraping_numOfTracks,"Scraping_numOfTracks/I");
    myEvent->Branch("Scraping_fractionOfGoodTracks",&Scraping_fractionOfGoodTracks,"Scraping_fractionOfGoodTracks/F");
  }

 if(runPileUp_){
    myEvent->Branch("npuVertices",&npuVertices,"npuVertices/I");
    myEvent->Branch("ootnpuVertices",&ootnpuVertices,"ootnpuVertices/I");
  }
	
  if (runtracks_){
    myEvent->Branch("Track_n",&Track_n,"Track_n/I");
    myEvent->Branch("Track_px",trk_px,"trk_px[Track_n]/F");
    myEvent->Branch("Track_py",trk_py,"trk_py[Track_n]/F");
    myEvent->Branch("Track_pz",trk_pz,"trk_pz[Track_n]/F");
    myEvent->Branch("Track_vx",trk_vx,"trk_vx[Track_n]/F");
    myEvent->Branch("Track_vy",trk_vy,"trk_vy[Track_n]/F");
    myEvent->Branch("Track_vz",trk_vz,"trk_vz[Track_n]/F");
    myEvent->Branch("Track_pt",trk_pt,"trk_pt[Track_n]/F");
    myEvent->Branch("Track_eta",trk_eta,"trk_eta[Track_n]/F");
    myEvent->Branch("Track_phi",trk_phi,"trk_phi[Track_n]/F");
  }
  if (runjets_){
    myEvent->Branch("Jet_n",&Jet_n,"Jet_n/I");
    myEvent->Branch("Jet_px",jet_px,"jet_px[Jet_n]/F");
    myEvent->Branch("Jet_py",jet_py,"jet_py[Jet_n]/F");
    myEvent->Branch("Jet_E",jet_E,"jet_E[Jet_n]/F");
    myEvent->Branch("Jet_pz",jet_pz,"jet_pz[Jet_n]/F");
    myEvent->Branch("Jet_vx",jet_vx,"jet_vx[Jet_n]/F");
    myEvent->Branch("Jet_vy",jet_vy,"jet_vy[Jet_n]/F");
    myEvent->Branch("Jet_vz",jet_vz,"jet_vz[Jet_n]/F");
    myEvent->Branch("Jet_pt",jet_pt,"jet_pt[Jet_n]/F");
    myEvent->Branch("Jet_eta",jet_eta,"jet_eta[Jet_n]/F");
    myEvent->Branch("Jet_phi",jet_phi,"jet_phi[Jet_n]/F");
    myEvent->Branch("Jet_emEnergyFraction",jet_emEnergyFraction,"jet_emEnergyFraction[Jet_n]/F");
    myEvent->Branch("Jet_energyFractionHadronic",jet_energyFractionHadronic,"jet_energyFractionHadronic[Jet_n]/F");
    myEvent->Branch("Jet_hitsInN90",jet_hitsInN90,"jet_hitsInN90[Jet_n]/I");
    myEvent->Branch("Jet_n90Hits",jet_n90Hits,"jet_n90Hits[Jet_n]/I");
    myEvent->Branch("Jet_nTowers",jet_nTowers,"jet_nTowers[Jet_n]/I");
    myEvent->Branch("Jet_fHPD",jet_fHPD,"jet_fHPD[Jet_n]/F");
    myEvent->Branch("Jet_fRBX",jet_fRBX,"jet_fRBX[Jet_n]/F");
    myEvent->Branch("Jet_RHF",jet_RHF,"jet_RHF[Jet_n]/F");
    myEvent->Branch("Jet_jecUncer",jet_jecUncer,"jet_jecUncer[Jet_n]/F");
    myEvent->Branch("Jet_jecCorr",jet_jecCorr,"jet_jecCorr[Jet_n]/F");
   //uncorrected jet infor
    myEvent->Branch("ucJet_px",ucjet_px,"ucjet_px[Jet_n]/F");
    myEvent->Branch("ucJet_py",ucjet_py,"ucjet_py[Jet_n]/F");
    myEvent->Branch("ucJet_E",ucjet_E,"ucjet_E[Jet_n]/F");
    myEvent->Branch("ucJet_pz",ucjet_pz,"ucjet_pz[Jet_n]/F");
    myEvent->Branch("ucJet_pt",ucjet_pt,"ucjet_pt[Jet_n]/F");
    myEvent->Branch("ucJet_eta",ucjet_eta,"ucjet_eta[Jet_n]/F");
    myEvent->Branch("ucJet_phi",ucjet_phi,"ucjet_phi[Jet_n]/F");
  }


  if (runpfjets_){
    myEvent->Branch("pfJet_n",&pfJet_n,"pfJet_n/I");
    myEvent->Branch("pfJet_px",pfjet_px,"pfjet_px[pfJet_n]/F");
    myEvent->Branch("pfJet_py",pfjet_py,"pfjet_py[pfJet_n]/F");
    myEvent->Branch("pfJet_E",pfjet_E,"pfjet_E[pfJet_n]/F");
    myEvent->Branch("pfJet_pz",pfjet_pz,"pfjet_pz[pfJet_n]/F");
    myEvent->Branch("pfJet_vx",pfjet_vx,"pfjet_vx[pfJet_n]/F");
    myEvent->Branch("pfJet_vy",pfjet_vy,"pfjet_vy[pfJet_n]/F");
    myEvent->Branch("pfJet_vz",pfjet_vz,"pfjet_vz[pfJet_n]/F");
    myEvent->Branch("pfJet_pt",pfjet_pt,"pfjet_pt[pfJet_n]/F");
    myEvent->Branch("pfJet_eta",pfjet_eta,"pfjet_eta[pfJet_n]/F");
    myEvent->Branch("pfJet_phi",pfjet_phi,"pfjet_phi[pfJet_n]/F");
    /*
    myEvent->Branch("pfJet_emEnergyFraction",pfjet_emEnergyFraction,"pfjet_emEnergyFraction[pfJet_n]/F");
    myEvent->Branch("pfJet_energyFractionHadronic",pfjet_energyFractionHadronic,"pfjet_energyFractionHadronic[pfJet_n]/F");
    myEvent->Branch("pfJet_hitsInN90",pfjet_hitsInN90,"pfjet_hitsInN90[pfJet_n]/I");
    myEvent->Branch("pfJet_n90Hits",pfjet_n90Hits,"pfjet_n90Hits[pfJet_n]/I");
    myEvent->Branch("pfJet_nTowers",pfjet_nTowers,"pfjet_nTowers[pfJet_n]/I");
    myEvent->Branch("pfJet_fHPD",pfjet_fHPD,"pfjet_fHPD[pfJet_n]/F");
    myEvent->Branch("pfJet_fRBX",pfjet_fRBX,"pfjet_fRBX[pfJet_n]/F");
    myEvent->Branch("pfJet_RHF",pfjet_RHF,"pfjet_RHF[pfJet_n]/F");
  */
    myEvent->Branch("pfJet_jecUncer",pfjet_jecUncer,"pfjet_jecUncer[pfJet_n]/F");
    myEvent->Branch("pfJet_jecCorr",pfjet_jecCorr,"pfjet_jecCorr[pfJet_n]/F");
    //uncorrected jet info
    myEvent->Branch("ucpfJet_px",ucpfjet_px,"ucpfjet_px[pfJet_n]/F");
    myEvent->Branch("ucpfJet_py",ucpfjet_py,"ucpfjet_py[pfJet_n]/F");
    myEvent->Branch("ucpfJet_E",ucpfjet_E,"ucpfjet_E[pfJet_n]/F");
    myEvent->Branch("ucpfJet_pz",ucpfjet_pz,"ucpfjet_pz[pfJet_n]/F");
    myEvent->Branch("ucpfJet_pt",ucpfjet_pt,"ucpfjet_pt[pfJet_n]/F");
    myEvent->Branch("ucpfJet_eta",ucpfjet_eta,"ucpfjet_eta[pfJet_n]/F");
    myEvent->Branch("ucpfJet_phi",ucpfjet_phi,"ucpfjet_phi[pfJet_n]/F");
}




  if(rungenjets_){
    myEvent->Branch("genJet_n",&genJet_n,"genJet_n/I");
    myEvent->Branch("genJet_px",genjet_px,"genjet_px[genJet_n]/F");
    myEvent->Branch("genJet_py",genjet_py,"genjet_py[genJet_n]/F");
    myEvent->Branch("genJet_E",genjet_E,"genjet_E[genJet_n]/F");
    myEvent->Branch("genJet_pz",genjet_pz,"genjet_pz[genJet_n]/F");
    myEvent->Branch("genJet_pt",genjet_pt,"genjet_pt[genJet_n]/F");
    myEvent->Branch("genJet_eta",genjet_eta,"genjet_eta[genJet_n]/F");
    myEvent->Branch("genJet_phi",genjet_phi,"genjet_phi[genJet_n]/F");
  }



  
  if (runelectrons_){
    myEvent->Branch("Electron_n",&Electron_n,"Electron_n/I");
    myEvent->Branch("Electron_px",electron_px,"electron_px[Electron_n]/F");
    myEvent->Branch("Electron_py",electron_py,"electron_py[Electron_n]/F");
    myEvent->Branch("Electron_pz",electron_pz,"electron_pz[Electron_n]/F");
    myEvent->Branch("Electron_vx",electron_vx,"electron_vx[Electron_n]/F");
    myEvent->Branch("Electron_vy",electron_vy,"electron_vy[Electron_n]/F");
    myEvent->Branch("Electron_vz",electron_vz,"electron_vz[Electron_n]/F");
    myEvent->Branch("Electron_pt",electron_pt,"electron_pt[Electron_n]/F");
    myEvent->Branch("Electron_eta",electron_eta,"electron_eta[Electron_n]/F");
    myEvent->Branch("Electron_phi",electron_phi,"electron_phi[Electron_n]/F");
    myEvent->Branch("Electron_energy",electron_energy,"electron_energy[Electron_n]/F");
    myEvent->Branch("Electron_charge",electron_charge,"electron_charge[Electron_n]/F");
    myEvent->Branch("Electron_trkIso",electron_trkIso,"electron_trkIso[Electron_n]/F");   
    myEvent->Branch("Electron_ecalIso",electron_ecalIso,"electron_ecalIso[Electron_n]/F");   
    myEvent->Branch("Electron_hcalIso",electron_hcalIso,"electron_hcalIso[Electron_n]/F");   
    myEvent->Branch("Electron_SigmaIetaIeta",electron_SigmaIetaIeta,"electron_SigmaIetaIeta[Electron_n]/F");   
    myEvent->Branch("Electron_dEtaIn",electron_dEtaIn,"electron_dEtaIn[Electron_n]/F");   
    myEvent->Branch("Electron_dPhiIn",electron_dPhiIn,"electron_dPhiIn[Electron_n]/F");   
    myEvent->Branch("Electron_HoE",electron_HoE,"electron_HoE[Electron_n]/F");   
    myEvent->Branch("Electron_sc_energy",electron_sc_energy,"electron_sc_energy[Electron_n]/F");
    myEvent->Branch("Electron_sc_eta",electron_sc_eta,"electron_sc_eta[Electron_n]/F");
    myEvent->Branch("Electron_sc_phi",electron_sc_phi,"electron_sc_phi[Electron_n]/F");
  }
  
  if (runmuons_){
    myEvent->Branch("Muon_n",&Muon_n,"Muon_n/I");
    myEvent->Branch("Muon_px",muon_px,"muon_px[Muon_n]/F");
    myEvent->Branch("Muon_py",muon_py,"muon_py[Muon_n]/F");
    myEvent->Branch("Muon_pz",muon_pz,"muon_pz[Muon_n]/F");
    myEvent->Branch("Muon_vx",muon_vx,"muon_vx[Muon_n]/F");
    myEvent->Branch("Muon_vy",muon_vy,"muon_vy[Muon_n]/F");
    myEvent->Branch("Muon_vz",muon_vz,"muon_vz[Muon_n]/F");
    myEvent->Branch("Muon_pt",muon_pt,"muon_pt[Muon_n]/F");
    myEvent->Branch("Muon_eta",muon_eta,"muon_eta[Muon_n]/F");
    myEvent->Branch("Muon_phi",muon_phi,"muon_phi[Muon_n]/F");
    myEvent->Branch("Muon_energy",muon_energy,"muon_energy[Muon_n]/F");
    myEvent->Branch("Muon_charge",muon_charge,"muon_charge[Muon_n]/F");
    myEvent->Branch("Muon_isGlobalMuon",muon_isGlobalMuon,"muon_isGlobalMuon[Muon_n]/O");
    myEvent->Branch("Muon_isTrackerMuon",muon_isTrackerMuon,"muon_isTrackerMuon[Muon_n]/O");
    myEvent->Branch("Muon_isStandAloneMuon",muon_isStandAloneMuon,"muon_isStandAloneMuon[Muon_n]/O");
    myEvent->Branch("Muon_InnerTrack_isNonnull",muon_InnerTrack_isNonnull,"muon_InnerTrack_isNonnull[Muon_n]/O");
    myEvent->Branch("Muon_OuterTrack_isNonnull",muon_OuterTrack_isNonnull,"muon_OuterTrack_isNonnull[Muon_n]/O");

    myEvent->Branch("Muon_OuterTrack_InnerPoint_x",muon_OuterTrack_InnerPoint_x,"muon_OuterTrack_InnerPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_InnerPoint_y",muon_OuterTrack_InnerPoint_y,"muon_OuterTrack_InnerPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_InnerPoint_z",muon_OuterTrack_InnerPoint_z,"muon_OuterTrack_InnerPoint_z[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_InnerPoint_px",muon_OuterTrack_InnerPoint_px,"muon_OuterTrack_InnerPoint_px[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_InnerPoint_py",muon_OuterTrack_InnerPoint_py,"muon_OuterTrack_InnerPoint_py[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_InnerPoint_pz",muon_OuterTrack_InnerPoint_pz,"muon_OuterTrack_InnerPoint_pz[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_x",muon_OuterTrack_OuterPoint_x,"muon_OuterTrack_OuterPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_y",muon_OuterTrack_OuterPoint_y,"muon_OuterTrack_OuterPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_z",muon_OuterTrack_OuterPoint_z,"muon_OuterTrack_OuterPoint_z[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_px",muon_OuterTrack_OuterPoint_px,"muon_OuterTrack_OuterPoint_px[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_py",muon_OuterTrack_OuterPoint_py,"muon_OuterTrack_OuterPoint_py[Muon_n]/F");
    myEvent->Branch("Muon_OuterTrack_OuterPoint_pz",muon_OuterTrack_OuterPoint_pz,"muon_OuterTrack_OuterPoint_pz[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_x",muon_InnerTrack_InnerPoint_x,"muon_InnerTrack_InnerPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_y",muon_InnerTrack_InnerPoint_y,"muon_InnerTrack_InnerPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_z",muon_InnerTrack_InnerPoint_z,"muon_InnerTrack_InnerPoint_z[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_px",muon_InnerTrack_InnerPoint_px,"muon_InnerTrack_InnerPoint_px[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_py",muon_InnerTrack_InnerPoint_py,"muon_InnerTrack_InnerPoint_py[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_InnerPoint_pz",muon_InnerTrack_InnerPoint_pz,"muon_InnerTrack_InnerPoint_pz[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_x",muon_InnerTrack_OuterPoint_x,"muon_InnerTrack_OuterPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_y",muon_InnerTrack_OuterPoint_y,"muon_InnerTrack_OuterPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_z",muon_InnerTrack_OuterPoint_z,"muon_InnerTrack_OuterPoint_z[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_px",muon_InnerTrack_OuterPoint_px,"muon_InnerTrack_OuterPoint_px[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_py",muon_InnerTrack_OuterPoint_py,"muon_InnerTrack_OuterPoint_py[Muon_n]/F");
    myEvent->Branch("Muon_InnerTrack_OuterPoint_pz",muon_InnerTrack_OuterPoint_pz,"muon_InnerTrack_OuterPoint_pz[Muon_n]/F");  
  
      
  if(isAOD_){
    myEvent->Branch("Muon_OuterPoint_x",muon_OuterPoint_x,"muon_OuterPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_OuterPoint_y",muon_OuterPoint_y,"muon_OuterPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_OuterPoint_z",muon_OuterPoint_z,"muon_OuterPoint_z[Muon_n]/F");
     //for Global,Tracker muon 
    myEvent->Branch("Muon_InnerPoint_x",muon_InnerPoint_x,"muon_InnerPoint_x[Muon_n]/F");
    myEvent->Branch("Muon_InnerPoint_y",muon_InnerPoint_y,"muon_InnerPoint_y[Muon_n]/F");
    myEvent->Branch("Muon_InnerPoint_z",muon_InnerPoint_z,"muon_InnerPoint_z[Muon_n]/F");

    }//isAOD_
  }
  
  if (runcosmicmuons_){
    myEvent->Branch("CosmicMuon_n",&CosmicMuon_n,"CosmicMuon_n/I");
    myEvent->Branch("CosmicMuon_px",cosmicmuon_px,"cosmicmuon_px[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_py",cosmicmuon_py,"cosmicmuon_py[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_pz",cosmicmuon_pz,"cosmicmuon_pz[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_pt",cosmicmuon_pt,"cosmicmuon_pt[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_eta",cosmicmuon_eta,"cosmicmuon_eta[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_phi",cosmicmuon_phi,"cosmicmuon_phi[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_energy",cosmicmuon_energy,"cosmicmuon_energy[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_charge",cosmicmuon_charge,"cosmicmuon_charge[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_isGlobalMuon",cosmicmuon_isGlobalMuon,"cosmicmuon_isGlobalMuon[CosmicMuon_n]/O");
    myEvent->Branch("CosmicMuon_isTrackerMuon",cosmicmuon_isTrackerMuon,"cosmicmuon_isTrackerMuon[CosmicMuon_n]/O");
    myEvent->Branch("CosmicMuon_isStandAloneMuon",cosmicmuon_isStandAloneMuon,"cosmicmuon_isStandAloneMuon[CosmicMuon_n]/O");
    myEvent->Branch("CosmicMuon_InnerTrack_isNonnull",cosmicmuon_InnerTrack_isNonnull,"cosmicmuon_InnerTrack_isNonnull[CosmicMuon_n]/O");
    myEvent->Branch("CosmicMuon_OuterTrack_isNonnull",cosmicmuon_OuterTrack_isNonnull,"cosmicmuon_OuterTrack_isNonnull[CosmicMuon_n]/O");
    
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_x",cosmicmuon_OuterTrack_InnerPoint_x,"cosmicmuon_OuterTrack_InnerPoint_x[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_y",cosmicmuon_OuterTrack_InnerPoint_y,"cosmicmuon_OuterTrack_InnerPoint_y[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_z",cosmicmuon_OuterTrack_InnerPoint_z,"cosmicmuon_OuterTrack_InnerPoint_z[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_px",cosmicmuon_OuterTrack_InnerPoint_px,"cosmicmuon_OuterTrack_InnerPoint_px[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_py",cosmicmuon_OuterTrack_InnerPoint_py,"cosmicmuon_OuterTrack_InnerPoint_py[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_InnerPoint_pz",cosmicmuon_OuterTrack_InnerPoint_pz,"cosmicmuon_OuterTrack_InnerPoint_pz[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_x",cosmicmuon_OuterTrack_OuterPoint_x,"cosmicmuon_OuterTrack_OuterPoint_x[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_y",cosmicmuon_OuterTrack_OuterPoint_y,"cosmicmuon_OuterTrack_OuterPoint_y[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_z",cosmicmuon_OuterTrack_OuterPoint_z,"cosmicmuon_OuterTrack_OuterPoint_z[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_px",cosmicmuon_OuterTrack_OuterPoint_px,"cosmicmuon_OuterTrack_OuterPoint_px[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_py",cosmicmuon_OuterTrack_OuterPoint_py,"cosmicmuon_OuterTrack_OuterPoint_py[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterTrack_OuterPoint_pz",cosmicmuon_OuterTrack_OuterPoint_pz,"cosmicmuon_OuterTrack_OuterPoint_pz[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_x",cosmicmuon_InnerTrack_InnerPoint_x,"cosmicmuon_InnerTrack_InnerPoint_x[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_y",cosmicmuon_InnerTrack_InnerPoint_y,"cosmicmuon_InnerTrack_InnerPoint_y[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_z",cosmicmuon_InnerTrack_InnerPoint_z,"cosmicmuon_InnerTrack_InnerPoint_z[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_px",cosmicmuon_InnerTrack_InnerPoint_px,"cosmicmuon_InnerTrack_InnerPoint_px[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_py",cosmicmuon_InnerTrack_InnerPoint_py,"cosmicmuon_InnerTrack_InnerPoint_py[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_InnerPoint_pz",cosmicmuon_InnerTrack_InnerPoint_pz,"cosmicmuon_InnerTrack_InnerPoint_pz[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_x",cosmicmuon_InnerTrack_OuterPoint_x,"cosmicmuon_InnerTrack_OuterPoint_x[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_y",cosmicmuon_InnerTrack_OuterPoint_y,"cosmicmuon_InnerTrack_OuterPoint_y[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_z",cosmicmuon_InnerTrack_OuterPoint_z,"cosmicmuon_InnerTrack_OuterPoint_z[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_px",cosmicmuon_InnerTrack_OuterPoint_px,"cosmicmuon_InnerTrack_OuterPoint_px[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_py",cosmicmuon_InnerTrack_OuterPoint_py,"cosmicmuon_InnerTrack_OuterPoint_py[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_InnerTrack_OuterPoint_pz",cosmicmuon_InnerTrack_OuterPoint_pz,"cosmicmuon_InnerTrack_OuterPoint_pz[CosmicMuon_n]/F");
   
  if(isAOD_){
    myEvent->Branch("CosmicMuon_OuterPoint_x",cosmicmuon_OuterPoint_x,"cosmicmuon_OuterPoint_x[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterPoint_y",cosmicmuon_OuterPoint_y,"cosmicmuon_OuterPoint_y[CosmicMuon_n]/F");
    myEvent->Branch("CosmicMuon_OuterPoint_z",cosmicmuon_OuterPoint_z,"cosmicmuon_OuterPoint_z[CosmicMuon_n]/F");
    }//AOD_
  }
  
  if (runtaus_){
    myEvent->Branch("Tau_n",&Tau_n,"Tau_n/I");
    myEvent->Branch("Tau_px",tau_px,"tau_px[Tau_n]/F");
    myEvent->Branch("Tau_py",tau_py,"tau_py[Tau_n]/F");
    myEvent->Branch("Tau_pz",tau_pz,"tau_pz[Tau_n]/F");
    myEvent->Branch("Tau_vx",tau_vx,"tau_vx[Tau_n]/F");
    myEvent->Branch("Tau_vy",tau_vy,"tau_vy[Tau_n]/F");
    myEvent->Branch("Tau_vz",tau_vz,"tau_vz[Tau_n]/F");
    myEvent->Branch("Tau_pt",tau_pt,"tau_pt[Tau_n]/F");
    myEvent->Branch("Tau_eta",tau_eta,"tau_eta[Tau_n]/F");
    myEvent->Branch("Tau_phi",tau_phi,"tau_phi[Tau_n]/F");
    myEvent->Branch("Tau_energy",tau_energy,"tau_energy[Tau_n]/F");
    myEvent->Branch("Tau_charge",tau_charge,"tau_charge[Tau_n]/F");
  }
  
  if( rungenParticleCandidates_ ){
    myEvent->Branch("gen_pthat",&gen_pthat,"gen_pthat/F");
	  
    //genlevel information from photons
    myEvent->Branch("ngenphotons",&ngenphotons,"ngenphotons/I");
    myEvent->Branch("gen_photonpt",gen_pho_pt,"gen_pho_pt[ngenphotons]/F");
    myEvent->Branch("gen_photoneta",gen_pho_eta,"gen_pho_eta[ngenphotons]/F");
    myEvent->Branch("gen_photonphi",gen_pho_phi,"gen_pho_phi[ngenphotons]/F");
    myEvent->Branch("gen_photonpx",gen_pho_px,"gen_pho_px[ngenphotons]/F");
    myEvent->Branch("gen_photonpy",gen_pho_py,"gen_pho_py[ngenphotons]/F");
    myEvent->Branch("gen_photonpz",gen_pho_pz,"gen_pho_pz[ngenphotons]/F");
    myEvent->Branch("gen_photonE",gen_pho_E,"gen_pho_E[ngenphotons]/F");
    myEvent->Branch("gen_photonstatus",gen_pho_status,"gen_pho_status[ngenphotons]/I");
    myEvent->Branch("gen_photonMotherID",gen_pho_motherID,"gen_pho_motherID[ngenphotons]/I");
    myEvent->Branch("gen_photonMotherPt",gen_pho_motherPt,"gen_pho_motherPt[ngenphotons]/F");
    myEvent->Branch("gen_photonMotherEta",gen_pho_motherEta,"gen_pho_motherEta[ngenphotons]/F");
    myEvent->Branch("gen_photonMotherPhi",gen_pho_motherPhi,"gen_pho_motherPhi[ngenphotons]/F");
    myEvent->Branch("gen_photonMotherStatus",gen_pho_motherStatus,"gen_pho_motherStatus[ngenphotons]/I");
    myEvent->Branch("gen_photonGrandmotherID",gen_pho_GrandmotherID,"gen_pho_GrandmotherID[ngenphotons]/I");
    myEvent->Branch("gen_photonGrandmotherPt",gen_pho_GrandmotherPt,"gen_pho_GrandmotherPt[ngenphotons]/F");
    myEvent->Branch("gen_photonGrandmotherEta",gen_pho_GrandmotherEta,"gen_pho_GrandmotherEta[ngenphotons]/F");
    myEvent->Branch("gen_photonGrandmotherPhi",gen_pho_GrandmotherPhi,"gen_pho_GrandmotherPhi[ngenphotons]/F");
    myEvent->Branch("gen_photonGrandmotherStatus",gen_pho_GrandmotherStatus,"gen_pho_GrandmotherStatus[ngenphotons]/I");
    
    
    myEvent->Branch("nhardphotons",&nhardphotons,"nhardphotons/I");
    myEvent->Branch("gen_hardphotonpt",gen_Hpho_pt,"gen_Hpho_pt[nhardphotons]/F");
    myEvent->Branch("gen_hardphotoneta",gen_Hpho_eta,"gen_Hpho_eta[nhardphotons]/F");
    myEvent->Branch("gen_hardphotonphi",gen_Hpho_phi,"gen_Hpho_phi[nhardphotons]/F");
    myEvent->Branch("gen_hardphotonpx",gen_Hpho_px,"gen_Hpho_px[nhardphotons]/F");
    myEvent->Branch("gen_hardphotonpy",gen_Hpho_py,"gen_Hpho_py[nhardphotons]/F");
    myEvent->Branch("gen_hardphotonpz",gen_Hpho_pz,"gen_Hpho_pz[nhardphotons]/F");
    myEvent->Branch("gen_hardphotonE",gen_Hpho_E,"gen_Hpho_E[nhardphotons]/F");
    
    //gen level graviton info
    myEvent->Branch("gen_gravitonpt",&gen_graviton_pt,"gen_graviton_pt/F");
    myEvent->Branch("gen_gravitoneta",&gen_graviton_eta,"gen_graviton_eta/F");
    myEvent->Branch("gen_gravitonphi",&gen_graviton_phi,"gen_graviton_phi/F");
    myEvent->Branch("gen_gravitonpx",&gen_graviton_px,"gen_graviton_px/F");
    myEvent->Branch("gen_gravitonpy",&gen_graviton_py,"gen_graviton_py/F");
    myEvent->Branch("gen_gravitonpz",&gen_graviton_pz,"gen_graviton_pz/F");
    myEvent->Branch("gen_gravitonE",&gen_graviton_E,"gen_graviton_E/F");
    
    //genlevel tree info of W+/W- 
    //genlevel tree information of Wdaughter
    myEvent->Branch("gen_Wdaughterpt",gen_Wdaughter_pt,"gen_Wdaughter_pt[2]/F");
    myEvent->Branch("gen_Wdaughtereta",gen_Wdaughter_eta,"gen_Wdaughter_eta[2]/F");
    myEvent->Branch("gen_Wdaughterphi",gen_Wdaughter_phi,"gen_Wdaughter_phi[2]/F");
    myEvent->Branch("gen_Wdaughterpx",gen_Wdaughter_px,"gen_Wdaughter_px[2]/F");
    myEvent->Branch("gen_Wdaughterpy",gen_Wdaughter_py,"gen_Wdaughter_py[2]/F");
    myEvent->Branch("gen_Wdaughterpz",gen_Wdaughter_pz,"gen_Wdaughter_pz[2]/F");
    myEvent->Branch("gen_WdaughterE",gen_Wdaughter_E,"gen_Wdaughter_E[2]/F");
    myEvent->Branch("gen_Wdaughter_charge",gen_Wdaughter_charge,"gen_Wdaughter_charge[2]/I");
    myEvent->Branch("gen_WdaughterID",gen_Wdaughter_ID,"gen_Wdaughter_ID[2]/I");
    
    //genlevel tree information of W
    myEvent->Branch("gen_Wbosonpt",&gen_Wboson_pt,"gen_Wboson_pt/F");
    myEvent->Branch("gen_Wbosoneta",&gen_Wboson_eta,"gen_Wboson_eta/F");
    myEvent->Branch("gen_Wbosonphi",&gen_Wboson_phi,"gen_Wboson_phi/F");
    myEvent->Branch("gen_Wbosonpx",&gen_Wboson_px,"gen_Wboson_px/F");
    myEvent->Branch("gen_Wbosonpy",&gen_Wboson_py,"gen_Wboson_py/F");
    myEvent->Branch("gen_Wbosonpz",&gen_Wboson_pz,"gen_Wboson_pz/F");
    myEvent->Branch("gen_WbosonE",&gen_Wboson_E,"gen_Wboson_E/F");
    myEvent->Branch("gen_Wbosoncharge",&gen_Wboson_charge,"gen_Wboson_charge/I");
    myEvent->Branch("gen_WbosonID",&gen_Wboson_ID,"gen_Wboson_ID/I");
    
    //genlevel tree information of Zdaughter
    myEvent->Branch("gen_Zdaughterpt",gen_Zdaughter_pt,"gen_Zdaughter_pt[2]/F");
    myEvent->Branch("gen_Zdaughtereta",gen_Zdaughter_eta,"gen_Zdaughter_eta[2]/F");
    myEvent->Branch("gen_Zdaughterphi",gen_Zdaughter_phi,"gen_Zdaughter_phi[2]/F");
    myEvent->Branch("gen_Zdaughterpx",gen_Zdaughter_px,"gen_Zdaughter_px[2]/F");
    myEvent->Branch("gen_Zdaughterpy",gen_Zdaughter_py,"gen_Zdaughter_py[2]/F");
    myEvent->Branch("gen_Zdaughterpz",gen_Zdaughter_pz,"gen_Zdaughter_pz[2]/F");
    myEvent->Branch("gen_ZdaughterE",gen_Zdaughter_E,"gen_Zdaughter_E[2]/F");
    myEvent->Branch("gen_Zdaughter_charge",gen_Zdaughter_charge,"gen_Zdaughter_charge[2]/I");
    myEvent->Branch("gen_ZdaughterID",gen_Zdaughter_ID,"gen_Zdaughter_ID[2]/I");
    
    //genlevel tree information of Z
    myEvent->Branch("gen_Zbosonpt",&gen_Zboson_pt,"gen_Zboson_pt/F");
    myEvent->Branch("gen_Zbosoneta",&gen_Zboson_eta,"gen_Zboson_eta/F");
    myEvent->Branch("gen_Zbosonphi",&gen_Zboson_phi,"gen_Zboson_phi/F");
    myEvent->Branch("gen_Zbosonpx",&gen_Zboson_px,"gen_Zboson_px/F");
    myEvent->Branch("gen_Zbosonpy",&gen_Zboson_py,"gen_Zboson_py/F");
    myEvent->Branch("gen_Zbosonpz",&gen_Zboson_pz,"gen_Zboson_pz/F");
    myEvent->Branch("gen_ZbosonE",&gen_Zboson_E,"gen_Zboson_E/F");
    
    myEvent->Branch("is_signal_event",&is_signal_event,"is_signal_event/O");
    myEvent->Branch("is_Z_event",&is_Z_event,"is_Z_event/O");
    myEvent->Branch("is_W_event",&is_W_event,"is_W_event/O");
    myEvent->Branch("is_Znunu_event",&is_Znunu_event,"is_Znunu_event/O");
    myEvent->Branch("is_Zelec_event",&is_Zelec_event,"is_Zelec_event/O");
    myEvent->Branch("is_Zmu_event",&is_Zmu_event,"is_Zmu_event/O");      
    myEvent->Branch("is_Ztau_event",&is_Ztau_event,"is_Ztau_event/O");   
    myEvent->Branch("is_Welec_event",&is_Welec_event,"is_Welec_event/O");
    myEvent->Branch("is_Wmu_event",&is_Wmu_event,"is_Wmu_event/O");      
    myEvent->Branch("is_Wtau_event",&is_Wtau_event,"is_Wtau_event/O");   
    myEvent->Branch("is_SingleHardPhoton_event",&is_SingleHardPhoton_event,"is_SingleHardPhoton_event/O");
    myEvent->Branch("is_diphoton_event",&is_diphoton_event,"is_diphoton_event/O");
    myEvent->Branch("is_isr_photon_event",&is_isr_photon_event,"is_isr_photon_event/O");
	  
    myEvent->Branch("n_signal_events",&n_signal_events,"n_signal_events/I");
    myEvent->Branch("n_Z_events",&n_Z_events,"n_Z_events/I");
    myEvent->Branch("n_W_events",&n_W_events,"n_W_events/I");
    myEvent->Branch("n_Znunu_events",&n_Znunu_events,"n_Znunu_events/I"); 
    myEvent->Branch("n_Zelec_events",&n_Zelec_events,"n_Zelec_events/I"); 
    myEvent->Branch("n_Zmu_events",&n_Zmu_events,"n_Zmu_events/I");       
    myEvent->Branch("n_Ztau_events",&n_Ztau_events,"n_Ztau_events/I");    
    myEvent->Branch("n_Welec_events",&n_Welec_events,"n_Welec_events/I"); 
    myEvent->Branch("n_Wmu_events",&n_Wmu_events,"n_Wmu_events/I");       
    myEvent->Branch("n_Wtau_events",&n_Wtau_events,"n_Wtau_events/I");    
    myEvent->Branch("n_SingleHardPhoton_events",&n_SingleHardPhoton_events,"n_SingleHardPhoton_events/I");
    myEvent->Branch("n_diphoton_events",&n_diphoton_events,"n_diphoton_events/I");
    
    //genlevel tree information of mu daughter
    myEvent->Branch("gen_MuonID",gen_Muon_ID,"gen_Muon_ID[3]/F");
    myEvent->Branch("gen_MuonStatus",gen_Muon_Status,"gen_Muon_Status[3]/F");
    myEvent->Branch("gen_MuonPt",gen_Muon_Pt,"gen_Muon_Pt[3]/F");
    myEvent->Branch("gen_MuonDaughterpt",gen_MuonDaughter_pt,"gen_MuonDaughter_pt[3]/F");
    myEvent->Branch("gen_MuonDaughtereta",gen_MuonDaughter_eta,"gen_MuonDaughter_eta[3]/F");
    myEvent->Branch("gen_MuonDaughterphi",gen_MuonDaughter_phi,"gen_MuonDaughter_phi[3]/F");
    myEvent->Branch("gen_MuonDaughterpx",gen_MuonDaughter_px,"gen_MuonDaughter_px[3]/F");
    myEvent->Branch("gen_MuonDaughterpy",gen_MuonDaughter_py,"gen_MuonDaughter_py[3]/F");
    myEvent->Branch("gen_MuonDaughterpz",gen_MuonDaughter_pz,"gen_MuonDaughter_pz[3]/F");
    myEvent->Branch("gen_MuonDaughterE",gen_MuonDaughter_E,"gen_MuonDaughter_E[3]/F");
    myEvent->Branch("gen_MuonDaughterCharge",gen_MuonDaughter_charge,"gen_MuonDaughter_charge[3]/I");
    myEvent->Branch("gen_MuonDaughterStatus",gen_MuonDaughter_status,"gen_MuonDaughter_status[3]/I");
    myEvent->Branch("gen_MuonDaughterID",gen_MuonDaughter_ID,"gen_MuonDaughter_ID[3]/I");
    
    //genlevel tree information of tau daughter
    myEvent->Branch("gen_tauID",gen_tau_ID,"gen_tau_ID[3]/F");
    myEvent->Branch("gen_tauStatus",gen_tau_Status,"gen_tau_Status[3]/F");
    myEvent->Branch("gen_tauPt",gen_tau_Pt,"gen_tau_Pt[3]/F");
    myEvent->Branch("gen_tauDaughterpt",gen_tauDaughter_pt,"gen_tauDaughter_pt[3]/F");
    myEvent->Branch("gen_tauDaughtereta",gen_tauDaughter_eta,"gen_tauDaughter_eta[3]/F");
    myEvent->Branch("gen_tauDaughterphi",gen_tauDaughter_phi,"gen_tauDaughter_phi[3]/F");
    myEvent->Branch("gen_tauDaughterpx",gen_tauDaughter_px,"gen_tauDaughter_px[3]/F");
    myEvent->Branch("gen_tauDaughterpy",gen_tauDaughter_py,"gen_tauDaughter_py[3]/F");
    myEvent->Branch("gen_tauDaughterpz",gen_tauDaughter_pz,"gen_tauDaughter_pz[3]/F");
    myEvent->Branch("gen_tauDaughterE",gen_tauDaughter_E,"gen_tauDaughter_E[3]/F");
    myEvent->Branch("gen_tauDaughterCharge",gen_tauDaughter_charge,"gen_tauDaughter_charge[3]/I");
    myEvent->Branch("gen_tauDaughterStatus",gen_tauDaughter_status,"gen_tauDaughter_status[3]/I");
    myEvent->Branch("gen_tauDaughterID",gen_tauDaughter_ID,"gen_tauDaughter_ID[3]/I");
    
  }//end of if( rungenParticleCandidates_ )
  
  if (runphotons_){
    //uncorrected photon information
    myEvent->Branch("Photon_n",&Photon_n,"Photon_n/I");
    myEvent->Branch("Photon_E",pho_E,"pho_E[Photon_n]/F");
    myEvent->Branch("Photon_pt",pho_pt,"pho_pt[Photon_n]/F");
    myEvent->Branch("Photon_eta",pho_eta,"pho_eta[Photon_n]/F");
    myEvent->Branch("Photon_phi",pho_phi,"pho_phi[Photon_n]/F");
    myEvent->Branch("Photon_theta",pho_theta,"pho_theta[Photon_n]/F");
    myEvent->Branch("Photon_et",pho_et,"pho_et[Photon_n]/F");
    myEvent->Branch("Photon_swissCross",pho_swissCross,"pho_swissCross[Photon_n]/F");
    myEvent->Branch("Photonr9",pho_r9,"pho_r9[Photon_n]/F");
    myEvent->Branch("Photon_e1x5",pho_e1x5,"pho_e1x5[Photon_n]/F");
    myEvent->Branch("Photon_e2x5",pho_e2x5,"pho_e2x5[Photon_n]/F");
    myEvent->Branch("Photon_e3x3",pho_e3x3,"pho_e3x3[Photon_n]/F");
    myEvent->Branch("Photon_e5x5",pho_e5x5,"pho_e5x5[Photon_n]/F");
    myEvent->Branch("Photon_r1x5",pho_r1x5,"pho_erx5[Photon_n]/F");
    myEvent->Branch("Photon_r2x5",pho_r2x5,"pho_erx5[Photon_n]/F");
    myEvent->Branch("Photon_maxEnergyXtal",pho_maxEnergyXtal,"pho_maxEnergyXtal[Photon_n]/F");
    myEvent->Branch("Photon_SigmaEtaEta",pho_SigmaEtaEta,"pho_SigmaEtaEta[Photon_n]/F");
    myEvent->Branch("Photon_SigmaIetaIeta",pho_SigmaIetaIeta,"pho_SigmaIetaIeta[Photon_n]/F");
    myEvent->Branch("Photon_SigmaEtaPhi",pho_SigmaEtaPhi,"pho_SigmaEtaPhi[Photon_n]/F");
    myEvent->Branch("Photon_SigmaIetaIphi",pho_SigmaIetaIphi,"pho_SigmaIetaIphi[Photon_n]/F");
    myEvent->Branch("Photon_SigmaPhiPhi",pho_SigmaPhiPhi,"pho_SigmaPhiPhi[Photon_n]/F");
    myEvent->Branch("Photon_SigmaIphiIphi",pho_SigmaIphiIphi,"pho_SigmaIphiIphi[Photon_n]/F");
    myEvent->Branch("Photon_Roundness",pho_roundness,"pho_roundness[Photon_n]/F");
    myEvent->Branch("Photon_Angle",pho_angle,"pho_angle[Photon_n]/F");
    myEvent->Branch("Photon_ecalRecHitSumEtConeDR03",pho_ecalRecHitSumEtConeDR03,"pho_ecalRecHitSumEtConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_hcalTowerSumEtConeDR03",pho_hcalTowerSumEtConeDR03,"pho_hcalTowerSumEtConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_trkSumPtSolidConeDR03",pho_trkSumPtSolidConeDR03,"pho_trkSumPtSolidConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_trkSumPtHollowConeDR03",pho_trkSumPtHollowConeDR03,"pho_trkSumPtHollowConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_nTrkSolidConeDR03",pho_nTrkSolidConeDR03,"pho_nTrkSolidConeDR03[Photon_n]/I");
    myEvent->Branch("Photon_nTrkHollowConeDR03",pho_nTrkHollowConeDR03,"pho_nTrkHollowConeDR03[Photon_n]/I");
    myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR03",pho_hcalDepth1TowerSumEtConeDR03,"pho_hcalDepth1TowerSumEtConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR03",pho_hcalDepth2TowerSumEtConeDR03,"pho_hcalDepth2TowerSumEtConeDR03[Photon_n]/F");
    myEvent->Branch("Photon_ecalRecHitSumEtConeDR04",pho_ecalRecHitSumEtConeDR04,"pho_ecalRecHitSumEtConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_hcalTowerSumEtConeDR04",pho_hcalTowerSumEtConeDR04,"pho_hcalTowerSumEtConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_trkSumPtSolidConeDR04",pho_trkSumPtSolidConeDR04,"pho_trkSumPtSolidConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_trkSumPtHollowConeDR04",pho_trkSumPtHollowConeDR04,"pho_trkSumPtHollowConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_nTrkSolidConeDR04",pho_nTrkSolidConeDR04,"pho_nTrkSolidConeDR04[Photon_n]/I");
    myEvent->Branch("Photon_nTrkHollowConeDR04",pho_nTrkHollowConeDR04,"pho_nTrkHollowConeDR04[Photon_n]/I");
    myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR04",pho_hcalDepth1TowerSumEtConeDR04,"pho_hcalDepth1TowerSumEtConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR04",pho_hcalDepth2TowerSumEtConeDR04,"pho_hcalDepth2TowerSumEtConeDR04[Photon_n]/F");
    myEvent->Branch("Photon_hasPixelSeed",pho_hasPixelSeed,"pho_hasPixelSeed[Photon_n]/O"); 
    myEvent->Branch("Photon_isEB",pho_isEB,"pho_isEB[Photon_n]/O");
    myEvent->Branch("Photon_isEE",pho_isEE,"pho_isEE[Photon_n]/O");
    myEvent->Branch("Photon_isEBGap",pho_isEBGap,"pho_isEBGap[Photon_n]/O");
    myEvent->Branch("Photon_isEEGap",pho_isEEGap,"pho_isEEGap[Photon_n]/O");
    myEvent->Branch("Photon_isEBEEGap",pho_isEBEEGap,"pho_isEBEEGap[Photon_n]/O");
    myEvent->Branch("Photon_e2e9",pho_e2e9,"pho_e2e9[Photon_n]/F");
    
    myEvent->Branch("Photon_HoE",pho_HoE,"pho_HoE[Photon_n]/F");
    myEvent->Branch("Photon_px",pho_px,"pho_px[Photon_n]/F");
    myEvent->Branch("Photon_py",pho_py,"pho_py[Photon_n]/F");
    myEvent->Branch("Photon_pz",pho_pz,"pho_pz[Photon_n]/F");
    myEvent->Branch("Photon_vx",pho_vx,"pho_vx[Photon_n]/F");
    myEvent->Branch("Photon_vy",pho_vy,"pho_vy[Photon_n]/F");
    myEvent->Branch("Photon_vz",pho_vz,"pho_vz[Photon_n]/F");
    myEvent->Branch("Photon_no_of_basic_clusters",pho_size,"pho_size[Photon_n]/I");
    
    myEvent->Branch("Photon_sc_energy",pho_sc_energy,"pho_sc_energy[Photon_n]/F");
    myEvent->Branch("Photon_sc_eta",pho_sc_eta,"pho_sc_eta[Photon_n]/F");
    myEvent->Branch("Photon_sc_phi",pho_sc_phi,"pho_sc_phi[Photon_n]/F");
    myEvent->Branch("Photon_sc_x",pho_sc_x,"pho_sc_x[Photon_n]/F");
    myEvent->Branch("Photon_sc_y",pho_sc_y,"pho_sc_y[Photon_n]/F");
    myEvent->Branch("Photon_sc_z",pho_sc_z,"pho_sc_z[Photon_n]/F");

    myEvent->Branch("Photon_etaWidth",pho_sc_etaWidth,"pho_sc_etaWidth[Photon_n]/F");
    myEvent->Branch("Photon_phiWidth",pho_sc_phiWidth,"pho_sc_phiWidth[Photon_n]/F");
    myEvent->Branch("Photon_sc_et",pho_sc_et,"pho_sc_et[Photon_n]/F");

    myEvent->Branch("matchphotonE",matchpho_E,"matchpho_E[Photon_n]/F");
    myEvent->Branch("matchphotonpt",matchpho_pt,"matchpho_pt[Photon_n]/F");
    myEvent->Branch("matchphotoneta",matchpho_eta,"matchpho_eta[Photon_n]/F");
    myEvent->Branch("matchphotonphi",matchpho_phi,"matchpho_phi[Photon_n]/F");
    myEvent->Branch("matchphotonpx",matchpho_px,"matchpho_px[Photon_n]/F");
    myEvent->Branch("matchphotonpy",matchpho_py,"matchpho_py[Photon_n]/F");
    myEvent->Branch("matchphotonpz",matchpho_pz,"matchpho_pz[Photon_n]/F");
    myEvent->Branch("ismatchedphoton",ismatchedpho,"ismatchedpho[Photon_n]/O");
    
    myEvent->Branch("Photon_hasConvTrk",pho_hasConvTrk,"pho_hasConvTrk[Photon_n]/O"); 
    myEvent->Branch("Photon_ntracks",pho_nTracks,"pho_nTracks[Photon_n]/I");
    myEvent->Branch("Photon_isconverted",pho_isConverted,"pho_isConverted[Photon_n]/O");
    myEvent->Branch("Photon_pairInvmass",pho_pairInvariantMass,"pho_pairInvariantMass[Photon_n]/F");
    myEvent->Branch("Photon_pairCotThetaSeperation",pho_pairCotThetaSeparation,"pho_pairCotThetaSeparation[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumX",pho_pairMomentum_x,"pho_pairMomentum_x[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumY",pho_pairMomentum_y,"pho_pairMomentum_y[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumZ",pho_pairMomentum_z,"pho_pairMomentum_z[Photon_n]/F");
    myEvent->Branch("Photon_EoverP",pho_EoverP,"pho_EoverP[Photon_n]/F");
    myEvent->Branch("Photon_ConvVx",pho_conv_vx,"pho_conv_vx[Photon_n]/F");
    myEvent->Branch("Photon_ConvVy",pho_conv_vy,"pho_conv_vy[Photon_n]/F");
    myEvent->Branch("Photon_ConvVz",pho_conv_vz,"pho_conv_vz[Photon_n]/F");
    myEvent->Branch("Photon_ZOfPrimaryVertex",pho_zOfPrimaryVertex,"pho_zOfPrimaryVertex[Photon_n]/F");
    myEvent->Branch("Photon_distOfMinimumApproach",pho_distOfMinimumApproach,"pho_distOfMinimumApproach[Photon_n]/F");
    myEvent->Branch("Photon_dPhiTracksAtVtx",pho_dPhiTracksAtVtx,"pho_dPhiTracksAtVtx[Photon_n]/F");
    myEvent->Branch("Photon_dPhiTracksAtEcal",pho_dPhiTracksAtEcal,"pho_dPhiTracksAtEcal[Photon_n]/F");
    myEvent->Branch("Photon_dEtaTracksAtEcal",pho_dEtaTracksAtEcal,"pho_dEtaTracksAtEcal[Photon_n]/F");
    if(runrechit_){
      myEvent->Branch("Photon_ncrys",ncrysPhoton,"ncrysPhoton[Photon_n]/I");
      myEvent->Branch("Photon_timing_xtal",pho_timing_xtal,"pho_timing_xtal[Photon_n][100]/F");
      myEvent->Branch("Photon_timingavg_xtal",pho_timingavg_xtal,"pho_timingavg_xtal[Photon_n]/F");
      myEvent->Branch("Photon_energy_xtal",pho_energy_xtal,"pho_energy_xtal[Photon_n][100]/F");
      myEvent->Branch("Photon_ieta_xtalEB",pho_ieta_xtalEB,"pho_ieta_xtalEB[Photon_n][100]/I");
      myEvent->Branch("Photon_iphi_xtalEB",pho_iphi_xtalEB,"pho_iphi_xtalEB[Photon_n][100]/I");
      myEvent->Branch("Photon_s9",pho_s9,"pho_s9[Photon_n]/F");
    }
    
    if(runHErechit_){
      myEvent->Branch("HERecHit_subset_n",&HERecHit_subset_n,"HERecHit_subset_n/I");
      myEvent->Branch("HERecHit_subset_detid",HERecHit_subset_detid,"HERecHit_subset_detid[HERecHit_subset_n]/i");
      myEvent->Branch("HERecHit_subset_energy",HERecHit_subset_energy,"HERecHit_subset_energy[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_time",HERecHit_subset_time,"HERecHit_subset_time[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_depth",HERecHit_subset_depth,"HERecHit_subset_depth[HERecHit_subset_n]/I");
      myEvent->Branch("HERecHit_subset_phi",HERecHit_subset_phi,"HERecHit_subset_phi[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_eta",HERecHit_subset_eta,"HERecHit_subset_eta[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_x",HERecHit_subset_x,"HERecHit_subset_x[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_y",HERecHit_subset_y,"HERecHit_subset_y[HERecHit_subset_n]/F");
      myEvent->Branch("HERecHit_subset_z",HERecHit_subset_z,"HERecHit_subset_z[HERecHit_subset_n]/F");
    }
    
  }//end of if (runphotons_)

  if(runrechit_){
    myEvent->Branch("EBRecHit_size",&EBRecHit_size,"EBRecHit_size/I");
    myEvent->Branch("EBRecHit_eta",EBRecHit_eta,"EBRecHit_eta[EBRecHit_size]/F");
    myEvent->Branch("EBRecHit_phi",EBRecHit_phi,"EBRecHit_phi[EBRecHit_size]/F");
    myEvent->Branch("EBRecHit_ieta",EBRecHit_ieta,"EBRecHit_ieta[EBRecHit_size]/I");
    myEvent->Branch("EBRecHit_iphi",EBRecHit_iphi,"EBRecHit_iphi[EBRecHit_size]/I");
    myEvent->Branch("EBRecHit_e",EBRecHit_e,"EBRecHit_e[EBRecHit_size]/F");
    myEvent->Branch("EBRecHit_et",EBRecHit_et,"EBRecHit_et[EBRecHit_size]/F");
    myEvent->Branch("EBRecHit_flag",EBRecHit_flag,"EBRecHit_flag[EBRecHit_size]/I");
    myEvent->Branch("EBRecHit_time",EBRecHit_time,"EBRecHit_time[EBRecHit_size]/F");

    myEvent->Branch("EERecHit_size",&EERecHit_size,"EERecHit_size/I");
    myEvent->Branch("EERecHit_eta",EERecHit_eta,"EERecHit_eta[EERecHit_size]/F");
    myEvent->Branch("EERecHit_phi",EERecHit_phi,"EERecHit_phi[EERecHit_size]/F");
    myEvent->Branch("EERecHit_ieta",EERecHit_ieta,"EERecHit_ieta[EERecHit_size]/I");
    myEvent->Branch("EERecHit_iphi",EERecHit_iphi,"EERecHit_iphi[EERecHit_size]/I");
    myEvent->Branch("EERecHit_e",EERecHit_e,"EERecHit_e[EERecHit_size]/F");
    myEvent->Branch("EERecHit_et",EERecHit_et,"EERecHit_et[EERecHit_size]/F");
    myEvent->Branch("EERecHit_flag",EERecHit_flag,"EERecHit_flag[EERecHit_size]/I");
    myEvent->Branch("EERecHit_time",EERecHit_time,"EERecHit_time[EERecHit_size]/F");
  }//end of if(runFroMIP_)
	
  if(runCSCseg_){
    myEvent->Branch("CSCseg_n", &CSCseg_n, "CSCseg_n/I");
    myEvent->Branch("CSCseg_time", CSCseg_time, "CSCseg_time[CSCseg_n]/F");
    myEvent->Branch("CSCseg_x", CSCseg_x, "CSCseg_x[CSCseg_n]/F");
    myEvent->Branch("CSCseg_y", CSCseg_y, "CSCseg_y[CSCseg_n]/F");
    myEvent->Branch("CSCseg_z", CSCseg_z, "CSCseg_z[CSCseg_n]/F");
    myEvent->Branch("CSCseg_phi", CSCseg_phi, "CSCseg_phi[CSCseg_n]/F");
    myEvent->Branch("CSCseg_DirectionX", CSCseg_DirectionX, "CSCseg_DirectionX[CSCseg_n]/F");
    myEvent->Branch("CSCseg_DirectionY", CSCseg_DirectionY, "CSCseg_DirectionY[CSCseg_n]/F");
    myEvent->Branch("CSCseg_DirectionZ", CSCseg_DirectionZ, "CSCseg_DirectionZ[CSCseg_n]/F");
  }// end of runCSCeg_

         
  if(runBeamHaloSummary_){                                                                                                       
   myEvent->Branch("isBeamHaloGlobalLoosePass",&isBeamHaloGlobalLoosePass,"isBeamHaloGlobalLoosePass/O");
   myEvent->Branch("isBeamHaloGlobalTightPass",&isBeamHaloGlobalTightPass,"isBeamHaloGloablTightPass/O");
   myEvent->Branch("isBeamHaloHcalLoosePass",&isBeamHaloHcalLoosePass,"isBeamHaloHcalLoosePass/O");
   myEvent->Branch("isBeamHaloHcalTightPass",&isBeamHaloHcalTightPass,"isBeamHaloHcalTightPass/O");
   myEvent->Branch("isBeamHaloCSCLoosePass",&isBeamHaloCSCLoosePass,"isBeamHaloCSCLoosePass/O");
   myEvent->Branch("isBeamHaloCSCTightPass",&isBeamHaloCSCTightPass,"isBeamHaloCSCTightPass/O");
   myEvent->Branch("isBeamHaloEcalLoosePass",&isBeamHaloEcalLoosePass,"isBeamHaloEcalLoosePass/O");
   myEvent->Branch("isBeamHaloEcalTightPass",&isBeamHaloEcalTightPass,"isBeamHaloEcalTightPass/O");
   myEvent->Branch("isBeamHaloIDTightPass",&isBeamHaloIDTightPass,"isBeamHaloIDTightPass/O");
   myEvent->Branch("isBeamHaloIDLoosePass",&isBeamHaloIDLoosePass,"isBeamHaloIDLoosePass/O");
   myEvent->Branch("isSmellsLikeHalo_Tag",&isSmellsLikeHalo_Tag, "isSmellsLikeHalo_Tag/O");
   myEvent->Branch("isLooseHalo_Tag",&isLooseHalo_Tag, "isLooseHalo_Tag/O");
   myEvent->Branch("isTightHalo_Tag",&isTightHalo_Tag, "isTightHalo_Tag/O");
   myEvent->Branch("isExtremeTightHalo_Tag",&isExtremeTightHalo_Tag, "isExtremeTightHalo_Tag/O");
  }//if(runBeamHaloSummary_)


  if(runRPChit_){
    myEvent->Branch("RPChit_n", &RPChit_n, "RPChit_n/I");
    myEvent->Branch("RPChit_x", RPChit_x, "RPChit_x[RPChit_n]/F");
    myEvent->Branch("RPChit_y", RPChit_y, "RPChit_y[RPChit_n]/F");
    myEvent->Branch("RPChit_z", RPChit_z, "RPChit_z[RPChit_n]/F");
    myEvent->Branch("RPChit_BunchX", RPChit_BunchX, "RPChit_BunchX[RPChit_n]/I");
  }//end of runRPChit_

  if(runmet_){
    //Calomet variables
    myEvent->Branch("CaloMetSigma",&CaloMetSig,"CaloMetSig/F");
    myEvent->Branch("CaloMetEz",&CaloMetEz,"CaloMetEz/F");
    myEvent->Branch("CaloEtFractionHadronic",&CaloEtFractionHadronic,"CaloEtFractionHadronic/F");
    myEvent->Branch("CaloEmEtFraction",&CaloEmEtFraction,"CaloEmEtFraction/F");
    myEvent->Branch("CaloHadEtInHB",&CaloHadEtInHB,"CaloHadEtInHB/F");
    myEvent->Branch("CaloHadEtInHE",&CaloHadEtInHE,"CaloHadEtInHE/F");
    myEvent->Branch("CaloHadEtInHO",&CaloHadEtInHO,"CaloHadEtInHO/F");
    myEvent->Branch("CaloHadEtInHF",&CaloHadEtInHF,"CaloHadEtInHF/F");
    myEvent->Branch("CaloEmEtInEB",&CaloEmEtInEB,"CaloEmEtInEB/F");
    myEvent->Branch("CaloEmEtInEE",&CaloEmEtInEE,"CaloEmEtInEE/F");
    myEvent->Branch("CaloEmEtInHF",&CaloEmEtInHF,"CaloEmEtInHF/F");
    myEvent->Branch("CaloMaxEtInEmTowers",&CaloMaxEtInEmTowers,"CaloMaxEtInEmTowers/F");
    myEvent->Branch("CaloMaxEtInHadTowers",&CaloMaxEtInHadTowers,"CaloMaxEtInHadTowers/F");
    myEvent->Branch("CaloMetPt",CaloMetPt,"CaloMetPt[6]/F");
    myEvent->Branch("CaloMetPx",CaloMetPx,"CaloMetPx[6]/F");
    myEvent->Branch("CaloMetPy",CaloMetPy,"CaloMetPy[6]/F");
    myEvent->Branch("CaloMetPhi",CaloMetPhi,"CaloMetPhi[6]/F");
    myEvent->Branch("CaloMetSumEt",CaloMetSumEt,"CaloMetSumEt[6]/F");

    if(rungenmet_){
      myEvent->Branch("genMetPt",&genMetPt,"genMetPt/F");
      myEvent->Branch("genMetPx",&genMetPx,"genMetPx/F");
      myEvent->Branch("genMetPy",&genMetPy,"genMetPy/F");
      myEvent->Branch("genMetPhi",&genMetPhi,"genMetPhi/F");
      myEvent->Branch("genMetSumEt",&genMetSumEt,"genMetSumEt/F");
    }
  }//end of if(runmet)
  if(runmet_&& runphotons_)
    myEvent->Branch("Delta_phi",&Delta_phi,"Delta_phi/F");
  if(runmet_ && runphotons_ && rungenmet_)
    myEvent->Branch("Delta_phiGEN",&Delta_phiGEN,"Delta_phiGEN/F");
  
  if(runPFmet_){
    myEvent->Branch("PFMetPt",PFMetPt,"PFMetPt[6]/F");
    myEvent->Branch("PFMetPx",PFMetPx,"PFMetPx[6]/F");
    myEvent->Branch("PFMetPy",PFMetPy,"PFMetPy[6]/F");
    myEvent->Branch("PFMetPhi",PFMetPhi,"PFMetPhi[6]/F");
    myEvent->Branch("PFMetSumEt",PFMetSumEt,"PFMetSumEt[6]/F");
  }//end of if(runmet)
  if(runPFmet_&& runphotons_)
    myEvent->Branch("Delta_phiPF",&Delta_phiPF,"Delta_phiPF/F");
  
  
  if(runTCmet_){
    myEvent->Branch("TCMetPt",TCMetPt,"TCMetPt[6]/F");
    myEvent->Branch("TCMetPx",TCMetPx,"TCMetPx[6]/F");
    myEvent->Branch("TCMetPy",TCMetPy,"TCMetPy[6]/F");
    myEvent->Branch("TCMetPhi",TCMetPhi,"TCMetPhi[6]/F");
    myEvent->Branch("TCMetSumEt",TCMetSumEt,"TCMetSumEt[6]/F");
  }//end of if(runmet)
  if(runTCmet_&& runphotons_)
    myEvent->Branch("Delta_phiTC",&Delta_phiTC,"Delta_phiTC/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {
  f->WriteTObject(myEvent);
  delete myEvent;
}

double Analyzer::GetE2OverE9( const DetId id, const EcalRecHitCollection & recHits)
{ ///////////start calculating e2/e9
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240
  // compute e2overe9
  //   | | | |
  //   +-+-+-+
  //   | |1|2|
  //   +-+-+-+
  //   | | | |
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
//   rechit 1 must have E_t > recHitEtThreshold
//   rechit 2 must have E_t > recHitEtThreshold2
//   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
//   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
//   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0


float recHitEtThreshold = 10.0; 
float recHitEtThreshold2 = 1.0;
bool avoidIeta85=false;
bool KillSecondHit=true;

if ( id.subdetId() == EcalBarrel ) {
  EBDetId ebId( id );
  // avoid recHits at |eta|=85 where one side of the neighbours is missing
  if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;
  // select recHits with Et above recHitEtThreshold
  float e1 = recHitE( id, recHits );
  float ete1=recHitApproxEt( id, recHits );
  // check that rechit E_t is above threshold
  if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;
  if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;
  
  float e2=-1;
  float ete2=0;
  float s9 = 0;
  // coordinates of 2nd hit relative to central hit
  int e2eta=0;
  int e2phi=0;

  // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
  for ( int deta = -1; deta <= +1; ++deta ) {
    for ( int dphi = -1; dphi <= +1; ++dphi ) {
      // compute 3x3 energy 
      float etmp=recHitE( id, recHits, deta, dphi );
      s9 += etmp;
      EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
      float eapproxet=recHitApproxEt( idtmp, recHits );
      // remember 2nd highest energy deposit (above threshold) in 3x3 array
      if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	e2=etmp;
	ete2=eapproxet;
	e2eta=deta;
	e2phi=dphi;
      }
    }
  }

  if ( e1 == 0 )  return 0;
  // return 0 if 2nd hit is below threshold
  if ( e2 == -1 ) return 0;
  // compute e2/e9 centered around 1st hit
  float e2nd=e1+e2;
  float e2e9=0;

  if (s9!=0) e2e9=e2nd/s9;
  // if central hit has higher energy than 2nd hit
  //  return e2/e9 if 1st hit is above E_t threshold
  if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;
  // if second hit has higher energy than 1st hit
  if ( e2 > e1 ) {
    // return 0 if user does not want to flag 2nd hit, or
    // hits are below E_t thresholds - note here we
    // now assume the 2nd hit to be the leading hit.
    
    if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
      return 0;
    }
    else {
      // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2 
      float s92nd=0;
      float e2nd_prime=0;
      int e2prime_eta=0;
      int e2prime_phi=0;

      EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

      for ( int deta = -1; deta <= +1; ++deta ) {
	for ( int dphi = -1; dphi <= +1; ++dphi ) {

	  // compute 3x3 energy
	  float etmp=recHitE( secondid, recHits, deta, dphi );
	  s92nd += etmp;

	  if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	    e2nd_prime=etmp;
	    e2prime_eta=deta;
	    e2prime_phi=dphi;
	  }
	}
      }
      // if highest energy hit around E2 is not the same as the input hit, return 0;
      if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	{
	  return 0;
	}
      // compute E2/E9 around second hit 
      float e2e9_2=0;
      if (s92nd!=0) e2e9_2=e2nd/s92nd;
      //   return the value of E2/E9 calculated around 2nd hit
      return e2e9_2;
    }
  }
 } else if ( id.subdetId() == EcalEndcap ) {
  // only used for EB at the moment 
  return 0;
 }
return 0;
}
//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
