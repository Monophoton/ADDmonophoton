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
// $Id: Analyzer.cc,v 1.28 2010/10/07 17:14:59 miceli Exp $
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

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

#include <map>
#include <string>

#include "ADDmonophoton/Analyzer/interface/Analyzer.h"

#include <iostream>
using namespace std;
using namespace ROOT::Math::VectorUtil ;

//utility function prototypes
float correct_phi(float phi);
float Theta(float eta);
float Pl(float P,float Pt);

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
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  PFmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PFmetTag")),
  TCmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("TCmetTag")),
  phoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag")),
  rechitBLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitBTag")),
  rechitELabel_(iConfig.getUntrackedParameter<edm::InputTag>("rechitETag")),
  hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults")),
  Tracks_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks")),
  Vertices_(iConfig.getUntrackedParameter<edm::InputTag>("Vertices")),
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
  runtaus_(iConfig.getUntrackedParameter<bool>("runtaus")),
  runHLT_(iConfig.getUntrackedParameter<bool>("runHLT")),
  runL1_(iConfig.getUntrackedParameter<bool>("runL1")),
  runscraping_(iConfig.getUntrackedParameter<bool>("runscraping")),
  runtracks_(iConfig.getUntrackedParameter<bool>("runtracks")),
  runrechit_(iConfig.getUntrackedParameter<bool>("runrechit")),
  runHErechit_(iConfig.getUntrackedParameter<bool>("runHErechit")),
  runvertex_(iConfig.getUntrackedParameter<bool>("runvertex")),
  init_(false)
{
  //now do what ever initialization is needed
  nevents =0;
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
  Electron_n       = 0; 
  Muon_n           = 0;
  CosmicMuon_n	 = 0;
  Tau_n            = 0;
}


Analyzer::~Analyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
  f->Close();
  delete f;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  bool debug_on = false;
  if(debug_on) cout<<"DEBUG: new event"<<endl;
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
	if(debug_on) cout<<"DEBUG: start saving stuff"<<endl;
  nevents++;
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
    std::vector<genPho>            mygenphoton_container;
    mygenphoton_container.clear();
    genPho genphoton;
    int ii =0;
    for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++) {
      //getting information from hard scattered Graviton 
      if(debug_on) cout<<"start gen particle collection event:"<<nevents<<" particle:"<<ii<<endl;
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
	//if(debug_on) cout<<" getting information from Z now"<<endl;
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
	if(debug_on) cout<<" W motherID:" << genparticle->mother()->pdgId()<<endl; 
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
	if(debug_on) cout<<"W pt:"<<genparticle->pt()<<endl;
	if(debug_on) cout<<"W daughters:"<<endl;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  if(abs(daughter->pdgId())==11) {is_Welec_event=true; n_Welec_events++;}
	  if(abs(daughter->pdgId())==13) {is_Wmu_event=true  ; n_Wmu_events++  ;}
	  if(abs(daughter->pdgId())==15) {is_Wtau_event=true ; n_Wtau_events++ ;}  
	  if(debug_on) cout<<"ID, Status,Pt:"<<abs(daughter->pdgId())<<"   "<<daughter->status()<<"   "<<daughter->pt()<<endl;
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
	if(debug_on) cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	int daughters   = genparticle->numberOfDaughters();
	int iDaughter=0;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  if(debug_on) cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
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
	if(debug_on) cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	int daughters   = genparticle->numberOfDaughters();
	int iDaughter=0;
	for(int i = 0;i<daughters;i++){
	  const reco::Candidate *daughter   = genparticle->daughter(i);
	  //if(debug_on) cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
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
      if(debug_on) cout<<"DEBUG 6"<<endl;
      if(debug_on) cout<<"getting gen photon information now"<<endl;
      //getting information from all photons
      //if (genparticle->pdgId()==22 && genparticle->status()==1 && genparticle->pt()>5.)
      if (genparticle->pdgId()==22){ 
	//doing it this way, as I want to sort it in pt after filling everything in the container
        if(genparticle->mother()!=NULL){
	  const reco::Candidate *mom = genparticle->mother();
	  genphoton.motherID = mom->pdgId();
	  genphoton.motherStatus = mom->status();
	  if(debug_on) cout<<"DEBUG mom pdgid is "<<mom->pdgId()<<" status is "<<mom->status()<<" pt "<<mom->pt()<< " e "<<mom->energy() <<endl;
	  genphoton.motherPt = mom->pt();
	  genphoton.motherEta = mom->eta();
	  genphoton.motherPhi = correct_phi(mom->phi());
	  if(genparticle->status()==1 && mom->pdgId()==2212 && mom->mother()==NULL)
	    is_isr_photon_event = true;
	  if(mom->mother()!=NULL){	  
	    if(debug_on) cout<<"DEBUG gramma is "<<mom->mother()->pdgId()<<endl;
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
	if(debug_on) cout<<"DEBUG done with pho"<<endl;
      }//end of if (genparticle->pdgId()==22 && genparticle->status()==1)
      if(debug_on) cout<<"DEBUG before checking mothers event:"<< nevents<<" particle:"<<ii<<endl;  
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
      if(debug_on) cout<<"DEBUG end loop over gen particles"<<endl;
      ii++;
    }//end of for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++)
    if(debug_on) cout<<"DEBUG: finish gen"<<endl;	
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
	//std::if(debug_on) cout<<"genphoton motherID:"<<mygenphoton_container[x].motherID<<endl;
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
	  if(debug_on) cout<<"[x]:"<<x<< "pho_pt:" << gen_pho_pt[x]<<" pho_status:"<< gen_pho_status[x]<<" motherID:" << gen_pho_motherID[x]<<" motherStatus: "<< gen_pho_motherStatus[x]<< " GrandmotherID: "<< gen_pho_GrandmotherID[x]<< " GrandmotherStatus:"<< gen_pho_GrandmotherStatus[x]<< endl;
        */
	// if(debug_on) cout<<"got the photon info right"<<endl;
      }//end of for loop
    }//end of if((mygenphoton_container.size()!=0)
    //if(debug_on) cout<<"mygenphoton_container loop ended"<<std::endl; 
  }//end of if(rungenParticleCandidates_)
  if(debug_on) cout<<"DEBUG: finish saving gen"<<endl;
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
     //Redone by AA, new interface in 3_6 for trigger information
     Handle<TriggerResults> HLTR;
     iEvent.getByLabel(hlTriggerResults_,HLTR);
     if (HLTR.isValid()) {
       //Get indexes of trigger in this event
       const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*HLTR);
       Int_t idx1 = triggerNames_.triggerIndex("HLT_MET50");
       Int_t idx2 = triggerNames_.triggerIndex("HLT_MET75");
       Int_t idx3 = triggerNames_.triggerIndex("HLT_Photon15_L1R");
       Int_t idx4 = triggerNames_.triggerIndex("HLT_Photon25_L1R");
       Int_t idx5 = triggerNames_.triggerIndex("HLT_DoubleEle10_SW_L1R");
       Int_t idx6 = triggerNames_.triggerIndex("HLC_DoubleMu3");
       Int_t idx7 = triggerNames_.triggerIndex("HLT_Photon20_L1R");
       Int_t idx8 = triggerNames_.triggerIndex("HLT_Photon20_Cleaned_L1R");
       Int_t idx9 = triggerNames_.triggerIndex("HLT_Photon30_L1R");
       Int_t idx9_1 = triggerNames_.triggerIndex("HLT_Photon30_L1R_8E29");
       Int_t idx9_2 = triggerNames_.triggerIndex("HLT_Photon30_L1R_1E31");
       Int_t idx10 = triggerNames_.triggerIndex("HLT_Photon30_Cleaned_L1R");
       Int_t idx11 = triggerNames_.triggerIndex("HLT_Photon30_Isol_EBOnly_Cleaned_L1R");
       Int_t idx12 = triggerNames_.triggerIndex("HLT_Photon35_Isol_Cleaned_L1R");
       Int_t idx13 = triggerNames_.triggerIndex("HLT_Photon50_Cleaned_L1R");
       Int_t idx14 = triggerNames_.triggerIndex("HLT_Photon70_NoHE_Cleaned_L1R");
       Int_t idx15 = triggerNames_.triggerIndex("HLT_Photon100_NoHE_Cleaned_L1R");
       Int_t idx16 = triggerNames_.triggerIndex("HLT_DoublePhoton17_L1R");//new
       Int_t idx17 = triggerNames_.triggerIndex("HLT_DoublePhoton5_CEP_L1R");
       Int_t idx18 = triggerNames_.triggerIndex("HLT_Photon100_NoHE_Cleaned_L1R_v1");
       Int_t idx19 = triggerNames_.triggerIndex("HLT_Photon10_Cleaned_L1R");
       Int_t idx20 = triggerNames_.triggerIndex("HLT_Photon15_Cleaned_L1R");
       Int_t idx21 = triggerNames_.triggerIndex("HLT_Photon17_SC17HE_L1R_v1");
       Int_t idx22 = triggerNames_.triggerIndex("HLT_Photon20_NoHE_L1R");
       Int_t idx23 = triggerNames_.triggerIndex("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1");
       Int_t idx24 = triggerNames_.triggerIndex("HLT_Photon35_Isol_Cleaned_L1R_v1");
       Int_t idx25 = triggerNames_.triggerIndex("HLT_Photon50_Cleaned_L1R_v1");
       Int_t idx26 = triggerNames_.triggerIndex("HLT_Photon50_NoHE_L1R");
       Int_t idx27 = triggerNames_.triggerIndex("HLT_Photon70_NoHE_Cleaned_L1R_v1");

       Int_t hsize = Int_t(HLTR->size());
       //There's a double check here: first check that the array index is in range, then check to see if it fired.
       //if it is out of range, then obviously it wasn't found, and didn't pass (typically means name wasn't in the HLT)

       //initialize HLT triggers to false
       HLT_MET50_event = false;
       HLT_MET75_event = false;
       HLT_Photon15_event = false;
       HLT_Photon25_event = false;
       HLT_DoubleEle10_event = false;
       HLT_DoubleMu3_event = false;
       HLT_Photon20_event = false;
       HLT_Photon20_Cleaned_event = false;
       HLT_Photon30_event = false;
       HLT_Photon30_Cleaned_event =false;//idx10
       HLT_Photon30_Isol_EBOnly_Cleaned_event=false;//idx11
       HLT_Photon35_Isol_Cleaned_event=false;//idx12
       HLT_Photon50_Cleaned_event=false;//idx13
       HLT_Photon70_NoHE_Cleaned_event=false;//idx14
       HLT_Photon100_NoHE_Cleaned_event=false;//idx15
       HLT_Photon30_L1R_8E29_event=false;
       HLT_Photon30_L1R_1E31_event=false;
       HLT_DoublePhoton17_L1R_event = false;//idx16
       HLT_DoublePhoton5_CEP_L1R_event = false;//idx17
       HLT_Photon100_NoHE_Cleaned_L1R_v1_event = false;//idx18
       HLT_Photon10_Cleaned_L1R_event = false;//idx19
       HLT_Photon15_Cleaned_L1R_event = false;//idx20
       HLT_Photon17_SC17HE_L1R_v1_event = false;//idx21
       HLT_Photon20_NoHE_L1R_event = false;//idx22
       HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event = false;//idx23
       HLT_Photon35_Isol_Cleaned_L1R_v1_event = false;//idx24
       HLT_Photon50_Cleaned_L1R_v1_event = false;//idx25
       HLT_Photon50_NoHE_L1R_event = false;//idx26
       HLT_Photon70_NoHE_Cleaned_L1R_v1_event = false;//idx27       

       if (idx1 < hsize)
	 if (HLTR->accept(idx1))
	   HLT_MET50_event = true;
       if (idx2 < hsize)
	 if (HLTR->accept(idx2))
	   HLT_MET75_event = true;
       if (idx3 < hsize)
	 if (HLTR->accept(idx3))
	   HLT_Photon15_event = true;
       if (idx4 < hsize)
	 if (HLTR->accept(idx4))
	   HLT_Photon25_event = true;
       if (idx5 < hsize)
	 if (HLTR->accept(idx5))
	   HLT_DoubleEle10_event = true;
       if (idx6 < hsize)
	 if (HLTR->accept(idx6))
	   HLT_DoubleMu3_event = true;
       if (idx7 < hsize)
	 if (HLTR->accept(idx7))
	   HLT_Photon20_event = true;
       if (idx8 < hsize)
	 if (HLTR->accept(idx8))
	   HLT_Photon20_Cleaned_event = true;
       if (idx9 < hsize)
	 if (HLTR->accept(idx9))
	   HLT_Photon30_event=true;
       if (idx9_1 < hsize)
	 if (HLTR->accept(idx9_1))
	   HLT_Photon30_L1R_8E29_event=true;
       if (idx9_2 < hsize)
	 if (HLTR->accept(idx9_2))
	   HLT_Photon30_L1R_1E31_event=true;
       if (idx10 < hsize)
	 if (HLTR->accept(idx10))
	   HLT_Photon30_Cleaned_event=true;
       if (idx11 < hsize)
	 if (HLTR->accept(idx11))
	   HLT_Photon30_Isol_EBOnly_Cleaned_event=true;
       if (idx12 < hsize)
	 if (HLTR->accept(idx12))
	   HLT_Photon35_Isol_Cleaned_event=true;
       if (idx13 < hsize)
	 if (HLTR->accept(idx13))
	   HLT_Photon50_Cleaned_event=true;
       if (idx14 < hsize)
	 if (HLTR->accept(idx14))
	   HLT_Photon70_NoHE_Cleaned_event=true;
       if (idx15 < hsize)
	 if (HLTR->accept(idx15))
	   HLT_Photon100_NoHE_Cleaned_event=true;
       if (idx16 < hsize)
         if (HLTR->accept(idx16))
           HLT_DoublePhoton17_L1R_event=true;
       if (idx17 < hsize)
         if (HLTR->accept(idx17))
           HLT_DoublePhoton5_CEP_L1R_event=true;
       if (idx18 < hsize)
         if (HLTR->accept(idx18))
           HLT_Photon100_NoHE_Cleaned_L1R_v1_event=true;
       if (idx19 < hsize)
         if (HLTR->accept(idx19))
           HLT_Photon10_Cleaned_L1R_event=true;
       if (idx20 < hsize)
         if (HLTR->accept(idx20))
           HLT_Photon15_Cleaned_L1R_event=true;
       if (idx21 < hsize)
         if (HLTR->accept(idx21))
           HLT_Photon17_SC17HE_L1R_v1_event=true;
       if (idx22 < hsize)
         if (HLTR->accept(idx22))
           HLT_Photon20_NoHE_L1R_event=true;
       if (idx23 < hsize)
         if (HLTR->accept(idx23))
           HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event=true;
       if (idx24 < hsize)
         if (HLTR->accept(idx24))
           HLT_Photon35_Isol_Cleaned_L1R_v1_event=true;
       if (idx25 < hsize)
         if (HLTR->accept(idx25))
           HLT_Photon50_Cleaned_L1R_v1_event=true;
       if (idx26 < hsize)
         if (HLTR->accept(idx26))
           HLT_Photon50_NoHE_L1R_event=true;
       if (idx27 < hsize)
         if (HLTR->accept(idx27))
           HLT_Photon70_NoHE_Cleaned_L1R_v1_event=true;

     }
   }
   
   
   if(runvertex_){
     Handle<reco::VertexCollection> recVtxs;
     iEvent.getByLabel(Vertices_, recVtxs);
     vector<Vertex> my_vertices;
	   
     my_vertices.clear();
     for(reco::VertexCollection::const_iterator v=recVtxs->begin();v!=recVtxs->end(); ++v){
       my_vertices.push_back(*v);
     }
     Vertex_n = 0;  
     for (unsigned int y = 0; y <  my_vertices.size();y++){
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
    iEvent.getByLabel("generalTracks",tkRef);    
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
	 Track_iter != tracks->end();++Track_iter) {
       if(Track_iter->pt()>5.){
	 myTrack_container.push_back(*Track_iter);
       }
     }
     
     
     if(myTrack_container.size()>1)
       std::sort(myTrack_container.begin(),myTrack_container.end(),PtSortCriterium3());
     Track_n = 0;
     for(unsigned int x=0;x < myTrack_container.size();x++){
       trk_pt[x]  = myTrack_container[x].pt();
       trk_px[x]  = myTrack_container[x].px();
       trk_py[x]  = myTrack_container[x].py();
       trk_pz[x]  = myTrack_container[x].pz();
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
     for(unsigned int x=0;x < mymuon_container.size();x++){
       muon_pt[x]  = mymuon_container[x].pt();
       muon_energy[x]  = mymuon_container[x].energy();
       muon_px[x]  = mymuon_container[x].px();
       muon_py[x]  = mymuon_container[x].py();
       muon_pz[x]  = mymuon_container[x].pz();
       muon_phi[x] = correct_phi(mymuon_container[x].phi());
       muon_eta[x] = mymuon_container[x].eta();
       muon_charge[x] = mymuon_container[x].charge();
       
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
  		
       muon_InnerTrack_isNonnull[x] =   mymuon_container[x].innerTrack().isNonnull();
       muon_OuterTrack_isNonnull[x] =   mymuon_container[x].outerTrack().isNonnull();
       
			
       if(mymuon_container[x].innerTrack().isNonnull()){
	 muon_InnerTrack_InnerPoint_x[x] = mymuon_container[x].innerTrack()->innerPosition().x();
	 muon_InnerTrack_InnerPoint_y[x] = mymuon_container[x].innerTrack()->innerPosition().y();
	 muon_InnerTrack_InnerPoint_z[x] = mymuon_container[x].innerTrack()->innerPosition().z();
	 muon_InnerTrack_InnerPoint_px[x] = mymuon_container[x].innerTrack()->innerMomentum().x();
	 muon_InnerTrack_InnerPoint_py[x] = mymuon_container[x].innerTrack()->innerMomentum().y();
	 muon_InnerTrack_InnerPoint_pz[x] = mymuon_container[x].innerTrack()->innerMomentum().z();
	 muon_InnerTrack_OuterPoint_x[x] = mymuon_container[x].innerTrack()->innerPosition().x();
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
	 muon_OuterTrack_OuterPoint_x[x] = mymuon_container[x].outerTrack()->innerPosition().x();
	 muon_OuterTrack_OuterPoint_y[x] = mymuon_container[x].outerTrack()->outerPosition().y();
	 muon_OuterTrack_OuterPoint_z[x] = mymuon_container[x].outerTrack()->outerPosition().z();
	 muon_OuterTrack_OuterPoint_px[x] = mymuon_container[x].outerTrack()->outerMomentum().x();
	 muon_OuterTrack_OuterPoint_py[x] = mymuon_container[x].outerTrack()->outerMomentum().y();
	 muon_OuterTrack_OuterPoint_pz[x] = mymuon_container[x].outerTrack()->outerMomentum().z();
       }
       Muon_n++;
     }//end of for loop
   }// if runmuons_
   
   if(runcosmicmuons_){
     //cosmic muon
     edm::Handle<reco::MuonCollection>cosmicMuonHandle;
     iEvent.getByLabel("muonsFromCosmics",cosmicMuonHandle);
     const reco::MuonCollection & cosmicmuons = *cosmicMuonHandle;
     vector <reco::Muon> mycosmicmuon_container;
     
     for(reco::MuonCollection::const_iterator cosmuon = cosmicmuons.begin(); cosmuon!=cosmicmuons.end(); ++cosmuon){
       //if(muon->isGlobalMuon())
       mycosmicmuon_container.push_back(*cosmuon);
     }
     CosmicMuon_n = 0;
     for(unsigned int x=0;x < mycosmicmuon_container.size();x++){
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
       
       cosmicmuon_InnerTrack_isNonnull[x] =   mycosmicmuon_container[x].innerTrack().isNonnull();
       cosmicmuon_OuterTrack_isNonnull[x] =   mycosmicmuon_container[x].outerTrack().isNonnull();
       
       if(mycosmicmuon_container[x].innerTrack().isNonnull()){
	 cosmicmuon_InnerTrack_InnerPoint_x[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().x();
	 cosmicmuon_InnerTrack_InnerPoint_y[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().y();
	 cosmicmuon_InnerTrack_InnerPoint_z[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().z();
	 cosmicmuon_InnerTrack_InnerPoint_px[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().x();
	 cosmicmuon_InnerTrack_InnerPoint_py[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().y();
	 cosmicmuon_InnerTrack_InnerPoint_pz[x] = mycosmicmuon_container[x].innerTrack()->innerMomentum().z();
	 cosmicmuon_InnerTrack_OuterPoint_x[x] = mycosmicmuon_container[x].innerTrack()->innerPosition().x();
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
	 cosmicmuon_OuterTrack_OuterPoint_x[x] = mycosmicmuon_container[x].outerTrack()->innerPosition().x();
	 cosmicmuon_OuterTrack_OuterPoint_y[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().y();
	 cosmicmuon_OuterTrack_OuterPoint_z[x] = mycosmicmuon_container[x].outerTrack()->outerPosition().z();
	 cosmicmuon_OuterTrack_OuterPoint_px[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().x();
	 cosmicmuon_OuterTrack_OuterPoint_py[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().y();
	 cosmicmuon_OuterTrack_OuterPoint_pz[x] = mycosmicmuon_container[x].outerTrack()->outerMomentum().z();
       }
       CosmicMuon_n++;
     }//end of for loop
   }//if runcosmicmuons_
	
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
       for(unsigned int x=0; x < myphoton_container.size();x++){
	 pho_E[x]                     =  myphoton_container[x].energy();
	 pho_pt[x]                    =  myphoton_container[x].pt();
	 pho_px[x]                    =  myphoton_container[x].px();
	 pho_py[x]                    =  myphoton_container[x].py();
	 pho_pz[x]                    =  myphoton_container[x].pz();
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
	 ismatchedpho[x]                         =  myphoton_container[x].genParticleRef().isNonnull();
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
	     pho_vertex_x[x]                   = aConv->conversionVertex().x();
	     pho_vertex_y[x]                   = aConv->conversionVertex().y();
	     pho_vertex_z[x]                   = aConv->conversionVertex().z();
	     pho_EoverP[x]                     = aConv->EoverP();
	     pho_zOfPrimaryVertex[x]           = aConv->zOfPrimaryVertexFromTracks();
	     pho_distOfMinimumApproach[x]      = aConv->distOfMinimumApproach();
	     pho_dPhiTracksAtVtx[x]            = aConv->dPhiTracksAtVtx();
	     pho_dPhiTracksAtEcal[x]           = aConv->dPhiTracksAtEcal();
	     pho_dEtaTracksAtEcal[x]           = aConv->dEtaTracksAtEcal();
	   }//end of if ( aConv->conversionVertex().isValid() )
	   else{
	     pho_nTracks[x]                    = 9999;
	     pho_isConverted[x]                = false;
	     pho_pairInvariantMass[x]          = -99.;
	     pho_pairCotThetaSeparation[x]     = -99.;
	     pho_pairMomentum_x[x]             = -99.;
	     pho_pairMomentum_y[x]             = -99.;
	     pho_pairMomentum_z[x]             = -99.;
	     pho_vertex_x[x]                   = -99.;
	     pho_vertex_y[x]                   = -99.;
	     pho_vertex_z[x]                   = -99.;
	     pho_EoverP[x]                     = -99.;
	     pho_zOfPrimaryVertex[x]           = -99.;
	     pho_distOfMinimumApproach[x]      = -99.;
	     pho_dPhiTracksAtVtx[x]            = -99.;
	     pho_dPhiTracksAtEcal[x]           = -99.;
	     pho_dEtaTracksAtEcal[x]           = -99.;
	   }//end of else
	 }//end of for (unsigned int iConv=0; iConv<conversions.size(); iConv++)
	 
	 if(runHErechit_ && pho_isEB[x]){
	   //Store HE rechits
	   edm::Handle<HBHERecHitCollection> hcalRecHitHandle;
	   iEvent.getByLabel(edm::InputTag("hbhereco"), hcalRecHitHandle);
	   const HBHERecHitCollection *hbhe =  hcalRecHitHandle.product();
	   
	   for(HBHERecHitCollection::const_iterator hh = hbhe->begin(); hh != hbhe->end() && HERecHit_subset_n<10000; hh++){
	     HcalDetId id(hh->detid());
	     if (id.subdet()==2){
	       edm::ESHandle<CaloGeometry> geoHandle;	       
	       iSetup.get<CaloGeometryRecord>().get(geoHandle);
	       const CaloGeometry* caloGeom = geoHandle.product();
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
	 }
	 Photon_n++;
       }//end of for loop over x
     }//if(myphoton_container.size!=0) 
     //to get the photon hit information from every crystal of SC
     if(runrechit_){ 
       Handle<EcalRecHitCollection> Brechit;//barrel
       Handle<EcalRecHitCollection> Erechit;//endcap
       iEvent.getByLabel(rechitBLabel_,Brechit);
       iEvent.getByLabel(rechitELabel_,Erechit);
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
	   for (unsigned int y =0; y < crystalinfo_container.size();y++){
	     pho_timing_xtal[x][y]         = crystalinfo_container[y].time;
	     pho_energy_xtal[x][y]         = crystalinfo_container[y].energy;
	     pho_ieta_xtalEB[x][y]           = crystalinfo_container[y].ieta;
	     pho_iphi_xtalEB[x][y]           = crystalinfo_container[y].iphi;
	   }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++
	   if(myphoton_container[x].isEB()){
	     std::vector<float> showershapes_barrel = EcalClusterTools::roundnessBarrelSuperClusters(*(myphoton_container[x].superCluster()),*barrelRecHits,0);
	     pho_roundness[x]    = (float)showershapes_barrel[0];
	     pho_angle[x]        = (float)showershapes_barrel[1];
	     pho_s9[x]           = pho_energy_xtal[x][0]/pho_e3x3[x];
	     pho_swissCross[x]   = EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), 
							   &(*barrelRecHits), &(*topology))+ 
	                           EcalClusterTools::eBottom( *(myphoton_container[x].superCluster()->seed()), 
							   &(*barrelRecHits), &(*topology))+ 
	                           EcalClusterTools::eLeft( *(myphoton_container[x].superCluster()->seed()), 
							   &(*barrelRecHits), &(*topology)) + 
	                           EcalClusterTools::eRight( *(myphoton_container[x].superCluster()->seed()), 
							   &(*barrelRecHits), &(*topology));
	     if(debug_on && 1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) 
	       cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl; 
	   }//end of if(myphoton_container[x].isEB())
	   else{ 
	     pho_roundness[x]   = -99.;
	     pho_angle[x]       = -99.;
	     pho_s9[x]          = pho_energy_xtal[x][0]/pho_e3x3[x];
	     pho_swissCross[x]   = EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), 
							   &(*endcapRecHits), &(*topology))+ 
	                           EcalClusterTools::eBottom( *(myphoton_container[x].superCluster()->seed()), 
							   &(*endcapRecHits), &(*topology)) + 
	                           EcalClusterTools::eLeft( *(myphoton_container[x].superCluster()->seed()), 
							   &(*endcapRecHits), &(*topology)) + 
	                           EcalClusterTools::eRight( *(myphoton_container[x].superCluster()->seed()), 
							   &(*endcapRecHits), &(*topology));
	     if(debug_on && 1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) {
	       cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm." << endl;
	       cout<<"This would be weird since there aren't spikes in the endcap of ECAL"<<endl; 
	     }
	   }//end of else (if !EB)
	 }//end of for loop over x
       }//if(myphoton_container.size!=0) 
     }//if(runrechit_)
     
   }//if(runphotons_)  
   if(runmet_){
     edm::Handle<edm::View<pat::MET> > metHandle;
     iEvent.getByLabel(metLabel_,metHandle);
     //const edm::View<pat::MET> & mets = *metHandle;
     edm::View<pat::MET>::const_iterator met;
     for ( met = metHandle->begin(); met != metHandle->end(); met++){
       //calomet variables
       CaloMetSig                            = met->mEtSig();
       CaloMetEt                             = met->pt();
       CaloMetEx                             = met->px();
       CaloMetEy                             = met->py();
       CaloMetPhi                            = correct_phi(met->phi());
       CaloEtFractionHadronic                = met->etFractionHadronic();
       CaloEmEtFraction                      = met->emEtFraction();
       CaloHadEtInHB                         = met->hadEtInHB();
       CaloHadEtInHO                         = met->hadEtInHO();
       CaloHadEtInHE                         = met->hadEtInHE();
       CaloHadEtInHF                         = met->hadEtInHF();
       CaloEmEtInEB                          = met->emEtInEB();
       CaloEmEtInEE                          = met->emEtInEE();
       CaloEmEtInHF                          = met->emEtInHF();
       CaloMetEz                             = met->e_longitudinal();
       CaloMaxEtInEmTowers                   = met->maxEtInEmTowers();
       CaloSumEt                             = met->sumEt();
       CaloMaxEtInHadTowers                  = met->maxEtInHadTowers();
       if(runphotons_==1)
	 if (myphoton_container.size()!=0)
	   Delta_phi                         = fabs(reco::deltaPhi(correct_phi(met->phi()),correct_phi(myphoton_container[0].phi())));
       if(rungenmet_){
	 const reco::GenMET *genMet = met->genMET();
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
       PFMetPt     = metsPF[0].et();
       PFMetPhi    = correct_phi(metsPF[0].phi());
       PFMetSumEt  = metsPF[0].sumEt();
       PFMetPx     = metsPF[0].px();
       PFMetPy     = metsPF[0].py();
       if(runphotons_==1)
	 if (myphoton_container.size()!=0)
	   Delta_phiPF  = fabs(reco::deltaPhi(PFMetPhi,correct_phi(myphoton_container[0].phi())));
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
       TCMetPt     = metsTC[0].et();
       TCMetPhi    = correct_phi(metsTC[0].phi());
       TCMetSumEt  = metsTC[0].sumEt();
       TCMetPx     = metsTC[0].px();
       TCMetPy     = metsTC[0].py();
       if(runphotons_==1)
	 if (myphoton_container.size()!=0)
	   Delta_phiTC  = fabs(reco::deltaPhi(TCMetPhi,correct_phi(myphoton_container[0].phi()))); 
     }
     else{
       LogWarning("METEventSelector") << "No Met results for InputTag " ;
       return;
     } 
     
   }//end of if(runTCmet_)
   if(runjets_){
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
     if(myjet_container.size()!=0){
       for(unsigned int x=0;x < myjet_container.size();x++){
	 jet_pt[x]  = myjet_container[x].pt();
	 jet_px[x]  = myjet_container[x].px();
	 jet_py[x]  = myjet_container[x].py();
	 jet_pz[x]  = myjet_container[x].pz();
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
	 Jet_n++;
       }//end of for loop
     }
   }
   if(runelectrons_){
     edm::Handle<edm::View<pat::Electron> > electronHandle;
     iEvent.getByLabel(eleLabel_,electronHandle);
     vector<pat::Electron> myelectron_container;
     
     const edm::View<pat::Electron> & electrons = *electronHandle;   // const ... &, we don't make a copy of it!
     for(edm::View<pat::Electron>::const_iterator electron = electrons.begin(); electron!=electrons.end(); ++electron){
       myelectron_container.push_back(*electron);
     }
     Electron_n = 0;
     for(unsigned int x=0;x < myelectron_container.size();x++){
       electron_pt[x]  = myelectron_container[x].pt();
       electron_energy[x]  = myelectron_container[x].energy();
       electron_px[x]  = myelectron_container[x].px();
       electron_py[x]  = myelectron_container[x].py();
       electron_pz[x]  = myelectron_container[x].pz();
       electron_phi[x] = correct_phi(myelectron_container[x].phi());
       electron_eta[x] = myelectron_container[x].eta();
       electron_charge[x] = myelectron_container[x].charge();
       electron_trkIso[x] = myelectron_container[x].trackIso();
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
     for(unsigned int x=0;x < mytau_container.size();x++){
       tau_pt[x]  = mytau_container[x].pt();
       tau_energy[x]  = mytau_container[x].energy();
       tau_px[x]  = mytau_container[x].px();
       tau_py[x]  = mytau_container[x].py();
       tau_pz[x]  = mytau_container[x].pz();
       tau_phi[x] = correct_phi(mytau_container[x].phi());
       tau_eta[x] = mytau_container[x].eta();
       tau_charge[x] = mytau_container[x].charge();
       Tau_n++;
     }//end of for loop
   }
   if (Photon_n>0)
     myEvent->Fill();
   if(debug_on) cout<<"DEBUG: analyze loop done"<<endl;
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob(){
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
  if(runHLT_){
    myEvent->Branch("HLT_MET50_event",&HLT_MET50_event,"HLT_MET50_event/O");
    myEvent->Branch("HLT_MET75_event",&HLT_MET75_event,"HLT_MET75_event/O");
    myEvent->Branch("HLT_Photon15_event",&HLT_Photon15_event,"HLT_Photon15_event/O");
    myEvent->Branch("HLT_Photon25_event",&HLT_Photon25_event,"HLT_Photon25_event/O");
    myEvent->Branch("HLT_DoubleEle10_event",&HLT_DoubleEle10_event,"HLT_DoubleEle10_event/O");
    myEvent->Branch("HLT_DoubleMu3_event",&HLT_DoubleMu3_event,"HLT_DoubleMu3_event/O");
    myEvent->Branch("HLT_Photon20_event", &HLT_Photon20_event,"HLT_Photon20_event/O");
    myEvent->Branch("HLT_Photon20_Cleaned_event",&HLT_Photon20_Cleaned_event,"HLT_Photon20_Cleaned_event/O");
    myEvent->Branch("HLT_Photon30_event", &HLT_Photon30_event,"HLT_Photon30_event/O");  
    myEvent->Branch("HLT_Photon30_L1R_8E29_event", &HLT_Photon30_L1R_8E29_event, "HLT_Photon30_L1R_8E29_event/O");
    myEvent->Branch("HLT_Photon30_L1R_1E31_event", &HLT_Photon30_L1R_1E31_event, "HLT_Photon30_L1R_1E31_event/O");
    myEvent->Branch("HLT_Photon30_Cleaned_event", &HLT_Photon30_Cleaned_event,"HLT_Photon30_Cleaned_event/O");
    myEvent->Branch("HLT_Photon30_Isol_EBOnly_Cleaned", &HLT_Photon30_Isol_EBOnly_Cleaned_event,"HLT_Photon30_Isol_EBOnly_event/O");
    myEvent->Branch("HLT_Photon35_Isol_Cleaned", &HLT_Photon35_Isol_Cleaned_event,"HLT_Photon35_Isol_Cleaned_event/O");
    myEvent->Branch("HLT_Photon50_Cleaned_event", &HLT_Photon50_Cleaned_event, "HLT_Photon50_Cleaned_event/O");
    myEvent->Branch("HLT_Photon70_NoHE_Cleaned_event",&HLT_Photon70_NoHE_Cleaned_event, "HLT_Photon70_NoHE_Cleaned_event/O");
    myEvent->Branch("HLT_Photon100_NoHE_Cleaned_event", &HLT_Photon100_NoHE_Cleaned_event,"HLT_Photon100_NoHE_Cleaned_event/O");
    myEvent->Branch("HLT_DoublePhoton17_L1R_event", & HLT_DoublePhoton17_L1R_event,"HLT_DoublePhoton17_L1R_event/O");
    myEvent->Branch("HLT_DoublePhoton5_CEP_L1R_event", & HLT_DoublePhoton5_CEP_L1R_event,"HLT_DoublePhoton5_CEP_L1R_event/O");
    myEvent->Branch("HLT_Photon100_NoHE_Cleaned_L1R_v1_event", & HLT_Photon100_NoHE_Cleaned_L1R_v1_event,"HLT_Photon100_NoHE_Cleaned_L1R_v1_event/O");
    myEvent->Branch("HLT_Photon10_Cleaned_L1R_event", & HLT_Photon10_Cleaned_L1R_event,"HLT_Photon10_Cleaned_L1R_event/O");
    myEvent->Branch("HLT_Photon15_Cleaned_L1R_event", & HLT_Photon15_Cleaned_L1R_event,"HLT_Photon15_Cleaned_L1R_event/O");
    myEvent->Branch("HLT_Photon17_SC17HE_L1R_v1_event", & HLT_Photon17_SC17HE_L1R_v1_event,"HLT_Photon17_SC17HE_L1R_v1_event/O");
    myEvent->Branch("HLT_Photon20_NoHE_L1R_event", & HLT_Photon20_NoHE_L1R_event,"HLT_Photon20_NoHE_L1R_event/O");
    myEvent->Branch("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event", & HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event,"HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1_event/O");
    myEvent->Branch("HLT_Photon35_Isol_Cleaned_L1R_v1_event", & HLT_Photon35_Isol_Cleaned_L1R_v1_event,"HLT_Photon35_Isol_Cleaned_L1R_v1_event/O");
    myEvent->Branch("HLT_Photon50_Cleaned_L1R_v1_event", & HLT_Photon50_Cleaned_L1R_v1_event,"HLT_Photon50_Cleaned_L1R_v1_event/O");
    myEvent->Branch("HLT_Photon50_NoHE_L1R_event", & HLT_Photon50_NoHE_L1R_event,"HLT_Photon50_NoHE_L1R_event/O");
    myEvent->Branch("HLT_Photon70_NoHE_Cleaned_L1R_v1_event", & HLT_Photon70_NoHE_Cleaned_L1R_v1_event,"HLT_Photon70_NoHE_Cleaned_L1R_v1_event/O");
  }
  
  if(runvertex_){
    myEvent->Branch("Vertex_n",&Vertex_n,"Vertex_n/I");
    myEvent->Branch("Vertex_x",vx,"vx[Vertex_n]/F");
    myEvent->Branch("Vertex_y",vy,"vy[Vertex_n]/F");
    myEvent->Branch("Vertex_z",vz,"vz[Vertex_n]/F");
    myEvent->Branch("Vertex_tracksize",vtracksize,"vtracksize[Vertex_n]/F");
    myEvent->Branch("Vertex_ndof",vndof,"vndof[Vertex_n]/F");
    myEvent->Branch("Vertex_chi2",chi2,"chi2[Vertex_n]/F");
    myEvent->Branch("Vertex_d0",v_d0,"v_d0[Vertex_n]/F");
    myEvent->Branch("Vertex_isFake",v_isFake,"v_isFake[Vertex_n]/O");
  }
  
  if(runscraping_){
    myEvent->Branch("Scraping_isScrapingEvent",&Scraping_isScrapingEvent,"Scraping_isScrapingEvent/O");
    myEvent->Branch("Scraping_numOfTracks",&Scraping_numOfTracks,"Scraping_numOfTracks/I");
    myEvent->Branch("Scraping_fractionOfGoodTracks",&Scraping_fractionOfGoodTracks,"Scraping_fractionOfGoodTracks/F");
  }
	
  if (runtracks_){
    myEvent->Branch("Track_n",&Track_n,"Track_n/I");
    myEvent->Branch("Track_px",trk_px,"trk_px[Track_n]/F");
    myEvent->Branch("Track_py",trk_py,"trk_py[Track_n]/F");
    myEvent->Branch("Track_pz",trk_pz,"trk_pz[Track_n]/F");
    myEvent->Branch("Track_pt",trk_pt,"trk_pt[Track_n]/F");
    myEvent->Branch("Track_eta",trk_eta,"trk_eta[Track_n]/F");
    myEvent->Branch("Track_phi",trk_phi,"trk_phi[Track_n]/F");
  }
  if (runjets_){
    myEvent->Branch("Jet_n",&Jet_n,"Jet_n/I");
    myEvent->Branch("Jet_px",jet_px,"jet_px[Jet_n]/F");
    myEvent->Branch("Jet_py",jet_py,"jet_py[Jet_n]/F");
    myEvent->Branch("Jet_pz",jet_pz,"jet_pz[Jet_n]/F");
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
  }
  
  if (runelectrons_){
    myEvent->Branch("Electron_n",&Electron_n,"Electron_n/I");
    myEvent->Branch("Electron_px",electron_px,"electron_px[Electron_n]/F");
    myEvent->Branch("Electron_py",electron_py,"electron_py[Electron_n]/F");
    myEvent->Branch("Electron_pz",electron_pz,"electron_pz[Electron_n]/F");
    myEvent->Branch("Electron_pt",electron_pt,"electron_pt[Electron_n]/F");
    myEvent->Branch("Electron_eta",electron_eta,"electron_eta[Electron_n]/F");
    myEvent->Branch("Electron_phi",electron_phi,"electron_phi[Electron_n]/F");
    myEvent->Branch("Electron_energy",electron_energy,"electron_energy[Electron_n]/F");
    myEvent->Branch("Electron_charge",electron_charge,"electron_charge[Electron_n]/F");
    myEvent->Branch("Electron_trkIso",electron_trkIso,"electron_trkIso[Electron_n]/F");   
  }
  
  if (runmuons_){
    myEvent->Branch("Muon_n",&Muon_n,"Muon_n/I");
    myEvent->Branch("Muon_px",muon_px,"muon_px[Muon_n]/F");
    myEvent->Branch("Muon_py",muon_py,"muon_py[Muon_n]/F");
    myEvent->Branch("Muon_pz",muon_pz,"muon_pz[Muon_n]/F");
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
    
  }
  
  if (runtaus_){
    myEvent->Branch("Tau_n",&Tau_n,"Tau_n/I");
    myEvent->Branch("Tau_px",tau_px,"tau_px[Tau_n]/F");
    myEvent->Branch("Tau_py",tau_py,"tau_py[Tau_n]/F");
    myEvent->Branch("Tau_pz",tau_pz,"tau_pz[Tau_n]/F");
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
    myEvent->Branch("is_Ztu_event",&is_Ztau_event,"is_Ztau_event/O");   
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
    
    myEvent->Branch("Photon_HoE",pho_HoE,"pho_HoE[Photon_n]/F");
    myEvent->Branch("Photon_px",pho_px,"pho_px[Photon_n]/F");
    myEvent->Branch("Photon_py",pho_py,"pho_py[Photon_n]/F");
    myEvent->Branch("Photon_pz",pho_pz,"pho_pz[Photon_n]/F");
    myEvent->Branch("Photon_no_of_basic_clusters",pho_size,"pho_size[Photon_n]/I");
    
    myEvent->Branch("Photon_sc_energy",pho_sc_energy,"pho_sc_energy[Photon_n]/F");
    myEvent->Branch("Photon_sc_eta",pho_sc_eta,"pho_sc_eta[Photon_n]/F");
    myEvent->Branch("Photon_sc_phi",pho_sc_phi,"pho_sc_phi[Photon_n]/F");
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
    
    myEvent->Branch("Photon_ntracks",pho_nTracks,"pho_nTracks[Photon_n]/I");
    myEvent->Branch("Photon_isconverted",pho_isConverted,"pho_isConverted[Photon_n]/O");
    myEvent->Branch("Photon_pairInvmass",pho_pairInvariantMass,"pho_pairInvariantMass[Photon_n]/F");
    myEvent->Branch("Photon_pairCotThetaSeperation",pho_pairCotThetaSeparation,"pho_pairCotThetaSeparation[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumX",pho_pairMomentum_x,"pho_pairMomentum_x[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumY",pho_pairMomentum_y,"pho_pairMomentum_y[Photon_n]/F");
    myEvent->Branch("Photon_pairmomentumZ",pho_pairMomentum_z,"pho_pairMomentum_z[Photon_n]/F");
    myEvent->Branch("Photon_EoverP",pho_EoverP,"pho_EoverP[Photon_n]/F");
    myEvent->Branch("Photon_vertexX",pho_vertex_x,"pho_vertex_x[Photon_n]/F");
    myEvent->Branch("Photon_vertexY",pho_vertex_y,"pho_vertex_y[Photon_n]/F");
    myEvent->Branch("Photon_vertexZ",pho_vertex_z,"pho_vertex_z[Photon_n]/F");
    myEvent->Branch("Photon_ZOfPrimaryVertex",pho_zOfPrimaryVertex,"pho_zOfPrimaryVertex[Photon_n]/F");
    myEvent->Branch("Photon_distOfMinimumApproach",pho_distOfMinimumApproach,"pho_distOfMinimumApproach[Photon_n]/F");
    myEvent->Branch("Photon_dPhiTracksAtVtx",pho_dPhiTracksAtVtx,"pho_dPhiTracksAtVtx[Photon_n]/F");
    myEvent->Branch("Photon_dPhiTracksAtEcal",pho_dPhiTracksAtEcal,"pho_dPhiTracksAtEcal[Photon_n]/F");
    myEvent->Branch("Photon_dEtaTracksAtVtx",pho_dEtaTracksAtEcal,"pho_dEtaTracksAtEcal[Photon_n]/F");
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
	
  if(runmet_){
    //Calomet variables
    myEvent->Branch("CaloMetSigma",&CaloMetSig,"CaloMetSig/F");
    //myEvent->Branch("CaloMetCorr",&CaloMetCorr,"CaloMetCorr/F");
    myEvent->Branch("CaloMetEt",&CaloMetEt,"CaloMetEt/F");
    myEvent->Branch("CaloMetEx",&CaloMetEx,"CaloMetEx/F");
    myEvent->Branch("CaloMetEy",&CaloMetEy,"CaloMetEy/F");
    myEvent->Branch("CaloMetEz",&CaloMetEz,"CaloMetEz/F");
    myEvent->Branch("CaloMetPhi",&CaloMetPhi,"CaloMetPhi/F");
    myEvent->Branch("CaloMetSumEt",&CaloSumEt,"CaloMetSumEt/F");
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
    myEvent->Branch("PFMetPt",&PFMetPt,"PFMetPt/F");
    myEvent->Branch("PFMetPx",&PFMetPx,"PFMetPx/F");
    myEvent->Branch("PFMetPy",&PFMetPy,"PFMetPy/F");
    myEvent->Branch("PFMetPhi",&PFMetPhi,"PFMetPhi/F");
    myEvent->Branch("PFMetSumEt",&PFMetSumEt,"PFMetSumEt/F");
  }//end of if(runmet)
  if(runPFmet_&& runphotons_)
    myEvent->Branch("Delta_phiPF",&Delta_phiPF,"Delta_phiPF/F");
  
  
  if(runTCmet_){
    myEvent->Branch("TCMetPt",&TCMetPt,"TCMetPt/F");
    myEvent->Branch("TCMetPx",&TCMetPx,"TCMetPx/F");
    myEvent->Branch("TCMetPy",&TCMetPy,"TCMetPy/F");
    myEvent->Branch("TCMetPhi",&TCMetPhi,"TCMetPhi/F");
    myEvent->Branch("TCMetSumEt",&TCMetSumEt,"TCMetSumEt/F");
  }//end of if(runmet)
  if(runTCmet_&& runphotons_)
    myEvent->Branch("Delta_phiTC",&Delta_phiTC,"Delta_phiTC/F");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {
  f->WriteTObject(myEvent);
  delete myEvent;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
