// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc Analysis/Analyzer/src/Analyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Sandhya Jain
//         Created:  Fri Apr 17 11:00:06 CEST 2009
// $Id: Analyzer.cc,v 1.3 2010/02/10 10:33:27 sandhya Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
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
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"//(phi1,phi2) does phi1-phi2
#include "DataFormats/Math/interface/deltaR.h"//(eta1,phi1,eta2,phi2)

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

#include <map>
#include <string>

#include "ADDmonophoton/Analyzer/interface/Analyzer.h"


using namespace std;
using namespace ROOT::Math::VectorUtil ;

//utility function prototypes
double correct_phi(double phi);
double Theta(double eta);
double Pl(double P,double Pt);

double Analyzer::rookFractionBarrelCalculator( const reco::SuperCluster &superCluster ,const EcalRecHitCollection &recHits){
  double rookFraction = 0.; // between 0 and 1

  // get recHit crystal IDs for this superCluster
  std::vector< std::pair<DetId, float> > myHitsPair = superCluster.hitsAndFractions();
  //make sure hits are in barrel ecal!
  bool isHitEcalBarrel = false;
  if ((myHitsPair[0].first).det() == DetId::Ecal && (myHitsPair[0].first).subdetId() == EcalBarrel ){ 
    isHitEcalBarrel = true;
  }
  if (isHitEcalBarrel == false){
    cout << "this superCluster is not in Barrel Ecal! rookFractionBarrelCalculator is returning nonsense value 5.0"<<endl;
    rookFraction = 5.0;
    return rookFraction;
  }
  std::vector<DetId> usedCrystals;
  for(unsigned int i=0; i< myHitsPair.size(); i++){
    usedCrystals.push_back(myHitsPair[i].first);
  }
  // get seed energy and position
  float seedEnergy = -20.;
  int seedIPhi = 500;
  int seedIEta = 500;
  for(unsigned int i=0; i<usedCrystals.size(); i++){
    //get pointer to recHit object
    EcalRecHitCollection::const_iterator myRH = recHits.find(usedCrystals[i]);
    EBDetId EBdetIdi( myRH->detid() );
    if(myRH->energy() > seedEnergy){ 
      seedEnergy = myRH->energy();
      seedIPhi   = EBdetIdi.iphi();
      seedIEta   = EBdetIdi.ieta();
    }// if energy is larger than seed E
  }// loop over clustered crystals
  if (seedIEta<0) seedIEta++; // account for no ieta = 0
  // select adjacent crystal with the most E (not diagonal!! hence the name ROOK)
  float adjacentEnergy = -20.;
  for(EcalRecHitCollection::const_iterator rh = recHits.begin(); rh != recHits.end(); rh++){
    EBDetId EBdetIdi( rh->detid() );
    int stampIPhi = EBdetIdi.iphi();
    int stampIEta = EBdetIdi.ieta();
    if (stampIEta<0) stampIEta++; // account for no ieta = 0
    int deltaIEta = abs(stampIEta - seedIEta);
    int deltaIPhi = abs(stampIPhi - seedIPhi);
    if (deltaIPhi > 180) deltaIPhi = 360 - deltaIPhi; // account for phi wrap around
    if( (deltaIEta==1 && deltaIPhi==0) || (deltaIEta==0 && deltaIPhi==1) ){
      if( rh->energy() > adjacentEnergy ){
	adjacentEnergy = rh->energy();
      }//if energy is greatest adjacent
    }
  }// loop over Ecal rec Hit collection
  rookFraction = adjacentEnergy/seedEnergy;
  return rookFraction;
}



//
// class decleration
//
class genPho{
public:
  genPho(){};  ~genPho(){};
  double pt,px,py,pz,eta,phi,E,motherPt,motherEta,motherPhi , GrandmotherPt,GrandmotherEta,GrandmotherPhi;
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
  runjets_(iConfig.getUntrackedParameter<bool>("runjets")),
  runtaus_(iConfig.getUntrackedParameter<bool>("runtaus")),
  runHLT_(iConfig.getUntrackedParameter<bool>("runHLT")),
  runL1_(iConfig.getUntrackedParameter<bool>("runL1")),
  runtracks_(iConfig.getUntrackedParameter<bool>("runtracks")),
  runrechit_(iConfig.getUntrackedParameter<bool>("runrechit")),
  runvertex_(iConfig.getUntrackedParameter<bool>("runvertex")),
  init_(false)
{
   //now do what ever initialization is needed
  nevents =0;
  n_signal_events = 0; n_Z_events   = 0; n_W_events    = 0 ;
  n_Zelec_events =  0; n_Zmu_events = 0; n_Ztau_events = 0; n_Znunu_events = 0;
  n_Welec_events =  0; n_Wmu_events = 0; n_Wtau_events = 0;
  n_diphoton_events=0;  n_SingleHardPhoton_events = 0;
  ngenphotons  = 0; //for every event, how many stable photons 
  nhardphotons = 0; //for every event, how many hard photons 
  Vertex_n         = 0;
  Track_n          = 0;
  Photon_n         = 0;
  Jet_n            = 0;
  Electron_n       = 0; 
  Muon_n           = 0;
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
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   
   RunNumber   = iEvent.id().run();
   EventNumber = iEvent.id().event();
   cout<<"RunNumber:"<<RunNumber<<"     Event:"<< EventNumber<<endl;
   nevents++;
   //cout<<"Event:"<<nevents<<endl;
   //getting handle to generator level information
   if( rungenParticleCandidates_ )
     {
       ngenphotons  = 0;
       nhardphotons = 0;
       is_signal_event = 0; is_Z_event   = 0; is_W_event    = 0; 
       is_Zelec_event  = 0; is_Zmu_event = 0; is_Ztau_event = 0; is_Znunu_event =0  ;
       is_Welec_event  = 0; is_Wmu_event = 0; is_Wtau_event = 0;
       is_SingleHardPhoton_event=0;  is_diphoton_event=0;
      
       Handle<GenParticleCollection> genParticles;
       iEvent.getByLabel("genParticles", genParticles); 
       std::vector<genPho>            mygenphoton_container;
       mygenphoton_container.clear();
       genPho genphoton;
       int ii =0;
       for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++) 
	 {
	   //getting information from hard scattered Graviton 
	   if (genparticle->pdgId()==39 && genparticle->status()==3) 
	     { 
	       is_signal_event = 1;
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
	 if (genparticle->pdgId()==23 && genparticle->status()==3) 
	   { 
	     is_Z_event = 1;
	     n_Z_events++;
	     //cout<"getting information from Z now"<<endl;
	     gen_Zboson_pt  = genparticle->pt();
	     gen_Zboson_px  = genparticle->px();
	     gen_Zboson_py  = genparticle->py();
	     gen_Zboson_pz  = genparticle->pz();
	     gen_Zboson_phi = correct_phi(genparticle->phi());
	     gen_Zboson_eta = genparticle->eta();
	     gen_Zboson_E   = genparticle->energy();
	     int daughters          = genparticle->numberOfDaughters();
	     int iDaughter=0;
	     for(int i = 0;i<daughters;i++)
	       {
		 const reco::Candidate *daughter   = genparticle->daughter(i);
		 //cout<<"genparticle->daughter(i)"<<genparticle->daughter(i)<<std::endl;
		 if(abs(daughter->pdgId())==12||abs(daughter->pdgId())==14||abs(daughter->pdgId())==16){is_Znunu_event=1; n_Znunu_events++;}
		 if(daughter->pdgId()==11) { is_Zelec_event=1; n_Zelec_events++;}
		 if(daughter->pdgId()==13) { is_Zmu_event=1  ; n_Zmu_events++  ;}
		 if(daughter->pdgId()==15) { is_Ztau_event=1 ; n_Ztau_events++  ;}
		 if(daughter->pdgId()!=23) 
		   {
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
	   is_W_event = 1;
	   cout<<" W motherID:" << genparticle->mother()->pdgId()<<endl; 
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
	   //cout<<"W pt:"<<genparticle->pt()<<endl;
	   //cout<<"W daughters:"<<endl;
	   for(int i = 0;i<daughters;i++)
	     {
	       const reco::Candidate *daughter   = genparticle->daughter(i);
	       if(abs(daughter->pdgId())==11) {is_Welec_event=1; n_Welec_events++;}
	       if(abs(daughter->pdgId())==13) {is_Wmu_event=1  ; n_Wmu_events++  ;}
	       if(abs(daughter->pdgId())==15) {is_Wtau_event=1 ; n_Wtau_events++ ;}  
	       cout<<"ID, Status,Pt:"<<abs(daughter->pdgId())<<"   "<<daughter->status()<<"   "<<daughter->pt()<<endl;
	       //getting leptons decaying from W
	       if(abs(daughter->pdgId())!=24) 
		 {
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
	 if( abs(genparticle->pdgId())==13 )
	   {
	     //cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	     int daughters   = genparticle->numberOfDaughters();
	     int iDaughter=0;
	     for(int i = 0;i<daughters;i++)
	       {
		 const reco::Candidate *daughter   = genparticle->daughter(i);
		 //cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
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
	 if( abs(genparticle->pdgId())==15 )
	   {
	     //cout<<"parent ID, Status, Pt:"<<abs(genparticle->pdgId())<<"   "<<genparticle->status()<<"  "<< genparticle->pt()<<endl;
	     int daughters   = genparticle->numberOfDaughters();
	     int iDaughter=0;
	     for(int i = 0;i<daughters;i++)
	       {
		 const reco::Candidate *daughter   = genparticle->daughter(i);
		 //cout<<"daughterID, status,Pt:"<<abs(daughter->pdgId())<<"   " <<daughter->status()<<"  "<< daughter->pt()<<endl;
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
	 
	 //cout<<"getting gen photon information now"<<endl;
	 //getting information from all photons
	 //if (genparticle->pdgId()==22 && genparticle->status()==1 && genparticle->pt()>5.)
	 if (genparticle->pdgId()==22)
	   { 
	     //doing it this way, as I want to sort it in pt after filling everything in the container
	     const reco::Candidate *mom = genparticle->mother();
	     genphoton.motherID = mom->pdgId();
	     genphoton.motherStatus = mom->status();
	     genphoton.motherPt = mom->pt();
	     genphoton.motherEta = mom->eta();
	     genphoton.motherPhi = correct_phi(mom->phi());
             genphoton.GrandmotherID = mom->mother()->pdgId();
             genphoton.GrandmotherStatus = mom->mother()->status();
             genphoton.GrandmotherPt = mom->mother()->pt();
             genphoton.GrandmotherEta = mom->mother()->eta();
             genphoton.GrandmotherPhi = correct_phi(mom->mother()->phi());
	     genphoton.pt  = genparticle->pt();
	     genphoton.px  = genparticle->px();
	     genphoton.py  = genparticle->py();
	     genphoton.pz  = genparticle->pz();
	     genphoton.phi = correct_phi(genparticle->phi());
	     genphoton.eta = genparticle->eta();
	     genphoton.E   = genparticle->energy();
	     genphoton.status = genparticle->status();
	     mygenphoton_container.push_back(genphoton); ngenphotons++;
	   }//end of if (genparticle->pdgId()==22 && genparticle->status()==1)
 
	 if (genparticle->pdgId()==22 && genparticle->status()==3) 
	   {
	     gen_Hpho_pt[nhardphotons]  = genparticle->pt();
	     gen_Hpho_px[nhardphotons]  = genparticle->px();
	     gen_Hpho_py[nhardphotons]  = genparticle->py();
	     gen_Hpho_pz[nhardphotons]  = genparticle->pz();
	     gen_Hpho_phi[nhardphotons] = correct_phi(genparticle->phi());
	     gen_Hpho_eta[nhardphotons] = genparticle->eta();
	     gen_Hpho_E[nhardphotons]   = genparticle->energy();
	     nhardphotons++;
	   }//end of if (genparticle->pdgId()==22 && genparticle->status()==3)

	 ii++;
       }//end of for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++)
     
       if (nhardphotons==1) 
	 {
	   n_SingleHardPhoton_events++;
	   is_SingleHardPhoton_event=1;
	 }
       if (nhardphotons==2)
	 {
	   n_diphoton_events++; 
	   is_diphoton_event=1;
	 }

       if(mygenphoton_container.size()!=0)
	 {
	   std::sort(mygenphoton_container.begin(),mygenphoton_container.end(),PtSortCriterium());
	   for(unsigned int x=0;x < mygenphoton_container.size(); x++)
             {
	       //std::cout<<"genphoton motherID:"<<mygenphoton_container[x].motherID<<endl;
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
	       if( gen_pho_pt[x] > 100. && (abs(gen_pho_motherID[x])<=6||abs(gen_pho_motherID[x])==11||abs(gen_pho_motherID[x])==9|| abs(gen_pho_motherID[x])==21 )) 
		 cout<<"[x]:"<<x<< "pho_pt:" << gen_pho_pt[x]<<" pho_status:"<< gen_pho_status[x]<<" motherID:" << gen_pho_motherID[x]<<" motherStatus: "<< gen_pho_motherStatus[x]<< " GrandmotherID: "<< gen_pho_GrandmotherID[x]<< " GrandmotherStatus:"<< gen_pho_GrandmotherStatus[x]<< endl;
	       // cout<<"got the photon info right"<<endl;
	     }//end of for loop
	 }//end of if((mygenphoton_container.size()!=0)
       //std::cout<<"mygenphoton_container loop ended"<<std::endl; 
     }//end of if(rungenParticleCandidates_)
   
   ///// L1
   if(runL1_)
     {
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

   /*
   /////HLT
   if(runHLT_==1){
     Handle<TriggerResults> HLTR;
     iEvent.getByLabel(hlTriggerResults_,HLTR);
     if (!init_) {
       init_=true;
       triggerNames_.init(*HLTR);
       hlNames_=triggerNames_.triggerNames();
     }

     // decision for each HL algorithm
     const unsigned int n(hlNames_.size());
     
     for(unsigned int i = 0; i<n ;i++)
       {
	 //cout<<hlNames_[i]<<" :"<<HLTR->accept(i)<<endl;
	 HLT_chosen[ hlNames_[i]]= HLTR->accept(i);
	 if(hlNames_[i]=="HLT_MET50")
	   HLT_MET50_event       = HLTR->accept(i);
	 if(hlNames_[i]=="HLT_MET75")
	   HLT_MET75_event       = HLTR->accept(i);
	 if(hlNames_[i]=="HLT_Photon15_L1R")
	   HLT_Photon15_event    = HLTR->accept(i);
	 if(hlNames_[i]=="HLT_Photon25_L1R")
	   HLT_Photon25_event    = HLTR->accept(i);
	 if(hlNames_[i]=="HLT_DoubleEle10_SW_L1R")
	   HLT_DoubleEle10_event = HLTR->accept(i);
	 if(hlNames_[i]=="HLT_DoubleMu3")
	   HLT_DoubleMu3_event   = HLTR->accept(i);

       }
   }
   */
   if(runvertex_){
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(Vertices_, recVtxs);
   vector<Vertex> my_vertices;
   Vertex_n = recVtxs->size();
   my_vertices.clear();
   for(reco::VertexCollection::const_iterator v=recVtxs->begin();v!=recVtxs->end(); ++v){
     my_vertices.push_back(*v);
   }
   for (unsigned int y = 0; y <  my_vertices.size();y++)
     {
       vx[y]     = my_vertices[y].x();
       vy[y]     = my_vertices[y].y();
       vz[y]     = my_vertices[y].z();
       chi2[y]   = my_vertices[y].chi2();
       vtracksize[y] = my_vertices[y].tracksSize();
       vndof[y] = my_vertices[y].ndof();
       v_isFake[y] = my_vertices[y].isFake();
       v_d0[y] = my_vertices[y].position().rho();
     }
   }

   if(runtracks_){
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(Tracks_,tracks);
   //cout<<"Track_n - all "<<tracks->size()<<endl;
   std::vector<reco::Track>  myTrack_container;
   myTrack_container.clear();
   /*
   for(reco::TrackCollection::const_iterator Track_iter = tracks->begin();
       Track_iter != tracks->end();++Track_iter) {
      if(Track_iter->pt()>5.){
        myTrack_container.push_back(*Track_iter);
        //cout<<"nTracks abv 5 GeV: "<<myTrack_container.size()<<endl;
      }
   }
   */
   //cout<<"nTracks abv 5 GeV after the loop:"<<myTrack_container.size()<<endl;
   //change it later to myTrack_container.size() when PAS is done 
   Track_n = tracks->size();
   /*
   if(myTrack_container.size()>1)
     std::sort(myTrack_container.begin(),myTrack_container.end(),PtSortCriterium3());
   for(unsigned int x=0;x < myTrack_container.size();x++)
     {
       trk_pt[x]  = myTrack_container[x].pt();
       trk_px[x]  = myTrack_container[x].px();
       trk_py[x]  = myTrack_container[x].py();
       trk_pz[x]  = myTrack_container[x].pz();
       trk_phi[x] = correct_phi(myTrack_container[x].phi());
       trk_eta[x] = myTrack_container[x].eta();
       //cout<<"pt: "<< trk_pt[x] << endl;
     }//end of for loop
*/
     }//end of if(runtracks_)
   
      
   std::vector<pat::Photon> myphoton_container;
   myphoton_container.clear();
   if(runphotons_)
     {
       
       Handle<reco::PhotonCollection> phoHandle;
       iEvent.getByLabel("photons", phoHandle);
       reco::PhotonCollection::const_iterator photon;
       /*
       edm::Handle<edm::View<pat::Photon> > allLayer1Photons;
       iEvent.getByLabel("allLayer1Photons", allLayer1Photons);
       edm::Handle<edm::View<pat::Photon> > selectedlayer1Photons;
       iEvent.getByLabel("selectedLayer1Photons", selectedlayer1Photons);
       
       edm::Handle<edm::View<pat::Photon> > phoHandle;
       iEvent.getByLabel(phoLabel_,phoHandle);
       //const edm::View<pat::Photon> & photons = *phoHandle;
       edm::View<pat::Photon>::const_iterator photon;
       */
       Photon_n = phoHandle->size();
       //cout<<"total photons reconstructed at RECO level:"<<recphotons->size()<<endl;
       //cout<<"Step 1 : number of all layer1 Photons:"<<allLayer1Photons->size()<<endl;
       //cout<<"Step 2 : number of selected layer1 Photons:"<<selectedlayer1Photons->size()<<endl;
       //cout<<"Step 3 : number of cleaned layer1 Photons:"<<phoHandle->size()<<endl;
       //if(allLayer1Photons->size()!=recphotons->size())cout<<"strange! allLayer1Photons!=recphotons "<<endl;
       //if(allLayer1Photons->size()!=selectedlayer1Photons->size())cout<<"strange! allLayer1Photons!= selectedlayer1Photons "<<endl;
       //if(selectedlayer1Photons->size()< phoHandle->size())cout<<"strange! selectedlayer1Photons->size() < cleanLayer1Photons"<<endl;

       //cout<<"photon container size:"<<Photon_n<<endl;
       
       for(photon = phoHandle->begin();photon!=phoHandle->end();++photon){
	 myphoton_container.push_back(*photon) ;
       }
       if(myphoton_container.size()!=0)
	 {
	   for(unsigned int x=0; x < myphoton_container.size();x++)
	     {
	       //cout<<"photon pt:"<<myphoton_container[x].pt()<<"  photon eta:"<<myphoton_container[x].eta()<<"  photon phi:"<<correct_phi(myphoton_container[x].phi())<<endl; 
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

	       if(myphoton_container[x].genParticleRef().isNonnull())
		 {
		   matchpho_E[x]                =  myphoton_container[x].genPhoton()->energy();
		   matchpho_pt[x]               =  myphoton_container[x].genPhoton()->pt();
		   matchpho_eta[x]              =  myphoton_container[x].genPhoton()->eta();
		   matchpho_phi[x]              =  correct_phi(myphoton_container[x].genPhoton()->phi());
		   matchpho_px[x]               =  myphoton_container[x].genPhoton()->px();
		   matchpho_py[x]               =  myphoton_container[x].genPhoton()->py();
		   matchpho_pz[x]               =  myphoton_container[x].genPhoton()->pz();
		 }
	       else
		 {
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
	       //cout<<"size of conversion vector:"<<conversions.size()<<endl;
	       for (unsigned int iConv=0; iConv<conversions.size(); iConv++) {
		 reco::ConversionRef aConv=conversions[iConv];
		 //cout<<"ntracks:"<<aConv->nTracks()<<endl;
		 //cout<<"isConverted:"<<aConv->isConverted()<<endl;
		 if ( aConv->nTracks() <2 ) continue; 
		 if ( aConv->conversionVertex().isValid() )
		   {
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
		 else
		   {
		     pho_nTracks[x]                    = 9999;
		     pho_isConverted[x]                = -99;
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
	     }//end of for loop over x
	 }//if(myphoton_container.size!=0) 
       //to get the photon hit information from every crystal of SC
       if(runrechit_)
	 { 
	   Handle<EcalRecHitCollection> Brechit;//barrel
	   Handle<EcalRecHitCollection> Erechit;//endcap
	   iEvent.getByLabel(rechitBLabel_,Brechit);
	   iEvent.getByLabel(rechitELabel_,Erechit);
	   const EcalRecHitCollection* barrelRecHits= Brechit.product();
	   const EcalRecHitCollection* endcapRecHits= Erechit.product();
	   edm::ESHandle<CaloTopology> pTopology;
	   iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
	   const CaloTopology *topology = theCaloTopo_.product();
	   if(myphoton_container.size()!=0)
	     {
	       for(unsigned int x=0; x < myphoton_container.size();x++)
		 {
		   
		   std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = myphoton_container[x].superCluster()->hitsAndFractions();
		   std::vector<CrystalInfo> crystalinfo_container;
		   crystalinfo_container.clear();
		   CrystalInfo crystal;
		   double timing_avg =0.0;
		   int ncrys   = 0;
		   ncrysPhoton[x]= 0;
		   vector< std::pair<DetId, float> >::const_iterator detitr;
		   for(detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr)
		     {
		       if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) 
			 {
			   EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
			   EcalRecHitCollection::const_iterator thishit;
			   if ( j!= Brechit->end())  thishit = j;
			   if ( j== Brechit->end())
			     {
			       std::cout<<"thishit not matched "<<std::endl;
			       continue;
			     }
			   //std::cout<<"thishit matched "<<std::endl;
			   EBDetId detId  = (EBDetId)((*detitr).first);
			   crystal.rawId  = thishit->id().rawId();
			   crystal.energy = thishit->energy();
			   crystal.time   = thishit->time();
			   crystal.ieta   = detId.ieta();
			   crystal.iphi   = detId.iphi();
			   //std::cout<<"thishit energy,time: "<<crystal.energy<<"   "<< crystal.time<<std::endl;
			   //calculate timing avg
			   if(crystal.energy > 0.1)
			     {
			       timing_avg  = timing_avg + crystal.time;
			       ncrys++;
			     } 

			 }//end of if ((*detitr).det() == DetId::Ecal && (*detitr).subdetId() == EcalBarrel)
		       else 
			   {
			     EcalRecHitCollection::const_iterator j= Erechit->find(((*detitr).first));
			     EcalRecHitCollection::const_iterator thishit;
			     if ( j!= Erechit->end())  thishit = j;
			     if ( j== Erechit->end())
			       {
				 std::cout<<"thishit not matched "<<std::endl;
				 continue;
			       }
			     //std::cout<<"thishit matched "<<std::endl;
			     EEDetId detId  = (EEDetId)((*detitr).first);
			     crystal.energy = thishit->energy();
			     crystal.time   = thishit->time();
			     crystal.rawId  = 999;
			     crystal.ieta   = -99;
                             crystal.iphi   = -99;
			     //std::cout<<"thishit energy,time: "<<crystal.energy<<"   "<< crystal.time<<std::endl;
			   if(crystal.energy > 0.1)
			     {
			       timing_avg  = timing_avg + crystal.time;
			       ncrys++;
			     } 
			   }
		       crystalinfo_container.push_back(crystal);  
		     }
		   std::sort(crystalinfo_container.begin(),crystalinfo_container.end(),EnergySortCriterium());
		   if (ncrys !=0) timing_avg = timing_avg/(double)ncrys;
		   else timing_avg = -99.;
		   //cout<<" total hits for this photon:"<<crystalinfo_container.size()<<endl;
		   //cout<<" total hits contibuting for timing(energy > 0.1):"<<ncrys<<endl;
		   ncrysPhoton[x] = crystalinfo_container.size(); 
		   cout<<"ncrysPhoton:"<<ncrysPhoton[x]<<endl;
		   pho_timingavg_xtal[x]      = timing_avg;
		   for (unsigned int y =0; y < 100.;y++)
		     {
		       pho_timing_xtal[x][y]         = -99.;
		       pho_energy_xtal[x][y]         = -99.;
		       pho_ieta_xtalEB[x][y]           = -99;
		       pho_iphi_xtalEB[x][y]           = -99;
		     }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++)
		   for (unsigned int y =0; y < crystalinfo_container.size();y++)
		     {
		       pho_timing_xtal[x][y]         = crystalinfo_container[y].time;
		       pho_energy_xtal[x][y]         = crystalinfo_container[y].energy;
		       pho_ieta_xtalEB[x][y]           = crystalinfo_container[y].ieta;
		       pho_iphi_xtalEB[x][y]           = crystalinfo_container[y].iphi;
		     }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++
		   if(myphoton_container[x].isEB())
		     {
		       std::vector<float> showershapes_barrel = EcalClusterTools::roundnessBarrelSuperClusters(*(myphoton_container[x].superCluster()),*barrelRecHits,0);
		       //cout<<"roundness for barrel photon:"<<showershapes_barrel[0]<<endl;
		       //cout<<"angle for barrel photon:"<<showershapes_barrel[1]<<endl;
		       pho_roundness[x]    = (double)showershapes_barrel[0];
		       pho_angle[x]        = (double)showershapes_barrel[1];
		       pho_s9[x]           = pho_energy_xtal[x][0]/pho_e3x3[x];
		       pho_rookFraction[x] = rookFractionBarrelCalculator(*( myphoton_container[x].superCluster() ), *barrelRecHits);
		       pho_swissCross[x]   = EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology))+ EcalClusterTools::eBottom( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eLeft( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology)) + EcalClusterTools::eRight( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology));
		       cout<<"etop: "<<EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl; 
		       cout<<"etop: "<<EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), &(*barrelRecHits), &(*topology))<<endl; 
		      if(1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl; 
 		     }//end of if(myphoton_container[x].isEB())
		    else{ 
		      pho_roundness[x]   = -99.;
		      pho_angle[x]       = -99.;
		      pho_s9[x]          = pho_energy_xtal[x][0]/pho_e3x3[x];
		      pho_rookFraction[x]= -99.;
		      pho_swissCross[x]   = EcalClusterTools::eTop( *(myphoton_container[x].superCluster()->seed()), &(*endcapRecHits), &(*topology))+ EcalClusterTools::eBottom( *(myphoton_container[x].superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eLeft( *(myphoton_container[x].superCluster()->seed()), &(*endcapRecHits), &(*topology)) + EcalClusterTools::eRight( *(myphoton_container[x].superCluster()->seed()), &(*endcapRecHits), &(*topology));
		      if(1-pho_swissCross[x]/pho_maxEnergyXtal[x] > 0.95) cout<<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm."<<endl; 
		    }//end of else
		   cout<<"rook fraction:"<<pho_rookFraction[x]<<endl;
		   cout<<"s9:"<<pho_s9[x]<<endl;
		   cout<<"pt,eta:"<<pho_pt[x]<<"  "<< pho_eta[x]<<endl;
		 }//end of for loop over x
	     }//if(myphoton_container.size!=0) 
	 }//if(runrechit_)
       //cout<<"got photon variables"<<endl; 
     }//if(runphotons_)  

   if(runmet_)
     {
       edm::Handle<edm::View<pat::MET> > metHandle;
       iEvent.getByLabel(metLabel_,metHandle);
       //const edm::View<pat::MET> & mets = *metHandle;
       edm::View<pat::MET>::const_iterator met;
       for ( met = metHandle->begin(); met != metHandle->end(); met++)
	 {
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
	   /*
	   genMetPhi    = correct_phi(genMet->phi());
	   genMetSumEt  = genMet->sumEt();
	   genMetPx     = genMet->px();
	   genMetPy     = genMet->py();
	   if(runphotons_==1)
	     if (myphoton_container.size()!=0)
	     Delta_phiGEN                        = fabs(reco::deltaPhi(correct_phi(genMet->phi()),correct_phi(myphoton_container[0].phi())));*/
	   }
	 }
     }
   if(runPFmet_)
     {
       edm::Handle<edm::View<pat::MET> > metPFHandle;
       iEvent.getByLabel(PFmetLabel_,metPFHandle);
       const edm::View<pat::MET> & metsPF = *metPFHandle;
       if ( metPFHandle.isValid() )
	 {
	   PFMetPt     = metsPF[0].et();
	   PFMetPhi    = correct_phi(metsPF[0].phi());
	   PFMetSumEt  = metsPF[0].sumEt();
	   PFMetPx     = metsPF[0].px();
	   PFMetPy     = metsPF[0].py();
	   if(runphotons_==1)
	     if (myphoton_container.size()!=0)
			 Delta_phiPF  = fabs(reco::deltaPhi(PFMetPhi,correct_phi(myphoton_container[0].phi())));
	 }
       else
	 {
	   LogWarning("METEventSelector") << "No Met results for InputTag " ;
	   return;
	 }
     }
   if(runTCmet_)
     {
       edm::Handle<edm::View<pat::MET> > metTCHandle;
       iEvent.getByLabel(TCmetLabel_,metTCHandle);
       const edm::View<pat::MET> & metsTC = *metTCHandle;
       if ( metTCHandle.isValid() )
	 {
	   TCMetPt     = metsTC[0].et();
	   TCMetPhi    = metsTC[0].phi();
	   TCMetSumEt  = metsTC[0].sumEt();
	   TCMetPx     = metsTC[0].px();
	   TCMetPy     = metsTC[0].py();
	   if(runphotons_==1)
             if (myphoton_container.size()!=0)
				 Delta_phiTC  = fabs(reco::deltaPhi(TCMetPhi,correct_phi(myphoton_container[0].phi()))); 
	 }
       else
	 {
	   LogWarning("METEventSelector") << "No Met results for InputTag " ;
	   return;
	 } 
       
     }//end of if(runTCmet_)
   
   if(runjets_)
     {
       edm::Handle<edm::View<pat::Jet> > jetHandle;
       iEvent.getByLabel(jetLabel_,jetHandle);
       const edm::View<pat::Jet> & jets = *jetHandle;
       Jet_n=jetHandle->size();
       size_t njetscounter=0;
       std::vector<pat::Jet>  myjet_container;
       myjet_container.clear();
       for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
	 if(jet_iter->pt()>50)  njetscounter++;
	 //const pat::Jet& myJet = jet_iter->correctedJet(*corr);
	 myjet_container.push_back(*jet_iter);
       }
       if(myjet_container.size()!=0){
	 for(unsigned int x=0;x < myjet_container.size();x++)
	   {
	     jet_pt[x]  = myjet_container[x].pt();
	     jet_px[x]  = myjet_container[x].px();
	     jet_py[x]  = myjet_container[x].py();
	     jet_pz[x]  = myjet_container[x].pz();
	     jet_phi[x] = correct_phi(myjet_container[x].phi());
	     jet_eta[x] = myjet_container[x].eta();
	     jet_emEnergyFraction[x]= myjet_container[x].emEnergyFraction();
	     jet_energyFractionHadronic[x] = myjet_container[x].energyFractionHadronic();
	   }//end of for loop
       }
     }
   
   if(runmuons_)
     {
       edm::Handle<edm::View<pat::Muon> > muonHandle;
       iEvent.getByLabel(muoLabel_,muonHandle);
       vector <pat::Muon> mymuon_container;
       Muon_n = muonHandle->size();
       const edm::View<pat::Muon> & muons = *muonHandle;   // const ... &, we don't make a copy of it!
       for(edm::View<pat::Muon>::const_iterator muon = muons.begin(); muon!=muons.end(); ++muon){
         mymuon_container.push_back(*muon);
       }
       for(unsigned int x=0;x < mymuon_container.size();x++)
	 {
	   muon_pt[x]  = mymuon_container[x].pt();
	   muon_energy[x]  = mymuon_container[x].energy();
	   muon_px[x]  = mymuon_container[x].px();
	   muon_py[x]  = mymuon_container[x].py();
	   muon_pz[x]  = mymuon_container[x].pz();
	   muon_phi[x] = correct_phi(mymuon_container[x].phi());
	   muon_eta[x] = mymuon_container[x].eta();
	   muon_charge[x] = mymuon_container[x].charge();
	 }//end of for loop
     }

  
   if(runelectrons_)
     {
       edm::Handle<edm::View<pat::Electron> > electronHandle;
       iEvent.getByLabel(eleLabel_,electronHandle);
       vector<pat::Electron> myelectron_container;
       Electron_n = electronHandle->size();
       const edm::View<pat::Electron> & electrons = *electronHandle;   // const ... &, we don't make a copy of it!
       for(edm::View<pat::Electron>::const_iterator electron = electrons.begin(); electron!=electrons.end(); ++electron){
         myelectron_container.push_back(*electron);
       }
       for(unsigned int x=0;x < myelectron_container.size();x++)
	 {
	   electron_pt[x]  = myelectron_container[x].pt();
	   electron_energy[x]  = myelectron_container[x].energy();
	   electron_px[x]  = myelectron_container[x].px();
	   electron_py[x]  = myelectron_container[x].py();
	   electron_pz[x]  = myelectron_container[x].pz();
	   electron_phi[x] = correct_phi(myelectron_container[x].phi());
	   electron_eta[x] = myelectron_container[x].eta();
	   electron_charge[x] = myelectron_container[x].charge();
	   electron_trkIso[x] = myelectron_container[x].trackIso();  
	 }//end of for loop
       // histocontainer_["nelectrons"]->Fill(electrons.size());
     }
   
   if(runtaus_)
     {
       edm::Handle<edm::View<pat::Tau> > tauHandle;
       iEvent.getByLabel(tauLabel_,tauHandle);
       vector <pat::Tau> mytau_container;
       Tau_n = tauHandle->size();
       const edm::View<pat::Tau> & taus = *tauHandle;   // const ... &, we don't make a copy of it!
       for(edm::View<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau){
         mytau_container.push_back(*tau);
       }
       for(unsigned int x=0;x < mytau_container.size();x++)
	 {
	   tau_pt[x]  = mytau_container[x].pt();
	   tau_energy[x]  = mytau_container[x].energy();
	   tau_px[x]  = mytau_container[x].px();
	   tau_py[x]  = mytau_container[x].py();
	   tau_pz[x]  = mytau_container[x].pz();
	   tau_phi[x] = correct_phi(mytau_container[x].phi());
	   tau_eta[x] = mytau_container[x].eta();
	   tau_charge[x] = mytau_container[x].charge();
	 }//end of for loop
     }
   //selecting events except the jet veto
   //   if(1-pho_swissCross[0]/pho_maxEnergyXtal[0]<0.95 && abs(pho_eta[0])<2.5  && pho_ecalRecHitSumEtConeDR04[0]<4.2+0.003*pho_pt[0] && pho_hcalTowerSumEtConeDR04[0]< 2.2+0.001*pho_pt[0]  && pho_trkSumPtHollowConeDR04[0] < 2+0.001*pho_pt[0] && pho_HoE[0]<0.05 && (pho_hasPixelSeed[0])==0 && HLT_Photon25_event==1 &&  Delta_phi >2.7 && trk_pt[0]<20. )
   myEvent->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob()
{
  f=new TFile(outFile_.c_str(),"RECREATE");
  //defining a tree here
  myEvent = new TTree("myEvent","a tree with histograms");
  myEvent->Branch("nevents",&nevents,"nevents/I");
  myEvent->Branch("run",&RunNumber,"RunNumber/I");
  myEvent->Branch("event",&EventNumber,"EventNumber/I");


  if(runHLT_)
    {
      myEvent->Branch("HLT_MET50_event",&HLT_MET50_event,"HLT_MET50_event/I");
      myEvent->Branch("HLT_MET75_event",&HLT_MET75_event,"HLT_MET75_event/I");
      myEvent->Branch("HLT_Photon15_event",&HLT_Photon15_event,"HLT_Photon15_event/I");
      myEvent->Branch("HLT_Photon25_event",&HLT_Photon25_event,"HLT_Photon25_event/I");
      myEvent->Branch("HLT_DoubleEle10_event",&HLT_DoubleEle10_event,"HLT_DoubleEle10_event/I");
      myEvent->Branch("HLT_DoubleMu3_event",&HLT_DoubleMu3_event,"HLT_DoubleMu3_event/I");
    }

  if(runvertex_)
    {
      myEvent->Branch("Vertex_n",&Vertex_n,"Vertex_n/I");
      myEvent->Branch("Vertex_x",vx,"vx[Vertex_n]/D");
      myEvent->Branch("Vertex_y",vy,"vy[Vertex_n]/D");
      myEvent->Branch("Vertex_z",vz,"vz[Vertex_n]/D");
      myEvent->Branch("Vertex_tracksize",vtracksize,"vtracksize[Vertex_n]/D");
      myEvent->Branch("Vertex_ndof",vndof,"vndof[Vertex_n]/D");
      myEvent->Branch("Vertex_chi2",chi2,"chi2[Vertex_n]/D");
      myEvent->Branch("Vertex_d0",v_d0,"v_d0[Vertex_n]/D");
      myEvent->Branch("Vertex_isFake",v_isFake,"v_isFake[Vertex_n]/D");
    }

  if (runtracks_)
    {
      myEvent->Branch("Track_n",&Track_n,"Track_n/I");
      myEvent->Branch("Track_px",trk_px,"trk_px[Track_n]/D");
      myEvent->Branch("Track_py",trk_py,"trk_py[Track_n]/D");
      myEvent->Branch("Track_pz",trk_pz,"trk_pz[Track_n]/D");
      myEvent->Branch("Track_pt",trk_pt,"trk_pt[Track_n]/D");
      myEvent->Branch("Track_eta",trk_eta,"trk_eta[Track_n]/D");
      myEvent->Branch("Track_phi",trk_phi,"trk_phi[Track_n]/D");
    }
  if (runjets_)
    {
      //for(vector<string>::iterator corr = JET_CORR.begin(); corr!=JET_CORR.end();++corr )
      //{
	  myEvent->Branch("Jet_n",&Jet_n,"Jet_n/I");
	  myEvent->Branch("Jet_px",jet_px,"jet_px[Jet_n]/D");
	  myEvent->Branch("Jet_py",jet_py,"jet_py[Jet_n]/D");
	  myEvent->Branch("Jet_pz",jet_pz,"jet_pz[Jet_n]/D");
	  myEvent->Branch("Jet_pt",jet_pt,"jet_pt[Jet_n]/D");
	  myEvent->Branch("Jet_eta",jet_eta,"jet_eta[Jet_n]/D");
	  myEvent->Branch("Jet_phi",jet_phi,"jet_phi[Jet_n]/D");
	  myEvent->Branch("Jet_emEnergyFraction",jet_emEnergyFraction,"jet_emEnergyFraction[Jet_n]/D");
	  myEvent->Branch("Jet_energyFractionHadronic",jet_energyFractionHadronic,"jet_energyFractionHadronic[Jet_n]/D");
	  //}
    }

  if (runelectrons_)
    {
      myEvent->Branch("Electron_n",&Electron_n,"Electron_n/I");
      myEvent->Branch("Electron_px",electron_px,"electron_px[Electron_n]/D");
      myEvent->Branch("Electron_py",electron_py,"electron_py[Electron_n]/D");
      myEvent->Branch("Electron_pz",electron_pz,"electron_pz[Electron_n]/D");
      myEvent->Branch("Electron_pt",electron_pt,"electron_pt[Electron_n]/D");
      myEvent->Branch("Electron_eta",electron_eta,"electron_eta[Electron_n]/D");
      myEvent->Branch("Electron_phi",electron_phi,"electron_phi[Electron_n]/D");
      myEvent->Branch("Electron_energy",electron_energy,"electron_energy[Electron_n]/D");
      myEvent->Branch("Electron_charge",electron_charge,"electron_charge[Electron_n]/D");
      myEvent->Branch("Electron_trkIso",electron_trkIso,"electron_trkIso[Electron_n]/D");
      
    }

  if (runmuons_)
    {
      myEvent->Branch("Muon_n",&Muon_n,"Muon_n/I");
      myEvent->Branch("Muon_px",muon_px,"muon_px[Muon_n]/D");
      myEvent->Branch("Muon_py",muon_py,"muon_py[Muon_n]/D");
      myEvent->Branch("Muon_pz",muon_pz,"muon_pz[Muon_n]/D");
      myEvent->Branch("Muon_pt",muon_pt,"muon_pt[Muon_n]/D");
      myEvent->Branch("Muon_eta",muon_eta,"muon_eta[Muon_n]/D");
      myEvent->Branch("Muon_phi",muon_phi,"muon_phi[Muon_n]/D");
      myEvent->Branch("Muon_energy",muon_energy,"muon_energy[Muon_n]/D");
      myEvent->Branch("Muon_charge",muon_charge,"muon_charge[Muon_n]/D");
    }

  if (runtaus_)
    {
      myEvent->Branch("Tau_n",&Tau_n,"Tau_n/I");
      myEvent->Branch("Tau_px",tau_px,"tau_px[Tau_n]/D");
      myEvent->Branch("Tau_py",tau_py,"tau_py[Tau_n]/D");
      myEvent->Branch("Tau_pz",tau_pz,"tau_pz[Tau_n]/D");
      myEvent->Branch("Tau_pt",tau_pt,"tau_pt[Tau_n]/D");
      myEvent->Branch("Tau_eta",tau_eta,"tau_eta[Tau_n]/D");
      myEvent->Branch("Tau_phi",tau_phi,"tau_phi[Tau_n]/D");
      myEvent->Branch("Tau_energy",tau_energy,"tau_energy[Tau_n]/D");
      myEvent->Branch("Tau_charge",tau_charge,"tau_charge[Tau_n]/D");
    }

  if( rungenParticleCandidates_ )
    {
      //genlevel information from photons
      myEvent->Branch("ngenphotons",&ngenphotons,"ngenphotons/I");
      myEvent->Branch("gen_photonpt",gen_pho_pt,"gen_pho_pt[ngenphotons]/D");
      myEvent->Branch("gen_photoneta",gen_pho_eta,"gen_pho_eta[ngenphotons]/D");
      myEvent->Branch("gen_photonphi",gen_pho_phi,"gen_pho_phi[ngenphotons]/D");
      myEvent->Branch("gen_photonpx",gen_pho_px,"gen_pho_px[ngenphotons]/D");
      myEvent->Branch("gen_photonpy",gen_pho_py,"gen_pho_py[ngenphotons]/D");
      myEvent->Branch("gen_photonpz",gen_pho_pz,"gen_pho_pz[ngenphotons]/D");
      myEvent->Branch("gen_photonE",gen_pho_E,"gen_pho_E[ngenphotons]/D");
      myEvent->Branch("gen_photonstatus",gen_pho_status,"gen_pho_status[ngenphotons]/I");
      myEvent->Branch("gen_photonMotherID",gen_pho_motherID,"gen_pho_motherID[ngenphotons]/I");
      myEvent->Branch("gen_photonMotherPt",gen_pho_motherPt,"gen_pho_motherPt[ngenphotons]/D");
      myEvent->Branch("gen_photonMotherEta",gen_pho_motherEta,"gen_pho_motherEta[ngenphotons]/D");
      myEvent->Branch("gen_photonMotherPhi",gen_pho_motherPhi,"gen_pho_motherPhi[ngenphotons]/D");
      myEvent->Branch("gen_photonMotherStatus",gen_pho_motherStatus,"gen_pho_motherStatus[ngenphotons]/I");
      myEvent->Branch("gen_photonGrandmotherID",gen_pho_GrandmotherID,"gen_pho_GrandmotherID[ngenphotons]/I");
      myEvent->Branch("gen_photonGrandmotherPt",gen_pho_GrandmotherPt,"gen_pho_GrandmotherPt[ngenphotons]/D");
      myEvent->Branch("gen_photonGrandmotherEta",gen_pho_GrandmotherEta,"gen_pho_GrandmotherEta[ngenphotons]/D");
      myEvent->Branch("gen_photonGrandmotherPhi",gen_pho_GrandmotherPhi,"gen_pho_GrandmotherPhi[ngenphotons]/D");
      myEvent->Branch("gen_photonGrandmotherStatus",gen_pho_GrandmotherStatus,"gen_pho_GrandmotherStatus[ngenphotons]/I");


      myEvent->Branch("nhardphotons",&nhardphotons,"nhardphotons/I");
      myEvent->Branch("gen_hardphotonpt",gen_Hpho_pt,"gen_Hpho_pt[nhardphotons]/D");
      myEvent->Branch("gen_hardphotoneta",gen_Hpho_eta,"gen_Hpho_eta[nhardphotons]/D");
      myEvent->Branch("gen_hardphotonphi",gen_Hpho_phi,"gen_Hpho_phi[nhardphotons]/D");
      myEvent->Branch("gen_hardphotonpx",gen_Hpho_px,"gen_Hpho_px[nhardphotons]/D");
      myEvent->Branch("gen_hardphotonpy",gen_Hpho_py,"gen_Hpho_py[nhardphotons]/D");
      myEvent->Branch("gen_hardphotonpz",gen_Hpho_pz,"gen_Hpho_pz[nhardphotons]/D");
      myEvent->Branch("gen_hardphotonE",gen_Hpho_E,"gen_Hpho_E[nhardphotons]/D");

      //gen level graviton info
      myEvent->Branch("gen_gravitonpt",&gen_graviton_pt,"gen_graviton_pt/D");
      myEvent->Branch("gen_gravitoneta",&gen_graviton_eta,"gen_graviton_eta/D");
      myEvent->Branch("gen_gravitonphi",&gen_graviton_phi,"gen_graviton_phi/D");
      myEvent->Branch("gen_gravitonpx",&gen_graviton_px,"gen_graviton_px/D");
      myEvent->Branch("gen_gravitonpy",&gen_graviton_py,"gen_graviton_py/D");
      myEvent->Branch("gen_gravitonpz",&gen_graviton_pz,"gen_graviton_pz/D");
      myEvent->Branch("gen_gravitonE",&gen_graviton_E,"gen_graviton_E/D");
    
      //genlevel tree info of W+/W- 
      //genlevel tree information of Wdaughter
      myEvent->Branch("gen_Wdaughterpt",gen_Wdaughter_pt,"gen_Wdaughter_pt[2]/D");
      myEvent->Branch("gen_Wdaughtereta",gen_Wdaughter_eta,"gen_Wdaughter_eta[2]/D");
      myEvent->Branch("gen_Wdaughterphi",gen_Wdaughter_phi,"gen_Wdaughter_phi[2]/D");
      myEvent->Branch("gen_Wdaughterpx",gen_Wdaughter_px,"gen_Wdaughter_px[2]/D");
      myEvent->Branch("gen_Wdaughterpy",gen_Wdaughter_py,"gen_Wdaughter_py[2]/D");
      myEvent->Branch("gen_Wdaughterpz",gen_Wdaughter_pz,"gen_Wdaughter_pz[2]/D");
      myEvent->Branch("gen_WdaughterE",gen_Wdaughter_E,"gen_Wdaughter_E[2]/D");
      myEvent->Branch("gen_Wdaughter_charge",gen_Wdaughter_charge,"gen_Wdaughter_charge[2]/I");
      myEvent->Branch("gen_WdaughterID",gen_Wdaughter_ID,"gen_Wdaughter_ID[2]/I");

      //genlevel tree information of W
      myEvent->Branch("gen_Wbosonpt",&gen_Wboson_pt,"gen_Wboson_pt/D");
      myEvent->Branch("gen_Wbosoneta",&gen_Wboson_eta,"gen_Wboson_eta/D");
      myEvent->Branch("gen_Wbosonphi",&gen_Wboson_phi,"gen_Wboson_phi/D");
      myEvent->Branch("gen_Wbosonpx",&gen_Wboson_px,"gen_Wboson_px/D");
      myEvent->Branch("gen_Wbosonpy",&gen_Wboson_py,"gen_Wboson_py/D");
      myEvent->Branch("gen_Wbosonpz",&gen_Wboson_pz,"gen_Wboson_pz/D");
      myEvent->Branch("gen_WbosonE",&gen_Wboson_E,"gen_Wboson_E/D");
      myEvent->Branch("gen_Wbosoncharge",&gen_Wboson_charge,"gen_Wboson_charge/I");
      myEvent->Branch("gen_WbosonID",&gen_Wboson_ID,"gen_Wboson_ID/I");
      
      //genlevel tree information of Zdaughter
      myEvent->Branch("gen_Zdaughterpt",gen_Zdaughter_pt,"gen_Zdaughter_pt[2]/D");
      myEvent->Branch("gen_Zdaughtereta",gen_Zdaughter_eta,"gen_Zdaughter_eta[2]/D");
      myEvent->Branch("gen_Zdaughterphi",gen_Zdaughter_phi,"gen_Zdaughter_phi[2]/D");
      myEvent->Branch("gen_Zdaughterpx",gen_Zdaughter_px,"gen_Zdaughter_px[2]/D");
      myEvent->Branch("gen_Zdaughterpy",gen_Zdaughter_py,"gen_Zdaughter_py[2]/D");
      myEvent->Branch("gen_Zdaughterpz",gen_Zdaughter_pz,"gen_Zdaughter_pz[2]/D");
      myEvent->Branch("gen_ZdaughterE",gen_Zdaughter_E,"gen_Zdaughter_E[2]/D");
      myEvent->Branch("gen_Zdaughter_charge",gen_Zdaughter_charge,"gen_Zdaughter_charge[2]/I");
      myEvent->Branch("gen_ZdaughterID",gen_Zdaughter_ID,"gen_Zdaughter_ID[2]/I");

      //genlevel tree information of Z
      myEvent->Branch("gen_Zbosonpt",&gen_Zboson_pt,"gen_Zboson_pt/D");
      myEvent->Branch("gen_Zbosoneta",&gen_Zboson_eta,"gen_Zboson_eta/D");
      myEvent->Branch("gen_Zbosonphi",&gen_Zboson_phi,"gen_Zboson_phi/D");
      myEvent->Branch("gen_Zbosonpx",&gen_Zboson_px,"gen_Zboson_px/D");
      myEvent->Branch("gen_Zbosonpy",&gen_Zboson_py,"gen_Zboson_py/D");
      myEvent->Branch("gen_Zbosonpz",&gen_Zboson_pz,"gen_Zboson_pz/D");
      myEvent->Branch("gen_ZbosonE",&gen_Zboson_E,"gen_Zboson_E/D");

      myEvent->Branch("is_signal_event",&is_signal_event,"is_signal_event/I");
      myEvent->Branch("is_Z_event",&is_Z_event,"is_Z_event/I");
      myEvent->Branch("is_W_event",&is_W_event,"is_W_event/I");
      myEvent->Branch("is_Znunu_event",&is_Znunu_event,"is_Znunu_event/I");
      myEvent->Branch("is_Zelec_event",&is_Zelec_event,"is_Zelec_event/I");
      myEvent->Branch("is_Zmu_event",&is_Zmu_event,"is_Zmu_event/I");      
      myEvent->Branch("is_Ztu_event",&is_Ztau_event,"is_Ztau_event/I");   
      myEvent->Branch("is_Welec_event",&is_Welec_event,"is_Welec_event/I");
      myEvent->Branch("is_Wmu_event",&is_Wmu_event,"is_Wmu_event/I");      
      myEvent->Branch("is_Wtau_event",&is_Wtau_event,"is_Wtau_event/I");   
      myEvent->Branch("is_SingleHardPhoton_event",&is_SingleHardPhoton_event,"is_SingleHardPhoton_event/I");
      myEvent->Branch("is_diphoton_event",&is_diphoton_event,"is_diphoton_event/I");

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
      myEvent->Branch("gen_MuonID",gen_Muon_ID,"gen_Muon_ID[3]/D");
      myEvent->Branch("gen_MuonStatus",gen_Muon_Status,"gen_Muon_Status[3]/D");
      myEvent->Branch("gen_MuonPt",gen_Muon_Pt,"gen_Muon_Pt[3]/D");
      myEvent->Branch("gen_MuonDaughterpt",gen_MuonDaughter_pt,"gen_MuonDaughter_pt[3]/D");
      myEvent->Branch("gen_MuonDaughtereta",gen_MuonDaughter_eta,"gen_MuonDaughter_eta[3]/D");
      myEvent->Branch("gen_MuonDaughterphi",gen_MuonDaughter_phi,"gen_MuonDaughter_phi[3]/D");
      myEvent->Branch("gen_MuonDaughterpx",gen_MuonDaughter_px,"gen_MuonDaughter_px[3]/D");
      myEvent->Branch("gen_MuonDaughterpy",gen_MuonDaughter_py,"gen_MuonDaughter_py[3]/D");
      myEvent->Branch("gen_MuonDaughterpz",gen_MuonDaughter_pz,"gen_MuonDaughter_pz[3]/D");
      myEvent->Branch("gen_MuonDaughterE",gen_MuonDaughter_E,"gen_MuonDaughter_E[3]/D");
      myEvent->Branch("gen_MuonDaughterCharge",gen_MuonDaughter_charge,"gen_MuonDaughter_charge[3]/I");
      myEvent->Branch("gen_MuonDaughterStatus",gen_MuonDaughter_status,"gen_MuonDaughter_status[3]/I");
      myEvent->Branch("gen_MuonDaughterID",gen_MuonDaughter_ID,"gen_MuonDaughter_ID[3]/I");

      //genlevel tree information of tau daughter
      myEvent->Branch("gen_tauID",gen_tau_ID,"gen_tau_ID[3]/D");
      myEvent->Branch("gen_tauStatus",gen_tau_Status,"gen_tau_Status[3]/D");
      myEvent->Branch("gen_tauPt",gen_tau_Pt,"gen_tau_Pt[3]/D");
      myEvent->Branch("gen_tauDaughterpt",gen_tauDaughter_pt,"gen_tauDaughter_pt[3]/D");
      myEvent->Branch("gen_tauDaughtereta",gen_tauDaughter_eta,"gen_tauDaughter_eta[3]/D");
      myEvent->Branch("gen_tauDaughterphi",gen_tauDaughter_phi,"gen_tauDaughter_phi[3]/D");
      myEvent->Branch("gen_tauDaughterpx",gen_tauDaughter_px,"gen_tauDaughter_px[3]/D");
      myEvent->Branch("gen_tauDaughterpy",gen_tauDaughter_py,"gen_tauDaughter_py[3]/D");
      myEvent->Branch("gen_tauDaughterpz",gen_tauDaughter_pz,"gen_tauDaughter_pz[3]/D");
      myEvent->Branch("gen_tauDaughterE",gen_tauDaughter_E,"gen_tauDaughter_E[3]/D");
      myEvent->Branch("gen_tauDaughterCharge",gen_tauDaughter_charge,"gen_tauDaughter_charge[3]/I");
      myEvent->Branch("gen_tauDaughterStatus",gen_tauDaughter_status,"gen_tauDaughter_status[3]/I");
      myEvent->Branch("gen_tauDaughterID",gen_tauDaughter_ID,"gen_tauDaughter_ID[3]/I");

    }//end of if( rungenParticleCandidates_ )

  if (runphotons_)
    {
      //uncorrected photon information
      myEvent->Branch("Photon_n",&Photon_n,"Photon_n/I");
      myEvent->Branch("Photon_E",pho_E,"pho_E[Photon_n]/D");
      myEvent->Branch("Photon_pt",pho_pt,"pho_pt[Photon_n]/D");
      myEvent->Branch("Photon_eta",pho_eta,"pho_eta[Photon_n]/D");
      myEvent->Branch("Photon_phi",pho_phi,"pho_phi[Photon_n]/D");
      myEvent->Branch("Photon_theta",pho_theta,"pho_theta[Photon_n]/D");
      myEvent->Branch("Photon_et",pho_et,"pho_et[Photon_n]/D");
      myEvent->Branch("Photon_swissCross",pho_swissCross,"pho_swissCross[Photon_n]/D");
      myEvent->Branch("Photonr9",pho_r9,"pho_r9[Photon_n]/D");
      myEvent->Branch("Photon_e1x5",pho_e1x5,"pho_e1x5[Photon_n]/D");
      myEvent->Branch("Photon_e2x5",pho_e2x5,"pho_e2x5[Photon_n]/D");
      myEvent->Branch("Photon_e3x3",pho_e3x3,"pho_e3x3[Photon_n]/D");
      myEvent->Branch("Photon_e5x5",pho_e5x5,"pho_e5x5[Photon_n]/D");
      myEvent->Branch("Photon_r1x5",pho_r1x5,"pho_erx5[Photon_n]/D");
      myEvent->Branch("Photon_r2x5",pho_r2x5,"pho_erx5[Photon_n]/D");
      myEvent->Branch("Photon_maxEnergyXtal",pho_maxEnergyXtal,"pho_maxEnergyXtal[Photon_n]/D");
      myEvent->Branch("Photon_SigmaEtaEta",pho_SigmaEtaEta,"pho_SigmaEtaEta[Photon_n]/D");
      myEvent->Branch("Photon_SigmaIetaIeta",pho_SigmaIetaIeta,"pho_SigmaIetaIeta[Photon_n]/D");
      myEvent->Branch("Photon_Roundness",pho_roundness,"pho_roundness[Photon_n]/D");
      myEvent->Branch("Photon_Angle",pho_angle,"pho_angle[Photon_n]/D");
      myEvent->Branch("Photon_ecalRecHitSumEtConeDR03",pho_ecalRecHitSumEtConeDR03,"pho_ecalRecHitSumEtConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_hcalTowerSumEtConeDR03",pho_hcalTowerSumEtConeDR03,"pho_hcalTowerSumEtConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_trkSumPtSolidConeDR03",pho_trkSumPtSolidConeDR03,"pho_trkSumPtSolidConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_trkSumPtHollowConeDR03",pho_trkSumPtHollowConeDR03,"pho_trkSumPtHollowConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_nTrkSolidConeDR03",pho_nTrkSolidConeDR03,"pho_nTrkSolidConeDR03[Photon_n]/I");
      myEvent->Branch("Photon_nTrkHollowConeDR03",pho_nTrkHollowConeDR03,"pho_nTrkHollowConeDR03[Photon_n]/I");
      myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR03",pho_hcalDepth1TowerSumEtConeDR03,"pho_hcalDepth1TowerSumEtConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR03",pho_hcalDepth2TowerSumEtConeDR03,"pho_hcalDepth2TowerSumEtConeDR03[Photon_n]/D");
      myEvent->Branch("Photon_ecalRecHitSumEtConeDR04",pho_ecalRecHitSumEtConeDR04,"pho_ecalRecHitSumEtConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_hcalTowerSumEtConeDR04",pho_hcalTowerSumEtConeDR04,"pho_hcalTowerSumEtConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_trkSumPtSolidConeDR04",pho_trkSumPtSolidConeDR04,"pho_trkSumPtSolidConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_trkSumPtHollowConeDR04",pho_trkSumPtHollowConeDR04,"pho_trkSumPtHollowConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_nTrkSolidConeDR04",pho_nTrkSolidConeDR04,"pho_nTrkSolidConeDR04[Photon_n]/I");
      myEvent->Branch("Photon_nTrkHollowConeDR04",pho_nTrkHollowConeDR04,"pho_nTrkHollowConeDR04[Photon_n]/I");
      myEvent->Branch("Photon_hcalDepth1TowerSumEtConeDR04",pho_hcalDepth1TowerSumEtConeDR04,"pho_hcalDepth1TowerSumEtConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_hcalDepth2TowerSumEtConeDR04",pho_hcalDepth2TowerSumEtConeDR04,"pho_hcalDepth2TowerSumEtConeDR04[Photon_n]/D");
      myEvent->Branch("Photon_hasPixelSeed",pho_hasPixelSeed,"pho_hasPixelSeed[Photon_n]/I"); 
      myEvent->Branch("Photon_isEB",pho_isEB,"pho_isEB[Photon_n]/D");
      myEvent->Branch("Photon_isEE",pho_isEE,"pho_isEE[Photon_n]/D");
      myEvent->Branch("Photon_isEBGap",pho_isEBGap,"pho_isEBGap[Photon_n]/D");
      myEvent->Branch("Photon_isEEGap",pho_isEEGap,"pho_isEEGap[Photon_n]/D");
      myEvent->Branch("Photon_isEBEEGap",pho_isEBEEGap,"pho_isEBEEGap[Photon_n]/D");

      myEvent->Branch("Photon_HoE",pho_HoE,"pho_HoE[Photon_n]/D");
      myEvent->Branch("Photon_px",pho_px,"pho_px[Photon_n]/D");
      myEvent->Branch("Photon_py",pho_py,"pho_py[Photon_n]/D");
      myEvent->Branch("Photon_pz",pho_pz,"pho_pz[Photon_n]/D");
      myEvent->Branch("Photon_no_of_basic_clusters",pho_size,"pho_size[Photon_n]/I");
      
      myEvent->Branch("Photon_sc_energy",pho_sc_energy,"pho_sc_energy[Photon_n]/D");
      myEvent->Branch("Photon_sc_eta",pho_sc_eta,"pho_sc_eta[Photon_n]/D");
      myEvent->Branch("Photon_sc_phi",pho_sc_phi,"pho_sc_phi[Photon_n]/D");
      myEvent->Branch("Photon_etaWidth",pho_sc_etaWidth,"pho_sc_etaWidth[Photon_n]/D");
      myEvent->Branch("Photon_phiWidth",pho_sc_phiWidth,"pho_sc_phiWidth[Photon_n]/D");
      myEvent->Branch("Photon_sc_et",pho_sc_et,"pho_sc_et[Photon_n]/D");
      /*
      myEvent->Branch("Photon_seedTime",pho_seedTime,"pho_seedTime[Photon_n]/D");
      myEvent->Branch("Photon_seedOutOfTimeChi2",pho_seedOutOfTimeChi2,"pho_seedOutOfTimeChi2[Photon_n]/D");
      myEvent->Branch("Photon_seedChi2",pho_seedChi2,"pho_seedChi2[Photon_n]/D");
      myEvent->Branch("Photon_seedRecoFlag",pho_seedRecoFlag,"pho_seedRecoFlag[Photon_n]/D");
      myEvent->Branch("Photon_seedSeverity",pho_seedSeverity,"pho_seedSeverity[Photon_n]/D");
      */
      myEvent->Branch("matchphotonE",matchpho_E,"matchpho_E[Photon_n]/D");
      myEvent->Branch("matchphotonpt",matchpho_pt,"matchpho_pt[Photon_n]/D");
      myEvent->Branch("matchphotoneta",matchpho_eta,"matchpho_eta[Photon_n]/D");
      myEvent->Branch("matchphotonphi",matchpho_phi,"matchpho_phi[Photon_n]/D");
      myEvent->Branch("matchphotonpx",matchpho_px,"matchpho_px[Photon_n]/D");
      myEvent->Branch("matchphotonpy",matchpho_py,"matchpho_py[Photon_n]/D");
      myEvent->Branch("matchphotonpz",matchpho_pz,"matchpho_pz[Photon_n]/D");
      myEvent->Branch("ismatchedphoton",ismatchedpho,"ismatchedpho[Photon_n]/I");
            
      myEvent->Branch("Photon_ntracks",pho_nTracks,"pho_nTracks[Photon_n]/I");
      myEvent->Branch("Photon_isconverted",pho_isConverted,"pho_isConverted[Photon_n]/I");
      myEvent->Branch("Photon_pairInvmass",pho_pairInvariantMass,"pho_pairInvariantMass[Photon_n]/D");
      myEvent->Branch("Photon_pairCotThetaSeperation",pho_pairCotThetaSeparation,"pho_pairCotThetaSeparation[Photon_n]/D");
      myEvent->Branch("Photon_pairmomentumX",pho_pairMomentum_x,"pho_pairMomentum_x[Photon_n]/D");
      myEvent->Branch("Photon_pairmomentumY",pho_pairMomentum_y,"pho_pairMomentum_y[Photon_n]/D");
      myEvent->Branch("Photon_pairmomentumZ",pho_pairMomentum_z,"pho_pairMomentum_z[Photon_n]/D");
      myEvent->Branch("Photon_EoverP",pho_EoverP,"pho_EoverP[Photon_n]/D");
      myEvent->Branch("Photon_vertexX",pho_vertex_x,"pho_vertex_x[Photon_n]/D");
      myEvent->Branch("Photon_vertexY",pho_vertex_y,"pho_vertex_y[Photon_n]/D");
      myEvent->Branch("Photon_vertexZ",pho_vertex_z,"pho_vertex_z[Photon_n]/D");
      myEvent->Branch("Photon_ZOfPrimaryVertex",pho_zOfPrimaryVertex,"pho_zOfPrimaryVertex[Photon_n]/D");
      myEvent->Branch("Photon_distOfMinimumApproach",pho_distOfMinimumApproach,"pho_distOfMinimumApproach[Photon_n]/D");
      myEvent->Branch("Photon_dPhiTracksAtVtx",pho_dPhiTracksAtVtx,"pho_dPhiTracksAtVtx[Photon_n]/D");
      myEvent->Branch("Photon_dPhiTracksAtEcal",pho_dPhiTracksAtEcal,"pho_dPhiTracksAtEcal[Photon_n]/D");
      myEvent->Branch("Photon_dEtaTracksAtVtx",pho_dEtaTracksAtEcal,"pho_dEtaTracksAtEcal[Photon_n]/D");
      if(runrechit_){
      myEvent->Branch("Photon_ncrys",ncrysPhoton,"ncrysPhoton[Photon_n]/I");
      myEvent->Branch("Photon_timing_xtal",pho_timing_xtal,"pho_timing_xtal[Photon_n][100]/D");
      myEvent->Branch("Photon_timingavg_xtal",pho_timingavg_xtal,"pho_timingavg_xtal[Photon_n]/D");
      myEvent->Branch("Photon_energy_xtal",pho_energy_xtal,"pho_energy_xtal[Photon_n][100]/D");
      myEvent->Branch("Photon_ieta_xtalEB",pho_ieta_xtalEB,"pho_ieta_xtalEB[Photon_n][100]/I");
      myEvent->Branch("Photon_iphi_xtalEB",pho_iphi_xtalEB,"pho_iphi_xtalEB[Photon_n][100]/I");
      myEvent->Branch("Photon_rookFraction",pho_rookFraction,"pho_rookFraction[Photon_n]/D");
      myEvent->Branch("Photon_s9",pho_s9,"pho_s9[Photon_n]/D");
      }
    }//end of if (runphotons_)
  
  if(runmet_)
    {
      //Calomet variables
      myEvent->Branch("CaloMetSigma",&CaloMetSig,"CaloMetSig/D");
      //myEvent->Branch("CaloMetCorr",&CaloMetCorr,"CaloMetCorr/D");
      myEvent->Branch("CaloMetEt",&CaloMetEt,"CaloMetEt/D");
      myEvent->Branch("CaloMetEx",&CaloMetEx,"CaloMetEx/D");
      myEvent->Branch("CaloMetEy",&CaloMetEy,"CaloMetEy/D");
      myEvent->Branch("CaloMetEz",&CaloMetEz,"CaloMetEz/D");
      myEvent->Branch("CaloMetPhi",&CaloMetPhi,"CaloMetPhi/D");
      myEvent->Branch("CaloMetSumEt",&CaloSumEt,"CaloMetSumEt/D");
      myEvent->Branch("CaloEtFractionHadronic",&CaloEtFractionHadronic,"CaloEtFractionHadronic/D");
      myEvent->Branch("CaloEmEtFraction",&CaloEmEtFraction,"CaloEmEtFraction/D");
      myEvent->Branch("CaloHadEtInHB",&CaloHadEtInHB,"CaloHadEtInHB/D");
      myEvent->Branch("CaloHadEtInHE",&CaloHadEtInHE,"CaloHadEtInHE/D");
      myEvent->Branch("CaloHadEtInHO",&CaloHadEtInHO,"CaloHadEtInHO/D");
      myEvent->Branch("CaloHadEtInHF",&CaloHadEtInHF,"CaloHadEtInHF/D");
      myEvent->Branch("CaloEmEtInEB",&CaloEmEtInEB,"CaloEmEtInEB/D");
      myEvent->Branch("CaloEmEtInEE",&CaloEmEtInEE,"CaloEmEtInEE/D");
      myEvent->Branch("CaloEmEtInHF",&CaloEmEtInHF,"CaloEmEtInHF/D");
      myEvent->Branch("CaloMaxEtInEmTowers",&CaloMaxEtInEmTowers,"CaloMaxEtInEmTowers/D");
      myEvent->Branch("CaloMaxEtInHadTowers",&CaloMaxEtInHadTowers,"CaloMaxEtInHadTowers/D");

      if(rungenmet_)
	{
	  myEvent->Branch("genMetPt",&genMetPt,"genMetPt/D");
	  myEvent->Branch("genMetPx",&genMetPx,"genMetPx/D");
	  myEvent->Branch("genMetPy",&genMetPy,"genMetPy/D");
	  myEvent->Branch("genMetPhi",&genMetPhi,"genMetPhi/D");
	  myEvent->Branch("genMetSumEt",&genMetSumEt,"genMetSumEt/D");
	}
    }//end of if(runmet)
  if(runmet_&& runphotons_)
    myEvent->Branch("Delta_phi",&Delta_phi,"Delta_phi/D");
  if(runmet_ && runphotons_ && rungenmet_)
    myEvent->Branch("Delta_phiGEN",&Delta_phiGEN,"Delta_phiGEN/D");

  if(runPFmet_)
    {
      myEvent->Branch("PFMetPt",&PFMetPt,"PFMetPt/D");
      myEvent->Branch("PFMetPx",&PFMetPx,"PFMetPx/D");
      myEvent->Branch("PFMetPy",&PFMetPy,"PFMetPy/D");
      myEvent->Branch("PFMetPhi",&PFMetPhi,"PFMetPhi/D");
      myEvent->Branch("PFMetSumEt",&PFMetSumEt,"PFMetSumEt/D");
    }//end of if(runmet)
  if(runPFmet_&& runphotons_)
    myEvent->Branch("Delta_phiPF",&Delta_phiPF,"Delta_phiPF/D");


  if(runTCmet_)
    {
      myEvent->Branch("TCMetPt",&TCMetPt,"TCMetPt/D");
      myEvent->Branch("TCMetPx",&TCMetPx,"TCMetPx/D");
      myEvent->Branch("TCMetPy",&TCMetPy,"TCMetPy/D");
      myEvent->Branch("TCMetPhi",&TCMetPhi,"TCMetPhi/D");
      myEvent->Branch("TCMetSumEt",&TCMetSumEt,"TCMetSumEt/D");
    }//end of if(runmet)
  if(runTCmet_&& runphotons_)
    myEvent->Branch("Delta_phiTC",&Delta_phiTC,"Delta_phiTC/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob() {
  //cout<<"no of diphoton events:"<<n_diphoton_events<<endl;
  //cout<<"no of gamJet events:"<<n_gamJet_events<<endl;
  f->WriteTObject(myEvent);
  delete myEvent;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
