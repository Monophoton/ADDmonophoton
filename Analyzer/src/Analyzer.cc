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
// $Id$
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
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

#include <map>
#include <string>

#include "/afs/cern.ch/user/s/sandhya/scratch0/CMSSW_3_1_2/src/Analysis/Analyzer/interface/Analyzer.h"

using namespace std;
using namespace ROOT::Math::VectorUtil ;

//utility function prototypes
double deltaphi(double phi1, double phi2);
double correct_phi(double phi);
double delta_R(double phi,double eta);
double Theta(double eta);
double Pl(double P,double Pt);


//
// class decleration
//

class PtSortCriterium{
public:
  bool operator() (reco::GenParticle p1,reco::GenParticle p2){
    return p1.pt() > p2.pt();
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
  runelectrons_(iConfig.getUntrackedParameter<bool>("runelectrons")),
  runmuons_(iConfig.getUntrackedParameter<bool>("runmuons")),
  runjets_(iConfig.getUntrackedParameter<bool>("runjets")),
  runtaus_(iConfig.getUntrackedParameter<bool>("runtaus")),
  runHLT_(iConfig.getUntrackedParameter<bool>("runHLT")),
  runtracks_(iConfig.getUntrackedParameter<bool>("runtracks")),
  runrechit_(iConfig.getUntrackedParameter<bool>("runrechit")),
  runvertex_(iConfig.getUntrackedParameter<bool>("runvertex")),
  init_(false)
{
   //now do what ever initialization is needed
  ngenphotons=0;
  nevents =0;
  n_diphoton_events=0;
  n_gamJet_events=0;
  JET_CORR.push_back("RAW");
  JET_CORR.push_back("OFF");
  JET_CORR.push_back("REL");
  JET_CORR.push_back("ABS");
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

   nevents++;
   cout<<"Event:"<<nevents<<endl;
   is_Z_event=0;
   is_nunu_event =0;
   is_elec_event =0;
   is_mu_event =0;
   is_tau_event =0;
   is_gamjet_event=0;
   is_diphoton_event=0;
  
   ///// L1
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
	 HLT_chosen[ hlNames_[i]]=HLTR->accept(i);
	 if(hlNames_[i]=="HLT_MET50")
	   is_HLT_MET50_event=HLTR->accept(i);
	 if(hlNames_[i]=="HLT_MET75")
	   is_HLT_MET75_event=HLTR->accept(i);
	 if(hlNames_[i]=="HLT_Photon25_L1R")
	   is_HLT_Photon25_event=HLTR->accept(i);
       }
   }

   if(runvertex_){
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(Vertices_, recVtxs);
   vector<Vertex> my_vertices;
   for(reco::VertexCollection::const_iterator v=recVtxs->begin();v!=recVtxs->end(); ++v){
     my_vertices.push_back(*v);
   }
   for (unsigned int y = 0; y <  my_vertices.size();y++)
     {
       vx[y]     = my_vertices[y].x();
       vy[y]     = my_vertices[y].y();
       vz[y]     = my_vertices[y].z();
       chi2[y]   = my_vertices[y].chi2();
     }
   }

   if(runtracks_){
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(Tracks_,tracks);
   nTracks=tracks->size();
   std::vector<reco::Track>  myTrack_container;
   myTrack_container.clear();
   for(reco::TrackCollection::const_iterator Track_iter = tracks->begin();
       Track_iter != tracks->end();++Track_iter) {
     myTrack_container.push_back(*Track_iter);
     
   }
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
     }//end of for loop
   }
         
   //getting handle to generator level information
   if( rungenParticleCandidates_ )
     {
       ngenphotons = 0;
       Handle<GenParticleCollection> genParticles;
       iEvent.getByLabel("genParticles", genParticles); 
       std::vector<reco::GenParticle>  mygenphoton_container;
       mygenphoton_container.clear();
       std::vector<reco::GenParticle> mygenZboson_container;
       mygenZboson_container.clear();
       std::vector<reco::GenParticle>  mygenhardphoton_container;
       mygenhardphoton_container.clear();
       int ii =0;
       for (GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++) {
	 //getting information from Z
	 if (genparticle->pdgId()==23 && genparticle->status()==3)  mygenZboson_container.push_back(*genparticle);
	 //getting information from photons
	 if (genparticle->pdgId()==22 && genparticle->status()==1)
	   {
	     mygenphoton_container.push_back(*genparticle); ngenphotons++;
	   }
	 if (genparticle->pdgId()==22 && genparticle->status()==3)  mygenhardphoton_container.push_back(*genparticle);
	 ii++;
       }//end of for loop

       if (mygenhardphoton_container.size()==1) 
	 {
	   n_gamJet_events++;
	   is_gamjet_event=1;
	 }
       if (mygenhardphoton_container.size()==2)
	 {
	   n_diphoton_events++; 
	   is_diphoton_event=1;
	 }
       std::cout<<"genphoton's container's size:"<<mygenphoton_container.size()<<std::endl;
       if(mygenphoton_container.size()!=0)
	 {
	   std::sort(mygenphoton_container.begin(),mygenphoton_container.end(),PtSortCriterium());
	  
	   for(unsigned int x=0;x < mygenphoton_container.size(); x++)
	     {
	       gen_pho_pt[x]  = mygenphoton_container[x].pt();
	       gen_pho_px[x]  = mygenphoton_container[x].px();
	       gen_pho_py[x]  = mygenphoton_container[x].py();
	       gen_pho_pz[x]  = mygenphoton_container[x].pz();
	       gen_pho_phi[x] = correct_phi(mygenphoton_container[x].phi());
	       gen_pho_eta[x] = mygenphoton_container[x].eta();
	       gen_pho_E[x]   = mygenphoton_container[x].energy();
	     }//end of for loop
	 }//end of if((mygenphoton_container.size()!=0)

       if(mygenhardphoton_container.size()!=0)
         {
	   for (unsigned int i =0 ; i < mygenhardphoton_container.size() ; i++)
	     {
	       gen_Hpho_pt[i]  = mygenhardphoton_container[i].pt();
	       gen_Hpho_px[i]  = mygenhardphoton_container[i].px();
	       gen_Hpho_py[i]  = mygenhardphoton_container[i].py();
	       gen_Hpho_pz[i]  = mygenhardphoton_container[i].pz();
	       gen_Hpho_phi[i] = correct_phi(mygenhardphoton_container[i].phi());
	       gen_Hpho_eta[i] = mygenhardphoton_container[i].eta();
	       gen_Hpho_E[i]   = mygenhardphoton_container[i].energy();
	     }
         }//end of if((mygenphoton_container.size()!=0))

       
       if(mygenZboson_container.size()!=0){
       //std::cout<<"found a Z"<<std::endl;
       is_Z_event=1;
       //std::cout<<"status of Z :"<<gen_status<<std::endl;
       int daughters          = mygenZboson_container[0].numberOfDaughters();
       //std::cout<<"no of Z daughters:"<<daughters<<std::endl;
       for(int i = 0;i<daughters;i++)
	 {
	   const reco::Candidate *daughter   = mygenZboson_container[0].daughter(i);
	   //std::cout<<"daughter's ID:"<<daughter->pdgId()<<std::endl;
	   //std::cout<<"daughter's status:"<<daughter->status()<<std::endl;
	   if((daughter->pdgId()==12)||(daughter->pdgId()==-12)||(daughter->pdgId()==14)||(daughter->pdgId()==-14)||(daughter->pdgId()==16)||(daughter->pdgId()==-16)) is_nunu_event=1;
	   if(daughter->pdgId()==11) is_elec_event=1;
	   if(daughter->pdgId()==13) is_mu_event=1;
	   if(daughter->pdgId()==15) is_tau_event=1;
	   
	   //getting leptons decaying from Z
	   int iDaughter=0;
	   if(daughter->pdgId()!=23) 
	     {
	       gen_Zdaughter_pt[iDaughter]  = daughter->pt();
	       gen_Zdaughter_px[iDaughter]  = daughter->px();
	       gen_Zdaughter_py[iDaughter]  = daughter->py();
	       gen_Zdaughter_pz[iDaughter]  = daughter->pz();
	       gen_Zdaughter_phi[iDaughter] = correct_phi(daughter->phi());
	       gen_Zdaughter_eta[iDaughter] = daughter->eta();
	       gen_Zdaughter_E[iDaughter]   = daughter->energy();
	       iDaughter++;
	     }
       
	   //getting Z
	   if(daughter->pdgId()==23) {
	     gen_Zboson_pt  = daughter->pt();
	     gen_Zboson_px  = daughter->px();
	     gen_Zboson_py  = daughter->py();
	     gen_Zboson_pz  = daughter->pz();
	     gen_Zboson_phi = correct_phi(daughter->phi());
	     gen_Zboson_eta = daughter->eta();
	     gen_Zboson_E   = daughter->energy();
	   }
	 }//end of for loop for daughters
       }//end of if(mygenZbosoncontainer.size!=0)
     }//end of if(rungenParticleCandidates_)
   std::vector<pat::Photon> myphoton_container;
   myphoton_container.clear();
   if(runphotons_)
     {
       edm::Handle<edm::View<pat::Photon> > phoHandle;
       iEvent.getByLabel(phoLabel_,phoHandle);
       //const edm::View<pat::Photon> & photons = *phoHandle;
       edm::View<pat::Photon>::const_iterator photon;
       nrecPhotons = phoHandle->size();
       cout<<"photon container size:"<<nrecPhotons<<endl;
       for(photon = phoHandle->begin();photon!=phoHandle->end();++photon){
	 myphoton_container.push_back(*photon) ;
       }
       if(myphoton_container.size()!=0)
	 {
	   for(unsigned int x=0; x < myphoton_container.size();x++)
	     {
	       cout<<"photon pt:"<<myphoton_container[x].pt()<<"  photon eta:"<<myphoton_container[x].eta()<<"  photon phi:"<<correct_phi(myphoton_container[x].phi())<<endl; 

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
	       pho_size[x]                  =  myphoton_container[x].superCluster()->clustersSize();
	       pho_sc_energy[x]             =  myphoton_container[x].superCluster()->energy();
	       pho_sc_eta[x]                =  myphoton_container[x].superCluster()->eta();
	       pho_sc_phi[x]                =  correct_phi(myphoton_container[x].superCluster()->phi());
	       pho_sc_etaWidth[x]           =  myphoton_container[x].superCluster()->etaWidth();
	       pho_sc_phiWidth[x]           =  myphoton_container[x].superCluster()->phiWidth();
	       pho_ecalRecHitIso[x]         =  myphoton_container[x].ecalRecHitSumEtConeDR03();
	       pho_hcalRecHitIso[x]         =  myphoton_container[x].hcalTowerSumEtConeDR03();
	       pho_HollowTrackConeIso[x]    =  myphoton_container[x].trkSumPtHollowConeDR03();
	       pho_HoE[x]                   =  myphoton_container[x].hadronicOverEm();              

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
	       ismatchedpho[x]              =  myphoton_container[x].genParticleRef().isNonnull();
	       /*
	       reco::ConversionRefVector ConvRef   = myphoton_container[x].conversions();
	       
	       if(myphoton_container[x].hasConversionTracks()==1 )
		 {
		   
		   pho_nTracks[x]                    = ConvRef[0]->nTracks();
		   pho_isConverted[x]                = ConvRef[0]->isConverted();
		   pho_pairInvariantMass[x]          = ConvRef[0]->pairInvariantMass();
		   pho_pairCotThetaSeparation[x]     = ConvRef[0]->pairCotThetaSeparation();
		   pho_pairMomentum_x[x]             = ConvRef[0]->pairMomentum().x();
		   pho_pairMomentum_y[x]             = ConvRef[0]->pairMomentum().y();
		   pho_pairMomentum_z[x]             = ConvRef[0]->pairMomentum().z();
		   pho_vertex_x[x]                   = ConvRef[0]->conversionVertex().x();
		   pho_vertex_y[x]                   = ConvRef[0]->conversionVertex().y();
		   pho_vertex_z[x]                   = ConvRef[0]->conversionVertex().z();
		   pho_EoverP[x]                     = ConvRef[0]->EoverP();
		   pho_zOfPrimaryVertex[x] = ConvRef[0]->zOfPrimaryVertexFromTracks();
		 }
	       else
		 {
		   pho_nTracks[x]                    = -99;
		   pho_isConverted[x]                =   0;
		   pho_pairInvariantMass[x]          = -99.;
		   pho_pairCotThetaSeparation[x]     = -99.;
		   pho_pairMomentum[x]               = -99.;
		   pho_pairMomentum_x[x]             = -99.;
		   pho_pairMomentum_y[x]             = -99.;
		   pho_pairMomentum_z[x]             = -99.;
		   pho_vertexPosition_x[x]           = -99.;
		   pho_vertexPosition_y[x]           = -99.;
		   pho_vertexPosition_z[x]           = -99.;
		   pho_EoverP[x]                     = -99.;
		   pho_zOfPrimaryVertex[x]           = -99.;
		 }
	       */
	       //to get the photon hit information from every crystal of SC
	       if(runrechit_)
		 {
		   std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = myphoton_container[x].superCluster()->hitsAndFractions();
		   Handle<EcalRecHitCollection> Brechit;//barrel                                                                            
		   Handle<EcalRecHitCollection> Erechit;//endcap                                                                            
		   iEvent.getByLabel(rechitBLabel_,Brechit);
		   iEvent.getByLabel(rechitELabel_,Erechit); 
		   std::vector<CrystalInfo> crystalinfo_container;
		   crystalinfo_container.clear();
		   CrystalInfo crystal;
		   double timing_avg =0.0;
		   int ncrys = 0;
		   vector< std::pair<DetId, float> >::const_iterator detitr;
		   
		   for(detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr)
		     {
		       if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) 
			 {
			   EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
			   EcalRecHitCollection::const_iterator thishit;
			   if (j!= Brechit->end())  thishit = j;
			   if ( j== Brechit->end())
			     {
			       std::cout<<"thishit not matched "<<std::endl;
			       continue;
			     }
			   EBDetId detId  = (EBDetId)((*detitr).first);
			   crystal.rawId  = thishit->id().rawId();
			   crystal.energy = thishit->energy();
			   crystal.time   = thishit->time();
			   crystal.ieta   = detId.ieta();
			   crystal.iphi   = detId.iphi();
			   //calculate timing avg
			   if(crystal.energy > 0.1)
			     {
			       timing_avg  = timing_avg + crystal.time;
			       ncrys++;
			     } 
			 }//end of if ((*detitr).det() == DetId::Ecal && (*detitr).subdetId() == EcalBarrel)
		       else 
			 {
			   crystal.rawId  = 999;
			   crystal.energy = -99;
			   crystal.time   = -99;
			   crystal.ieta   = -99;
			   crystal.iphi   = -99;
			 }
		       crystalinfo_container.push_back(crystal);  
		     }
		   std::sort(crystalinfo_container.begin(),crystalinfo_container.end(),EnergySortCriterium());
		   if (ncrys !=0) timing_avg = timing_avg/(double)ncrys;
		   else timing_avg = -99.;
		   cout<<" total hits for this photon:"<<crystalinfo_container.size()<<endl;
		   cout<<" total hits contibuting for timing:"<<ncrys<<endl;
		   pho_timingavg_xtalEB[x]      = timing_avg;
		   for (unsigned int y =0; y < crystalinfo_container.size();y++)
		     {
		       pho_timing_xtalEB[x][y]         = crystalinfo_container[y].time;
		       pho_energy_xtalEB[x][y]         = crystalinfo_container[y].energy;
		       pho_ieta_xtalEB[x][y]           = crystalinfo_container[y].ieta;
		       pho_iphi_xtalEB[x][y]           = crystalinfo_container[y].iphi;
		     }
		 }//if(runrechit_)
	     }//end of for loop over x
	 }//if(myphoton_container.size!=0) 
       cout<<"got photon variables"<<endl; 
     }
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
	       Delta_phi                             = deltaphi(correct_phi(met->phi()),correct_phi(myphoton_container[0].phi()));
	 }
     }
   
   if(runjets_)
     {
       edm::Handle<edm::View<pat::Jet> > jetHandle;
       iEvent.getByLabel(jetLabel_,jetHandle);
       const edm::View<pat::Jet> & jets = *jetHandle;
       number_of_jets=jetHandle->size();
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
       nmuons = muonHandle->size();
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
       nelectrons = electronHandle->size();
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
       ntaus = tauHandle->size();
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
   myEvent->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob(const edm::EventSetup&)
{
  f=new TFile(outFile_.c_str(),"RECREATE");
  //defining a tree here
  myEvent = new TTree("myEvent","a tree with histograms");
  myEvent->Branch("nevents",&nevents,"nevents/I");
  
  if(runHLT_)
    {
      myEvent->Branch("is_HLT_MET50_event",&is_HLT_MET50_event,"is_HLT_MET50_event/I");
      myEvent->Branch("is_HLT_MET75_event",&is_HLT_MET75_event,"is_HLT_MET75_event/I");
      myEvent->Branch("is_HLT_Photon25_event",&is_HLT_Photon25_event,"is_HLT_Photon25_event/I");
    }

  if(runvertex_)
    {
      myEvent->Branch("vx",vx,"vx[100]/F");
      myEvent->Branch("vy",vy,"vy[100]/F");
      myEvent->Branch("vz",vz,"vz[100]/F");
      myEvent->Branch("chi2",chi2,"chi2[100]/F");
    }

  if (runtracks_)
    {
      myEvent->Branch("nTracks",&nTracks,"nTracks/I");
      myEvent->Branch("trackpx",trk_px,"trk_px[250]/D");
      myEvent->Branch("trackpy",trk_py,"trk_py[250]/D");
      myEvent->Branch("trackpz",trk_pz,"trk_pz[250]/D");
      myEvent->Branch("trackpt",trk_pt,"trk_pt[250]/D");
      myEvent->Branch("tracketa",trk_eta,"trk_eta[250]/D");
      myEvent->Branch("trackphi",trk_phi,"trk_phi[250]/D");
    }
  if (runjets_)
    {
      //for(vector<string>::iterator corr = JET_CORR.begin(); corr!=JET_CORR.end();++corr )
      //{
	  myEvent->Branch("number_of_jets",&number_of_jets,"number_of_jets/I");
	  myEvent->Branch("jetpx",jet_px,"jet_px[250]/D");
	  myEvent->Branch("jetpy",jet_py,"jet_py[250]/D");
	  myEvent->Branch("jetpz",jet_pz,"jet_pz[250]/D");
	  myEvent->Branch("jetpt",jet_pt,"jet_pt[250]/D");
	  myEvent->Branch("jeteta",jet_eta,"jet_eta[250]/D");
	  myEvent->Branch("jetphi",jet_phi,"jet_phi[250]/D");
	  myEvent->Branch("jetemEnergyFraction",jet_emEnergyFraction,"jet_emEnergyFraction[250]/D");
	  myEvent->Branch("jetenergyFractionHadronic",jet_energyFractionHadronic,"jet_energyFractionHadronic[250]/D");
	  //}
    }

  if (runelectrons_)
    {
      myEvent->Branch("number_of_electrons",&nelectrons,"nelectrons/I");
      myEvent->Branch("electronpx",electron_px,"electron_px[250]/D");
      myEvent->Branch("electronpy",electron_py,"electron_py[250]/D");
      myEvent->Branch("electronpz",electron_pz,"electron_pz[250]/D");
      myEvent->Branch("electronpt",electron_pt,"electron_pt[250]/D");
      myEvent->Branch("electroneta",electron_eta,"electron_eta[250]/D");
      myEvent->Branch("electronphi",electron_phi,"electron_phi[250]/D");
      myEvent->Branch("electronenergy",electron_energy,"electron_energy[250]/D");
      myEvent->Branch("electroncharge",electron_charge,"electron_charge[250]/D");
      myEvent->Branch("electrontrkIso",electron_trkIso,"electron_trkIso[250]/D");
      
    }

  if (runmuons_)
    {
      myEvent->Branch("number_of_muons",&nmuons,"nmuons/I");
      myEvent->Branch("muonpx",muon_px,"muon_px[250]/D");
      myEvent->Branch("muonpy",muon_py,"muon_py[250]/D");
      myEvent->Branch("muonpz",muon_pz,"muon_pz[250]/D");
      myEvent->Branch("muonpt",muon_pt,"muon_pt[250]/D");
      myEvent->Branch("muoneta",muon_eta,"muon_eta[250]/D");
      myEvent->Branch("muonphi",muon_phi,"muon_phi[250]/D");
      myEvent->Branch("muonenergy",muon_energy,"muon_energy[250]/D");
      myEvent->Branch("muoncharge",muon_charge,"muon_charge[250]/D");
    }

  if (runtaus_)
    {
      myEvent->Branch("number_of_taus",&ntaus,"ntaus/I");
      myEvent->Branch("taupx",tau_px,"tau_px[250]/D");
      myEvent->Branch("taupy",tau_py,"tau_py[250]/D");
      myEvent->Branch("taupz",tau_pz,"tau_pz[250]/D");
      myEvent->Branch("taupt",tau_pt,"tau_pt[250]/D");
      myEvent->Branch("taueta",tau_eta,"tau_eta[250]/D");
      myEvent->Branch("tauphi",tau_phi,"tau_phi[250]/D");
      myEvent->Branch("tauenergy",tau_energy,"tau_energy[250]/D");
      myEvent->Branch("taucharge",tau_charge,"tau_charge[250]/D");
    }

  if( rungenParticleCandidates_ )
    {
      //genlevel information
      myEvent->Branch("ngenphotons",&ngenphotons,"ngenphotons/I");
      myEvent->Branch("gen_photonpt",gen_pho_pt,"gen_pho_pt[250]/D");
      myEvent->Branch("gen_photoneta",gen_pho_eta,"gen_pho_eta[250]/D");
      myEvent->Branch("gen_photonphi",gen_pho_phi,"gen_pho_phi[250]/D");
      myEvent->Branch("gen_photonpx",gen_pho_px,"gen_pho_px[250]/D");
      myEvent->Branch("gen_photonpy",gen_pho_py,"gen_pho_py[250]/D");
      myEvent->Branch("gen_photonpz",gen_pho_pz,"gen_pho_pz[250]/D");
      myEvent->Branch("gen_photonE",gen_pho_E,"gen_pho_E[250]/D");

      myEvent->Branch("gen_hardphotonpt",gen_Hpho_pt,"gen_Hpho_pt[2]/D");
      myEvent->Branch("gen_hardphotoneta",gen_Hpho_eta,"gen_Hpho_eta[2]/D");
      myEvent->Branch("gen_hardphotonphi",gen_Hpho_phi,"gen_Hpho_phi[2]/D");
      myEvent->Branch("gen_hardphotonpx",gen_Hpho_px,"gen_Hpho_px[2]/D");
      myEvent->Branch("gen_hardphotonpy",gen_Hpho_py,"gen_Hpho_py[2]/D");
      myEvent->Branch("gen_hardphotonpz",gen_Hpho_pz,"gen_Hpho_pz[2]/D");
      myEvent->Branch("gen_hardphotonE",gen_Hpho_E,"gen_Hpho_E[2]/D");
      myEvent->Branch("gen_hardphotonSerialNumber",gen_Hpho_SerialNumber,"gen_Hpho_SerialNumber[2]/I");
      myEvent->Branch("is_gamjet_event",&is_gamjet_event,"is_gamjet_event/I");
      myEvent->Branch("is_diphoton_event",&is_diphoton_event,"is_diphoton_event/I");


      //genlevel tree information of Zdaughter
      myEvent->Branch("gen_Zdaughterpt",gen_Zdaughter_pt,"gen_Zdaughter_pt[2]/D");
      myEvent->Branch("gen_Zdaughtereta",gen_Zdaughter_eta,"gen_Zdaughter_eta[2]/D");
      myEvent->Branch("gen_Zdaughterphi",gen_Zdaughter_phi,"gen_Zdaughter_phi[2]/D");
      myEvent->Branch("gen_Zdaughterpx",gen_Zdaughter_px,"gen_Zdaughter_px[2]/D");
      myEvent->Branch("gen_Zdaughterpy",gen_Zdaughter_py,"gen_Zdaughter_py[2]/D");
      myEvent->Branch("gen_Zdaughterpz",gen_Zdaughter_pz,"gen_Zdaughter_pz[2]/D");
      myEvent->Branch("gen_ZdaughterE",gen_Zdaughter_E,"gen_Zdaughter_E[2]/D");

      //genlevel tree information of Z
      myEvent->Branch("is_nunu_event",&is_nunu_event,"is_nunu_event/I");
      myEvent->Branch("is_elec_event",&is_elec_event,"is_elec_event/I");
      myEvent->Branch("is_mu_event",&is_mu_event,"is_mu_event/I");
      myEvent->Branch("is_tau_event",&is_tau_event,"is_tau_event/I");
      myEvent->Branch("is_Z_event",&is_Z_event,"is_Z_event/I");
      myEvent->Branch("gen_Zbosonpt",&gen_Zboson_pt,"gen_Zboson_pt/D");
      myEvent->Branch("gen_Zbosoneta",&gen_Zboson_eta,"gen_Zboson_eta/D");
      myEvent->Branch("gen_Zbosonphi",&gen_Zboson_phi,"gen_Zboson_phi/D");
      myEvent->Branch("gen_Zbosonpx",&gen_Zboson_px,"gen_Zboson_px/D");
      myEvent->Branch("gen_Zbosonpy",&gen_Zboson_py,"gen_Zboson_py/D");
      myEvent->Branch("gen_Zbosonpz",&gen_Zboson_pz,"gen_Zboson_pz/D");
      myEvent->Branch("gen_ZbosonE",&gen_Zboson_E,"gen_Zboson_E/D");
    }//end of if( rungenParticleCandidates_ )

  if (runphotons_)
    {
      //uncorrected photon information
      myEvent->Branch("nrecPhotons",&nrecPhotons,"nrecPhotons/I");
      myEvent->Branch("photonE",pho_E,"pho_E[250]/D");
      myEvent->Branch("photonpt",pho_pt,"pho_pt[250]/D");
      myEvent->Branch("photoneta",pho_eta,"pho_eta[250]/D");
      myEvent->Branch("photonphi",pho_phi,"pho_phi[250]/D");
      myEvent->Branch("photonr9",pho_r9,"pho_r9[250]/D");
      myEvent->Branch("photon_theta",pho_theta,"pho_theta[250]/D");
      myEvent->Branch("photon_et",pho_et,"pho_et[250]/D");
      myEvent->Branch("photon_ecalRecHitIso",pho_ecalRecHitIso,"pho_ecalRecHitIso[250]/D");
      myEvent->Branch("photon_hcalRecHitIso",pho_hcalRecHitIso,"pho_hcalRecHitIso[250]/D");
      myEvent->Branch("photon_HollowTrackConeIso",pho_HollowTrackConeIso,"pho_HollowTrackConeIso[250]/D");
      myEvent->Branch("photon_HoE",pho_HoE,"pho_HoE[250]/D");
      myEvent->Branch("photonpx",pho_px,"pho_px[250]/D");
      myEvent->Branch("photonpy",pho_py,"pho_py[250]/D");
      myEvent->Branch("photonpz",pho_pz,"pho_pz[250]/D");
      myEvent->Branch("photon_no_of_basic_clusters",pho_size,"pho_size[250]/I");
      
      myEvent->Branch("photon_sc_energy",pho_sc_energy,"pho_sc_energy[250]/D");
      myEvent->Branch("photon_sc_eta",pho_sc_eta,"pho_sc_eta[250]/D");
      myEvent->Branch("photon_sc_phi",pho_sc_phi,"pho_sc_phi[250]/D");
      myEvent->Branch("photon_etaWidth",pho_sc_etaWidth,"pho_sc_etaWidth[250]/D");
      myEvent->Branch("photon_phiWidth",pho_sc_phiWidth,"pho_sc_phiWidth[250]/D");
      myEvent->Branch("photon_sc_et",pho_sc_et,"pho_sc_et[250]/D");
      
      myEvent->Branch("matchphotonE",matchpho_E,"matchpho_E[250]/D");
      myEvent->Branch("matchphotonpt",matchpho_pt,"matchpho_pt[250]/D");
      myEvent->Branch("matchphotoneta",matchpho_eta,"matchpho_eta[250]/D");
      myEvent->Branch("matchphotonphi",matchpho_phi,"matchpho_phi[250]/D");
      myEvent->Branch("matchphotonpx",matchpho_px,"matchpho_px[250]/D");
      myEvent->Branch("matchphotonpy",matchpho_py,"matchpho_py[250]/D");
      myEvent->Branch("matchphotonpz",matchpho_pz,"matchpho_pz[250]/D");
      myEvent->Branch("ismatchedphoton",ismatchedpho,"ismatchedpho[250]/I");
           
      /*
      myEvent->Branch("photon_ntracks",pho_nTracks,"pho_nTracks[250]/I");
      myEvent->Branch("photon_isconverted",pho_isConverted,"pho_isConverted[250]/I");
      myEvent->Branch("photon_pairInvmass",pho_pairInvariantMass,"pho_pairInvariantMass[250]/D");
      myEvent->Branch("photon_pairCotThetaSeperation",pho_pairCotThetaSeparation,"pho_pairCotThetaSeparatio\
n[250]/D");
      myEvent->Branch("photon_pairmomentumX",pho_pairMomentum_x,"pho_pairMomentum_x[250]/D");
      myEvent->Branch("photon_pairmomentumY",pho_pairMomentum_y,"pho_pairMomentum_y[250]/D");
      myEvent->Branch("photon_pairmomentumZ",pho_pairMomentum_z,"pho_pairMomentum_z[250]/D");
      myEvent->Branch("photon_EoverP",pho_EoverP,"pho_EoverP[250]/D");
      myEvent->Branch("photon_vertexX",pho_vertex_x,"pho_vertex_x[250]/D");
      myEvent->Branch("photon_vertexY",pho_vertex_y,"pho_vertex_y[250]/D");
      myEvent->Branch("photon_vertexZ",pho_vertex_z,"pho_vertex_z[250]/D");
      myEvent->Branch("photon_ZOfPrimaryVertex",pho_zOfPrimaryVertex,"pho_zOfPrimaryVertex[250]/D");
      */

      myEvent->Branch("photon_timing_xtalEB",pho_timing_xtalEB,"pho_timing_xtalEB[250][100]/D");
      myEvent->Branch("photon_timingavg_xtalEB",pho_timingavg_xtalEB,"pho_timingavg_xtalEB[250]/D");
      myEvent->Branch("photon_energy_xtalEB",pho_energy_xtalEB,"pho_energy_xtalEB[250][100]/D");
      myEvent->Branch("photon_ieta_xtalEB",pho_ieta_xtalEB,"pho_ieta_xtalEB[250][100]/D");
      myEvent->Branch("photon_iphi_xtalEB",pho_iphi_xtalEB,"pho_iphi_xtalEB[250][100]/D");
    }//end of if (runphotons_)
  
  if(runmet_)
    {
      //Calomet variables
      myEvent->Branch("CaloMetSigma",&CaloMetSig,"CaloMetSig/D");
      myEvent->Branch("CaloMetCorr",&CaloMetCorr,"CaloMetCorr/D");
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
    }//end of if(runmet)
  if(runmet_&& runphotons_)
    myEvent->Branch("Delta_phi",&Delta_phi,"Delta_phi/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob() {
  cout<<"no of diphoton events:"<<n_diphoton_events<<endl;
  cout<<"no of gamJet events:"<<n_gamJet_events<<endl;
  f->WriteTObject(myEvent);
  delete myEvent;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
