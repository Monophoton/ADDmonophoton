#include "ADDmonophoton/Skimmer/interface/MonoPhotonSkimmer.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include <vector>
#include <TMath.h>
#include <iostream>
using namespace std;

MonoPhotonSkimmer::MonoPhotonSkimmer(const edm::ParameterSet& iConfig){

    //now do what ever initialization is needed
  _phoTag = iConfig.getParameter<edm::InputTag>("phoTag");
   _selectEE = iConfig.getParameter<bool>("selectEE");
 
  _ecalisoOffsetEB = iConfig.getParameter<double>("ecalisoOffsetEB");
  _ecalisoSlopeEB = iConfig.getParameter<double>("ecalisoSlopeEB");

  _hcalisoOffsetEB = iConfig.getParameter<double>("hcalisoOffsetEB");
  _hcalisoSlopeEB  = iConfig.getParameter<double>("hcalisoSlopeEB");

  _hadoveremEB = iConfig.getParameter<double>("hadoveremEB");
  _minPhoEtEB = iConfig.getParameter<double>("minPhoEtEB");
  _scHighEtThreshEB = iConfig.getParameter<double>("scHighEtThreshEB");
                          
  _ecalisoOffsetEE = iConfig.getParameter<double>("ecalisoOffsetEE");
  _ecalisoSlopeEE = iConfig.getParameter<double>("ecalisoSlopeEE");

  _hcalisoOffsetEE = iConfig.getParameter<double>("hcalisoOffsetEE");
  _hcalisoSlopeEE  = iConfig.getParameter<double>("hcalisoSlopeEE");

  _hadoveremEE = iConfig.getParameter<double>("hadoveremEE");
  _minPhoEtEE = iConfig.getParameter<double>("minPhoEtEE");
  _scHighEtThreshEE = iConfig.getParameter<double>("scHighEtThreshEE");

  
}
 

MonoPhotonSkimmer::~MonoPhotonSkimmer()
{
  
  // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
 
}
 
 
 //
 // member functions
 //
 
bool MonoPhotonSkimmer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
 
  Handle<reco::PhotonCollection> photonColl;
  iEvent.getByLabel(_phoTag, photonColl);
  const reco::PhotonCollection *photons = photonColl.product();  
  //Iterate over photon collection.
  std::vector<reco::Photon> PreselPhotons;
  reco::PhotonCollection::const_iterator pho;
  for (pho = (*photons).begin(); pho!= (*photons).end(); pho++){  
    if (!pho->isEB() && !_selectEE) continue;
    
    double ecalisocut = 0;
    double hcalisocut = 0;
    double hadoverem  = 0;
    double minphoet   = 0;
    double hiEtThresh = 0;
    if (pho->isEB()){
      ecalisocut = _ecalisoOffsetEB + _ecalisoSlopeEB * pho->pt();
      hcalisocut = _hcalisoOffsetEB + _hcalisoSlopeEB * pho->pt();
      hadoverem  = _hadoveremEB;
      minphoet   = _minPhoEtEB;
      hiEtThresh = _scHighEtThreshEB;
    }
    else{
      ecalisocut = _ecalisoOffsetEE + _ecalisoSlopeEE * pho->pt();
      hcalisocut = _hcalisoOffsetEE + _hcalisoSlopeEE * pho->pt();
      hadoverem  = _hadoveremEE;
      minphoet   = _minPhoEtEE;
      hiEtThresh = _scHighEtThreshEE;
    }
    
    if (pho->pt() > hiEtThresh){ 
      PreselPhotons.push_back(*pho);
    }
    else{
      if (pho->ecalRecHitSumEtConeDR04() < ecalisocut
	  && pho->hcalTowerSumEtConeDR04() < hcalisocut
	  && pho->hadronicOverEm() < hadoverem
	  && pho->pt() > minphoet
	  )
	PreselPhotons.push_back(*pho);
    }
  }//Loop over pat::Photons 
  if (PreselPhotons.size() > 0 ) return true;

  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void MonoPhotonSkimmer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MonoPhotonSkimmer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoPhotonSkimmer);
