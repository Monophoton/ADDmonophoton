// -*- C++ -*-
//
// Package:    ADDmononphoton/Skimmer
// Class:      MonoPhotonSkimmer
// 
/*class MonoPhotonSkimmer
 
Description: EDFilter used to preselect monophoton events
             
 
         Implementation:
// Original Author:  Andrew Askew
//         Created:  Fri Sep 3 07:51:10 EST 2010
// $Id: MonoPhotonSkimmer.h,v $
//
//
*/
// system include files

#include <memory>
 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <string>
 
//
// class declaration
//
 
using namespace edm;
 
class MonoPhotonSkimmer : public edm::EDFilter {
 public:
  explicit MonoPhotonSkimmer(const edm::ParameterSet&);
  ~MonoPhotonSkimmer();
 
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  edm::InputTag _phoTag;
  bool _selectEE;         //Do you want to select EE photons?
                            //True enables this.

  double _ecalisoOffsetEB; //Photon Preselection has linearized cuts.
  double _ecalisoSlopeEB;  //slope * photonpt + offset is the isolation
                            //threshold.  This is ECAL EB.

  double _hcalisoOffsetEB; //Linearized cut on HCAL towers, EB.
  double _hcalisoSlopeEB;

  double _hadoveremEB;     //Flat selection cut on HadOverEM.

  double _minPhoEtEB;      //Minimum Photon ET threshold, EB.

  double _scHighEtThreshEB;//Photon preselction is turned OFF
                            //for photons above a certain ET.  This is 
                            //that threshold.

  double _ecalisoOffsetEE; //As above, but separately set for EE.
  double _ecalisoSlopeEE;
  double _hcalisoOffsetEE;
  double _hcalisoSlopeEE;
  double _hadoveremEE;
  double _minPhoEtEE;
  double _scHighEtThreshEE;
    
};
