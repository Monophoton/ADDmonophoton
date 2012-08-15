#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "TMath.h"

double correct_phi(double phi){
	return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}
double Theta(double eta){
  double theta = 2. * atan(exp(-eta));
  return theta;
}
double Pl(double P, double Pt){
  double pl = sqrt(pow(P,2)-pow(Pt,2));
  return pl;
}

float correct_phi(float phi){
        return (phi >= 0 ? phi : (2*TMath::Pi() + phi));
}
float Theta(float eta){
  float theta = 2. * atan(exp(-eta));
  return theta;
}
float Pl(float P, float Pt){
  float pl = sqrt(pow(P,2)-pow(Pt,2));
  return pl;
}



//used for E2E9 value calculations
float recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}


float recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY
  
  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );
  
  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}



float recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits )
{
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}




float recHitE( const DetId id, const EcalRecHitCollection &recHits, bool useTimingInfo ){
  if ( id.rawId() == 0 ) return 0;

 //These values are taken from:RecoLocalCalo/EcalRecAlgos/python/ecalCleaningAlgo.py
 double e4e1Threshold_barrel_  = 0.080;
 double e4e1Threshold_endcap_  = 0.300;
 double ignoreOutOfTimeThresh_ = 1e9;
 
 
   float threshold = e4e1Threshold_barrel_;

   if ( id.subdetId() == EcalEndcap) threshold = e4e1Threshold_endcap_; 
 
   EcalRecHitCollection::const_iterator it = recHits.find( id );
   if ( it != recHits.end() ){
     float ene= (*it).energy();
 
     // ignore out of time in EB when making e4e1 if so configured
     if (useTimingInfo){
        if (id.subdetId()==EcalBarrel &&
          it->checkFlag(EcalRecHit::kOutOfTime) 
           && ene>ignoreOutOfTimeThresh_) return 0;
     }
 
     // ignore hits below threshold
     if (ene < threshold) return 0;
 
     // else return the energy of this hit
     return ene;
   }
   return 0;
}

/// four neighbours in the swiss cross around id
const std::vector<DetId> neighbours(const DetId& id){                                                                                                                              
        std::vector<DetId> ret;
    
   if ( id.subdetId() == EcalBarrel) {
    
     ret.push_back( EBDetId::offsetBy( id,  1, 0 ));
     ret.push_back( EBDetId::offsetBy( id, -1, 0 ));
     ret.push_back( EBDetId::offsetBy( id,  0, 1 ));
     ret.push_back( EBDetId::offsetBy( id,  0,-1 ));
   }
   // nobody understands what polymorphism is for, sgrunt !
   else  if (id.subdetId() == EcalEndcap) {
     ret.push_back( EEDetId::offsetBy( id,  1, 0 ));
     ret.push_back( EEDetId::offsetBy( id, -1, 0 ));
     ret.push_back( EEDetId::offsetBy( id,  0, 1 ));
     ret.push_back( EEDetId::offsetBy( id,  0,-1 ));
    }
    
    
   return ret;
    
 }  

