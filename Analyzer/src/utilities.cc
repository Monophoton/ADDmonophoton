#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
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
