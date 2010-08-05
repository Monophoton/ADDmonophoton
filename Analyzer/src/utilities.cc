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
