#include "TMath.h"

double pi=acos(-1.0);

double deltaphi(double phi1, double phi2){
 double delta_phi=phi1-phi2;
 if(delta_phi<0) delta_phi=-delta_phi;
 if(delta_phi>=(2*pi-delta_phi)) delta_phi= 2.0*pi-delta_phi;
 return delta_phi;
}
double correct_phi(double phi){
  return (phi >= 0 ? phi : (2*pi + phi));
}
double delta_R(double phi,double eta){
  double deltar=pow((pow(phi,2)+pow(eta,2)),0.5);
  return deltar;
}
double Theta(double eta){
  double theta = 2. * atan(exp(-eta));
  return theta;
}
double Pl(double P, double Pt){
  double pl = sqrt(pow(P,2)-pow(Pt,2));
  return pl;
}
