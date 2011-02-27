//////////////////////////////////////////////////////////
// Original Author: Tia Miceli
// Wed Aug 18 15:35:09 2010
// Purpose:
// Provide common tools for removing non collision backgrounds.
//
// Important notice. It is assumed that the rechit arrays are
// those belonging to the photon, and that they are energy
// sorted. Element 0 should be the seed/highest-energy hit.
//
//////////////////////////////////////////////////////////

#ifndef NonCollisionBG_h
#define NonCollisionBG_h

#include <TROOT.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"

using namespace std;
using namespace ROOT;

class NonCollisionBG {

  public :
  NonCollisionBG();
  virtual ~NonCollisionBG();
  
  //public variables for constructor/destructor/FisherDCosmic
  TFile* file_COSMIC;
  TFile* file_PHOTON;
  
  //Cosmic functions
  float sigmaRCosmic(float Photon_sc_x, float Photon_sc_y, float Photon_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n);
  float fisherCosmic3(float Photon_Roundness, float Photon_Angle, bool Photon_phHasConversionTracks);
  float fisherCosmic2(float Photon_Roundness, float Photon_Angle);
   
  //Halo functions
  bool isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime=false);
  bool isTrackHalo(float photonSCphi, int nCosMu, float CosTrackX[], float CosTracksY[]);
  bool isCSCHalo(float photonSCphi, int nAllCSCSegments, float CSCSegmentX[], float CSCSegmentY[], float CSCSegmentTime[], bool useTime=false);
  float deltaPhiCSCHalo(float photonSCphi,float CSCSegmentX,float CSCSegmentY);
  
  //Spike functions
  float simpleBarrelE2E9(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit = 100); //assume first entry is seed! (highest energy)
  bool isE2E9Spike(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit = 100);
  bool isTimeSpike(vector<float> &tRH);
   
  private:
  //Cosmic functions
  vector<pair <bool,TVector3> > getEBIntersections( const TVector3& X,  const TVector3& P);
  double h(double z, double& length);
  double OuterEBRho(double zprimeprime, int iterations = 100 );
  double tDCA(const TVector3& X, const TVector3& P);
  vector <double> outerIntersectionTs(const TVector3& X, const TVector3& P, double timeDCA);
  vector <double> innerIntersectionTs(double R, const TVector3& X, const TVector3& P);
  float prob(TH1F* h, float value);
  float prob2(TH1F* h, bool Photon_phHasConversionTracks);
    
  //Halo functions
  float absDeltaPhi(float phi1, float phi2);
  float fixPhi(float phi); //makes range [-pi,pi] to [0,2pi]
    
  //Spike functions
    
  ClassDef(NonCollisionBG,0)
};

NonCollisionBG::NonCollisionBG(){
  file_PHOTON = TFile::Open("photon_MC.root","READ");
  file_COSMIC = TFile::Open("data_cosmics.root","READ");
}

NonCollisionBG::~NonCollisionBG(){
  file_PHOTON->Close();
  file_COSMIC->Close();
}

#endif 

#ifdef __MAKECINT__
#pragma link C++ class NonCollisionBG+;
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NonCollisionBG+;
#endif


