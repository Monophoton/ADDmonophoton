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
#include <iostream>
using namespace std;

class NonCollisionBG {

    public :
    NonCollisionBG();
    virtual ~NonCollisionBG();
    
    //Cosmic functions
    
    
    //Halo functions
    bool isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime=false);
    bool isTrackHalo(float photonSCphi, int nCosMu, float CosTrackX[], float CosTracksY[]);
    
    //Spike functions
    float simpleBarrelE2E9(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit = 100); //assume first entry is seed! (highest energy)
    bool isE2E9Spike(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit = 100);
    bool isTimeSpike(vector<float> &tRH);
    
    
    
    
    
    
    private:
    //Cosmic functions
    
    
    //Halo functions
    float absDeltaPhi(float phi1, float phi2);
    float fixPhi(float phi); //makes range [-pi,pi] to [0,2pi]
    
    //Spike functions
    
    
    
    ClassDef(NonCollisionBG,0)
};

NonCollisionBG::NonCollisionBG(){
}

NonCollisionBG::~NonCollisionBG(){
    //Delete stuff here made with 'new'!
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


