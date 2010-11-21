//////////////////////////////////////////////////////////
// Original Author: Tia Miceli
// Wed Aug 18 15:35:09 2010
// Purpose:
// Provide common tools for removing non collision backgrounds.
//////////////////////////////////////////////////////////

#ifndef NonCollisionBG_cxx
#define NonCollisionBG_cxx
#include "NonCollisionBG.h"

ClassImp(NonCollisionBG)
//////////////////////////////////////////////////////////
// 
// Cosmic functions
// 
//////////////////////////////////////////////////////////

//shruti's functions




//////////////////////////////////////////////////////////
// 
// Halo functions
// 
//////////////////////////////////////////////////////////

// Example useage NonCollisionBG::isHEHalo(Photon_sc_phi[x], HERecHit_subset_n, HERecHit_subset_x, HERecHit_subset_y, HERecHit_subset_energy, HERecHit_subset_time, false)
// useTime is defualt false (specified in .h file)
bool NonCollisionBG::isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime){

    float HERecHitEnergy_MIN_CUT = 1.;
    float deltaPhi_he_MAX_CUT = 0.2;
    float rho_he_MIN_CUT = 115.;
    float rho_he_MAX_CUT = 130.;
    float HERecHitTime_MAX_CUT = 0.;
    
    for(int hehit = 0; hehit < nAllHERecHits; hehit++){
	if (HERecHitEnergy[hehit]>HERecHitEnergy_MIN_CUT){
	    bool deltaPhiHE_isHalo = false;
	    bool rhoHE_isHalo = false;
	    bool timeHE_isHalo = false;
	    
	    float HERecHitPhi = TMath::ATan2(HERecHitY[hehit],HERecHitX[hehit]);
	    float deltaPhi_he = absDeltaPhi(fixPhi(photonSCphi),fixPhi(HERecHitPhi));

	    float rho_he = sqrt( HERecHitX[hehit]*HERecHitX[hehit] + HERecHitY[hehit]*HERecHitY[hehit]);
	    
	    if(deltaPhi_he<deltaPhi_he_MAX_CUT) deltaPhiHE_isHalo = true;
	    if(rho_he_MIN_CUT < rho_he && rho_he < rho_he_MAX_CUT) rhoHE_isHalo = true;
	    if(HERecHitTime[hehit]<HERecHitTime_MAX_CUT) timeHE_isHalo = true;
	    
	    if(useTime){
		if(deltaPhiHE_isHalo && rhoHE_isHalo && timeHE_isHalo) return true;
	    }else{ // don't useTime
		if(deltaPhiHE_isHalo && rhoHE_isHalo) return true;
	    }
	}// if HE hit meets energy requirements
	
    }// loop over HE hits
    return false;
}

//Example usage isTrackHalo(Photon_sc_phi[x], CosmicMuon_n, CosmicMuon_OuterTrack_InnerPoint_x, CosmicMuon_OuterTrack_InnerPoint_y)
// where Photon_sc_phi [0,2pi]
bool NonCollisionBG::isTrackHalo(float photonSCphi, int nCosMu, float CosTrackX[], float CosTrackY[]){
    float deltaPhi_cosTrack_MAX_CUT = 0.2;
    float rho_cosTrack_MIN_CUT = 115.;
    float rho_cosTrack_MAX_CUT = 170.;
    
    for(int CosTrack = 0; CosTrack < nCosMu; CosTrack++){
	bool deltaPhiCosTrack0_isHalo = false;
	bool rhoCosTrack0_isHalo = false;
	
	float CosTrackPhi = TMath::ATan2(CosTrackY[CosTrack],CosTrackX[CosTrack]);
	float deltaPhi0 = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CosTrackPhi));
	
	float rho0 = sqrt(CosTrackX[CosTrack]*CosTrackX[CosTrack] + CosTrackY[CosTrack]*CosTrackY[CosTrack]);
	
	if(deltaPhi0 <= deltaPhi_cosTrack_MAX_CUT) deltaPhiCosTrack0_isHalo = true;
	if(rho_cosTrack_MIN_CUT < rho0 && rho0 < rho_cosTrack_MAX_CUT) rhoCosTrack0_isHalo = true;
	if( deltaPhiCosTrack0_isHalo && rhoCosTrack0_isHalo) return true;
    }// loop over CosTracks
    
    return false;
}

//before using this, make sure phi1 and phi2 are defined on the same range
float NonCollisionBG::absDeltaPhi(float phi1, float phi2){
    float dPhi = fabs(phi1 - phi2);
    if(dPhi > TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
    return dPhi;
}

float NonCollisionBG::fixPhi(float phi){
    if (phi < 0.) phi+=2.0*TMath::Pi();
    return phi;
}

//////////////////////////////////////////////////////////
// 
// Spike functions
// 
//////////////////////////////////////////////////////////

// Example usage NonCollisionBG::simpleBarrelE2E9(Photon_ncrys[x], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH)
// arrayLimit is default to 100 because that's what is so for our ntuples (specified in .h file)
float NonCollisionBG::simpleBarrelE2E9(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit){
    float e9 = 0;
    float NeighborHiE = -10.;
    
    for(int cry=0;cry<arrayLimit && cry<nRH;cry++){
    
	//only consider barrel hits because only they have ieta & iphi !
	if(abs(ietaRH[cry])>85) continue;
	
	//get rid of no ieta=0 problem :P
	int ieta_pho_cry = ietaRH[cry];
	if(ieta_pho_cry<0) ieta_pho_cry++;
	
	int ieta_pho_cry_MaxE = ietaRH[0];
	if(ieta_pho_cry_MaxE<0) ieta_pho_cry_MaxE++;
	
	int delIEta = ieta_pho_cry_MaxE - ieta_pho_cry;
	
	//take care of phi wrapping
	int delIPhi = iphiRH[0] - iphiRH[cry];
	if(delIPhi==-359) delIPhi= 1;
	if(delIPhi== 359) delIPhi=-1;
	
	if(abs(delIPhi)<=1 && abs(delIEta)<=1){
	    e9+=eRH[cry];
	    if(eRH[cry] > NeighborHiE && cry!=0) //find highest neighbor to seed (that is not the seed itself)
		NeighborHiE = eRH[cry];
	}
    }
    return (NeighborHiE+eRH[0])/e9;
}

// Example usage NonCollisionBG::isE2E9Spike(Photon_ncrys[x], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH)
// arrayLimit is default to 100 because that's what is so for our ntuples (specified in .h file)
bool NonCollisionBG::isE2E9Spike(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit){
    if(simpleBarrelE2E9(nRH,ietaRH,iphiRH,eRH,arrayLimit)>0.95)
	return true;
    else
	return false;
}

// Example usage NonCollisionBG::isTimeSpike(Photon_ncrys[x],thisPho_tRH)
bool NonCollisionBG::isTimeSpike(vector<float> &tRH){
    if(tRH[0]<-3.5)
	return true;
    else
	return false;
}

#endif 
