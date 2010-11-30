//////////////////////////////////////////////////////////
// Original Author: Tia Miceli
// Wed Aug 18 15:35:09 2010
// Purpose:
// This test file makes data structures that mimick the
// monophoton ntuples. And shows how to instantiate the
// NonCollisionBG object and use it's members.
//////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TSystem.h"
#include "NonCollisionBG.C"
#include <iostream>
#include "TMath.h"
#include <TH2.h>
#include "TFile.h"


using namespace std;



void test(){
    NonCollisionBG *bg=new NonCollisionBG();
    //==================================initializing fake data to look like monophoton ntuples

    int HERecHit_subset_n = 3;
    float HERecHit_subset_x[3]={120.,0.,sqrt(120.)};
    float HERecHit_subset_y[3]={0.,120.,sqrt(120.)};
    float HERecHit_subset_energy[3]={2,3.,4.};
    float HERecHit_subset_time[3]={-1.,1.,-1.};
    float Photon_sc_x[3] = {60., 59., 10.};
    float Photon_sc_y[3] = {30.,-50.,-40.};
    float Photon_sc_z[3]= {45., -23., 30.};

    
    int CosmicMuon_n = 3;
    float CosmicMuon_OuterTrack_InnerPoint_x[4]={0.,   120. , sqrt(120.)};
    float CosmicMuon_OuterTrack_InnerPoint_y[4]={120., 0.   , sqrt(120.)};
    float CosmicMuon_OuterTrack_InnerPoint_z[4]={50.0, 25.0, 45.0};
    float CosmicMuon_OuterTrack_InnerPoint_px[4]={25.0, 0.0, 120.0};
    float CosmicMuon_OuterTrack_InnerPoint_py[4]={45.0,5.0,45.0};
    float CosmicMuon_OuterTrack_InnerPoint_pz[4]= {7.0,45.0,0.0};
    int Photon_ncrys[1]={9};
    int Photon_ieta_xtalEB[1][9];
    int Photon_iphi_xtalEB[1][9];
    float Photon_energy_xtal[1][9];
    float Photon_timing_xtal[1][9];
    
    Photon_iphi_xtalEB[0][4] = 359;		Photon_iphi_xtalEB[0][1] = 360;		Photon_iphi_xtalEB[0][2] = 1;
    Photon_ieta_xtalEB[0][4] = -2;		Photon_ieta_xtalEB[0][1] = -2;		Photon_ieta_xtalEB[0][2] = -2;
    Photon_energy_xtal[0][4] = 0.;		Photon_energy_xtal[0][1] = 0.;		Photon_energy_xtal[0][2] = 0.;
    Photon_timing_xtal[0][4] = 0.;		Photon_timing_xtal[0][1] = 0.;		Photon_timing_xtal[0][2] = 0.;
    
    
    Photon_iphi_xtalEB[0][3] = 359;		Photon_iphi_xtalEB[0][0] = 360;		Photon_iphi_xtalEB[0][5] = 1;
    Photon_ieta_xtalEB[0][3] = -1;		Photon_ieta_xtalEB[0][0] = -1;		Photon_ieta_xtalEB[0][5] = -1;
    Photon_energy_xtal[0][3] = 0.;		Photon_energy_xtal[0][0] = 94.;		Photon_energy_xtal[0][5] = 2.;
    Photon_timing_xtal[0][3] = 0.;		Photon_timing_xtal[0][0] = -20.;	Photon_timing_xtal[0][5] = 0.;
    
    
    Photon_iphi_xtalEB[0][6] = 359;		Photon_iphi_xtalEB[0][7] = 360;		Photon_iphi_xtalEB[0][8] = 1;
    Photon_ieta_xtalEB[0][6] = 1;		Photon_ieta_xtalEB[0][7] = 1;		Photon_ieta_xtalEB[0][8] = 1;
    Photon_energy_xtal[0][6] = 4.;		Photon_energy_xtal[0][7] = 0.;		Photon_energy_xtal[0][8] = 0.;
    Photon_timing_xtal[0][6] = 0.;		Photon_timing_xtal[0][7] = 0.;		Photon_timing_xtal[0][8] = 0.;
    //======fake parameters to check cosmics===========
    int pho = 0; // this is the photon index, as if we were looping over a bunch of photons
    float Photon_sc_phi[3]={0.,0.,0.};
    int i = 2;
    float Photon_phHasConversionTracks[3]= {0, 1, 1};
    float Photon_Roundness[3] = {0.9, 0.1, 0.8};
    float Photon_Angle[3]= {0.9,0.3,1.0};
    //for i =0 { its more like signal event, for i = 1, bkg like event}
    
    //===================================== initialized fake data
    
    
    
    // Cosmic functions
    float sigR = bg->sigmaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n);
    cout << "sigR = " << sigR << endl;
    
    float fd2 = bg->fisherCosmic2(Photon_Roundness[i],Photon_Angle[i]);
    cout<<"fd2 = "<<fd2<<endl;
    
    float fd3 = bg->fisherCosmic3(Photon_Roundness[i],Photon_Angle[i],Photon_phHasConversionTracks[i]);
    cout<<"fd3 = "<<fd3<<endl;
    
    // Halo functions
    bool myIsHEHalo = bg->isHEHalo(Photon_sc_phi[pho], HERecHit_subset_n, HERecHit_subset_x, HERecHit_subset_y, HERecHit_subset_energy, HERecHit_subset_time);
    cout<<"myIsHEHalo = "<<myIsHEHalo<<endl;
    
    bool myIsTrackHalo = bg->isTrackHalo(Photon_sc_phi[pho], CosmicMuon_n, CosmicMuon_OuterTrack_InnerPoint_x, CosmicMuon_OuterTrack_InnerPoint_y);
    cout<<"myIsTrackHalo = "<<myIsTrackHalo<<endl;
    
    // Spike functions
    //YOU MUST DO THIS TO USE THESE SPIKE FUNCTIONS
    vector<int> thisPho_ietaRH;
    vector<int> thisPho_iphiRH;
    vector<float> thisPho_eRH;
    vector<float> thisPho_tRH;
    for(int cry=0; cry<Photon_ncrys[pho];cry++){
	thisPho_ietaRH.push_back(Photon_ieta_xtalEB[pho][cry]);
	thisPho_iphiRH.push_back(Photon_iphi_xtalEB[pho][cry]);
	thisPho_eRH.push_back(Photon_energy_xtal[pho][cry]);
	thisPho_tRH.push_back(Photon_timing_xtal[pho][cry]);
    }
    
    float mySimpleBarrelE2E9 = bg->simpleBarrelE2E9(Photon_ncrys[0], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH);
    cout<<"mySimpleBarrelE2E9 = "<<mySimpleBarrelE2E9<<endl;
    
    bool myIsE2E9Spike = bg->isE2E9Spike(Photon_ncrys[0], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH);
    cout<<"myIsE2E9Spike = "<<myIsE2E9Spike<<endl;
    
    bool myIsTimeSpike = bg->isTimeSpike(thisPho_tRH);
    cout<<"myIsTimeSpike = "<<myIsTimeSpike<<endl;

    
    
}
