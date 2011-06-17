//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  9 16:56:27 2011 by ROOT version 5.22/00d
// from TChain myEvent/
//////////////////////////////////////////////////////////

#ifndef makeHistogramsMonoPho_8_h
#define makeHistogramsMonoPho_8_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <vector>
#include "TVector3.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class makeHistogramsMonoPho_8 {
    public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    
    // Declaration of leaf types
   Int_t           nevents;
   UInt_t          run;
   UInt_t          event;
   UInt_t          luminosityBlock;
   UInt_t          beamCrossing;
   UInt_t          totalIntensityBeam1;
   UInt_t          totalIntensityBeam2;
   Float_t         avgInsDelLumi;
   Float_t         avgInsDelLumiErr;
   Float_t         avgInsRecLumi;
   Float_t         avgInsRecLumiErr;
   Int_t           ntriggers;
   vector<string>  *triggernames;
   vector<int>     *triggerprescales;
   vector<bool>    *ifTriggerpassed;
   Float_t         trobjpt[100][100][10];   //[ntriggers]
   Float_t         trobjeta[100][100][10];   //[ntriggers]
   Float_t         trobjphi[100][100][10];   //[ntriggers]
   Int_t           Vertex_n;
   Float_t         Vertex_x[200];   //[Vertex_n]
   Float_t         Vertex_y[200];   //[Vertex_n]
   Float_t         Vertex_z[200];   //[Vertex_n]
   Int_t           Vertex_tracksize[200];   //[Vertex_n]
   Int_t           Vertex_ndof[200];   //[Vertex_n]
   Float_t         Vertex_chi2[200];   //[Vertex_n]
   Float_t         Vertex_d0[200];   //[Vertex_n]
   Bool_t          Vertex_isFake[200];   //[Vertex_n]
   Bool_t          Scraping_isScrapingEvent;
   Int_t           Scraping_numOfTracks;
   Float_t         Scraping_fractionOfGoodTracks;
   Int_t           Track_n;
   Float_t         Track_px[400];   //[Track_n]
   Float_t         Track_py[400];   //[Track_n]
   Float_t         Track_pz[400];   //[Track_n]
   Float_t         Track_vx[400];   //[Track_n]
   Float_t         Track_vy[400];   //[Track_n]
   Float_t         Track_vz[400];   //[Track_n]
   Float_t         Track_pt[400];   //[Track_n]
   Float_t         Track_eta[400];   //[Track_n]
   Float_t         Track_phi[400];   //[Track_n]
   Int_t           Jet_n;
   Float_t         Jet_px[200];   //[Jet_n]
   Float_t         Jet_py[200];   //[Jet_n]
   Float_t         Jet_E[200];   //[Jet_n]
   Float_t         Jet_pz[200];   //[Jet_n]
   Float_t         Jet_vx[200];   //[Jet_n]
   Float_t         Jet_vy[200];   //[Jet_n]
   Float_t         Jet_vz[200];   //[Jet_n]
   Float_t         Jet_pt[200];   //[Jet_n]
   Float_t         Jet_eta[200];   //[Jet_n]
   Float_t         Jet_phi[200];   //[Jet_n]
   Float_t         Jet_emEnergyFraction[200];   //[Jet_n]
   Float_t         Jet_energyFractionHadronic[200];   //[Jet_n]
   Int_t           Jet_hitsInN90[200];   //[Jet_n]
   Int_t           Jet_n90Hits[200];   //[Jet_n]
   Int_t           Jet_nTowers[200];   //[Jet_n]
   Float_t         Jet_fHPD[200];   //[Jet_n]
   Float_t         Jet_fRBX[200];   //[Jet_n]
   Float_t         Jet_RHF[200];   //[Jet_n]
   Float_t         Jet_jecCorr[200];   //[Jet_n]
   Float_t         ucJet_px[200];   //[Jet_n]
   Float_t         ucJet_py[200];   //[Jet_n]
   Float_t         ucJet_E[200];   //[Jet_n]
   Float_t         ucJet_pz[200];   //[Jet_n]
   Float_t         ucJet_pt[200];   //[Jet_n]
   Float_t         ucJet_eta[200];   //[Jet_n]
   Float_t         ucJet_phi[200];   //[Jet_n]
   Int_t           pfJet_n;
   Float_t         pfJet_px[200];   //[pfJet_n]
   Float_t         pfJet_py[200];   //[pfJet_n]
   Float_t         pfJet_E[200];   //[pfJet_n]
   Float_t         pfJet_pz[200];   //[pfJet_n]
   Float_t         pfJet_vx[200];   //[pfJet_n]
   Float_t         pfJet_vy[200];   //[pfJet_n]
   Float_t         pfJet_vz[200];   //[pfJet_n]
   Float_t         pfJet_pt[200];   //[pfJet_n]
   Float_t         pfJet_eta[200];   //[pfJet_n]
   Float_t         pfJet_phi[200];   //[pfJet_n]
   Float_t         pfJet_jecCorr[200];   //[pfJet_n]
   Float_t         ucpfJet_px[200];   //[pfJet_n]
   Float_t         ucpfJet_py[200];   //[pfJet_n]
   Float_t         ucpfJet_E[200];   //[pfJet_n]
   Float_t         ucpfJet_pz[200];   //[pfJet_n]
   Float_t         ucpfJet_pt[200];   //[pfJet_n]
   Float_t         ucpfJet_eta[200];   //[pfJet_n]
   Float_t         ucpfJet_phi[200];   //[pfJet_n]
   Int_t           Electron_n;
   Float_t         Electron_px[200];   //[Electron_n]
   Float_t         Electron_py[200];   //[Electron_n]
   Float_t         Electron_pz[200];   //[Electron_n]
   Float_t         Electron_vx[200];   //[Electron_n]
   Float_t         Electron_vy[200];   //[Electron_n]
   Float_t         Electron_vz[200];   //[Electron_n]
   Float_t         Electron_pt[200];   //[Electron_n]
   Float_t         Electron_eta[200];   //[Electron_n]
   Float_t         Electron_phi[200];   //[Electron_n]
   Float_t         Electron_energy[200];   //[Electron_n]
   Float_t         Electron_charge[200];   //[Electron_n]
   Float_t         Electron_trkIso[200];   //[Electron_n]
   Float_t         Electron_ecalIso[200];   //[Electron_n]
   Float_t         Electron_hcalIso[200];   //[Electron_n]
   Float_t         Electron_SigmaIetaIeta[200];   //[Electron_n]
   Float_t         Electron_dEtaIn[200];   //[Electron_n]
   Float_t         Electron_dPhiIn[200];   //[Electron_n]
   Float_t         Electron_HoE[200];   //[Electron_n]
   Float_t         Electron_sc_energy[200];   //[Electron_n]
   Float_t         Electron_sc_eta[200];   //[Electron_n]
   Float_t         Electron_sc_phi[200];   //[Electron_n]
   Int_t           Muon_n;
   Float_t         Muon_px[200];   //[Muon_n]
   Float_t         Muon_py[200];   //[Muon_n]
   Float_t         Muon_pz[200];   //[Muon_n]
   Float_t         Muon_vx[200];   //[Muon_n]
   Float_t         Muon_vy[200];   //[Muon_n]
   Float_t         Muon_vz[200];   //[Muon_n]
   Float_t         Muon_pt[200];   //[Muon_n]
   Float_t         Muon_eta[200];   //[Muon_n]
   Float_t         Muon_phi[200];   //[Muon_n]
   Float_t         Muon_energy[200];   //[Muon_n]
   Float_t         Muon_charge[200];   //[Muon_n]
   Bool_t          Muon_isGlobalMuon[200];   //[Muon_n]
   Bool_t          Muon_isTrackerMuon[200];   //[Muon_n]
   Bool_t          Muon_isStandAloneMuon[200];   //[Muon_n]
   Bool_t          Muon_InnerTrack_isNonnull[200];   //[Muon_n]
   Bool_t          Muon_OuterTrack_isNonnull[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_x[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_y[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_z[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_px[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_py[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_InnerPoint_pz[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_x[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_y[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_z[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_px[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_py[200];   //[Muon_n]
   Float_t         Muon_OuterTrack_OuterPoint_pz[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_x[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_y[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_z[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_px[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_py[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_InnerPoint_pz[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_x[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_y[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_z[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_px[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_py[200];   //[Muon_n]
   Float_t         Muon_InnerTrack_OuterPoint_pz[200];   //[Muon_n]
   Float_t         Muon_OuterPoint_x[200];   //[Muon_n]
   Float_t         Muon_OuterPoint_y[200];   //[Muon_n]
   Float_t         Muon_OuterPoint_z[200];   //[Muon_n]
   Float_t         Muon_InnerPoint_x[200];   //[Muon_n]
   Float_t         Muon_InnerPoint_y[200];   //[Muon_n]
   Float_t         Muon_InnerPoint_z[200];   //[Muon_n]
   Int_t           CosmicMuon_n;
   Float_t         CosmicMuon_px[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_py[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_pz[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_pt[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_eta[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_phi[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_energy[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_charge[200];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isGlobalMuon[200];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isTrackerMuon[200];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_isStandAloneMuon[200];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_InnerTrack_isNonnull[200];   //[CosmicMuon_n]
   Bool_t          CosmicMuon_OuterTrack_isNonnull[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_x[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_y[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_z[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_px[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_py[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_x[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_y[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_z[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_px[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_py[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterTrack_OuterPoint_pz[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_x[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_y[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_z[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_px[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_py[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_x[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_y[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_z[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_px[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_py[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_InnerTrack_OuterPoint_pz[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterPoint_x[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterPoint_y[200];   //[CosmicMuon_n]
   Float_t         CosmicMuon_OuterPoint_z[200];   //[CosmicMuon_n]
   Int_t           Tau_n;
   Float_t         Tau_px[100];   //[Tau_n]
   Float_t         Tau_py[100];   //[Tau_n]
   Float_t         Tau_pz[100];   //[Tau_n]
   Float_t         Tau_vx[100];   //[Tau_n]
   Float_t         Tau_vy[100];   //[Tau_n]
   Float_t         Tau_vz[100];   //[Tau_n]
   Float_t         Tau_pt[100];   //[Tau_n]
   Float_t         Tau_eta[100];   //[Tau_n]
   Float_t         Tau_phi[100];   //[Tau_n]
   Float_t         Tau_energy[100];   //[Tau_n]
   Float_t         Tau_charge[100];   //[Tau_n]
   vector<string>  *genTauDecayMode1;
   Int_t           oneProng0Pi0[100];   //[Tau_n]
   Int_t           oneProng1Pi0[100];   //[Tau_n]
   Int_t           oneProng2Pi0[100];   //[Tau_n]
   Int_t           threeProng0Pi0[100];   //[Tau_n]
   Int_t           threeProng1Pi0[100];   //[Tau_n]
   Int_t           tauelectron[100];   //[Tau_n]
   Int_t           taumuon[100];   //[Tau_n]
   Int_t           nthreeProng1Pi0;
   Int_t           ntauelectron;
   Int_t           ntaumuon;
   Double_t        genHadTauPt[100];   //[Tau_n]
   Double_t        genHadTauEta[100];   //[Tau_n]
   Double_t        genHadTauPhi[100];   //[Tau_n]
   Int_t           nPions[100];   //[Tau_n]
   Int_t           PionPdgId[100][5];   //[Tau_n]
   Double_t        PionPt[100][5];   //[Tau_n]
   Double_t        PionEta[100][5];   //[Tau_n]
   Double_t        PionPhi[100][5];   //[Tau_n]
   Int_t           nPi0[100];   //[Tau_n]
   Int_t           Pi0PdgId[100][5];   //[Tau_n]
   Double_t        Pi0Pt[100][5];   //[Tau_n]
   Double_t        Pi0Eta[100][5];   //[Tau_n]
   Double_t        Pi0Phi[100][5];   //[Tau_n]
   Int_t           nPhotons[100];   //[Tau_n]
   Double_t        PhotonPt[100][5];   //[Tau_n]
   Double_t        PhotonEta[100][5];   //[Tau_n]
   Double_t        PhotonPhi[100][5];   //[Tau_n]
   Int_t           PhotonPdgId[100][5];   //[Tau_n]
   Int_t           Photon_n;
   Float_t         Photon_E[200];   //[Photon_n]
   Float_t         Photon_pt[200];   //[Photon_n]
   Float_t         Photon_eta[200];   //[Photon_n]
   Float_t         Photon_phi[200];   //[Photon_n]
   Float_t         Photon_theta[200];   //[Photon_n]
   Float_t         Photon_et[200];   //[Photon_n]
   Float_t         Photon_swissCross[200];   //[Photon_n]
   Float_t         Photonr9[200];
   Float_t         Photon_e6e2[200];   //[Photon_n]
   Float_t         Photon_e4e1[200];   //[Photon_n]
   Float_t         Photonr200[200];   //[Photon_n]
   Float_t         Photon_e1x5[200];   //[Photon_n]
   Float_t         Photon_e2x5[200];   //[Photon_n]
   Float_t         Photon_e3x3[200];   //[Photon_n]
   Float_t         Photon_e5x5[200];   //[Photon_n]
   Float_t         Photon_r1x5[200];   //[Photon_n]
   Float_t         Photon_r2x5[200];   //[Photon_n]
   Float_t         Photon_maxEnergyXtal[200];   //[Photon_n]
   Float_t         Photon_SigmaEtaEta[200];   //[Photon_n]
   Float_t         Photon_SigmaIetaIeta[200];   //[Photon_n]
   Float_t         Photon_SigmaEtaPhi[200];   //[Photon_n]
   Float_t         Photon_SigmaIetaIphi[200];   //[Photon_n]
   Float_t         Photon_SigmaPhiPhi[200];   //[Photon_n]
   Float_t         Photon_SigmaIphiIphi[200];   //[Photon_n]
   Float_t         Photon_Roundness[200];   //[Photon_n]
   Float_t         Photon_Angle[200];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR03[200];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR03[200];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR03[200];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR03[200];   //[Photon_n]
   Float_t         Photon_ecalRecHitSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalTowerSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_trkSumPtSolidConeDR04[200];   //[Photon_n]
   Float_t         Photon_trkSumPtHollowConeDR04[200];   //[Photon_n]
   Int_t           Photon_nTrkSolidConeDR04[200];   //[Photon_n]
   Int_t           Photon_nTrkHollowConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalDepth1TowerSumEtConeDR04[200];   //[Photon_n]
   Float_t         Photon_hcalDepth2TowerSumEtConeDR04[200];   //[Photon_n]
   Bool_t          Photon_hasPixelSeed[200];   //[Photon_n]
   Bool_t          Photon_isEB[200];   //[Photon_n]
   Bool_t          Photon_isEE[200];   //[Photon_n]
   Bool_t          Photon_isEBGap[200];   //[Photon_n]
   Bool_t          Photon_isEEGap[200];   //[Photon_n]
   Bool_t          Photon_isEBEEGap[200];   //[Photon_n]
   Float_t         Photon_e2e200[200];   //[Photon_n]
   Float_t         Photon_HoE[200];   //[Photon_n]
   Float_t         Photon_px[200];   //[Photon_n]
   Float_t         Photon_py[200];   //[Photon_n]
   Float_t         Photon_pz[200];   //[Photon_n]
   Float_t         Photon_vx[200];   //[Photon_n]
   Float_t         Photon_vy[200];   //[Photon_n]
   Float_t         Photon_vz[200];   //[Photon_n]
   Int_t           Photon_no_of_basic_clusters[200];   //[Photon_n]
   Float_t         Photon_sc_energy[200];   //[Photon_n]
   Float_t         Photon_sc_eta[200];   //[Photon_n]
   Float_t         Photon_sc_phi[200];   //[Photon_n]
   Float_t         Photon_sc_x[200];   //[Photon_n]
   Float_t         Photon_sc_y[200];   //[Photon_n]
   Float_t         Photon_sc_z[200];   //[Photon_n]
   Float_t         Photon_etaWidth[200];   //[Photon_n]
   Float_t         Photon_phiWidth[200];   //[Photon_n]
   Float_t         Photon_sc_et[200];   //[Photon_n]
   Float_t         matchphotonE[200];   //[Photon_n]
   Float_t         matchphotonpt[200];   //[Photon_n]
   Float_t         matchphotoneta[200];   //[Photon_n]
   Float_t         matchphotonphi[200];   //[Photon_n]
   Float_t         matchphotonpx[200];   //[Photon_n]
   Float_t         matchphotonpy[200];   //[Photon_n]
   Float_t         matchphotonpz[200];   //[Photon_n]
   Bool_t          ismatchedphoton[200];   //[Photon_n]
   Bool_t          Photon_hasConvTrk[200];   //[Photon_n]
   Int_t           Photon_ntracks[200];   //[Photon_n]
   Bool_t          Photon_isconverted[200];   //[Photon_n]
   Float_t         Photon_pairInvmass[200];   //[Photon_n]
   Float_t         Photon_pairCotThetaSeperation[200];   //[Photon_n]
   Float_t         Photon_pairmomentumX[200];   //[Photon_n]
   Float_t         Photon_pairmomentumY[200];   //[Photon_n]
   Float_t         Photon_pairmomentumZ[200];   //[Photon_n]
   Float_t         Photon_EoverP[200];   //[Photon_n]
   Float_t         Photon_ConvVx[200];   //[Photon_n]
   Float_t         Photon_ConvVy[200];   //[Photon_n]
   Float_t         Photon_ConvVz[200];   //[Photon_n]
   Float_t         Photon_ZOfPrimaryVertex[200];   //[Photon_n]
   Float_t         Photon_distOfMinimumApproach[200];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtVtx[200];   //[Photon_n]
   Float_t         Photon_dPhiTracksAtEcal[200];   //[Photon_n]
   Float_t         Photon_dEtaTracksAtEcal[200];   //[Photon_n]
   Int_t           Photon_ncrys[200];   //[Photon_n]
   Float_t         Photon_timing_xtal[200][100];   //[Photon_n]
   Float_t         Photon_timingavg_xtal[200];   //[Photon_n]
   Float_t         Photon_energy_xtal[200][100];   //[Photon_n]
   Int_t           Photon_ieta_xtalEB[200][100];   //[Photon_n]
   Int_t           Photon_iphi_xtalEB[200][100];   //[Photon_n]
   Int_t           Photon_recoFlag_xtalEB[200][100];   //[Photon_n]
   Float_t         Photon_timeError_xtal[200][100];   //[Photon_n]
   Float_t         Photon_s200[200];   //[Photon_n]
   Float_t         Photon_s9[200];
   Int_t           HERecHit_subset_n;
   UInt_t          HERecHit_subset_detid[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_energy[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_time[10000];   //[HERecHit_subset_n]
   Int_t           HERecHit_subset_depth[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_phi[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_eta[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_x[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_y[10000];   //[HERecHit_subset_n]
   Float_t         HERecHit_subset_z[10000];   //[HERecHit_subset_n]
   Int_t           EBRecHit_size;
   Float_t         EBRecHit_eta[10000];   //[EBRecHit_size]
   Float_t         EBRecHit_phi[10000];   //[EBRecHit_size]
   Int_t           EBRecHit_ieta[10000];   //[EBRecHit_size]
   Int_t           EBRecHit_iphi[10000];   //[EBRecHit_size]
   Float_t         EBRecHit_e[10000];   //[EBRecHit_size]
   Float_t         EBRecHit_et[10000];   //[EBRecHit_size]
   Int_t           EBRecHit_flag[10000];   //[EBRecHit_size]
   Float_t         EBRecHit_time[10000];   //[EBRecHit_size]
   Int_t           EERecHit_size;
   Float_t         EERecHit_eta[10000];   //[EERecHit_size]
   Float_t         EERecHit_phi[10000];   //[EERecHit_size]
   Int_t           EERecHit_ieta[10000];   //[EERecHit_size]
   Int_t           EERecHit_iphi[10000];   //[EERecHit_size]
   Float_t         EERecHit_e[10000];   //[EERecHit_size]
   Float_t         EERecHit_et[10000];   //[EERecHit_size]
   Int_t           EERecHit_flag[10000];   //[EERecHit_size]
   Float_t         EERecHit_time[10000];   //[EERecHit_size]
   Bool_t          isBeamHaloGlobalLoosePass;
   Bool_t          isBeamHaloGlobalTightPass;
   Bool_t          isBeamHaloHcalLoosePass;
   Bool_t          isBeamHaloHcalTightPass;
   Bool_t          isBeamHaloCSCLoosePass;
   Bool_t          isBeamHaloCSCTightPass;
   Bool_t          isBeamHaloEcalLoosePass;
   Bool_t          isBeamHaloEcalTightPass;
   Bool_t          isBeamHaloIDTightPass;
   Bool_t          isBeamHaloIDLoosePass;
   Bool_t          isSmellsLikeHalo_Tag;
   Bool_t          isLooseHalo_Tag;
   Bool_t          isTightHalo_Tag;
   Bool_t          isExtremeTightHalo_Tag;
   Float_t         CaloMetSigma;
   Float_t         CaloMetEz;
   Float_t         CaloEtFractionHadronic;
   Float_t         CaloEmEtFraction;
   Float_t         CaloHadEtInHB;
   Float_t         CaloHadEtInHE;
   Float_t         CaloHadEtInHO;
   Float_t         CaloHadEtInHF;
   Float_t         CaloEmEtInEB;
   Float_t         CaloEmEtInEE;
   Float_t         CaloEmEtInHF;
   Float_t         CaloMaxEtInEmTowers;
   Float_t         CaloMaxEtInHadTowers;
   Float_t         CaloMetPt[6];
   Float_t         CaloMetPx[6];
   Float_t         CaloMetPy[6];
   Float_t         CaloMetPhi[6];
   Float_t         CaloMetSumEt[6];
   Float_t         Delta_phi;
   Float_t         PFMetPt[6];
   Float_t         PFMetPx[6];
   Float_t         PFMetPy[6];
   Float_t         PFMetPhi[6];
   Float_t         PFMetSumEt[6];
   Float_t         Delta_phiPF;
   Float_t         TCMetPt[6];
   Float_t         TCMetPx[6];
   Float_t         TCMetPy[6];
   Float_t         TCMetPhi[6];
   Float_t         TCMetSumEt[6];
   Float_t         Delta_phiTC;
   Int_t           ucPhoton_n;
   Float_t         ucPhoton_E[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pt[200];   //[ucPhoton_n]
   Float_t         ucPhoton_eta[200];   //[ucPhoton_n]
   Float_t         ucPhoton_phi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_theta[200];   //[ucPhoton_n]
   Float_t         ucPhoton_et[200];   //[ucPhoton_n]
   Float_t         ucPhoton_swissCross[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e6e2[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e4e1[200];   //[ucPhoton_n]
   Float_t         ucPhotonr200[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e1x5[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e2x5[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e3x3[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e5x5[200];   //[ucPhoton_n]
   Float_t         ucPhoton_r1x5[200];   //[ucPhoton_n]
   Float_t         ucPhoton_r2x5[200];   //[ucPhoton_n]
   Float_t         ucPhoton_maxEnergyXtal[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaEtaEta[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaIetaIeta[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaEtaPhi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaIetaIphi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaPhiPhi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_SigmaIphiIphi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_Roundness[200];   //[ucPhoton_n]
   Float_t         ucPhoton_Angle[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ecalRecHitSumEtConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalTowerSumEtConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_trkSumPtSolidConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_trkSumPtHollowConeDR03[200];   //[ucPhoton_n]
   Int_t           ucPhoton_nTrkSolidConeDR03[200];   //[ucPhoton_n]
   Int_t           ucPhoton_nTrkHollowConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalDepth1TowerSumEtConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalDepth2TowerSumEtConeDR03[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ecalRecHitSumEtConeDR04[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalTowerSumEtConeDR04[200];   //[ucPhoton_n]
   Float_t         ucPhoton_trkSumPtSolidConeDR04[200];   //[ucPhoton_n]
   Float_t         ucPhoton_trkSumPtHollowConeDR04[200];   //[ucPhoton_n]
   Int_t           ucPhoton_nTrkSolidConeDR04[200];   //[ucPhoton_n]
   Int_t           ucPhoton_nTrkHollowConeDR04[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalDepth1TowerSumEtConeDR04[200];   //[ucPhoton_n]
   Float_t         ucPhoton_hcalDepth2TowerSumEtConeDR04[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_hasPixelSeed[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isEB[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isEE[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isEBGap[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isEEGap[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isEBEEGap[200];   //[ucPhoton_n]
   Float_t         ucPhoton_e2e200[200];   //[ucPhoton_n]
   Float_t         ucPhoton_HoE[200];   //[ucPhoton_n]
   Float_t         ucPhoton_px[200];   //[ucPhoton_n]
   Float_t         ucPhoton_py[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pz[200];   //[ucPhoton_n]
   Float_t         ucPhoton_vx[200];   //[ucPhoton_n]
   Float_t         ucPhoton_vy[200];   //[ucPhoton_n]
   Float_t         ucPhoton_vz[200];   //[ucPhoton_n]
   Int_t           ucPhoton_no_of_basic_clusters[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_energy[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_eta[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_phi[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_x[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_y[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_z[200];   //[ucPhoton_n]
   Float_t         ucPhoton_etaWidth[200];   //[ucPhoton_n]
   Float_t         ucPhoton_phiWidth[200];   //[ucPhoton_n]
   Float_t         ucPhoton_sc_et[200];   //[ucPhoton_n]
   Float_t         matchucphotonE[200];   //[ucPhoton_n]
   Float_t         matchucphotonpt[200];   //[ucPhoton_n]
   Float_t         matchucphotoneta[200];   //[ucPhoton_n]
   Float_t         matchucphotonphi[200];   //[ucPhoton_n]
   Float_t         matchucphotonpx[200];   //[ucPhoton_n]
   Float_t         matchucphotonpy[200];   //[ucPhoton_n]
   Float_t         matchucphotonpz[200];   //[ucPhoton_n]
   Bool_t          ismatcheducphoton[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_hasConvTrk[200];   //[ucPhoton_n]
   Int_t           ucPhoton_ntracks[200];   //[ucPhoton_n]
   Bool_t          ucPhoton_isconverted[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pairInvmass[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pairCotThetaSeperation[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pairmomentumX[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pairmomentumY[200];   //[ucPhoton_n]
   Float_t         ucPhoton_pairmomentumZ[200];   //[ucPhoton_n]
   Float_t         ucPhoton_EoverP[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ConvVx[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ConvVy[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ConvVz[200];   //[ucPhoton_n]
   Float_t         ucPhoton_ZOfPrimaryVertex[200];   //[ucPhoton_n]
   Float_t         ucPhoton_distOfMinimumApproach[200];   //[ucPhoton_n]
   Float_t         ucPhoton_dPhiTracksAtVtx[200];   //[ucPhoton_n]
   Float_t         ucPhoton_dPhiTracksAtEcal[200];   //[ucPhoton_n]
   Float_t         ucPhoton_dEtaTracksAtEcal[200];   //[ucPhoton_n]
   Int_t           ucPhoton_ncrys[200];   //[ucPhoton_n]
   Float_t         ucPhoton_timing_xtal[200][100];   //[ucPhoton_n]
   Float_t         ucPhoton_timingavg_xtal[200];   //[ucPhoton_n]
   Float_t         ucPhoton_energy_xtal[200][100];   //[ucPhoton_n]
   Int_t           ucPhoton_ieta_xtalEB[200][100];   //[ucPhoton_n]
   Int_t           ucPhoton_iphi_xtalEB[200][100];   //[ucPhoton_n]
   Int_t           ucPhoton_recoFlag_xtalEB[200][100];   //[ucPhoton_n]
   Float_t         ucPhoton_timeError_xtal[200][100];   //[ucPhoton_n]
   Float_t         ucPhoton_s200[200];   //[ucPhoton_n]
   Int_t           ucHERecHit_subset_n;
   UInt_t          ucHERecHit_subset_detid[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_energy[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_time[10000];   //[ucHERecHit_subset_n]
   Int_t           ucHERecHit_subset_depth[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_phi[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_eta[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_x[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_y[10000];   //[ucHERecHit_subset_n]
   Float_t         ucHERecHit_subset_z[10000];   //[ucHERecHit_subset_n]
   Int_t           CaloTower_n;
   Float_t         CaloTower_eta[5000];   //[CaloTower_n]
   Float_t         CaloTower_phi[5000];   //[CaloTower_n]
   Float_t         CaloTower_E[5000];   //[CaloTower_n]
   Float_t         CaloTower_Et[5000];   //[CaloTower_n]
   Float_t         CaloTower_emEnergy[5000];   //[CaloTower_n]
   Float_t         CaloTower_hadEnergy[5000];   //[CaloTower_n]
   Float_t         CaloTower_p[5000];   //[CaloTower_n]
   Float_t         CaloTower_EMEt[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadEt[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadPhi[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadEta[5000];   //[CaloTower_n]
   Float_t         CaloTower_EMPhi[5000];   //[CaloTower_n]
   Float_t         CaloTower_EMEta[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadX[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadY[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadZ[5000];   //[CaloTower_n]
   Float_t         CaloTower_HE_E[5000];   //[CaloTower_n]
   Float_t         CaloTower_HB_E[5000];   //[CaloTower_n]
   Float_t         CaloTower_EMTime[5000];   //[CaloTower_n]
   Float_t         CaloTower_HadTime[5000];   //[CaloTower_n]

   // List of branches
   TBranch        *b_nevents;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_BXNumber;   //!
   TBranch        *b_totalIntensityBeam1;   //!
   TBranch        *b_totalIntensityBeam2;   //!
   TBranch        *b_avgInsDelLumi;   //!
   TBranch        *b_avgInsDelLumiErr;   //!
   TBranch        *b_avgInsRecLumi;   //!
   TBranch        *b_avgInsRecLumiErr;   //!
   TBranch        *b_ntriggers;   //!
   TBranch        *b_triggernames;   //!
   TBranch        *b_triggerprescales;   //!
   TBranch        *b_ifTriggerpassed;   //!
   TBranch        *b_trobjpt;   //!
   TBranch        *b_trobjeta;   //!
   TBranch        *b_trobjphi;   //!
   TBranch        *b_Vertex_n;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_tracksize;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_d0;   //!
   TBranch        *b_Vertex_isFake;   //!
   TBranch        *b_Scraping_isScrapingEvent;   //!
   TBranch        *b_Scraping_numOfTracks;   //!
   TBranch        *b_Scraping_fractionOfGoodTracks;   //!
   TBranch        *b_Track_n;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_Jet_n;   //!
   TBranch        *b_Jet_px;   //!
   TBranch        *b_Jet_py;   //!
   TBranch        *b_Jet_E;   //!
   TBranch        *b_Jet_pz;   //!
   TBranch        *b_Jet_vx;   //!
   TBranch        *b_Jet_vy;   //!
   TBranch        *b_Jet_vz;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_emEnergyFraction;   //!
   TBranch        *b_Jet_energyFractionHadronic;   //!
   TBranch        *b_Jet_hitsInN90;   //!
   TBranch        *b_Jet_n90Hits;   //!
   TBranch        *b_Jet_nTowers;   //!
   TBranch        *b_Jet_fHPD;   //!
   TBranch        *b_Jet_fRBX;   //!
   TBranch        *b_Jet_RHF;   //!
   TBranch        *b_Jet_jecCorr;   //!
   TBranch        *b_ucJet_px;   //!
   TBranch        *b_ucJet_py;   //!
   TBranch        *b_ucJet_E;   //!
   TBranch        *b_ucJet_pz;   //!
   TBranch        *b_ucJet_pt;   //!
   TBranch        *b_ucJet_eta;   //!
   TBranch        *b_ucJet_phi;   //!
   TBranch        *b_pfJet_n;   //!
   TBranch        *b_pfJet_px;   //!
   TBranch        *b_pfJet_py;   //!
   TBranch        *b_pfJet_E;   //!
   TBranch        *b_pfJet_pz;   //!
   TBranch        *b_pfJet_vx;   //!
   TBranch        *b_pfJet_vy;   //!
   TBranch        *b_pfJet_vz;   //!
   TBranch        *b_pfJet_pt;   //!
   TBranch        *b_pfJet_eta;   //!
   TBranch        *b_pfJet_phi;   //!
   TBranch        *b_pfJet_jecCorr;   //!
   TBranch        *b_ucpfJet_px;   //!
   TBranch        *b_ucpfJet_py;   //!
   TBranch        *b_ucpfJet_E;   //!
   TBranch        *b_ucpfJet_pz;   //!
   TBranch        *b_ucpfJet_pt;   //!
   TBranch        *b_ucpfJet_eta;   //!
   TBranch        *b_ucpfJet_phi;   //!
   TBranch        *b_Electron_n;   //!
   TBranch        *b_Electron_px;   //!
   TBranch        *b_Electron_py;   //!
   TBranch        *b_Electron_pz;   //!
   TBranch        *b_Electron_vx;   //!
   TBranch        *b_Electron_vy;   //!
   TBranch        *b_Electron_vz;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_trkIso;   //!
   TBranch        *b_Electron_ecalIso;   //!
   TBranch        *b_Electron_hcalIso;   //!
   TBranch        *b_Electron_SigmaIetaIeta;   //!
   TBranch        *b_Electron_dEtaIn;   //!
   TBranch        *b_Electron_dPhiIn;   //!
   TBranch        *b_Electron_HoE;   //!
   TBranch        *b_Electron_sc_energy;   //!
   TBranch        *b_Electron_sc_eta;   //!
   TBranch        *b_Electron_sc_phi;   //!
   TBranch        *b_Muon_n;   //!
   TBranch        *b_Muon_px;   //!
   TBranch        *b_Muon_py;   //!
   TBranch        *b_Muon_pz;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_InnerTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterTrack_isNonnull;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_x;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_y;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_z;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_px;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_py;   //!
   TBranch        *b_Muon_OuterTrack_InnerPoint_pz;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_x;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_y;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_z;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_px;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_py;   //!
   TBranch        *b_Muon_OuterTrack_OuterPoint_pz;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_x;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_y;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_z;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_px;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_py;   //!
   TBranch        *b_Muon_InnerTrack_InnerPoint_pz;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_x;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_y;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_z;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_px;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_py;   //!
   TBranch        *b_Muon_InnerTrack_OuterPoint_pz;   //!
   TBranch        *b_Muon_OuterPoint_x;   //!
   TBranch        *b_Muon_OuterPoint_y;   //!
   TBranch        *b_Muon_OuterPoint_z;   //!
   TBranch        *b_Muon_InnerPoint_x;   //!
   TBranch        *b_Muon_InnerPoint_y;   //!
   TBranch        *b_Muon_InnerPoint_z;   //!
   TBranch        *b_CosmicMuon_n;   //!
   TBranch        *b_CosmicMuon_px;   //!
   TBranch        *b_CosmicMuon_py;   //!
   TBranch        *b_CosmicMuon_pz;   //!
   TBranch        *b_CosmicMuon_pt;   //!
   TBranch        *b_CosmicMuon_eta;   //!
   TBranch        *b_CosmicMuon_phi;   //!
   TBranch        *b_CosmicMuon_energy;   //!
   TBranch        *b_CosmicMuon_charge;   //!
   TBranch        *b_CosmicMuon_isGlobalMuon;   //!
   TBranch        *b_CosmicMuon_isTrackerMuon;   //!
   TBranch        *b_CosmicMuon_isStandAloneMuon;   //!
   TBranch        *b_CosmicMuon_InnerTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterTrack_isNonnull;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_z;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_px;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_py;   //!
   TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_pz;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_z;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_px;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_py;   //!
   TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_pz;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_x;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_y;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_z;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_px;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_py;   //!
   TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_pz;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_z;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_px;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_py;   //!
   TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_pz;   //!
   TBranch        *b_CosmicMuon_OuterPoint_x;   //!
   TBranch        *b_CosmicMuon_OuterPoint_y;   //!
   TBranch        *b_CosmicMuon_OuterPoint_z;   //!
   TBranch        *b_Tau_n;   //!
   TBranch        *b_Tau_px;   //!
   TBranch        *b_Tau_py;   //!
   TBranch        *b_Tau_pz;   //!
   TBranch        *b_Tau_vx;   //!
   TBranch        *b_Tau_vy;   //!
   TBranch        *b_Tau_vz;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_genTauDecayMode1;   //!
   TBranch        *b_oneProng0Pi0;   //!
   TBranch        *b_oneProng1Pi0;   //!
   TBranch        *b_oneProng2Pi0;   //!
   TBranch        *b_threeProng0Pi0;   //!
   TBranch        *b_threeProng1Pi0;   //!
   TBranch        *b_tauelectron;   //!
   TBranch        *b_taumuon;   //!
   TBranch        *b_nthreeProng1Pi0;   //!
   TBranch        *b_ntauelectron;   //!
   TBranch        *b_ntaumuon;   //!
   TBranch        *b_genHadTauPt;   //!
   TBranch        *b_genHadTauEta;   //!
   TBranch        *b_genHadTauPhi;   //!
   TBranch        *b_nPions;   //!
   TBranch        *b_PionPdgId;   //!
   TBranch        *b_PionPt;   //!
   TBranch        *b_PionEta;   //!
   TBranch        *b_PionPhi;   //!
   TBranch        *b_nPi0;   //!
   TBranch        *b_Pi0PdgId;   //!
   TBranch        *b_Pi0Pt;   //!
   TBranch        *b_Pi0Eta;   //!
   TBranch        *b_Pi0Phi;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonPdgId;   //!
   TBranch        *b_Photon_n;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_theta;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_swissCross;   //!
   TBranch        *b_Photon_e6e2;   //!
   TBranch        *b_Photon_e4e1;   //!
   TBranch        *b_Photonr9;   //!
   TBranch        *b_Photon_e1x5;   //!
   TBranch        *b_Photon_e2x5;   //!
   TBranch        *b_Photon_e3x3;   //!
   TBranch        *b_Photon_e5x5;   //!
   TBranch        *b_Photon_r1x5;   //!
   TBranch        *b_Photon_r2x5;   //!
   TBranch        *b_Photon_maxEnergyXtal;   //!
   TBranch        *b_Photon_SigmaEtaEta;   //!
   TBranch        *b_Photon_SigmaIetaIeta;   //!
   TBranch        *b_Photon_SigmaEtaPhi;   //!
   TBranch        *b_Photon_SigmaIetaIphi;   //!
   TBranch        *b_Photon_SigmaPhiPhi;   //!
   TBranch        *b_Photon_SigmaIphiIphi;   //!
   TBranch        *b_Photon_Roundness;   //!
   TBranch        *b_Photon_Angle;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_nTrkSolidConeDR03;   //!
   TBranch        *b_Photon_nTrkHollowConeDR03;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_nTrkSolidConeDR04;   //!
   TBranch        *b_Photon_nTrkHollowConeDR04;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_Photon_isEB;   //!
   TBranch        *b_Photon_isEE;   //!
   TBranch        *b_Photon_isEBGap;   //!
   TBranch        *b_Photon_isEEGap;   //!
   TBranch        *b_Photon_isEBEEGap;   //!
   TBranch        *b_Photon_e2e9;   //!
   TBranch        *b_Photon_HoE;   //!
   TBranch        *b_Photon_px;   //!
   TBranch        *b_Photon_py;   //!
   TBranch        *b_Photon_pz;   //!
   TBranch        *b_Photon_vx;   //!
   TBranch        *b_Photon_vy;   //!
   TBranch        *b_Photon_vz;   //!
   TBranch        *b_Photon_no_of_basic_clusters;   //!
   TBranch        *b_Photon_sc_energy;   //!
   TBranch        *b_Photon_sc_eta;   //!
   TBranch        *b_Photon_sc_phi;   //!
   TBranch        *b_Photon_sc_x;   //!
   TBranch        *b_Photon_sc_y;   //!
   TBranch        *b_Photon_sc_z;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_sc_et;   //!
   TBranch        *b_matchphotonE;   //!
   TBranch        *b_matchphotonpt;   //!
   TBranch        *b_matchphotoneta;   //!
   TBranch        *b_matchphotonphi;   //!
   TBranch        *b_matchphotonpx;   //!
   TBranch        *b_matchphotonpy;   //!
   TBranch        *b_matchphotonpz;   //!
   TBranch        *b_ismatchedphoton;   //!
   TBranch        *b_Photon_hasConvTrk;   //!
   TBranch        *b_Photon_ntracks;   //!
   TBranch        *b_Photon_isconverted;   //!
   TBranch        *b_Photon_pairInvmass;   //!
   TBranch        *b_Photon_pairCotThetaSeperation;   //!
   TBranch        *b_Photon_pairmomentumX;   //!
   TBranch        *b_Photon_pairmomentumY;   //!
   TBranch        *b_Photon_pairmomentumZ;   //!
   TBranch        *b_Photon_EoverP;   //!
   TBranch        *b_Photon_ConvVx;   //!
   TBranch        *b_Photon_ConvVy;   //!
   TBranch        *b_Photon_ConvVz;   //!
   TBranch        *b_Photon_ZOfPrimaryVertex;   //!
   TBranch        *b_Photon_distOfMinimumApproach;   //!
   TBranch        *b_Photon_dPhiTracksAtVtx;   //!
   TBranch        *b_Photon_dPhiTracksAtEcal;   //!
   TBranch        *b_Photon_dEtaTracksAtEcal;   //!
   TBranch        *b_Photon_ncrys;   //!
   TBranch        *b_Photon_timing_xtal;   //!
   TBranch        *b_Photon_timingavg_xtal;   //!
   TBranch        *b_Photon_energy_xtal;   //!
   TBranch        *b_Photon_ieta_xtalEB;   //!
   TBranch        *b_Photon_iphi_xtalEB;   //!
   TBranch        *b_Photon_recoFlag_xtalEB;   //!
   TBranch        *b_Photon_timeError_xtal;   //!
   TBranch        *b_Photon_s9;   //!
   TBranch        *b_HERecHit_subset_n;   //!
   TBranch        *b_HERecHit_subset_detid;   //!
   TBranch        *b_HERecHit_subset_energy;   //!
   TBranch        *b_HERecHit_subset_time;   //!
   TBranch        *b_HERecHit_subset_depth;   //!
   TBranch        *b_HERecHit_subset_phi;   //!
   TBranch        *b_HERecHit_subset_eta;   //!
   TBranch        *b_HERecHit_subset_x;   //!
   TBranch        *b_HERecHit_subset_y;   //!
   TBranch        *b_HERecHit_subset_z;   //!
   TBranch        *b_EBRecHit_size;   //!
   TBranch        *b_EBRecHit_eta;   //!
   TBranch        *b_EBRecHit_phi;   //!
   TBranch        *b_EBRecHit_ieta;   //!
   TBranch        *b_EBRecHit_iphi;   //!
   TBranch        *b_EBRecHit_e;   //!
   TBranch        *b_EBRecHit_et;   //!
   TBranch        *b_EBRecHit_flag;   //!
   TBranch        *b_EBRecHit_time;   //!
   TBranch        *b_EERecHit_size;   //!
   TBranch        *b_EERecHit_eta;   //!
   TBranch        *b_EERecHit_phi;   //!
   TBranch        *b_EERecHit_ieta;   //!
   TBranch        *b_EERecHit_iphi;   //!
   TBranch        *b_EERecHit_e;   //!
   TBranch        *b_EERecHit_et;   //!
   TBranch        *b_EERecHit_flag;   //!
   TBranch        *b_EERecHit_time;   //!
   TBranch        *b_isBeamHaloGlobalLoosePass;   //!
   TBranch        *b_isBeamHaloGloablTightPass;   //!
   TBranch        *b_isBeamHaloHcalLoosePass;   //!
   TBranch        *b_isBeamHaloHcalTightPass;   //!
   TBranch        *b_isBeamHaloCSCLoosePass;   //!
   TBranch        *b_isBeamHaloCSCTightPass;   //!
   TBranch        *b_isBeamHaloEcalLoosePass;   //!
   TBranch        *b_isBeamHaloEcalTightPass;   //!
   TBranch        *b_isBeamHaloIDTightPass;   //!
   TBranch        *b_isBeamHaloIDLoosePass;   //!
   TBranch        *b_isSmellsLikeHalo_Tag;   //!
   TBranch        *b_isLooseHalo_Tag;   //!
   TBranch        *b_isTightHalo_Tag;   //!
   TBranch        *b_isExtremeTightHalo_Tag;   //!
   TBranch        *b_CaloMetSig;   //!
   TBranch        *b_CaloMetEz;   //!
   TBranch        *b_CaloEtFractionHadronic;   //!
   TBranch        *b_CaloEmEtFraction;   //!
   TBranch        *b_CaloHadEtInHB;   //!
   TBranch        *b_CaloHadEtInHE;   //!
   TBranch        *b_CaloHadEtInHO;   //!
   TBranch        *b_CaloHadEtInHF;   //!
   TBranch        *b_CaloEmEtInEB;   //!
   TBranch        *b_CaloEmEtInEE;   //!
   TBranch        *b_CaloEmEtInHF;   //!
   TBranch        *b_CaloMaxEtInEmTowers;   //!
   TBranch        *b_CaloMaxEtInHadTowers;   //!
   TBranch        *b_CaloMetPt;   //!
   TBranch        *b_CaloMetPx;   //!
   TBranch        *b_CaloMetPy;   //!
   TBranch        *b_CaloMetPhi;   //!
   TBranch        *b_CaloMetSumEt;   //!
   TBranch        *b_Delta_phi;   //!
   TBranch        *b_PFMetPt;   //!
   TBranch        *b_PFMetPx;   //!
   TBranch        *b_PFMetPy;   //!
   TBranch        *b_PFMetPhi;   //!
   TBranch        *b_PFMetSumEt;   //!
   TBranch        *b_Delta_phiPF;   //!
   TBranch        *b_TCMetPt;   //!
   TBranch        *b_TCMetPx;   //!
   TBranch        *b_TCMetPy;   //!
   TBranch        *b_TCMetPhi;   //!
   TBranch        *b_TCMetSumEt;   //!
   TBranch        *b_Delta_phiTC;   //!
   TBranch        *b_ucPhoton_n;   //!
   TBranch        *b_ucPhoton_E;   //!
   TBranch        *b_ucPhoton_pt;   //!
   TBranch        *b_ucPhoton_eta;   //!
   TBranch        *b_ucPhoton_phi;   //!
   TBranch        *b_ucPhoton_theta;   //!
   TBranch        *b_ucPhoton_et;   //!
   TBranch        *b_ucPhoton_swissCross;   //!
   TBranch        *b_ucPhoton_e6e2;   //!
   TBranch        *b_ucPhoton_e4e1;   //!
   TBranch        *b_ucPhotonr9;   //!
   TBranch        *b_ucPhoton_e1x5;   //!
   TBranch        *b_ucPhoton_e2x5;   //!
   TBranch        *b_ucPhoton_e3x3;   //!
   TBranch        *b_ucPhoton_e5x5;   //!
   TBranch        *b_ucPhoton_r1x5;   //!
   TBranch        *b_ucPhoton_r2x5;   //!
   TBranch        *b_ucPhoton_maxEnergyXtal;   //!
   TBranch        *b_ucPhoton_SigmaEtaEta;   //!
   TBranch        *b_ucPhoton_SigmaIetaIeta;   //!
   TBranch        *b_ucPhoton_SigmaEtaPhi;   //!
   TBranch        *b_ucPhoton_SigmaIetaIphi;   //!
   TBranch        *b_ucPhoton_SigmaPhiPhi;   //!
   TBranch        *b_ucPhoton_SigmaIphiIphi;   //!
   TBranch        *b_ucPhoton_Roundness;   //!
   TBranch        *b_ucPhoton_Angle;   //!
   TBranch        *b_ucPhoton_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_ucPhoton_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_ucPhoton_trkSumPtSolidConeDR03;   //!
   TBranch        *b_ucPhoton_trkSumPtHollowConeDR03;   //!
   TBranch        *b_ucPhoton_nTrkSolidConeDR03;   //!
   TBranch        *b_ucPhoton_nTrkHollowConeDR03;   //!
   TBranch        *b_ucPhoton_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_ucPhoton_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_ucPhoton_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_ucPhoton_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_ucPhoton_trkSumPtSolidConeDR04;   //!
   TBranch        *b_ucPhoton_trkSumPtHollowConeDR04;   //!
   TBranch        *b_ucPhoton_nTrkSolidConeDR04;   //!
   TBranch        *b_ucPhoton_nTrkHollowConeDR04;   //!
   TBranch        *b_ucPhoton_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_ucPhoton_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_ucPhoton_hasPixelSeed;   //!
   TBranch        *b_ucPhoton_isEB;   //!
   TBranch        *b_ucPhoton_isEE;   //!
   TBranch        *b_ucPhoton_isEBGap;   //!
   TBranch        *b_ucPhoton_isEEGap;   //!
   TBranch        *b_ucPhoton_isEBEEGap;   //!
   TBranch        *b_ucPhoton_e2e9;   //!
   TBranch        *b_ucPhoton_HoE;   //!
   TBranch        *b_ucPhoton_px;   //!
   TBranch        *b_ucPhoton_py;   //!
   TBranch        *b_ucPhoton_pz;   //!
   TBranch        *b_ucPhoton_vx;   //!
   TBranch        *b_ucPhoton_vy;   //!
   TBranch        *b_ucPhoton_vz;   //!
   TBranch        *b_ucPhoton_no_of_basic_clusters;   //!
   TBranch        *b_ucPhoton_sc_energy;   //!
   TBranch        *b_ucPhoton_sc_eta;   //!
   TBranch        *b_ucPhoton_sc_phi;   //!
   TBranch        *b_ucPhoton_sc_x;   //!
   TBranch        *b_ucPhoton_sc_y;   //!
   TBranch        *b_ucPhoton_sc_z;   //!
   TBranch        *b_ucPhoton_etaWidth;   //!
   TBranch        *b_ucPhoton_phiWidth;   //!
   TBranch        *b_ucPhoton_sc_et;   //!
   TBranch        *b_matchucphotonE;   //!
   TBranch        *b_matchucphotonpt;   //!
   TBranch        *b_matchucphotoneta;   //!
   TBranch        *b_matchucphotonphi;   //!
   TBranch        *b_matchucphotonpx;   //!
   TBranch        *b_matchucphotonpy;   //!
   TBranch        *b_matchucphotonpz;   //!
   TBranch        *b_ismatcheducphoton;   //!
   TBranch        *b_ucPhoton_hasConvTrk;   //!
   TBranch        *b_ucPhoton_ntracks;   //!
   TBranch        *b_ucPhoton_isconverted;   //!
   TBranch        *b_ucPhoton_pairInvmass;   //!
   TBranch        *b_ucPhoton_pairCotThetaSeperation;   //!
   TBranch        *b_ucPhoton_pairmomentumX;   //!
   TBranch        *b_ucPhoton_pairmomentumY;   //!
   TBranch        *b_ucPhoton_pairmomentumZ;   //!
   TBranch        *b_ucPhoton_EoverP;   //!
   TBranch        *b_ucPhoton_ConvVx;   //!
   TBranch        *b_ucPhoton_ConvVy;   //!
   TBranch        *b_ucPhoton_ConvVz;   //!
   TBranch        *b_ucPhoton_ZOfPrimaryVertex;   //!
   TBranch        *b_ucPhoton_distOfMinimumApproach;   //!
   TBranch        *b_ucPhoton_dPhiTracksAtVtx;   //!
   TBranch        *b_ucPhoton_dPhiTracksAtEcal;   //!
   TBranch        *b_ucPhoton_dEtaTracksAtEcal;   //!
   TBranch        *b_ucPhoton_ncrys;   //!
   TBranch        *b_ucPhoton_timing_xtal;   //!
   TBranch        *b_ucPhoton_timingavg_xtal;   //!
   TBranch        *b_ucPhoton_energy_xtal;   //!
   TBranch        *b_ucPhoton_ieta_xtalEB;   //!
   TBranch        *b_ucPhoton_iphi_xtalEB;   //!
   TBranch        *b_ucPhoton_recoFlag_xtalEB;   //!
   TBranch        *b_ucPhoton_timeError_xtal;   //!
   TBranch        *b_ucPhoton_s9;   //!
   TBranch        *b_ucHERecHit_subset_n;   //!
   TBranch        *b_ucHERecHit_subset_detid;   //!
   TBranch        *b_ucHERecHit_subset_energy;   //!
   TBranch        *b_ucHERecHit_subset_time;   //!
   TBranch        *b_ucHERecHit_subset_depth;   //!
   TBranch        *b_ucHERecHit_subset_phi;   //!
   TBranch        *b_ucHERecHit_subset_eta;   //!
   TBranch        *b_ucHERecHit_subset_x;   //!
   TBranch        *b_ucHERecHit_subset_y;   //!
   TBranch        *b_ucHERecHit_subset_z;   //!
   TBranch        *b_CaloTower_n;   //!
   TBranch        *b_CaloTower_eta;   //!
   TBranch        *b_CaloTower_phi;   //!
   TBranch        *b_CaloTower_E;   //!
   TBranch        *b_CaloTower_Et;   //!
   TBranch        *b_CaloTower_emEnergy;   //!
   TBranch        *b_CaloTower_hadEnergy;   //!
   TBranch        *b_CaloTower_p;   //!
   TBranch        *b_CaloTower_EMEt;   //!
   TBranch        *b_CaloTower_HadEt;   //!
   TBranch        *b_CaloTower_HadPhi;   //!
   TBranch        *b_CaloTower_HadEta;   //!
   TBranch        *b_CaloTower_EMPhi;   //!
   TBranch        *b_CaloTower_EMEta;   //!
   TBranch        *b_CaloTower_HadX;   //!
   TBranch        *b_CaloTower_HadY;   //!
   TBranch        *b_CaloTower_HadZ;   //!
   TBranch        *b_CaloTower_HE_E;   //!
   TBranch        *b_CaloTower_HB_E;   //!
   TBranch        *b_CaloTower_EMTime;   //!
   TBranch        *b_CaloTower_HadTime;   //!
    
    makeHistogramsMonoPho_8(TTree *tree=0,bool applyScraping=true);
    virtual ~makeHistogramsMonoPho_8();
    //virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    Cut();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree,bool applyScraping);
    virtual void     Loop(char subset, string outName, bool applyTriggerCut=true, bool applyCosMuCut=true, bool applyScraping=true, bool onlyOneVertex=false);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    float deltaPhi(float phi1, float phi2);
    float deltaR(float eta1, float phi1, float eta2, float phi2);
    float fixPhi(float phi);
    bool eventHasGoodVertex(int v_n, int v_ndof[], float v_z[], float v_d0[], bool v_isFake[], bool onlyOne=false);
    bool passesPhotonID(float ncrys,float pt, float sc_eta, float HoE, float trkSumPtHollowConeDR04,
                        float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float sigmaIEIE);
    float invMass(float E, float px, float py, float pz, float E2, float px2, float py2, float pz2);
    float massT(float phi, float et, float metphi, float met);
    bool eventHasMu(int nCosMu, bool CosOuterTrack_isNonnull[], bool CosTrack_isSA[]);
    bool phoHasMu(int nMu, bool MuOuterTrack_isNonnull[], bool MuTrack_isSA[],float MuTrack_Pt[], float MuTrack_Eta[],
                    float MuTrack_Phi[], float MuTrack_Vx[], float MuTrack_Vy[], float MuTrack_Vz[], float phoPhi,
                    float phoEta, float phoVx, float phoVy, float phoVz, float muMinEt=-10., bool restrictEta=false, 
                    float maxEta=1.4442, bool restrictVertex=false, float maxVertexDist=10000.);
    float computeLICTD(int nCry, float cryTime[], float cryEnergy[]);
    int GetSeed(int Ncry, float cryE[]);
    float GetE2overE9(int Ncry, int ieta[], int iphi[], float e[]);
    bool triggerPassed(unsigned int Ntriggers, vector <string> myDesiredTrigger, vector <string> storedTriggers_name, vector <bool> storedTriggers_bool);
    float newMetPhi(float MetPt, float MetPhi, float PhoPt, float PhoPhi, float seedT);
    float newMetPt(float MetPt, float MetPhi, float PhoPt, float PhoPhi, float seedT);
    
    //Halo functions
    float absDeltaPhi(float phi1, float phi2);
    bool isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime=false);
    bool isTrackHalo(float photonSCphi, int nCosMu, bool CosOuterTrack_isNonnull[], float CosTrackX[], float CosTracksY[]);
    bool isCSCHalo(float photonSCphi, int nAllCSCSegments, float CSCSegmentX[], float CSCSegmentY[], float CSCSegmentTime[], bool useTime=false);
    float deltaPhiCSCHalo(float photonSCphi,float CSCSegmentX,float CSCSegmentY);
    
    //Cosmic functions
    //TFile* file_COSMIC;
    //TFile* file_PHOTON;
    float sigmaRCosmic(float Photon_sc_x, float Photon_sc_y, float Photon_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n);
    float deltaRCosmic(float Photon_sc_x, float Photon_sc_y, float Photon_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n);
    //float fisherCosmic3(float Photon_Roundness, float Photon_Angle, bool Photon_hasConvTrk);
    //float fisherCosmic2(float Photon_Roundness, float Photon_Angle);
    vector<pair <bool,TVector3> > getEBIntersections( const TVector3& X,  const TVector3& P);
    double h(double z, double& length);
    double OuterEBRho(double zprimeprime, int iterations = 5 );
    double tDCA(const TVector3& X, const TVector3& P);
    vector <double> outerIntersectionTs(const TVector3& X, const TVector3& P, double timeDCA);
    vector <double> innerIntersectionTs(double R, const TVector3& X, const TVector3& P);
    //float prob(TH1F* h, float value);
    //float prob2(TH1F* h, bool Photon_hasConvTrk);
    bool passesTightPhotonID(float pt, float HoE, float trkSumPtHollowConeDR04,
        float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float sigmaIEIE, float eta,
        float minEt=-10., bool restrictEta=false, float maxEta=1.4442);
    bool passesLoosePhotonID(float pt, float HoE, float trkSumPtHollowConeDR04, float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float eta, float minEt=-10., bool restrictEta=false, float maxEta=1.4442);
    bool passesEMID(float pt, float HoE, float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float eta, float minEt=-10., bool restrictEta=false, float maxEta=1.4442);
    bool hasObj(float phoPhi, float phoEta, float phoVx, float phoVy, float phoVz, int objN, float objEt[], float objEta[], float objPhi[], float objV[], float objVy[], float objVz[], float objMinEt=-10., bool restrictEta=false, float maxEta=1.4442, bool restrictVertex=false, float maxVertexDist=10000.);
    bool hasOtherPho(int phoID, int candIdx, int PhoN, float et[],float HoE[],float trkSumPtHollowConeDR04[],
        float ecalRecHitSumEtConeDR04[], float hcalTowerSumEtConeDR04[], float sigmaIEIE[], float eta[], float vx[], float vy[], float vz[],
        float minEt=-10., bool restrictEta=false, float maxEta=1.4442, bool restrictVertex=false, float maxVertexDist=10000.);
};

#endif

#ifdef makeHistogramsMonoPho_8_cxx
makeHistogramsMonoPho_8::makeHistogramsMonoPho_8(TTree *tree,bool applyScraping)
{
    //file_PHOTON = TFile::Open("photon_MC.root","READ");
    //file_COSMIC = TFile::Open("data_cosmics.root","READ");
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        
#ifdef SINGLE_TREE
        // The following code should be used if you want this class to access
        // a single tree instead of a chain
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
        if (!f) {
            f = new TFile("Memory Directory");
            f->cd("Rint:/");
        }
        tree = (TTree*)gDirectory->Get("myEvent");
        
#else // SINGLE_TREE
        
        // The following code should be used if you want this class to access a chain
        // of trees.
        TChain * chain = new TChain("myEvent","");
        chain->Add("/uscms/home/miceli/scratch/monophoton/subsets/F/fakeCand_skimOfNtuples_runB_6.root/myEvent");
        tree = chain;
#endif // SINGLE_TREE
        
    }
    Init(tree,applyScraping);
}

makeHistogramsMonoPho_8::~makeHistogramsMonoPho_8()
{
//    file_PHOTON->Close();
//    file_COSMIC->Close();
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t makeHistogramsMonoPho_8::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t makeHistogramsMonoPho_8::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
    TChain *chain = (TChain*)fChain;
    if (chain->GetTreeNumber() != fCurrent) {
        fCurrent = chain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void makeHistogramsMonoPho_8::Init(TTree *tree,bool applyScraping)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set object pointer
    triggernames = 0;
    triggerprescales = 0;
    ifTriggerpassed = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("nevents", &nevents, &b_nevents);
    fChain->SetBranchAddress("run", &run, &b_RunNumber);
    fChain->SetBranchAddress("event", &event, &b_EventNumber);
    fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_LumiNumber);
    fChain->SetBranchAddress("beamCrossing", &beamCrossing, &b_BXNumber);
    fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
    fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
    fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
    fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
    fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
    fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
    fChain->SetBranchAddress("ntriggers", &ntriggers, &b_ntriggers);
    fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
    fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
    fChain->SetBranchAddress("ifTriggerpassed", &ifTriggerpassed, &b_ifTriggerpassed);
    fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
    fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
    fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
    fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
    fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize, &b_Vertex_tracksize);
    fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof, &b_Vertex_ndof);
    fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2, &b_Vertex_chi2);
    fChain->SetBranchAddress("Vertex_d0", Vertex_d0, &b_Vertex_d0);
    fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake, &b_Vertex_isFake);
    if(applyScraping){
        fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent, &b_Scraping_isScrapingEvent);
        fChain->SetBranchAddress("Scraping_numOfTracks", &Scraping_numOfTracks, &b_Scraping_numOfTracks);
        fChain->SetBranchAddress("Scraping_fractionOfGoodTracks", &Scraping_fractionOfGoodTracks, &b_Scraping_fractionOfGoodTracks);
    }
    fChain->SetBranchAddress("Track_n", &Track_n, &b_Track_n);
    fChain->SetBranchAddress("Track_px", Track_px, &b_Track_px);
    fChain->SetBranchAddress("Track_py", Track_py, &b_Track_py);
    fChain->SetBranchAddress("Track_pz", Track_pz, &b_Track_pz);
    fChain->SetBranchAddress("Track_vx", Track_vx, &b_Track_vx);
    fChain->SetBranchAddress("Track_vy", Track_vy, &b_Track_vy);
    fChain->SetBranchAddress("Track_vz", Track_vz, &b_Track_vz);
    fChain->SetBranchAddress("Track_pt", Track_pt, &b_Track_pt);
    fChain->SetBranchAddress("Track_eta", Track_eta, &b_Track_eta);
    fChain->SetBranchAddress("Track_phi", Track_phi, &b_Track_phi);
    fChain->SetBranchAddress("Jet_n", &Jet_n, &b_Jet_n);
    fChain->SetBranchAddress("Jet_px", Jet_px, &b_Jet_px);
    fChain->SetBranchAddress("Jet_py", Jet_py, &b_Jet_py);
    fChain->SetBranchAddress("Jet_E", Jet_E, &b_Jet_E);
    fChain->SetBranchAddress("Jet_pz", Jet_pz, &b_Jet_pz);
    fChain->SetBranchAddress("Jet_vx", Jet_vx, &b_Jet_vx);
    fChain->SetBranchAddress("Jet_vy", Jet_vy, &b_Jet_vy);
    fChain->SetBranchAddress("Jet_vz", Jet_vz, &b_Jet_vz);
    fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
    fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
    fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
    fChain->SetBranchAddress("Jet_emEnergyFraction", Jet_emEnergyFraction, &b_Jet_emEnergyFraction);
    fChain->SetBranchAddress("Jet_energyFractionHadronic", Jet_energyFractionHadronic, &b_Jet_energyFractionHadronic);
    fChain->SetBranchAddress("Jet_hitsInN90", Jet_hitsInN90, &b_Jet_hitsInN90);
    fChain->SetBranchAddress("Jet_n90Hits", Jet_n90Hits, &b_Jet_n90Hits);
    fChain->SetBranchAddress("Jet_nTowers", Jet_nTowers, &b_Jet_nTowers);
    fChain->SetBranchAddress("Jet_fHPD", Jet_fHPD, &b_Jet_fHPD);
    fChain->SetBranchAddress("Jet_fRBX", Jet_fRBX, &b_Jet_fRBX);
    fChain->SetBranchAddress("Jet_RHF", Jet_RHF, &b_Jet_RHF);
    fChain->SetBranchAddress("Electron_n", &Electron_n, &b_Electron_n);
    fChain->SetBranchAddress("Electron_px", Electron_px, &b_Electron_px);
    fChain->SetBranchAddress("Electron_py", Electron_py, &b_Electron_py);
    fChain->SetBranchAddress("Electron_pz", Electron_pz, &b_Electron_pz);
    fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
    fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
    fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
    fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
    fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
    fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
    fChain->SetBranchAddress("Electron_energy", Electron_energy, &b_Electron_energy);
    fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
    fChain->SetBranchAddress("Electron_trkIso", Electron_trkIso, &b_Electron_trkIso);
    fChain->SetBranchAddress("Electron_ecalIso", Electron_ecalIso, &b_Electron_ecalIso);
    fChain->SetBranchAddress("Electron_hcalIso", Electron_hcalIso, &b_Electron_hcalIso);
    fChain->SetBranchAddress("Electron_SigmaIetaIeta", Electron_SigmaIetaIeta, &b_Electron_SigmaIetaIeta);
    fChain->SetBranchAddress("Electron_dEtaIn", Electron_dEtaIn, &b_Electron_dEtaIn);
    fChain->SetBranchAddress("Electron_dPhiIn", Electron_dPhiIn, &b_Electron_dPhiIn);
    fChain->SetBranchAddress("Electron_HoE", Electron_HoE, &b_Electron_HoE);
    fChain->SetBranchAddress("Electron_sc_energy", Electron_sc_energy, &b_Electron_sc_energy);
    fChain->SetBranchAddress("Electron_sc_eta", Electron_sc_eta, &b_Electron_sc_eta);
    fChain->SetBranchAddress("Electron_sc_phi", Electron_sc_phi, &b_Electron_sc_phi);
    fChain->SetBranchAddress("Muon_n", &Muon_n, &b_Muon_n);
    fChain->SetBranchAddress("Muon_px", Muon_px, &b_Muon_px);
    fChain->SetBranchAddress("Muon_py", Muon_py, &b_Muon_py);
    fChain->SetBranchAddress("Muon_pz", Muon_pz, &b_Muon_pz);
    fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
    fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
    fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
    fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
    fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
    fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
    fChain->SetBranchAddress("Muon_energy", Muon_energy, &b_Muon_energy);
    fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
    fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
    fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
    fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
    fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull, &b_Muon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull, &b_Muon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_x", Muon_OuterTrack_InnerPoint_x, &b_Muon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_y", Muon_OuterTrack_InnerPoint_y, &b_Muon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_z", Muon_OuterTrack_InnerPoint_z, &b_Muon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_px", Muon_OuterTrack_InnerPoint_px, &b_Muon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_py", Muon_OuterTrack_InnerPoint_py, &b_Muon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_pz", Muon_OuterTrack_InnerPoint_pz, &b_Muon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_x", Muon_OuterTrack_OuterPoint_x, &b_Muon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_y", Muon_OuterTrack_OuterPoint_y, &b_Muon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_z", Muon_OuterTrack_OuterPoint_z, &b_Muon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_px", Muon_OuterTrack_OuterPoint_px, &b_Muon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_py", Muon_OuterTrack_OuterPoint_py, &b_Muon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_pz", Muon_OuterTrack_OuterPoint_pz, &b_Muon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_x", Muon_InnerTrack_InnerPoint_x, &b_Muon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_y", Muon_InnerTrack_InnerPoint_y, &b_Muon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_z", Muon_InnerTrack_InnerPoint_z, &b_Muon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_px", Muon_InnerTrack_InnerPoint_px, &b_Muon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_py", Muon_InnerTrack_InnerPoint_py, &b_Muon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_pz", Muon_InnerTrack_InnerPoint_pz, &b_Muon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_x", Muon_InnerTrack_OuterPoint_x, &b_Muon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_y", Muon_InnerTrack_OuterPoint_y, &b_Muon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_z", Muon_InnerTrack_OuterPoint_z, &b_Muon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_px", Muon_InnerTrack_OuterPoint_px, &b_Muon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_py", Muon_InnerTrack_OuterPoint_py, &b_Muon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_pz", Muon_InnerTrack_OuterPoint_pz, &b_Muon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_n", &CosmicMuon_n, &b_CosmicMuon_n);
    fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px, &b_CosmicMuon_px);
    fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py, &b_CosmicMuon_py);
    fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz, &b_CosmicMuon_pz);
    fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt, &b_CosmicMuon_pt);
    fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta, &b_CosmicMuon_eta);
    fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi, &b_CosmicMuon_phi);
    fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy, &b_CosmicMuon_energy);
    fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge, &b_CosmicMuon_charge);
    fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon, &b_CosmicMuon_isGlobalMuon);
    fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon, &b_CosmicMuon_isTrackerMuon);
    fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon, &b_CosmicMuon_isStandAloneMuon);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull, &b_CosmicMuon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull, &b_CosmicMuon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_x", CosmicMuon_OuterTrack_InnerPoint_x, &b_CosmicMuon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_y", CosmicMuon_OuterTrack_InnerPoint_y, &b_CosmicMuon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_z", CosmicMuon_OuterTrack_InnerPoint_z, &b_CosmicMuon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_px", CosmicMuon_OuterTrack_InnerPoint_px, &b_CosmicMuon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_py", CosmicMuon_OuterTrack_InnerPoint_py, &b_CosmicMuon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_pz", CosmicMuon_OuterTrack_InnerPoint_pz, &b_CosmicMuon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_x", CosmicMuon_OuterTrack_OuterPoint_x, &b_CosmicMuon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_y", CosmicMuon_OuterTrack_OuterPoint_y, &b_CosmicMuon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_z", CosmicMuon_OuterTrack_OuterPoint_z, &b_CosmicMuon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_px", CosmicMuon_OuterTrack_OuterPoint_px, &b_CosmicMuon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_py", CosmicMuon_OuterTrack_OuterPoint_py, &b_CosmicMuon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_pz", CosmicMuon_OuterTrack_OuterPoint_pz, &b_CosmicMuon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_x", CosmicMuon_InnerTrack_InnerPoint_x, &b_CosmicMuon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_y", CosmicMuon_InnerTrack_InnerPoint_y, &b_CosmicMuon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_z", CosmicMuon_InnerTrack_InnerPoint_z, &b_CosmicMuon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_px", CosmicMuon_InnerTrack_InnerPoint_px, &b_CosmicMuon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_py", CosmicMuon_InnerTrack_InnerPoint_py, &b_CosmicMuon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_pz", CosmicMuon_InnerTrack_InnerPoint_pz, &b_CosmicMuon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_x", CosmicMuon_InnerTrack_OuterPoint_x, &b_CosmicMuon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_y", CosmicMuon_InnerTrack_OuterPoint_y, &b_CosmicMuon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_z", CosmicMuon_InnerTrack_OuterPoint_z, &b_CosmicMuon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_px", CosmicMuon_InnerTrack_OuterPoint_px, &b_CosmicMuon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_py", CosmicMuon_InnerTrack_OuterPoint_py, &b_CosmicMuon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_pz", CosmicMuon_InnerTrack_OuterPoint_pz, &b_CosmicMuon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Tau_n", &Tau_n, &b_Tau_n);
    fChain->SetBranchAddress("Tau_px", Tau_px, &b_Tau_px);
    fChain->SetBranchAddress("Tau_py", Tau_py, &b_Tau_py);
    fChain->SetBranchAddress("Tau_pz", Tau_pz, &b_Tau_pz);
    fChain->SetBranchAddress("Tau_vx", Tau_vx, &b_Tau_vx);
    fChain->SetBranchAddress("Tau_vy", Tau_vy, &b_Tau_vy);
    fChain->SetBranchAddress("Tau_vz", Tau_vz, &b_Tau_vz);
    fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
    fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
    fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
    fChain->SetBranchAddress("Tau_energy", Tau_energy, &b_Tau_energy);
    fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
    fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
    fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
    fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
    fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
    fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
    fChain->SetBranchAddress("Photon_theta", Photon_theta, &b_Photon_theta);
    fChain->SetBranchAddress("Photon_et", Photon_et, &b_Photon_et);
    fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross, &b_Photon_swissCross);
    fChain->SetBranchAddress("Photonr9", Photonr9, &b_Photonr9);
    fChain->SetBranchAddress("Photon_e1x5", Photon_e1x5, &b_Photon_e1x5);
    fChain->SetBranchAddress("Photon_e2x5", Photon_e2x5, &b_Photon_e2x5);
    fChain->SetBranchAddress("Photon_e3x3", Photon_e3x3, &b_Photon_e3x3);
    fChain->SetBranchAddress("Photon_e5x5", Photon_e5x5, &b_Photon_e5x5);
    fChain->SetBranchAddress("Photon_r1x5", Photon_r1x5, &b_Photon_r1x5);
    fChain->SetBranchAddress("Photon_r2x5", Photon_r2x5, &b_Photon_r2x5);
    fChain->SetBranchAddress("Photon_maxEnergyXtal", Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
    fChain->SetBranchAddress("Photon_SigmaEtaEta", Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
    fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
    fChain->SetBranchAddress("Photon_SigmaEtaPhi", Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
    fChain->SetBranchAddress("Photon_SigmaIetaIphi", Photon_SigmaIetaIphi, &b_Photon_SigmaIetaIphi);
    fChain->SetBranchAddress("Photon_SigmaPhiPhi", Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
    fChain->SetBranchAddress("Photon_SigmaIphiIphi", Photon_SigmaIphiIphi, &b_Photon_SigmaIphiIphi);
    fChain->SetBranchAddress("Photon_Roundness", Photon_Roundness, &b_Photon_Roundness);
    fChain->SetBranchAddress("Photon_Angle", Photon_Angle, &b_Photon_Angle);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
    fChain->SetBranchAddress("Photon_isEB", Photon_isEB, &b_Photon_isEB);
    fChain->SetBranchAddress("Photon_isEE", Photon_isEE, &b_Photon_isEE);
    fChain->SetBranchAddress("Photon_isEBGap", Photon_isEBGap, &b_Photon_isEBGap);
    fChain->SetBranchAddress("Photon_isEEGap", Photon_isEEGap, &b_Photon_isEEGap);
    fChain->SetBranchAddress("Photon_isEBEEGap", Photon_isEBEEGap, &b_Photon_isEBEEGap);
    fChain->SetBranchAddress("Photon_HoE", Photon_HoE, &b_Photon_HoE);
    fChain->SetBranchAddress("Photon_px", Photon_px, &b_Photon_px);
    fChain->SetBranchAddress("Photon_py", Photon_py, &b_Photon_py);
    fChain->SetBranchAddress("Photon_pz", Photon_pz, &b_Photon_pz);
    fChain->SetBranchAddress("Photon_vx", Photon_vx, &b_Photon_vx);
    fChain->SetBranchAddress("Photon_vy", Photon_vy, &b_Photon_vy);
    fChain->SetBranchAddress("Photon_vz", Photon_vz, &b_Photon_vz);
    fChain->SetBranchAddress("Photon_no_of_basic_clusters", Photon_no_of_basic_clusters, &b_Photon_no_of_basic_clusters);
    fChain->SetBranchAddress("Photon_sc_energy", Photon_sc_energy, &b_Photon_sc_energy);
    fChain->SetBranchAddress("Photon_sc_eta", Photon_sc_eta, &b_Photon_sc_eta);
    fChain->SetBranchAddress("Photon_sc_phi", Photon_sc_phi, &b_Photon_sc_phi);
    fChain->SetBranchAddress("Photon_sc_x", Photon_sc_x, &b_Photon_sc_x);
    fChain->SetBranchAddress("Photon_sc_y", Photon_sc_y, &b_Photon_sc_y);
    fChain->SetBranchAddress("Photon_sc_z", Photon_sc_z, &b_Photon_sc_z);
    fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
    fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
    fChain->SetBranchAddress("Photon_sc_et", Photon_sc_et, &b_Photon_sc_et);
    fChain->SetBranchAddress("matchphotonE", matchphotonE, &b_matchphotonE);
    fChain->SetBranchAddress("matchphotonpt", matchphotonpt, &b_matchphotonpt);
    fChain->SetBranchAddress("matchphotoneta", matchphotoneta, &b_matchphotoneta);
    fChain->SetBranchAddress("matchphotonphi", matchphotonphi, &b_matchphotonphi);
    fChain->SetBranchAddress("matchphotonpx", matchphotonpx, &b_matchphotonpx);
    fChain->SetBranchAddress("matchphotonpy", matchphotonpy, &b_matchphotonpy);
    fChain->SetBranchAddress("matchphotonpz", matchphotonpz, &b_matchphotonpz);
    fChain->SetBranchAddress("ismatchedphoton", ismatchedphoton, &b_ismatchedphoton);
    fChain->SetBranchAddress("Photon_hasConvTrk", Photon_hasConvTrk, &b_Photon_hasConvTrk);
    fChain->SetBranchAddress("Photon_ntracks", Photon_ntracks, &b_Photon_ntracks);
    fChain->SetBranchAddress("Photon_isconverted", Photon_isconverted, &b_Photon_isconverted);
    fChain->SetBranchAddress("Photon_pairInvmass", Photon_pairInvmass, &b_Photon_pairInvmass);
    fChain->SetBranchAddress("Photon_pairCotThetaSeperation", Photon_pairCotThetaSeperation, &b_Photon_pairCotThetaSeperation);
    fChain->SetBranchAddress("Photon_pairmomentumX", Photon_pairmomentumX, &b_Photon_pairmomentumX);
    fChain->SetBranchAddress("Photon_pairmomentumY", Photon_pairmomentumY, &b_Photon_pairmomentumY);
    fChain->SetBranchAddress("Photon_pairmomentumZ", Photon_pairmomentumZ, &b_Photon_pairmomentumZ);
    fChain->SetBranchAddress("Photon_EoverP", Photon_EoverP, &b_Photon_EoverP);
    fChain->SetBranchAddress("Photon_ConvVx", Photon_ConvVx, &b_Photon_ConvVx);
    fChain->SetBranchAddress("Photon_ConvVy", Photon_ConvVy, &b_Photon_ConvVy);
    fChain->SetBranchAddress("Photon_ConvVz", Photon_ConvVz, &b_Photon_ConvVz);
    fChain->SetBranchAddress("Photon_ZOfPrimaryVertex", Photon_ZOfPrimaryVertex, &b_Photon_ZOfPrimaryVertex);
    fChain->SetBranchAddress("Photon_distOfMinimumApproach", Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
    fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
    fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
    fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
    fChain->SetBranchAddress("Photon_ncrys", Photon_ncrys, &b_Photon_ncrys);
    fChain->SetBranchAddress("Photon_timing_xtal", Photon_timing_xtal, &b_Photon_timing_xtal);
    fChain->SetBranchAddress("Photon_timingavg_xtal", Photon_timingavg_xtal, &b_Photon_timingavg_xtal);
    fChain->SetBranchAddress("Photon_energy_xtal", Photon_energy_xtal, &b_Photon_energy_xtal);
    fChain->SetBranchAddress("Photon_ieta_xtalEB", Photon_ieta_xtalEB, &b_Photon_ieta_xtalEB);
    fChain->SetBranchAddress("Photon_iphi_xtalEB", Photon_iphi_xtalEB, &b_Photon_iphi_xtalEB);
    fChain->SetBranchAddress("Photon_s9", Photon_s9, &b_Photon_s9);
    fChain->SetBranchAddress("HERecHit_subset_n", &HERecHit_subset_n, &b_HERecHit_subset_n);
    fChain->SetBranchAddress("HERecHit_subset_detid", HERecHit_subset_detid, &b_HERecHit_subset_detid);
    fChain->SetBranchAddress("HERecHit_subset_energy", HERecHit_subset_energy, &b_HERecHit_subset_energy);
    fChain->SetBranchAddress("HERecHit_subset_time", HERecHit_subset_time, &b_HERecHit_subset_time);
    fChain->SetBranchAddress("HERecHit_subset_depth", HERecHit_subset_depth, &b_HERecHit_subset_depth);
    fChain->SetBranchAddress("HERecHit_subset_phi", HERecHit_subset_phi, &b_HERecHit_subset_phi);
    fChain->SetBranchAddress("HERecHit_subset_eta", HERecHit_subset_eta, &b_HERecHit_subset_eta);
    fChain->SetBranchAddress("HERecHit_subset_x", HERecHit_subset_x, &b_HERecHit_subset_x);
    fChain->SetBranchAddress("HERecHit_subset_y", HERecHit_subset_y, &b_HERecHit_subset_y);
    fChain->SetBranchAddress("HERecHit_subset_z", HERecHit_subset_z, &b_HERecHit_subset_z);
    /*
    fChain->SetBranchAddress("CSCseg_n", &CSCseg_n, &b_CSCseg_n);
    fChain->SetBranchAddress("CSCseg_time", CSCseg_time, &b_CSCseg_time);
    fChain->SetBranchAddress("CSCseg_x", CSCseg_x, &b_CSCseg_x);
    fChain->SetBranchAddress("CSCseg_y", CSCseg_y, &b_CSCseg_y);
    fChain->SetBranchAddress("CSCseg_z", CSCseg_z, &b_CSCseg_z);
    fChain->SetBranchAddress("CSCseg_phi", CSCseg_phi, &b_CSCseg_phi);
    fChain->SetBranchAddress("CSCseg_DirectionX", CSCseg_DirectionX, &b_CSCseg_DirectionX);
    fChain->SetBranchAddress("CSCseg_DirectionY", CSCseg_DirectionY, &b_CSCseg_DirectionY);
    fChain->SetBranchAddress("CSCseg_DirectionZ", CSCseg_DirectionZ, &b_CSCseg_DirectionZ);
    fChain->SetBranchAddress("RPChit_n", &RPChit_n, &b_RPChit_n);
    fChain->SetBranchAddress("RPChit_x", RPChit_x, &b_RPChit_x);
    fChain->SetBranchAddress("RPChit_y", RPChit_y, &b_RPChit_y);
    fChain->SetBranchAddress("RPChit_z", RPChit_z, &b_RPChit_z);
    fChain->SetBranchAddress("RPChit_BunchX", RPChit_BunchX, &b_RPChit_BunchX);
    */
    fChain->SetBranchAddress("CaloMetSigma", &CaloMetSigma, &b_CaloMetSig);
    fChain->SetBranchAddress("CaloMetEz", &CaloMetEz, &b_CaloMetEz);
    fChain->SetBranchAddress("CaloEtFractionHadronic", &CaloEtFractionHadronic, &b_CaloEtFractionHadronic);
    fChain->SetBranchAddress("CaloEmEtFraction", &CaloEmEtFraction, &b_CaloEmEtFraction);
    fChain->SetBranchAddress("CaloHadEtInHB", &CaloHadEtInHB, &b_CaloHadEtInHB);
    fChain->SetBranchAddress("CaloHadEtInHE", &CaloHadEtInHE, &b_CaloHadEtInHE);
    fChain->SetBranchAddress("CaloHadEtInHO", &CaloHadEtInHO, &b_CaloHadEtInHO);
    fChain->SetBranchAddress("CaloHadEtInHF", &CaloHadEtInHF, &b_CaloHadEtInHF);
    fChain->SetBranchAddress("CaloEmEtInEB", &CaloEmEtInEB, &b_CaloEmEtInEB);
    fChain->SetBranchAddress("CaloEmEtInEE", &CaloEmEtInEE, &b_CaloEmEtInEE);
    fChain->SetBranchAddress("CaloEmEtInHF", &CaloEmEtInHF, &b_CaloEmEtInHF);
    fChain->SetBranchAddress("CaloMaxEtInEmTowers", &CaloMaxEtInEmTowers, &b_CaloMaxEtInEmTowers);
    fChain->SetBranchAddress("CaloMaxEtInHadTowers", &CaloMaxEtInHadTowers, &b_CaloMaxEtInHadTowers);
    fChain->SetBranchAddress("CaloMetPt", CaloMetPt, &b_CaloMetPt);
    fChain->SetBranchAddress("CaloMetPx", CaloMetPx, &b_CaloMetPx);
    fChain->SetBranchAddress("CaloMetPy", CaloMetPy, &b_CaloMetPy);
    fChain->SetBranchAddress("CaloMetPhi", CaloMetPhi, &b_CaloMetPhi);
    fChain->SetBranchAddress("CaloMetSumEt", CaloMetSumEt, &b_CaloMetSumEt);
    fChain->SetBranchAddress("Delta_phi", &Delta_phi, &b_Delta_phi);
    fChain->SetBranchAddress("PFMetPt", PFMetPt, &b_PFMetPt);
    fChain->SetBranchAddress("PFMetPx", PFMetPx, &b_PFMetPx);
    fChain->SetBranchAddress("PFMetPy", PFMetPy, &b_PFMetPy);
    fChain->SetBranchAddress("PFMetPhi", PFMetPhi, &b_PFMetPhi);
    fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt, &b_PFMetSumEt);
    fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF, &b_Delta_phiPF);
    fChain->SetBranchAddress("TCMetPt", TCMetPt, &b_TCMetPt);
    fChain->SetBranchAddress("TCMetPx", TCMetPx, &b_TCMetPx);
    fChain->SetBranchAddress("TCMetPy", TCMetPy, &b_TCMetPy);
    fChain->SetBranchAddress("TCMetPhi", TCMetPhi, &b_TCMetPhi);
    fChain->SetBranchAddress("TCMetSumEt", TCMetSumEt, &b_TCMetSumEt);
    fChain->SetBranchAddress("Delta_phiTC", &Delta_phiTC, &b_Delta_phiTC);
    Notify();
}

Bool_t makeHistogramsMonoPho_8::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}

void makeHistogramsMonoPho_8::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t makeHistogramsMonoPho_8::Cut(){
    //Int_t makeHistogramsMonoPho_8::Cut(Long64_t entry){
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

float makeHistogramsMonoPho_8::deltaPhi(float phi1, float phi2){
    float result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2.*TMath::Pi();
    while (result <= -1.*TMath::Pi()) result += 2.*TMath::Pi();
    return result;
}

float makeHistogramsMonoPho_8::deltaR(float eta1, float phi1, float eta2, float phi2){
    float dPhi = deltaPhi(phi1,phi2);
    float dEta = eta1-eta2;
    return sqrt(dPhi*dPhi + dEta*dEta);
}

float makeHistogramsMonoPho_8::fixPhi(float phi){
    return (phi >= 0. ? phi : (2.*TMath::Pi() + phi));
}

//before using this, make sure phi1 and phi2 are defined on the same range
float makeHistogramsMonoPho_8::absDeltaPhi(float phi1, float phi2){
    float dPhi = fabs(phi1 - phi2);
    if(dPhi > TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
    return dPhi;
}

//////////////////////////////////////////////////////////
// 
// Halo functions
// 
//////////////////////////////////////////////////////////

// Example useage isHEHalo(Photon_sc_phi[x], HERecHit_subset_n, HERecHit_subset_x, HERecHit_subset_y, HERecHit_subset_energy, HERecHit_subset_time, false)
// useTime is defualt false (specified in .h file)
bool makeHistogramsMonoPho_8::isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime){
    
    float HERecHitEnergy_MIN_CUT = 1.;
    float deltaPhi_he_MAX_CUT = 0.2;
    float rho_he_MIN_CUT = 115.;
    float rho_he_MAX_CUT = 130.;
    float HERecHitTime_MAX_CUT = 0.;
    
    for(int hehit = 0; hehit < nAllHERecHits && hehit < 10000; hehit++){
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
bool makeHistogramsMonoPho_8::isTrackHalo(float photonSCphi, int nCosMu, bool CosOuterTrack_isNonnull[], float CosTrackX[], float CosTrackY[]){
    float deltaPhi_cosTrack_MAX_CUT = 0.2;
    float rho_cosTrack_MIN_CUT = 115.;
    float rho_cosTrack_MAX_CUT = 200.;
    
    for(int CosTrack = 0; CosTrack < nCosMu && CosTrack<200; CosTrack++){
        bool deltaPhiCosTrack0_isHalo = false;
        bool rhoCosTrack0_isHalo = false;
        if(!CosOuterTrack_isNonnull[CosTrack]) continue;
        float CosTrackPhi = TMath::ATan2(CosTrackY[CosTrack],CosTrackX[CosTrack]);
        float deltaPhi0 = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CosTrackPhi));
        
        float rho0 = sqrt(CosTrackX[CosTrack]*CosTrackX[CosTrack] + CosTrackY[CosTrack]*CosTrackY[CosTrack]);
        
        if(deltaPhi0 <= deltaPhi_cosTrack_MAX_CUT) deltaPhiCosTrack0_isHalo = true;
        if(rho_cosTrack_MIN_CUT < rho0 && rho0 < rho_cosTrack_MAX_CUT) rhoCosTrack0_isHalo = true;
        if( deltaPhiCosTrack0_isHalo && rhoCosTrack0_isHalo) return true;
    }// loop over CosTracks
    
    return false;
}

// Example useage isCSCHalo(Photon_sc_phi[x], CSCseg_n, CSCseg_x, CSCseg_y, CSCseg_time, false)
// useTime is defualt false (specified in .h file)
bool makeHistogramsMonoPho_8::isCSCHalo(float photonSCphi, int nAllCSCSegments, float CSCSegmentX[], float CSCSegmentY[], float CSCSegmentTime[], bool useTime){
    
    float deltaPhi_csc_MAX_CUT = 0.2;
    float CSCSegmentTime_MAX_CUT = 0.;
    
    for(int cscseg = 0; cscseg < nAllCSCSegments && cscseg<10000; cscseg++){
        bool deltaPhiCSC_isHalo = false;
        //bool rhoCSC_isHalo = false;
        bool timeCSC_isHalo = false;
        
        float CSCSegmentPhi = TMath::ATan2(CSCSegmentY[cscseg],CSCSegmentX[cscseg]);
        float deltaPhi_csc = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CSCSegmentPhi));
        
        //float rho_csc = sqrt( CSCSegmentX[cscseg]*CSCSegmentX[cscseg] + CSCSegmentY[cscseg]*CSCSegmentY[cscseg]);
        
        if(deltaPhi_csc<deltaPhi_csc_MAX_CUT) deltaPhiCSC_isHalo = true;
        if(CSCSegmentTime[cscseg]<CSCSegmentTime_MAX_CUT) timeCSC_isHalo = true;
        
        if(useTime){
            if(deltaPhiCSC_isHalo && timeCSC_isHalo) return true;
        }else{ // don't useTime
            if(deltaPhiCSC_isHalo) return true;
        }
    }// loop over CSC segments
    return false;
}

// Example useage deltaPhiCSCHalo(Photon_sc_phi[x], CSCseg_x[x], CSCseg_y[x])
float makeHistogramsMonoPho_8::deltaPhiCSCHalo(float photonSCphi, float CSCSegmentX, float CSCSegmentY){
    
    float CSCSegmentPhi = TMath::ATan2(CSCSegmentY,CSCSegmentX);
    float deltaPhi_csc = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CSCSegmentPhi));
    //cout<<"CSCSegmentY: "<<CSCSegmentY<<"CSCSegmentX: "<<CSCSegmentX<<"PhotonSCphi:"<< fixPhi(photonSCphi)<<"CSCSegmentPhi: "<<fixPhi(CSCSegmentPhi)<<"deltaPhi_csc: " <<deltaPhi_csc<<endl; 
    return deltaPhi_csc;
}

//////////////////////////////////////////////////////////
// 
// Cosmic functions
//
//////////////////////////////////////////////////////////
/*
// Example useage fisherCosmic3(Photon_Roundness[pho],Photon_Angle[pho],Photon_hasConvTrk[pho])
float makeHistogramsMonoPho_8::fisherCosmic3(float Photon_Roundness, float Photon_Angle, bool Photon_hasConvTrk){
    
    TH1F *sigRound          = (TH1F*)file_PHOTON->Get("h_Rho_MC");                                                                                                                     
    TH1F *sigAngle          = (TH1F*)file_PHOTON->Get("h_Angle_MC");
    TH1F *sighasconversion  = (TH1F*)file_PHOTON->Get("h_hasconversion_MC");                                                                                                  
    TH1F *bkgRound          = (TH1F*)file_COSMIC->Get("h_Rho_cos");
    TH1F *bkgAngle          = (TH1F*)file_COSMIC->Get("h_Angle_cos");
    TH1F *bkghasconversion  = (TH1F*)file_COSMIC->Get("h_hasconversion_cos");
    
    float sigProb1          = prob(sigRound, Photon_Roundness);
    float bkgProb1          = prob(bkgRound, Photon_Roundness);
    
    float sigProb2          = prob(sigAngle, Photon_Angle);
    float bkgProb2          = prob(bkgAngle, Photon_Angle);
    
    float sigProb3          = prob2(sighasconversion,Photon_hasConvTrk);
    float bkgProb3          = prob2(bkghasconversion,Photon_hasConvTrk);
    
    //cout<<"(signal like event if prob_Signal ~1 and bkg like event if prob_Signal ~0)"<<endl;
    return sigProb1*sigProb2*sigProb3/(sigProb1*sigProb2*sigProb3 + bkgProb1*bkgProb2*bkgProb3);
    
}

// Example useage fisherCosmic2(Photon_Roundness[pho],Photon_Angle[pho])
float makeHistogramsMonoPho_8::fisherCosmic2(float Photon_Roundness, float Photon_Angle){
    //======= calculation fisher-discriminant  ===
    TH1F *sigRound          = (TH1F*)file_PHOTON->Get("h_Rho_MC");
    TH1F *sigAngle          = (TH1F*)file_PHOTON->Get("h_Angle_MC");
    TH1F *bkgRound          = (TH1F*)file_COSMIC->Get("h_Rho_cos");
    TH1F *bkgAngle          = (TH1F*)file_COSMIC->Get("h_Angle_cos");
    
    float sigProb1          = prob(sigRound, Photon_Roundness);
    float bkgProb1          = prob(bkgRound, Photon_Roundness);
    
    float sigProb2          = prob(sigAngle, Photon_Angle);
    float bkgProb2          = prob(bkgAngle, Photon_Angle);
    
    //probability of finding signal or background using shower shape variables(signal like event if prob_2 ~1 and bkg like event if prob_2 ~0)
    return sigProb1*sigProb2/(sigProb1*sigProb2 + bkgProb1*bkgProb2);
}

float makeHistogramsMonoPho_8::prob(TH1F* h, float value){
    int iBin = h->FindBin(value);
    float norm = h->Integral(1, h->GetNbinsX());
    return float(h->GetBinContent(iBin))/norm;
}

float makeHistogramsMonoPho_8::prob2(TH1F* h, bool Photon_hasConvTrk){
    float value=0.;
    if(Photon_hasConvTrk){value=1.;}
    else{value=0.;}
    int iBin = h->FindBin(value);
    float norm = h->Integral(1, h->GetNbinsX());
    return float(h->GetBinContent(iBin))/norm;
}

*/
float makeHistogramsMonoPho_8::sigmaRCosmic(float Pho_sc_x, float Pho_sc_y, float Pho_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n){
    float dr_min_cos = 99999;
    float dx_min_cos = 10000;
    float dy_min_cos = 10000;
    float dz_min_cos = 10000;
    
    vector< pair<bool,TVector3> > solutions;
    
    for(int iMuon=0 ; iMuon <CosmicMuon_n && iMuon<200; iMuon++){   
        bool singlesol = false;
        bool bothsol = false;
        
        float averageSolution_x1 = 0.;
        float averageSolution_y1 = 0.;
        float averageSolution_z1 = 0.;
        float averageSolution_x2 = 0.;
        float averageSolution_y2 = 0.;
        float averageSolution_z2 = 0.;
        TVector3 point(CosmicMuon_OuterTrack_InnerPoint_x[iMuon],CosmicMuon_OuterTrack_InnerPoint_y[iMuon],CosmicMuon_OuterTrack_InnerPoint_z[iMuon]);		
        TVector3 direction(CosmicMuon_OuterTrack_InnerPoint_px[iMuon],CosmicMuon_OuterTrack_InnerPoint_py[iMuon],CosmicMuon_OuterTrack_InnerPoint_pz[iMuon]);
        solutions= getEBIntersections(point, direction);
        
        if(solutions.size()==4){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
            if(solutions[2].first==true && solutions[3].first==false){
                averageSolution_x2 = solutions[2].second.X();
                averageSolution_y2 = solutions[2].second.Y();
                averageSolution_z2 = solutions[2].second.Z();
            }else if(solutions[2].first==false && solutions[3].first==true){
                averageSolution_x2 = solutions[3].second.X();
                averageSolution_y2 = solutions[3].second.Y();
                averageSolution_z2 = solutions[3].second.Z();
            }else if(solutions[2].first==true && solutions[3].first==true){
                averageSolution_x2 = (solutions[2].second.X() + solutions[3].second.X())/2.;
                averageSolution_y2 = (solutions[2].second.Y() + solutions[3].second.Y())/2.;
                averageSolution_z2 = (solutions[2].second.Z() + solutions[3].second.Z())/2.;
            }
            
        }else if(solutions.size()==2){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
        }else if(solutions.size()==1){
            if(solutions[0].first == true){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }
        }
        //========distance of separation between first muon intersection in ECAL and SC position 
        
        //cout << "**dx1** "<< dx2_cos << " dy1: " << dy2_cos << "dz1: "<< dz2_cos << endl;
        
        //single solution means that the muon track only reported a DCA (solutions.size==1)
        // or that the muon grazed the ecal, never entering the tracking region
        // (solutions.size==2) where both of these solutions are intersections with the outer barrel only
        if(solutions.size()==2 || solutions.size()==1) singlesol = true;
        //both solutions means that the muon track passed through the tracking region so that
        // there is an intersection with the outer barrel and inner barrel on it's way in, and
        // there is an intersection with the inner barrel and outer barrel on it's way out hence
        // (solutions.size()==4)
        if(solutions.size()==4) bothsol = true;
        
        //initialize these
        float dr_cos = -1.;
        float dx_cos = 10000.;
        float dy_cos = 10000.;
        float dz_cos = 10000.;
        float dr2_cos = -1.;
        float dx2_cos = 10000.;
        float dy2_cos = 10000.;
        float dz2_cos = 10000.;
        
        if(bothsol == true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dx2_cos = Pho_sc_x - averageSolution_x2;
            dy2_cos = Pho_sc_y - averageSolution_y2;
            dz2_cos = Pho_sc_z - averageSolution_z2;
            dr_cos  = sqrt((dx_cos*dx_cos/4.4/4.4) + (dy_cos*dy_cos/4.1/4.1) + (dz_cos*dz_cos/4.6/4.6));
            dr2_cos  = sqrt((dx2_cos*dx2_cos/4.4/4.4) + (dy2_cos*dy2_cos/4.1/4.1) + (dz2_cos*dz2_cos/4.6/4.6));
            if( dr2_cos < dr_min_cos ){
                dr_min_cos = dr2_cos;
                dx_min_cos = dx2_cos;
                dy_min_cos = dy2_cos;
                dz_min_cos = dz2_cos;
            }
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
            
        }
        if (singlesol== true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dr_cos  = sqrt((dx_cos*dx_cos/5.2/5.2) + (dy_cos*dy_cos/4.8/4.8) + (dz_cos*dz_cos/9/9));
            
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
        }
        //cout << "mydrcos: "<< dr_min_cos << endl;
    }//loop muon
    
    return dr_min_cos; 
    
}

float makeHistogramsMonoPho_8::deltaRCosmic(float Pho_sc_x, float Pho_sc_y, float Pho_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n){
    float dr_min_cos = 99999;
    float dx_min_cos = 10000;
    float dy_min_cos = 10000;
    float dz_min_cos = 10000;
    
    vector< pair<bool,TVector3> > solutions;
    
    for(int iMuon=0 ; iMuon <CosmicMuon_n && iMuon<200; iMuon++){   
        bool singlesol = false;
        bool bothsol = false;
        
        float averageSolution_x1 = 0.;
        float averageSolution_y1 = 0.;
        float averageSolution_z1 = 0.;
        float averageSolution_x2 = 0.;
        float averageSolution_y2 = 0.;
        float averageSolution_z2 = 0.;
        TVector3 point(CosmicMuon_OuterTrack_InnerPoint_x[iMuon],CosmicMuon_OuterTrack_InnerPoint_y[iMuon],CosmicMuon_OuterTrack_InnerPoint_z[iMuon]);		
        TVector3 direction(CosmicMuon_OuterTrack_InnerPoint_px[iMuon],CosmicMuon_OuterTrack_InnerPoint_py[iMuon],CosmicMuon_OuterTrack_InnerPoint_pz[iMuon]);
        solutions= getEBIntersections(point, direction);
        
        if(solutions.size()==4){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
            if(solutions[2].first==true && solutions[3].first==false){
                averageSolution_x2 = solutions[2].second.X();
                averageSolution_y2 = solutions[2].second.Y();
                averageSolution_z2 = solutions[2].second.Z();
            }else if(solutions[2].first==false && solutions[3].first==true){
                averageSolution_x2 = solutions[3].second.X();
                averageSolution_y2 = solutions[3].second.Y();
                averageSolution_z2 = solutions[3].second.Z();
            }else if(solutions[2].first==true && solutions[3].first==true){
                averageSolution_x2 = (solutions[2].second.X() + solutions[3].second.X())/2.;
                averageSolution_y2 = (solutions[2].second.Y() + solutions[3].second.Y())/2.;
                averageSolution_z2 = (solutions[2].second.Z() + solutions[3].second.Z())/2.;
            }
            
        }else if(solutions.size()==2){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
        }else if(solutions.size()==1){
            if(solutions[0].first == true){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }
        }
        //========distance of separation between first muon intersection in ECAL and SC position 
        
        //cout << "**dx1** "<< dx2_cos << " dy1: " << dy2_cos << "dz1: "<< dz2_cos << endl;
        
        //single solution means that the muon track only reported a DCA (solutions.size==1)
        // or that the muon grazed the ecal, never entering the tracking region
        // (solutions.size==2) where both of these solutions are intersections with the outer barrel only
        if(solutions.size()==2 || solutions.size()==1) singlesol = true;
        //both solutions means that the muon track passed through the tracking region so that
        // there is an intersection with the outer barrel and inner barrel on it's way in, and
        // there is an intersection with the inner barrel and outer barrel on it's way out hence
        // (solutions.size()==4)
        if(solutions.size()==4) bothsol = true;
        
        //initialize these
        float dr_cos = -1.;
        float dx_cos = 10000.;
        float dy_cos = 10000.;
        float dz_cos = 10000.;
        float dr2_cos = -1.;
        float dx2_cos = 10000.;
        float dy2_cos = 10000.;
        float dz2_cos = 10000.;
        
        if(bothsol == true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dx2_cos = Pho_sc_x - averageSolution_x2;
            dy2_cos = Pho_sc_y - averageSolution_y2;
            dz2_cos = Pho_sc_z - averageSolution_z2;
            dr_cos  = sqrt((dx_cos*dx_cos) + (dy_cos*dy_cos) + (dz_cos*dz_cos));
            dr2_cos  = sqrt((dx2_cos*dx2_cos) + (dy2_cos*dy2_cos) + (dz2_cos*dz2_cos));
            if( dr2_cos < dr_min_cos ){
                dr_min_cos = dr2_cos;
                dx_min_cos = dx2_cos;
                dy_min_cos = dy2_cos;
                dz_min_cos = dz2_cos;
            }
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
            
        }
        if (singlesol== true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dr_cos  = sqrt((dx_cos*dx_cos) + (dy_cos*dy_cos) + (dz_cos*dz_cos));
            
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
        }
        //cout << "mydrcos: "<< dr_min_cos << endl;
    }//loop muon
    
    return dr_min_cos; 
    
}

vector<pair <bool,TVector3> > makeHistogramsMonoPho_8::getEBIntersections(const TVector3& X, const TVector3& P){
    
    //Ecal Barrel parameters
    float OuterEBZ   = 291.12;
    float InnerEBZ   = 270.89;
    //note OuterEBRho is defined as a function because it varies with z
    float InnerEBRho = 129.;
    
    double t_closestZ = tDCA(X,P);  
    TVector3 X_closestZ = P*t_closestZ + X;
    double OuterEBRho_closestZ = OuterEBRho( X_closestZ.Z() );
    
    if( fabs(X_closestZ.Perp())>OuterEBRho_closestZ){
        TVector3 PointOnEBtemp(0.,0.,0.);
        pair<bool,TVector3> PointOnEB(false,PointOnEBtemp);
        if( fabs(X_closestZ.Z())< OuterEBZ){
            PointOnEB.first = true;
            double phi = X_closestZ.Phi();
            double x = OuterEBRho_closestZ*sin(phi);
            double y = OuterEBRho_closestZ*cos(phi);
            TVector3 temp(x,y,X_closestZ.Z());
            PointOnEB.second = temp;
        }
        vector<pair<bool,TVector3> > vectorOfEBPoints;
        vectorOfEBPoints.push_back(PointOnEB);
        return vectorOfEBPoints;
        
    }   
    
    
    else if( fabs(InnerEBRho)<fabs(X_closestZ.Perp()) && fabs(X_closestZ.Perp())<fabs(OuterEBRho_closestZ)){
        vector<double> outerTs = outerIntersectionTs(X,P,t_closestZ);
        
        //intersection = X0 + t P_X, where X = <x, y, z>
        TVector3 intersectionOuterEB1 = X + P * outerTs[0];
        TVector3 intersectionOuterEB2 = X + P * outerTs[1];
        //cout << "intersectionOuterEB1.Z()" << intersectionOuterEB1.Z()<< "OuterEBZ"<<OuterEBZ<< endl;
        //cout << "intersectionOuterEB2.Z()" << intersectionOuterEB2.Z()<< "OuterEBZ"<<OuterEBZ<< endl;
        vector<pair<bool,TVector3> > vectorOfEBPoints;
        
        // mytest z of 1st outer solution
        
        if(fabs(intersectionOuterEB1.Z())<OuterEBZ){ 
            vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB1));
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
	
        if(fabs(intersectionOuterEB2.Z())<OuterEBZ){
            vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB2));
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
        
        return vectorOfEBPoints;
        
    }else{ // rho_closestZ must be inside the Inner radius of the EB
        // we will report 4 intersections 2 with Ecal inner surface,
        // and 2 with the outer surface.
        
        vector<double> innerTs = innerIntersectionTs(129.,X,P);
        vector<double> outerTs = outerIntersectionTs(X,P,t_closestZ);
        
        TVector3 intersectionInnerEB1 = X +P * innerTs[0];
        TVector3 intersectionInnerEB2 = X + P * innerTs[1];
        TVector3 intersectionOuterEB1 = X + P * outerTs[0];
        TVector3 intersectionOuterEB2 = X + P * outerTs[1];
        
        vector< pair<bool,TVector3> > vectorOfEBPoints;
        if(fabs(intersectionOuterEB1.Z())<OuterEBZ){
            vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB1));
            
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
        
        if(fabs(intersectionInnerEB1.Z())<InnerEBZ){
            vectorOfEBPoints.push_back(make_pair(true,intersectionInnerEB1));
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
        
        if(fabs(intersectionInnerEB2.Z())<InnerEBZ){
            vectorOfEBPoints.push_back(make_pair(true,intersectionInnerEB2));
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
        
        // mytest z of 2nd outer solution
        if(fabs(intersectionOuterEB2.Z())<OuterEBZ){
            vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB2));
            //cout << "*2nd outer solution" <<intersectionOuterEB2.Z() <<  "OuterEBZ " << OuterEBZ <<endl;
        }else{
            TVector3 zero(0.,0.,0.);
            vectorOfEBPoints.push_back(make_pair(false,zero));
        }
        
        return vectorOfEBPoints;
    }
    
}

double makeHistogramsMonoPho_8::h(double z, double& length){
    // all barrel crystals point toward the tangent of a 6.67 cm circle around the IP
    // so they are semi-projective
    double r = 6.67;
    
    // z is the distance from the IP to the z of the inner surface of the crystal.
    // zprime is the distance from the IP to the point where the tangent line described
    // above crosses the z-axis in terms of z.
    // note: 129 cm is the radius of the inner barrel surface
    double zprime = -2.*r*r*z + sqrt(4*r*r*r*r*z*z + 4*r*r*(z*z+129.*129.)*(129.*129.-r*r));
    zprime = zprime/2./(129.*129. - r*r);
    
    // height is the projection of the crystal length onto the xy-plane.
    // note 23 cm is the length of all of the crystals
    double height = 23.*r/zprime;
    //cout << "height " << height <<endl;
    // length is the projection of the crystal length onto the z-axis
    length = sqrt(23.*23. - height*height);
    
    return height;
    
}

// iterations = 5 by default in the .h file
double makeHistogramsMonoPho_8::OuterEBRho(double zprimeprime, int iterations ){
    // zprimeprime is just a funny name for the z location where you would
    // like to know the radius of the outer envelope of the barrel Ecal.
    // This is computed numerically below.
    double lengthOfTheCrystal = 23.0;
    
    double minDZ = 1000;
    double bestH = 1000;
    for(int i = 0; i < iterations; ++i) {
        // newZ is a trial Z point we are testing
        double newZ = zprimeprime - lengthOfTheCrystal*(1. - (1./iterations)*i);
        
        // get length and height of this crystal whose front face is at newZ
        double length = 0;
        double height = h(newZ, length); // length is passed by reference!
        // face z postion, newZ, we must vary newZ until
        // zprimeprime = newZ + length
        double dz = fabs(zprimeprime - (newZ + length));
        if( dz < minDZ ){
            minDZ = dz;
            bestH = height;
            ////////cout << "bestH" << bestH << endl;
            //       //////cout << "z\'\' =  " << zprimeprime 
            //    << "\t new z = " << newZ << "\t length = " << length 
            //    << "\t approx value (z + length) =  " << newZ + length 
            //    << "\t height = " << height
            //    << endl;
        }
    }
    
    return bestH+129.0; // return the radius of the outer EB surface
}

double makeHistogramsMonoPho_8::tDCA(const TVector3& X, const TVector3& P){
    
    // get position, and momentum components
    double x0 = X.X();
    double y0 = X.Y();
    double px = P.X();
    double py = P.Y();
    
    // solve for t of closest approach to z-axis:
    double t_closestZ = -(y0*py + x0*px)/(px*px + py*py);
    
    return t_closestZ;
}

vector <double> makeHistogramsMonoPho_8::outerIntersectionTs(const TVector3& X, const TVector3& P, double timeDCA){
    // get limits of t inside which the outer intersection solutions must exist
    double x0 = X.X();
    double y0 = X.Y();
    double z0 = X.Z();
    double px = P.X();
    double py = P.Y();
    double pz = P.Z();
    
    //vector <double> tLimits = innerIntersectionTs(155.,X,P);
    vector <double> tLimits = innerIntersectionTs(155.,X,P);
    // make steps of 1cm
    TVector3 startPoint = tLimits[0]*P + X;
    TVector3 endPoint   = tLimits[1]*P + X;
    
    double pathLength = sqrt( (endPoint.X()-startPoint.X())*(endPoint.X()-startPoint.X())
                             + (endPoint.Y()-startPoint.Y())*(endPoint.Y()-startPoint.Y())
                             + (endPoint.Z()-startPoint.Z())*(endPoint.Z()-startPoint.Z()) );
    
    double timeStep = (tLimits[1]-tLimits[0])/pathLength;
    double testTime = tLimits[0];
    double t1 = testTime;
    double smallestDiff1 = 1000.;
    double t2 = testTime;
    double smallestDiff2 = 1000.;
    
    while(1){
        testTime+=timeStep;
        if(testTime>tLimits[1]) break;
        // look at quadratic form for t described in innerIntersectionTs
        // here R is no longer a constant, but varies with z. We now scan
        // t to make this equation as close to 0 as possible.
        double testZ = pz*testTime + z0;
        //cout << "testZ " << testZ << "z0 "<< z0<< "y0" << y0 << "z0 "<< z0 <<endl;
        double testR = OuterEBRho(testZ);
        double testDiff = fabs(  testTime*testTime * (px*px + py*py)
                               + testTime          * 2.*(px*x0 + py*y0)
                               + x0*x0 + y0*y0 - testR*testR
                               );
        
        if(testTime<timeDCA){
            if(testDiff<smallestDiff1){
                smallestDiff1 = testDiff;
                
                t1 = testTime;
                
            }
        }else{ 
            if(testDiff<smallestDiff2){
                smallestDiff2 = testDiff;
                t2 = testTime;
                
            }
        }//else timeDCA
    }//while
    
    vector<double> t;
    t.push_back(t1);
    t.push_back(t2);
    return t;
}

vector <double> makeHistogramsMonoPho_8::innerIntersectionTs(double R, const TVector3& X, const TVector3& P){  
    double x0 = X.X();
    double y0 = X.Y();
    double px = P.X();
    double py = P.Y();
    
    // this b_over_2 = b/2 from eq 3.
    double b_over_2 = x0 * px + y0 * py ;
    
    // this thingUnderSqrt = (b/2)^2 - a c from eq 3.
    double thingUnderSqrt = b_over_2 * b_over_2 - (px*px + py*py) * (x0*x0 + y0*y0 - R*R) ;
    vector <double> t;
    t.push_back(   (- x0*px - y0*py - sqrt(thingUnderSqrt)) / (px*px + py*py)  );
    t.push_back(   (- x0*px - y0*py + sqrt(thingUnderSqrt)) / (px*px + py*py)  );
    return t;
    
}

//////////////////////////////////////
//
// Selection functions
//
//////////////////////////////////////
//select events with a good vertex
//http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DPGAnalysis/Skims/python/GoodVertex_cfg.py?revision=1.5&view=markup
//example useage eventHasGoodVertex( Vertex_n,  Vertex_ndof,  Vertex_z,  Vertex_d0, Vertex_isFake);
bool makeHistogramsMonoPho_8::eventHasGoodVertex(int v_n, int v_ndof[], float v_z[], float v_d0[], bool v_isFake[], bool onlyOne){
    int N_goodVerticies=0;
    for(int v=0; v<v_n && v<200 ;v++)
        if(v_ndof[v]>4 && fabs(v_z[v])<=24. && fabs(v_d0[v])<=2 && v_isFake[v]!=true) N_goodVerticies++;
    if(onlyOne && N_goodVerticies==1){
        return true;
    }else if(!onlyOne && N_goodVerticies>0){
        return true;
    }else{
        return false;
    }
    
}

//example useage passesPhotonID(Photon_ncrys[x], Photon_pt[x],Photon_sc_eta[x],Photon_HoE[x],Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x],Photon_SigmaIetaIeta[x])
bool makeHistogramsMonoPho_8::passesPhotonID( float ncrys,float pt, float sc_eta, float HoE, float trkSumPtHollowConeDR04,
                                        float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float sigmaIEIE){
    bool thisPassesPhotonID=false;
    if(ncrys<=0) return false;
    if(  pt>80.
        &&  fabs(sc_eta) <              1.4442
        &&  trkSumPtHollowConeDR04 <    2.0+0.001*pt
        &&  ecalRecHitSumEtConeDR04<    4.2+0.006*pt
        &&  hcalTowerSumEtConeDR04 <    2.2+0.0025*pt
        &&  HoE <		        0.05
        &&  sigmaIEIE < 0.013
        ) thisPassesPhotonID=true;
    return thisPassesPhotonID;

}

//example makeHistogramsMonoPho_8::invMass(Photon_E[i], Photon_px[i], Photon_py[i]), Photon_pz[i], Photon_E[j], Photon_px[j], Photon_py[j]), Photon_pz[j])
float makeHistogramsMonoPho_8::invMass(float E, float px, float py, float pz, float E2, float px2, float py2, float pz2){
    return sqrt( 2. * E  * E2
                -2. * px * px2
                -2. * py * py2
                -2. * pz * pz2
                );
}

//useage makeHistogramsMonoPho_8::massT(Photon_phi[i],Photon_et[i],TCMetPhi[0],TCMetPt[0])
float makeHistogramsMonoPho_8::massT(float phi, float et, float metphi, float met){
    float dPhi= deltaPhi(fixPhi(phi),fixPhi(metphi));
    return sqrt( 2.*et*met*(1.-cos(dPhi)) );
}

//useage makeHistogramsMonoPho_8::eventHasMu(CosmicMuon_n, CosmicMuon_OuterTrack_isNonnull, CosmicMuon_isStandAloneMuon)
bool makeHistogramsMonoPho_8::eventHasMu(int nCosMu, bool CosOuterTrack_isNonnull[], bool CosTrack_isSA[]){
    bool answer=false;
    for(int i=0; i<nCosMu && i<200 ; i++){
        if(CosOuterTrack_isNonnull[i] && CosTrack_isSA[i]) answer=true;
    }
    return answer;
}

//useage makeHistogramsMonoPho_8::phoHasMu(CosmicMuon_n, CosmicMuon_OuterTrack_isNonnull, CosmicMuon_isStandAloneMuon,Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx])
bool makeHistogramsMonoPho_8::phoHasMu(int nMu, bool MuOuterTrack_isNonnull[], bool MuTrack_isSA[],float MuTrack_Pt[], float MuTrack_Eta[],
                                float MuTrack_Phi[], float MuTrack_Vx[], float MuTrack_Vy[], float MuTrack_Vz[], float phoPhi,
                                float phoEta, float phoVx, float phoVy, float phoVz, float muMinEt, bool restrictEta, 
                                float maxEta, bool restrictVertex, float maxVertexDist){
    bool answer=false;
    for(int i=0; i<nMu && i<200; i++){
        float dR=deltaR(MuTrack_Eta[i],MuTrack_Phi[i],phoEta,phoPhi);
        if(dR<0.3) continue; //ah! this jet lines up with the photon!
        if(MuTrack_Pt[i]<muMinEt) continue;
        if(restrictEta && fabs(MuTrack_Eta[i])>maxEta) continue;
        if(!MuOuterTrack_isNonnull[i]) continue;
        if(!MuTrack_isSA[i]) continue;
        float d = sqrt( (phoVx-MuTrack_Vx[i])*(phoVx-MuTrack_Vx[i]) + (phoVy-MuTrack_Vy[i])*(phoVy-MuTrack_Vy[i]) + (phoVz-MuTrack_Vz[i])*(phoVz-MuTrack_Vz[i]));
        if(restrictVertex && d>maxVertexDist) continue;
        answer=true;
    }
    return answer;
}

//useage makeHistogramsMonoPho_8::computeLICTD(Photon_ncrys[x], Photon_timing_xtal[x], Photon_energy_xtal[x])
float makeHistogramsMonoPho_8::computeLICTD(int nCry, float cryTime[], float cryEnergy[]){
    int seed = GetSeed(nCry, cryEnergy);
    if(seed<0) return -99;
    Float_t LICTD =0;
    Int_t crysCrys=-1;
    Int_t crysThresh=0;
    for(int k=0;k<nCry&&k<100;++k){
        if (seed==k) continue;
        Float_t crysE = cryEnergy[k];
        if(crysE > 1.){
	    crysThresh++;
	    Float_t tdiff = cryTime[seed] - cryTime[k];
	    if(fabs(tdiff) > fabs(LICTD)){
                LICTD = tdiff;
                crysCrys=k;
	    }
        }
    }
    return LICTD;
}

//  seed finding function
//useage makeHistogramsMonoPho_8::GetSeed(Photon_ncrys[x], Photon_energy_xtal[x])
int makeHistogramsMonoPho_8::GetSeed(int Ncry, float cryE[]){
    int seedIdx=-1;
    float seedE=-99;
    for(int k=0;k<Ncry&&k<100;++k){
        if(cryE[k] > seedE){
            seedIdx = k;
            seedE=cryE[seedIdx];
        }
    }
    return seedIdx;  
}

//useage makeHistogramsMonoPho_8::GetE2overE9(Photon_ncrys[x], Photon_ieta_xtalEB[x], Photon_iphi_xtalEB[x], Photon_energy_xtal[x])
float makeHistogramsMonoPho_8::GetE2overE9(int Ncry, int ieta[], int iphi[], float e[]){
    int seed = GetSeed(Ncry, e);
    if(seed<0) return -1;
    if(abs(ieta[seed])>85) return -2;
    float e9=0;
    float secondEMax=-99;
    for(int cry=0;cry<100 && cry<Ncry;cry++){
        // skip the seed crystal
        if(cry==seed) continue;
        
        //only consider barrel hits
        if(abs(ieta[cry])>85) continue;
        
        //get rid of no ieta=0 problem :P
        int ieta_cry = ieta[cry];
        if(ieta_cry<0) ieta_cry++;
        int ieta_seed = ieta[seed];
        if(ieta_seed<0) ieta_seed++;
        
        //compute deltaIEta and deltaIPhi
        int delIEta = ieta_seed - ieta_cry;
        int delIPhi = iphi[seed] - iphi[cry];
        if(delIPhi==-359) delIPhi= 1;
        if(delIPhi==359)  delIPhi=-1;
        if(abs(delIPhi)<=1 && abs(delIEta)<=1){
            e9+=e[cry];
            if(e[cry]>secondEMax) secondEMax=e[cry];
        }
    }//loop over cry
    
        if(e9>0. && secondEMax>0.){
            return secondEMax/e9;
        }
        return -3;
}


//useage makeHistogramsMonoPho_8::triggerPassed(ntriggers,myDesiredTrigger,triggernames,ifTriggerpassed)
bool makeHistogramsMonoPho_8::triggerPassed(unsigned int Ntriggers, vector <string> myDesiredTrigger, vector <string> storedTriggers_name, vector <bool> storedTriggers_bool){
    bool kTriggerPassed=false;
    
    for (unsigned int t=0;t<Ntriggers && t<storedTriggers_name.size();++t){
        string trig = storedTriggers_name[t];//string trig = (*storedTriggers_name)[t];
        //cout<<storedTriggers_bool[t]<<endl;
         if ( storedTriggers_bool[t]==true){
          //cout << "trig: " << trig << endl;
          for (unsigned int jt=0;jt<myDesiredTrigger.size();++jt){
            if (trig.compare(myDesiredTrigger[jt])==0)
              kTriggerPassed=true;
          }    
        }    
      }    
      //cout <<"kTriggerPassed  "<<kTriggerPassed<<endl;
      return kTriggerPassed;
}

//example useage newMetPt(TCMetPt[0], TCMetPhi[0], Photon_pt[x], Photon_phi[x],seedT)
float makeHistogramsMonoPho_8::newMetPt(float MetPt, float MetPhi, float PhoPt, float PhoPhi, float seedT){
	float NEWMET = MetPt;
	if(fabs(seedT)>3.){
		float METX = MetPt*cos(MetPhi);
		float METY = MetPt*sin(MetPhi);
		METX = METX - PhoPt*cos(PhoPhi);
		METY = METY - PhoPt*sin(PhoPhi);
		NEWMET = sqrt(METX*METX + METY*METY);
	}
	return NEWMET;
}

//example useage newMetPhi(TCMetPt[0], TCMetPhi[0], Photon_pt[x], Photon_phi[x],seedT)
float makeHistogramsMonoPho_8::newMetPhi(float MetPt, float MetPhi, float PhoPt, float PhoPhi, float seedT){
	float NEWPHI = MetPhi;
	if(fabs(seedT)>3.){
		float METX = MetPt*cos(MetPhi);
		float METY = MetPt*sin(MetPhi);
		METX = METX - PhoPt*cos(PhoPhi);
		METY = METY - PhoPt*sin(PhoPhi);
		NEWPHI = atan2(METY,METX);
	}
        NEWPHI = fixPhi(NEWPHI);
	return NEWPHI;
}

// passesTightPhotonID(Photon_et[x],Photon_HoE[x],Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x])
bool makeHistogramsMonoPho_8::passesTightPhotonID(float pt, float HoE, float trkSumPtHollowConeDR04,
    float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float sigmaIEIE, float eta, 
    float minEt, bool restrictEta, float maxEta){
    bool thisPassesPhotonID=false;
    if(restrictEta && fabs(eta)>maxEta) return false;
    if(     pt > minEt
        &&  trkSumPtHollowConeDR04 <    2.0+0.001*pt
        &&  ecalRecHitSumEtConeDR04<    4.2+0.006*pt
        &&  hcalTowerSumEtConeDR04 <    2.2+0.0025*pt
        &&  HoE <		        0.05
        &&  sigmaIEIE < 0.013
        ) thisPassesPhotonID=true;
    return thisPassesPhotonID;

}

//passesLoosePhotonID(Photon_et[x],Photon_HoE[x],Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x])
bool makeHistogramsMonoPho_8::passesLoosePhotonID(float pt, float HoE, float trkSumPtHollowConeDR04, float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float eta, float minEt, bool restrictEta, float maxEta){
    bool thisPassesPhotonID=false;
    if(restrictEta && fabs(eta)>maxEta) return false;
    if(     pt > minEt
        &&  trkSumPtHollowConeDR04 <    3.5+0.001*pt
        &&  ecalRecHitSumEtConeDR04<    4.2+0.006*pt
        &&  hcalTowerSumEtConeDR04 <    2.2+0.0025*pt
        &&  HoE <		        0.05
        ) thisPassesPhotonID=true;
    return thisPassesPhotonID;

}

//passesEMID(Photon_et[x],Photon_HoE[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x])
bool makeHistogramsMonoPho_8::passesEMID(float pt, float HoE, float ecalRecHitSumEtConeDR04, float hcalTowerSumEtConeDR04, float eta, float minEt, bool restrictEta, float maxEta){
    bool thisPassesPhotonID=false;
    if(restrictEta && fabs(eta)>maxEta) return false;
    if(     pt > minEt
        &&  ecalRecHitSumEtConeDR04<    4.2+0.006*pt
        &&  hcalTowerSumEtConeDR04 <    2.2+0.0025*pt
        &&  HoE <		        0.05
        ) thisPassesPhotonID=true;
    return thisPassesPhotonID;

}

//hasObj(Photon_phi[leadPhoIdx], Photon_eta[leadPhoIdx],Photon_vx[leadPhoIdx],Photon_vy[leadPhoIdx],Photon_vz[leadPhoIdx],Track_n,Track_pt,Track_eta,Track_phi,Track_vx,Track_vy,Track_vz,5.)
bool makeHistogramsMonoPho_8::hasObj(float phoPhi, float phoEta, float phoVx, float phoVy, float phoVz,
                                int objN, float objEt[], float objEta[], float objPhi[],
                                float objVx[], float objVy[], float objVz[], float objMinEt,
                                bool restrictEta, float maxEta, bool restrictVertex, float maxVertexDist){
    bool hasObjAnswer = false;
    for(int i=0; i<objN && i<200;i++){
        float d = sqrt( (phoVx-objVx[i])*(phoVx-objVx[i]) + (phoVy-objVy[i])*(phoVy-objVy[i]) + (phoVz-objVz[i])*(phoVz-objVz[i]));
        if(restrictVertex && d>maxVertexDist) continue;
        if(restrictEta && fabs(objEta[i])>maxEta) continue;
        if(objEt[i]<objMinEt) continue;
        float dR=deltaR(objEta[i],objPhi[i],phoEta,phoPhi);
        if(dR<0.3) continue; //ah! this jet lines up with the photon!
        
        hasObjAnswer = true;
    }
    return hasObjAnswer;
}

//hasOtherPho(0,leadPhoIdx,Photon_n,Photon_et,Photon_HoE,Photon_trkSumPtHollowConeDR04,Photon_ecalRecHitSumEtConeDR04,Photon_hcalTowerSumEtConeDR04)
//first int indicates 0->tight ID, 1->loose ID, 2->EM ID
bool makeHistogramsMonoPho_8::hasOtherPho(int phoID, int candIdx, int PhoN, float et[],float HoE[],float trkSumPtHollowConeDR04[],
    float ecalRecHitSumEtConeDR04[], float hcalTowerSumEtConeDR04[], float sigmaIEIE[], float eta[], float vx[], float vy[], float vz[],
    float minEt, bool restrictEta, float maxEta, bool restrictVertex, float maxVertexDist){
    bool hasOtherPhoAnswer = false;
    for(int x=0; x<PhoN && x<200; x++){
        if(x==candIdx) continue;
        float d = sqrt( (vx[candIdx]-vx[x])*(vx[candIdx]-vx[x]) + (vy[candIdx]-vy[x])*(vy[candIdx]-vy[x]) + (vz[candIdx]-vz[x])*(vz[candIdx]-vz[x]));
        if(restrictVertex && d>maxVertexDist) continue;
        if(restrictEta && fabs(eta[x])>maxEta) continue;
        if(phoID==0 && passesTightPhotonID(Photon_et[x],Photon_HoE[x],Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],hcalTowerSumEtConeDR04[x], sigmaIEIE[x], eta[x], minEt, restrictEta, maxEta))
            hasOtherPhoAnswer = true;
        if(phoID==1 && passesLoosePhotonID(Photon_et[x],Photon_HoE[x],Photon_trkSumPtHollowConeDR04[x],Photon_ecalRecHitSumEtConeDR04[x],hcalTowerSumEtConeDR04[x], eta[x], minEt, restrictEta, maxEta))
            hasOtherPhoAnswer = true;
        if(phoID==2 && passesEMID(Photon_et[x],Photon_HoE[x],Photon_ecalRecHitSumEtConeDR04[x],Photon_hcalTowerSumEtConeDR04[x], eta[x], minEt, restrictEta, maxEta))
            hasOtherPhoAnswer = true;
    }
    return hasOtherPhoAnswer;
}
#endif // #ifdef makeHistogramsMonoPho_8_cxx
