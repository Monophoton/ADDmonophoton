/*
 g++ -c `root-config --cflags` zgSubsetCompiled.cpp
 g++ -o testCode `root-config --libs` zgSubsetCompiled.o
 ./testCode
*/

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <TDCacheFile.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <TH1F.h>
#include <TH2F.h>
using namespace std;

Float_t dRCalc(Float_t etaC, Float_t phiC, Float_t etaCrys, Float_t phiCrys, bool debug=false){
    //   Float_t ClusPhi = 0;
    //   if (phiC>0) ClusPhi = phiC;
    //   else ClusPhi = TMath::Pi()*2. + phiC;
    //   Float_t CrysPhi = 0;
    //   if (phiC>0) CrysPhi = phiCrys;
    //   else CrysPhi = TMath::Pi()*2. + phiCrys;
    Float_t ClusPhi = phiC;
    Float_t CrysPhi = phiCrys;
    
    Float_t dphi = fabs(CrysPhi -ClusPhi);
    Float_t tempdphi=dphi;
    if (dphi > TMath::Pi()) tempdphi = TMath::Pi()*2. - dphi;
    dphi = tempdphi;
    if (debug){
        cout << "CrysPhi: " << CrysPhi << endl;
        cout << "ClusPhi: " << ClusPhi << endl;
        cout << "dphi: " << dphi << endl;
        cout << "deta: " << fabs(etaC-etaCrys) << endl;
    }
    Float_t deta = fabs(etaC - etaCrys);
    Float_t dR = sqrt(deta*deta + dphi*dphi);
    return dR;
    
}

void WPlotter(void){
    TChain *fChain = new TChain("myEvent");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runA_0.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runA_1.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runA_2.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runA_3.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runA_4.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_10.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_11.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_12.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_13.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_14.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_15.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_16.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_17.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_18.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_19.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_20.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_21.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_22.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_23.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_24.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_25.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_26.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_28.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_29.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_30.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_31.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_32.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_5.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_6.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_7.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_8.root");
fChain->Add("/uscms/home/miceli/nobackup/monophoton/subsets/W/WenuCand_skimOfNtuples_runB_9.root");






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
    
    // Set object pointer
    triggernames = 0;
    triggerprescales = 0;
    ifTriggerpassed = 0;
    
    Int_t           Vertex_n;
    Float_t         Vertex_x[10];   //[Vertex_n]
    Float_t         Vertex_y[10];   //[Vertex_n]
    Float_t         Vertex_z[10];   //[Vertex_n]
    Int_t           Vertex_tracksize[10];   //[Vertex_n]
    Int_t           Vertex_ndof[10];   //[Vertex_n]
    Float_t         Vertex_chi2[10];   //[Vertex_n]
    Float_t         Vertex_d0[10];   //[Vertex_n]
    Bool_t          Vertex_isFake[10];   //[Vertex_n]
    Bool_t          Scraping_isScrapingEvent;
    Int_t           Scraping_numOfTracks;
    Float_t         Scraping_fractionOfGoodTracks;
    Int_t           Track_n;
    Float_t         Track_px[100];   //[Track_n]
    Float_t         Track_py[100];   //[Track_n]
    Float_t         Track_pz[100];   //[Track_n]
    Float_t         Track_vx[100];   //[Track_n]
    Float_t         Track_vy[100];   //[Track_n]
    Float_t         Track_vz[100];   //[Track_n]
    Float_t         Track_pt[100];   //[Track_n]
    Float_t         Track_eta[100];   //[Track_n]
    Float_t         Track_phi[100];   //[Track_n]
    Int_t           Jet_n;
    Float_t         Jet_px[100];   //[Jet_n]
    Float_t         Jet_py[100];   //[Jet_n]
    Float_t         Jet_E[100];   //[Jet_n]
    Float_t         Jet_pz[100];   //[Jet_n]
    Float_t         Jet_vx[100];   //[Jet_n]
    Float_t         Jet_vy[100];   //[Jet_n]
    Float_t         Jet_vz[100];   //[Jet_n]
    Float_t         Jet_pt[100];   //[Jet_n]
    Float_t         Jet_eta[100];   //[Jet_n]
    Float_t         Jet_phi[100];   //[Jet_n]
    Float_t         Jet_emEnergyFraction[100];   //[Jet_n]
    Float_t         Jet_energyFractionHadronic[100];   //[Jet_n]
    Int_t           Jet_hitsInN90[100];   //[Jet_n]
    Int_t           Jet_n90Hits[100];   //[Jet_n]
    Int_t           Jet_nTowers[100];   //[Jet_n]
    Float_t         Jet_fHPD[100];   //[Jet_n]
    Float_t         Jet_fRBX[100];   //[Jet_n]
    Float_t         Jet_RHF[100];   //[Jet_n]
    Int_t           Electron_n;
    Float_t         Electron_px[100];   //[Electron_n]
    Float_t         Electron_py[100];   //[Electron_n]
    Float_t         Electron_pz[100];   //[Electron_n]
    Float_t         Electron_vx[100];   //[Electron_n]
    Float_t         Electron_vy[100];   //[Electron_n]
    Float_t         Electron_vz[100];   //[Electron_n]
    Float_t         Electron_pt[100];   //[Electron_n]
    Float_t         Electron_eta[100];   //[Electron_n]
    Float_t         Electron_phi[100];   //[Electron_n]
    Float_t         Electron_energy[100];   //[Electron_n]
    Float_t         Electron_charge[100];   //[Electron_n]
    Float_t         Electron_trkIso[100];   //[Electron_n]
    Float_t         Electron_ecalIso[100];   //[Electron_n]
    Float_t         Electron_hcalIso[100];   //[Electron_n]
    Float_t         Electron_SigmaIetaIeta[100];   //[Electron_n]
    Float_t         Electron_dEtaIn[100];   //[Electron_n]
    Float_t         Electron_dPhiIn[100];   //[Electron_n]
    Float_t         Electron_HoE[100];   //[Electron_n]
    Float_t         Electron_sc_energy[100];   //[Electron_n]
    Float_t         Electron_sc_eta[100];   //[Electron_n]
    Float_t         Electron_sc_phi[100];   //[Electron_n]
    Int_t           Muon_n;
    Float_t         Muon_px[100];   //[Muon_n]
    Float_t         Muon_py[100];   //[Muon_n]
    Float_t         Muon_pz[100];   //[Muon_n]
    Float_t         Muon_vx[100];   //[Muon_n]
    Float_t         Muon_vy[100];   //[Muon_n]
    Float_t         Muon_vz[100];   //[Muon_n]
    Float_t         Muon_pt[100];   //[Muon_n]
    Float_t         Muon_eta[100];   //[Muon_n]
    Float_t         Muon_phi[100];   //[Muon_n]
    Float_t         Muon_energy[100];   //[Muon_n]
    Float_t         Muon_charge[100];   //[Muon_n]
    Bool_t          Muon_isGlobalMuon[100];   //[Muon_n]
    Bool_t          Muon_isTrackerMuon[100];   //[Muon_n]
    Bool_t          Muon_isStandAloneMuon[100];   //[Muon_n]
    Bool_t          Muon_InnerTrack_isNonnull[100];   //[Muon_n]
    Bool_t          Muon_OuterTrack_isNonnull[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_x[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_y[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_z[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_px[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_py[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_InnerPoint_pz[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_x[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_y[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_z[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_px[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_py[100];   //[Muon_n]
    Float_t         Muon_OuterTrack_OuterPoint_pz[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_x[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_y[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_z[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_px[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_py[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_InnerPoint_pz[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_x[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_y[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_z[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_px[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_py[100];   //[Muon_n]
    Float_t         Muon_InnerTrack_OuterPoint_pz[100];   //[Muon_n]
    Int_t           CosmicMuon_n;
    Float_t         CosmicMuon_px[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_py[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_pz[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_pt[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_eta[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_phi[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_energy[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_charge[100];   //[CosmicMuon_n]
    Bool_t          CosmicMuon_isGlobalMuon[100];   //[CosmicMuon_n]
    Bool_t          CosmicMuon_isTrackerMuon[100];   //[CosmicMuon_n]
    Bool_t          CosmicMuon_isStandAloneMuon[100];   //[CosmicMuon_n]
    Bool_t          CosmicMuon_InnerTrack_isNonnull[100];   //[CosmicMuon_n]
    Bool_t          CosmicMuon_OuterTrack_isNonnull[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_x[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_y[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_z[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_px[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_py[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_InnerPoint_pz[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_x[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_y[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_z[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_px[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_py[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_OuterTrack_OuterPoint_pz[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_x[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_y[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_z[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_px[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_py[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_InnerPoint_pz[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_x[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_y[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_z[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_px[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_py[100];   //[CosmicMuon_n]
    Float_t         CosmicMuon_InnerTrack_OuterPoint_pz[100];   //[CosmicMuon_n]
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
    Int_t           Photon_n;
    Float_t         Photon_E[100];   //[Photon_n]
    Float_t         Photon_pt[100];   //[Photon_n]
    Float_t         Photon_eta[100];   //[Photon_n]
    Float_t         Photon_phi[100];   //[Photon_n]
    Float_t         Photon_theta[100];   //[Photon_n]
    Float_t         Photon_et[100];   //[Photon_n]
    Float_t         Photon_swissCross[100];   //[Photon_n]
    Float_t         Photonr9[100];   //[Photon_n]
    Float_t         Photon_e1x5[100];   //[Photon_n]
    Float_t         Photon_e2x5[100];   //[Photon_n]
    Float_t         Photon_e3x3[100];   //[Photon_n]
    Float_t         Photon_e5x5[100];   //[Photon_n]
    Float_t         Photon_r1x5[100];   //[Photon_n]
    Float_t         Photon_r2x5[100];   //[Photon_n]
    Float_t         Photon_maxEnergyXtal[100];   //[Photon_n]
    Float_t         Photon_SigmaEtaEta[100];   //[Photon_n]
    Float_t         Photon_SigmaIetaIeta[100];   //[Photon_n]
    Float_t         Photon_SigmaEtaPhi[100];   //[Photon_n]
    Float_t         Photon_SigmaIetaIphi[100];   //[Photon_n]
    Float_t         Photon_SigmaPhiPhi[100];   //[Photon_n]
    Float_t         Photon_SigmaIphiIphi[100];   //[Photon_n]
    Float_t         Photon_Roundness[100];   //[Photon_n]
    Float_t         Photon_Angle[100];   //[Photon_n]
    Float_t         Photon_ecalRecHitSumEtConeDR03[100];   //[Photon_n]
    Float_t         Photon_hcalTowerSumEtConeDR03[100];   //[Photon_n]
    Float_t         Photon_trkSumPtSolidConeDR03[100];   //[Photon_n]
    Float_t         Photon_trkSumPtHollowConeDR03[100];   //[Photon_n]
    Int_t           Photon_nTrkSolidConeDR03[100];   //[Photon_n]
    Int_t           Photon_nTrkHollowConeDR03[100];   //[Photon_n]
    Float_t         Photon_hcalDepth1TowerSumEtConeDR03[100];   //[Photon_n]
    Float_t         Photon_hcalDepth2TowerSumEtConeDR03[100];   //[Photon_n]
    Float_t         Photon_ecalRecHitSumEtConeDR04[100];   //[Photon_n]
    Float_t         Photon_hcalTowerSumEtConeDR04[100];   //[Photon_n]
    Float_t         Photon_trkSumPtSolidConeDR04[100];   //[Photon_n]
    Float_t         Photon_trkSumPtHollowConeDR04[100];   //[Photon_n]
    Int_t           Photon_nTrkSolidConeDR04[100];   //[Photon_n]
    Int_t           Photon_nTrkHollowConeDR04[100];   //[Photon_n]
    Float_t         Photon_hcalDepth1TowerSumEtConeDR04[100];   //[Photon_n]
    Float_t         Photon_hcalDepth2TowerSumEtConeDR04[100];   //[Photon_n]
    Bool_t          Photon_hasPixelSeed[100];   //[Photon_n]
    Bool_t          Photon_isEB[100];   //[Photon_n]
    Bool_t          Photon_isEE[100];   //[Photon_n]
    Bool_t          Photon_isEBGap[100];   //[Photon_n]
    Bool_t          Photon_isEEGap[100];   //[Photon_n]
    Bool_t          Photon_isEBEEGap[100];   //[Photon_n]
    Float_t         Photon_HoE[100];   //[Photon_n]
    Float_t         Photon_px[100];   //[Photon_n]
    Float_t         Photon_py[100];   //[Photon_n]
    Float_t         Photon_pz[100];   //[Photon_n]
    Float_t         Photon_vx[100];   //[Photon_n]
    Float_t         Photon_vy[100];   //[Photon_n]
    Float_t         Photon_vz[100];   //[Photon_n]
    Int_t           Photon_no_of_basic_clusters[100];   //[Photon_n]
    Float_t         Photon_sc_energy[100];   //[Photon_n]
    Float_t         Photon_sc_eta[100];   //[Photon_n]
    Float_t         Photon_sc_phi[100];   //[Photon_n]
    Float_t         Photon_sc_x[100];   //[Photon_n]
    Float_t         Photon_sc_y[100];   //[Photon_n]
    Float_t         Photon_sc_z[100];   //[Photon_n]
    Float_t         Photon_etaWidth[100];   //[Photon_n]
    Float_t         Photon_phiWidth[100];   //[Photon_n]
    Float_t         Photon_sc_et[100];   //[Photon_n]
    Float_t         matchphotonE[100];   //[Photon_n]
    Float_t         matchphotonpt[100];   //[Photon_n]
    Float_t         matchphotoneta[100];   //[Photon_n]
    Float_t         matchphotonphi[100];   //[Photon_n]
    Float_t         matchphotonpx[100];   //[Photon_n]
    Float_t         matchphotonpy[100];   //[Photon_n]
    Float_t         matchphotonpz[100];   //[Photon_n]
    Bool_t          ismatchedphoton[100];   //[Photon_n]
    Bool_t          Photon_hasConvTrk[100];   //[Photon_n]
    Int_t           Photon_ntracks[100];   //[Photon_n]
    Bool_t          Photon_isconverted[100];   //[Photon_n]
    Float_t         Photon_pairInvmass[100];   //[Photon_n]
    Float_t         Photon_pairCotThetaSeperation[100];   //[Photon_n]
    Float_t         Photon_pairmomentumX[100];   //[Photon_n]
    Float_t         Photon_pairmomentumY[100];   //[Photon_n]
    Float_t         Photon_pairmomentumZ[100];   //[Photon_n]
    Float_t         Photon_EoverP[100];   //[Photon_n]
    Float_t         Photon_ConvVx[100];   //[Photon_n]
    Float_t         Photon_ConvVy[100];   //[Photon_n]
    Float_t         Photon_ConvVz[100];   //[Photon_n]
    Float_t         Photon_ZOfPrimaryVertex[100];   //[Photon_n]
    Float_t         Photon_distOfMinimumApproach[100];   //[Photon_n]
    Float_t         Photon_dPhiTracksAtVtx[100];   //[Photon_n]
    Float_t         Photon_dPhiTracksAtEcal[100];   //[Photon_n]
    Float_t         Photon_dEtaTracksAtEcal[100];   //[Photon_n]
    Int_t           Photon_ncrys[100];   //[Photon_n]
    Float_t         Photon_timing_xtal[100][100];   //[Photon_n]
    Float_t         Photon_timingavg_xtal[100];   //[Photon_n]
    Float_t         Photon_energy_xtal[100][100];   //[Photon_n]
    Int_t           Photon_ieta_xtalEB[100][100];   //[Photon_n]
    Int_t           Photon_iphi_xtalEB[100][100];   //[Photon_n]
    Float_t         Photon_s9[100];   //[Photon_n]
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
    Int_t           CSCseg_n;
    Float_t         CSCseg_time[10000];   //[CSCseg_n]
    Float_t         CSCseg_x[10000];   //[CSCseg_n]
    Float_t         CSCseg_y[10000];   //[CSCseg_n]
    Float_t         CSCseg_z[10000];   //[CSCseg_n]
    Float_t         CSCseg_phi[10000];   //[CSCseg_n]
    Float_t         CSCseg_DirectionX[10000];   //[CSCseg_n]
    Float_t         CSCseg_DirectionY[10000];   //[CSCseg_n]
    Float_t         CSCseg_DirectionZ[10000];   //[CSCseg_n]
    Int_t           RPChit_n;
    Float_t         RPChit_x[10000];   //[RPChit_n]
    Float_t         RPChit_y[10000];   //[RPChit_n]
    Float_t         RPChit_z[10000];   //[RPChit_n]
    Int_t           RPChit_BunchX[10000];   //[RPChit_n]
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
    
    fChain->SetBranchAddress("nevents", &nevents);
    fChain->SetBranchAddress("run", &run);
    fChain->SetBranchAddress("event", &event);
    fChain->SetBranchAddress("luminosityBlock", &luminosityBlock);
    fChain->SetBranchAddress("beamCrossing", &beamCrossing);
    fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1);
    fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2);
    fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi);
    fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr);
    fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi);
    fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr);
    fChain->SetBranchAddress("ntriggers", &ntriggers);
    fChain->SetBranchAddress("triggernames", &triggernames);
    fChain->SetBranchAddress("triggerprescales", &triggerprescales);
    fChain->SetBranchAddress("ifTriggerpassed", &ifTriggerpassed);
    fChain->SetBranchAddress("Vertex_n", &Vertex_n);
    fChain->SetBranchAddress("Vertex_x", Vertex_x);
    fChain->SetBranchAddress("Vertex_y", Vertex_y);
    fChain->SetBranchAddress("Vertex_z", Vertex_z);
    fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize);
    fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof);
    fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2);
    fChain->SetBranchAddress("Vertex_d0", Vertex_d0);
    fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake);
    fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent);
    fChain->SetBranchAddress("Scraping_numOfTracks", &Scraping_numOfTracks);
    fChain->SetBranchAddress("Scraping_fractionOfGoodTracks", &Scraping_fractionOfGoodTracks);
    fChain->SetBranchAddress("Track_n", &Track_n);
    fChain->SetBranchAddress("Track_px", Track_px);
    fChain->SetBranchAddress("Track_py", Track_py);
    fChain->SetBranchAddress("Track_pz", Track_pz);
    fChain->SetBranchAddress("Track_vx", Track_vx);
    fChain->SetBranchAddress("Track_vy", Track_vy);
    fChain->SetBranchAddress("Track_vz", Track_vz);
    fChain->SetBranchAddress("Track_pt", Track_pt);
    fChain->SetBranchAddress("Track_eta", Track_eta);
    fChain->SetBranchAddress("Track_phi", Track_phi);
    fChain->SetBranchAddress("Jet_n", &Jet_n);
    fChain->SetBranchAddress("Jet_px", Jet_px);
    fChain->SetBranchAddress("Jet_py", Jet_py);
    fChain->SetBranchAddress("Jet_E", Jet_E);
    fChain->SetBranchAddress("Jet_pz", Jet_pz);
    fChain->SetBranchAddress("Jet_vx", Jet_vx);
    fChain->SetBranchAddress("Jet_vy", Jet_vy);
    fChain->SetBranchAddress("Jet_vz", Jet_vz);
    fChain->SetBranchAddress("Jet_pt", Jet_pt);
    fChain->SetBranchAddress("Jet_eta", Jet_eta);
    fChain->SetBranchAddress("Jet_phi", Jet_phi);
    fChain->SetBranchAddress("Jet_emEnergyFraction", Jet_emEnergyFraction);
    fChain->SetBranchAddress("Jet_energyFractionHadronic", Jet_energyFractionHadronic);
    fChain->SetBranchAddress("Jet_hitsInN90", Jet_hitsInN90);
    fChain->SetBranchAddress("Jet_n90Hits", Jet_n90Hits);
    fChain->SetBranchAddress("Jet_nTowers", Jet_nTowers);
    fChain->SetBranchAddress("Jet_fHPD", Jet_fHPD);
    fChain->SetBranchAddress("Jet_fRBX", Jet_fRBX);
    fChain->SetBranchAddress("Jet_RHF", Jet_RHF);
    fChain->SetBranchAddress("Electron_n", &Electron_n);
    fChain->SetBranchAddress("Electron_px", Electron_px);
    fChain->SetBranchAddress("Electron_py", Electron_py);
    fChain->SetBranchAddress("Electron_pz", Electron_pz);
    fChain->SetBranchAddress("Electron_vx", Electron_vx);
    fChain->SetBranchAddress("Electron_vy", Electron_vy);
    fChain->SetBranchAddress("Electron_vz", Electron_vz);
    fChain->SetBranchAddress("Electron_pt", Electron_pt);
    fChain->SetBranchAddress("Electron_eta", Electron_eta);
    fChain->SetBranchAddress("Electron_phi", Electron_phi);
    fChain->SetBranchAddress("Electron_energy", Electron_energy);
    fChain->SetBranchAddress("Electron_charge", Electron_charge);
    fChain->SetBranchAddress("Electron_trkIso", Electron_trkIso);
    fChain->SetBranchAddress("Electron_ecalIso", Electron_ecalIso);
    fChain->SetBranchAddress("Electron_hcalIso", Electron_hcalIso);
    fChain->SetBranchAddress("Electron_SigmaIetaIeta", Electron_SigmaIetaIeta);
    fChain->SetBranchAddress("Electron_dEtaIn", Electron_dEtaIn);
    fChain->SetBranchAddress("Electron_dPhiIn", Electron_dPhiIn);
    fChain->SetBranchAddress("Electron_HoE", Electron_HoE);
    fChain->SetBranchAddress("Electron_sc_energy", Electron_sc_energy);
    fChain->SetBranchAddress("Electron_sc_eta", Electron_sc_eta);
    fChain->SetBranchAddress("Electron_sc_phi", Electron_sc_phi);
    fChain->SetBranchAddress("Muon_n", &Muon_n);
    fChain->SetBranchAddress("Muon_px", Muon_px);
    fChain->SetBranchAddress("Muon_py", Muon_py);
    fChain->SetBranchAddress("Muon_pz", Muon_pz);
    fChain->SetBranchAddress("Muon_vx", Muon_vx);
    fChain->SetBranchAddress("Muon_vy", Muon_vy);
    fChain->SetBranchAddress("Muon_vz", Muon_vz);
    fChain->SetBranchAddress("Muon_pt", Muon_pt);
    fChain->SetBranchAddress("Muon_eta", Muon_eta);
    fChain->SetBranchAddress("Muon_phi", Muon_phi);
    fChain->SetBranchAddress("Muon_energy", Muon_energy);
    fChain->SetBranchAddress("Muon_charge", Muon_charge);
    fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon);
    fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon);
    fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon);
    fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_x", Muon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_y", Muon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_z", Muon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_px", Muon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_py", Muon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_pz", Muon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_x", Muon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_y", Muon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_z", Muon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_px", Muon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_py", Muon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_pz", Muon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_x", Muon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_y", Muon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_z", Muon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_px", Muon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_py", Muon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_pz", Muon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_x", Muon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_y", Muon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_z", Muon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_px", Muon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_py", Muon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_pz", Muon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_n", &CosmicMuon_n);
    fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px);
    fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py);
    fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz);
    fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt);
    fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta);
    fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi);
    fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy);
    fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge);
    fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon);
    fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon);
    fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_x", CosmicMuon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_y", CosmicMuon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_z", CosmicMuon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_px", CosmicMuon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_py", CosmicMuon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_pz", CosmicMuon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_x", CosmicMuon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_y", CosmicMuon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_z", CosmicMuon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_px", CosmicMuon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_py", CosmicMuon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_pz", CosmicMuon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_x", CosmicMuon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_y", CosmicMuon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_z", CosmicMuon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_px", CosmicMuon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_py", CosmicMuon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_pz", CosmicMuon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_x", CosmicMuon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_y", CosmicMuon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_z", CosmicMuon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_px", CosmicMuon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_py", CosmicMuon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_pz", CosmicMuon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Tau_n", &Tau_n);
    fChain->SetBranchAddress("Tau_px", Tau_px);
    fChain->SetBranchAddress("Tau_py", Tau_py);
    fChain->SetBranchAddress("Tau_pz", Tau_pz);
    fChain->SetBranchAddress("Tau_vx", Tau_vx);
    fChain->SetBranchAddress("Tau_vy", Tau_vy);
    fChain->SetBranchAddress("Tau_vz", Tau_vz);
    fChain->SetBranchAddress("Tau_pt", Tau_pt);
    fChain->SetBranchAddress("Tau_eta", Tau_eta);
    fChain->SetBranchAddress("Tau_phi", Tau_phi);
    fChain->SetBranchAddress("Tau_energy", Tau_energy);
    fChain->SetBranchAddress("Tau_charge", Tau_charge);
    fChain->SetBranchAddress("Photon_n", &Photon_n);
    fChain->SetBranchAddress("Photon_E", Photon_E);
    fChain->SetBranchAddress("Photon_pt", Photon_pt);
    fChain->SetBranchAddress("Photon_eta", Photon_eta);
    fChain->SetBranchAddress("Photon_phi", Photon_phi);
    fChain->SetBranchAddress("Photon_theta", Photon_theta);
    fChain->SetBranchAddress("Photon_et", Photon_et);
    fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross);
    fChain->SetBranchAddress("Photonr9", Photonr9);
    fChain->SetBranchAddress("Photon_e1x5", Photon_e1x5);
    fChain->SetBranchAddress("Photon_e2x5", Photon_e2x5);
    fChain->SetBranchAddress("Photon_e3x3", Photon_e3x3);
    fChain->SetBranchAddress("Photon_e5x5", Photon_e5x5);
    fChain->SetBranchAddress("Photon_r1x5", Photon_r1x5);
    fChain->SetBranchAddress("Photon_r2x5", Photon_r2x5);
    fChain->SetBranchAddress("Photon_maxEnergyXtal", Photon_maxEnergyXtal);
    fChain->SetBranchAddress("Photon_SigmaEtaEta", Photon_SigmaEtaEta);
    fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta);
    fChain->SetBranchAddress("Photon_SigmaEtaPhi", Photon_SigmaEtaPhi);
    fChain->SetBranchAddress("Photon_SigmaIetaIphi", Photon_SigmaIetaIphi);
    fChain->SetBranchAddress("Photon_SigmaPhiPhi", Photon_SigmaPhiPhi);
    fChain->SetBranchAddress("Photon_SigmaIphiIphi", Photon_SigmaIphiIphi);
    fChain->SetBranchAddress("Photon_Roundness", Photon_Roundness);
    fChain->SetBranchAddress("Photon_Angle", Photon_Angle);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", Photon_ecalRecHitSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", Photon_hcalTowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", Photon_trkSumPtSolidConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", Photon_nTrkSolidConeDR03);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", Photon_nTrkHollowConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", Photon_hcalDepth1TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", Photon_hcalDepth2TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", Photon_ecalRecHitSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", Photon_hcalTowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", Photon_trkSumPtHollowConeDR04);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", Photon_nTrkSolidConeDR04);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", Photon_nTrkHollowConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", Photon_hcalDepth1TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", Photon_hcalDepth2TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed);
    fChain->SetBranchAddress("Photon_isEB", Photon_isEB);
    fChain->SetBranchAddress("Photon_isEE", Photon_isEE);
    fChain->SetBranchAddress("Photon_isEBGap", Photon_isEBGap);
    fChain->SetBranchAddress("Photon_isEEGap", Photon_isEEGap);
    fChain->SetBranchAddress("Photon_isEBEEGap", Photon_isEBEEGap);
    fChain->SetBranchAddress("Photon_HoE", Photon_HoE);
    fChain->SetBranchAddress("Photon_px", Photon_px);
    fChain->SetBranchAddress("Photon_py", Photon_py);
    fChain->SetBranchAddress("Photon_pz", Photon_pz);
    fChain->SetBranchAddress("Photon_vx", Photon_vx);
    fChain->SetBranchAddress("Photon_vy", Photon_vy);
    fChain->SetBranchAddress("Photon_vz", Photon_vz);
    fChain->SetBranchAddress("Photon_no_of_basic_clusters", Photon_no_of_basic_clusters);
    fChain->SetBranchAddress("Photon_sc_energy", Photon_sc_energy);
    fChain->SetBranchAddress("Photon_sc_eta", Photon_sc_eta);
    fChain->SetBranchAddress("Photon_sc_phi", Photon_sc_phi);
    fChain->SetBranchAddress("Photon_sc_x", Photon_sc_x);
    fChain->SetBranchAddress("Photon_sc_y", Photon_sc_y);
    fChain->SetBranchAddress("Photon_sc_z", Photon_sc_z);
    fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth);
    fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth);
    fChain->SetBranchAddress("Photon_sc_et", Photon_sc_et);
    fChain->SetBranchAddress("matchphotonE", matchphotonE);
    fChain->SetBranchAddress("matchphotonpt", matchphotonpt);
    fChain->SetBranchAddress("matchphotoneta", matchphotoneta);
    fChain->SetBranchAddress("matchphotonphi", matchphotonphi);
    fChain->SetBranchAddress("matchphotonpx", matchphotonpx);
    fChain->SetBranchAddress("matchphotonpy", matchphotonpy);
    fChain->SetBranchAddress("matchphotonpz", matchphotonpz);
    fChain->SetBranchAddress("ismatchedphoton", ismatchedphoton);
    fChain->SetBranchAddress("Photon_hasConvTrk", Photon_hasConvTrk);
    fChain->SetBranchAddress("Photon_ntracks", Photon_ntracks);
    fChain->SetBranchAddress("Photon_isconverted", Photon_isconverted);
    fChain->SetBranchAddress("Photon_pairInvmass", Photon_pairInvmass);
    fChain->SetBranchAddress("Photon_pairCotThetaSeperation", Photon_pairCotThetaSeperation);
    fChain->SetBranchAddress("Photon_pairmomentumX", Photon_pairmomentumX);
    fChain->SetBranchAddress("Photon_pairmomentumY", Photon_pairmomentumY);
    fChain->SetBranchAddress("Photon_pairmomentumZ", Photon_pairmomentumZ);
    fChain->SetBranchAddress("Photon_EoverP", Photon_EoverP);
    fChain->SetBranchAddress("Photon_ConvVx", Photon_ConvVx);
    fChain->SetBranchAddress("Photon_ConvVy", Photon_ConvVy);
    fChain->SetBranchAddress("Photon_ConvVz", Photon_ConvVz);
    fChain->SetBranchAddress("Photon_ZOfPrimaryVertex", Photon_ZOfPrimaryVertex);
    fChain->SetBranchAddress("Photon_distOfMinimumApproach", Photon_distOfMinimumApproach);
    fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", Photon_dPhiTracksAtVtx);
    fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", Photon_dPhiTracksAtEcal);
    fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", Photon_dEtaTracksAtEcal);
    fChain->SetBranchAddress("Photon_ncrys", Photon_ncrys);
    fChain->SetBranchAddress("Photon_timing_xtal", Photon_timing_xtal);
    fChain->SetBranchAddress("Photon_timingavg_xtal", Photon_timingavg_xtal);
    fChain->SetBranchAddress("Photon_energy_xtal", Photon_energy_xtal);
    fChain->SetBranchAddress("Photon_ieta_xtalEB", Photon_ieta_xtalEB);
    fChain->SetBranchAddress("Photon_iphi_xtalEB", Photon_iphi_xtalEB);
    fChain->SetBranchAddress("Photon_s9", Photon_s9);
    fChain->SetBranchAddress("HERecHit_subset_n", &HERecHit_subset_n);
    fChain->SetBranchAddress("HERecHit_subset_detid", HERecHit_subset_detid);
    fChain->SetBranchAddress("HERecHit_subset_energy", HERecHit_subset_energy);
    fChain->SetBranchAddress("HERecHit_subset_time", HERecHit_subset_time);
    fChain->SetBranchAddress("HERecHit_subset_depth", HERecHit_subset_depth);
    fChain->SetBranchAddress("HERecHit_subset_phi", HERecHit_subset_phi);
    fChain->SetBranchAddress("HERecHit_subset_eta", HERecHit_subset_eta);
    fChain->SetBranchAddress("HERecHit_subset_x", HERecHit_subset_x);
    fChain->SetBranchAddress("HERecHit_subset_y", HERecHit_subset_y);
    fChain->SetBranchAddress("HERecHit_subset_z", HERecHit_subset_z);
    fChain->SetBranchAddress("CSCseg_n", &CSCseg_n);
    fChain->SetBranchAddress("CSCseg_time", CSCseg_time);
    fChain->SetBranchAddress("CSCseg_x", CSCseg_x);
    fChain->SetBranchAddress("CSCseg_y", CSCseg_y);
    fChain->SetBranchAddress("CSCseg_z", CSCseg_z);
    fChain->SetBranchAddress("CSCseg_phi", CSCseg_phi);
    fChain->SetBranchAddress("CSCseg_DirectionX", CSCseg_DirectionX);
    fChain->SetBranchAddress("CSCseg_DirectionY", CSCseg_DirectionY);
    fChain->SetBranchAddress("CSCseg_DirectionZ", CSCseg_DirectionZ);
    fChain->SetBranchAddress("RPChit_n", &RPChit_n);
    fChain->SetBranchAddress("RPChit_x", RPChit_x);
    fChain->SetBranchAddress("RPChit_y", RPChit_y);
    fChain->SetBranchAddress("RPChit_z", RPChit_z);
    fChain->SetBranchAddress("RPChit_BunchX", RPChit_BunchX);
    fChain->SetBranchAddress("CaloMetSigma", &CaloMetSigma);
    fChain->SetBranchAddress("CaloMetEz", &CaloMetEz);
    fChain->SetBranchAddress("CaloEtFractionHadronic", &CaloEtFractionHadronic);
    fChain->SetBranchAddress("CaloEmEtFraction", &CaloEmEtFraction);
    fChain->SetBranchAddress("CaloHadEtInHB", &CaloHadEtInHB);
    fChain->SetBranchAddress("CaloHadEtInHE", &CaloHadEtInHE);
    fChain->SetBranchAddress("CaloHadEtInHO", &CaloHadEtInHO);
    fChain->SetBranchAddress("CaloHadEtInHF", &CaloHadEtInHF);
    fChain->SetBranchAddress("CaloEmEtInEB", &CaloEmEtInEB);
    fChain->SetBranchAddress("CaloEmEtInEE", &CaloEmEtInEE);
    fChain->SetBranchAddress("CaloEmEtInHF", &CaloEmEtInHF);
    fChain->SetBranchAddress("CaloMaxEtInEmTowers", &CaloMaxEtInEmTowers);
    fChain->SetBranchAddress("CaloMaxEtInHadTowers", &CaloMaxEtInHadTowers);
    fChain->SetBranchAddress("CaloMetPt", CaloMetPt);
    fChain->SetBranchAddress("CaloMetPx", CaloMetPx);
    fChain->SetBranchAddress("CaloMetPy", CaloMetPy);
    fChain->SetBranchAddress("CaloMetPhi", CaloMetPhi);
    fChain->SetBranchAddress("CaloMetSumEt", CaloMetSumEt);
    fChain->SetBranchAddress("Delta_phi", &Delta_phi);
    fChain->SetBranchAddress("PFMetPt", PFMetPt);
    fChain->SetBranchAddress("PFMetPx", PFMetPx);
    fChain->SetBranchAddress("PFMetPy", PFMetPy);
    fChain->SetBranchAddress("PFMetPhi", PFMetPhi);
    fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt);
    fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF);
    fChain->SetBranchAddress("TCMetPt", TCMetPt);
    fChain->SetBranchAddress("TCMetPx", TCMetPx);
    fChain->SetBranchAddress("TCMetPy", TCMetPy);
    fChain->SetBranchAddress("TCMetPhi", TCMetPhi);
    fChain->SetBranchAddress("TCMetSumEt", TCMetSumEt);
    fChain->SetBranchAddress("Delta_phiTC", &Delta_phiTC);
    
    
    TFile *aa = new TFile("WPlots.root","RECREATE");
    TH1F *SCEta = new TH1F("SCEta","Supercluster Eta",300,-1.5,1.5);
    TH2F *SCEtaVsTime = new TH2F("SCEtaVsTime","Supercluster Eta Vs. Time",300,-1.5,1.5,100,-25,25);
    TH1F *NAddCrys = new TH1F("NAddCrys","Number of additional crystals E>1 GeV",100,0,100);
    TH1F *SigIetaIeta = new TH1F("SigIetaIeta","Width",75,0,0.025);
    TH1F *Timing =new TH1F("Timing","Time",400,-50,50);
    TH1F *ClusterTimeDiff = new TH1F("ClusterTimeDiff","LICTD",100,-50,50);
    TH2F *WidthVsTime = new TH2F("WidthVsTime","Width Versus Time",100,-50,50,150,0,0.03);
    TH2F *WidthVsTimeDiff = new TH2F("WidthVsTimeDiff","Width Versus TimeDiff",100,-50,50,150,0,0.03);
    TH1F *PhotonET = new TH1F("PhotonET","Photon transverse energy",1000,0,1000);
    TH1F *SeedEH = new TH1F("SeedEH","Seed Energy",200,0,200);
    TH2F *TimeVsTime = new TH2F("TimeVsTime","Time vs. Time",100,-50,50,100,-50,50);
    TH1F *TcMET = new TH1F("TcMET","TcMET",200,0,100);
    Float_t nentries = fChain->GetEntries();
    cout << "nentries: " << nentries << endl;
    for (int entry=0; entry<nentries; ++entry){
      fChain->GetEntry(entry);
      if (entry%10000==0) cout << "entry: " << entry << " / " << nentries << endl;
      if (Scraping_isScrapingEvent) continue;
      //Check trigger bits:
      Bool_t kTriggerPassed=kFALSE;
      vector <string> TriggerList;
      string bl1("HLT_Photon30_Cleaned_L1R");
      TriggerList.push_back(bl1);
      string bl2("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1");
      TriggerList.push_back(bl2);
      string bl3("HLT_Photon35_Isol_Cleaned_L1R_v1");
      TriggerList.push_back(bl3);
      string bl4("HLT_Photon30_L1R");
      TriggerList.push_back(bl4);
      string bl5("HLT_MET80_v1");
      //TriggerList.push_back(bl5);

      for (int jet=0;jet<ntriggers && jet<triggernames->size();++jet){
	string trig = (*triggernames)[jet];
	if ( (*ifTriggerpassed)[jet]==1){
	  //cout << "trig: " << trig << endl;
	  for (int jt=0;jt<TriggerList.size();++jt){
	    if (trig.compare(TriggerList[jt])==0)
	      kTriggerPassed=kTRUE;
	  }
	}
      }
      if (!kTriggerPassed) continue;

      Bool_t WriteEvent = kFALSE;
      int NumCandidateElectrons=0;
      int NumCandidatePhotons=0;
      Int_t IEle[100];
      Int_t nEle=0;
      for (int pho=0;pho<Photon_n;++pho){
	if(//tight barrel 30 photon
	   fabs(Photon_sc_eta[pho])<1.4442
	   && Photon_pt[pho] > 30.
	   && Photon_ecalRecHitSumEtConeDR04[pho] < 4.2+0.006*Photon_pt[pho]
	   && Photon_hcalTowerSumEtConeDR04[pho] < 2.2+0.0025*Photon_pt[pho]
	   && Photon_HoE[pho] < 0.05
	   && Photon_trkSumPtHollowConeDR04[pho] < 2.0+0.001*Photon_pt[pho]
	   && Photon_hasPixelSeed[pho]==1//pixel req
	   //&& TCMetPt[0]>30
	   ){
	  IEle[nEle] = pho;
	  nEle++;
	}
      }
      // if (TCMetPt[0]<20) continue;
       if (nEle>0){


       }

	//}
      for (int j=0;j<nEle;++j){

	Int_t EleA = IEle[j];
	Bool_t hasCosTag=kFALSE;
	for (int kk=0;kk<CosmicMuon_n;++kk){
	  Float_t cosRho = sqrt( CosmicMuon_OuterTrack_InnerPoint_x[kk]*CosmicMuon_OuterTrack_InnerPoint_x[kk]
			 + CosmicMuon_OuterTrack_InnerPoint_y[kk]*CosmicMuon_OuterTrack_InnerPoint_y[kk]);
	  Float_t cosDphi = dRCalc(0,Photon_sc_phi[EleA],0,CosmicMuon_phi[kk]);
	  if (cosRho>115 && cosRho<170
	      && cosDphi<0.2)
	    hasCosTag=kTRUE;
	}
	if (hasCosTag) continue;

	Float_t SeedTime=-999;
	Float_t SeedE = -999;
	Int_t crysIdx=-1;
	for (int k=0;k<Photon_ncrys[EleA]&&k<100;++k){
	  Float_t crysE = Photon_energy_xtal[EleA][k];
	  if (crysE > SeedE){
	    SeedE = crysE;
	    SeedTime = Photon_timing_xtal[EleA][k];
	    crysIdx = k;
	  }
	}

	if (fabs(SeedTime)>3) continue;
	Float_t LICTD =0;
	Int_t crysCrys=-1;
	Int_t crysThresh=0;
	for (int k=0;k<Photon_ncrys[EleA]&&k<100;++k){
	  if (crysIdx==k) continue;
	  Float_t crysE = Photon_energy_xtal[EleA][k];
	  if (crysE > 1.){
	    crysThresh++;
	    Float_t tdiff = Photon_timing_xtal[EleA][crysIdx] -
	      Photon_timing_xtal[EleA][k];
	    if (fabs(tdiff) > fabs(LICTD)){
	      LICTD = tdiff;
	      crysCrys=k;
	    }
	  }
	}
	Float_t SecondTime = Photon_timing_xtal[EleA][crysCrys];
	ClusterTimeDiff->Fill(LICTD);

	if (//kTRUE
	    fabs(LICTD)<5
	    && Photon_SigmaIetaIeta[EleA]>0.005 
	    //&& Photon_SigmaIphiIphi[EleA]>0.001
	    //&& (fabs(SeedTime)<3||fabs(SecondTime)<3) 
	    ){
	  if (SeedTime < 5)
	    SCEtaVsTime->Fill(Photon_sc_eta[EleA], SeedTime);
	  else SCEtaVsTime->Fill(Photon_sc_eta[EleA], SeedTime-25);
	  // 	  if (SeedTime>-1){
	  // 	    cout << "Run: " << run << " Event: " << event << endl;
	  // 	    for (int k=0;k<Photon_ncrys[EleA]&&k<100;++k){
	  // 	      cout << "IEta: " << Photon_ieta_xtalEB[EleA][k];
	  // 	      cout << " IPhi: " << Photon_iphi_xtalEB[EleA][k];
	  // 	      cout << " Energy: " << Photon_energy_xtal[EleA][k];
	  // 	      cout << " Time: " << Photon_timing_xtal[EleA][k] << endl;
	  // 	    }
	  // 	  }
	  
	  
	  TimeVsTime->Fill(SeedTime, Photon_timing_xtal[EleA][crysCrys]);
	  NAddCrys->Fill(crysThresh);
	  SeedEH->Fill(SeedE);
	  SigIetaIeta->Fill(Photon_SigmaIetaIeta[EleA]);
	  Timing->Fill(SeedTime);
	  WidthVsTime->Fill(SeedTime, Photon_SigmaIetaIeta[EleA]);
	  WidthVsTimeDiff->Fill(LICTD, Photon_SigmaIetaIeta[EleA]);
	  PhotonET->Fill(Photon_pt[EleA]);
	  TcMET->Fill(TCMetPt[0]);
	  SCEta->Fill(Photon_sc_eta[EleA]);
	}
      }
    }
    SCEta->Write();
    SCEtaVsTime->Write();
    TcMET->Write();
    TimeVsTime->Write();
    NAddCrys->Write();
    ClusterTimeDiff->Write();
    SeedEH->Write();
    Timing->Write();
    SigIetaIeta->Write();
    WidthVsTime->Write();
    WidthVsTimeDiff->Write();
    PhotonET->Write();
    aa->Close();
}


